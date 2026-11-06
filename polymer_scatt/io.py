"""Auditable loaders for reduced 1-D scattering data; no unit guessing.

Arrays are returned exactly as parsed: masks, nonpositive intensities and NaNs
remain available to the analysis engine's explicit quality checks. ZIP members
are read in memory one at a time and are never extracted to the filesystem.
"""
from __future__ import annotations

import csv
import hashlib
import io
import re
import zipfile
from pathlib import Path, PurePosixPath
from typing import Iterator

import numpy as np


TEXT_SUFFIXES = {".dat", ".txt", ".csv", ".tsv"}
Q_UNITS = {"A^-1", "nm^-1", "unknown"}


class DataFormatError(ValueError):
    """Input is malformed, ambiguous, or exceeds a configured read budget."""


def file_hash(path: str | Path) -> str:
    """SHA-256 of the exact source bytes (including the archive, if applicable)."""
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _natural_key(value: str):
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", value)]


def _declared_unit(text: str) -> str:
    compact = re.sub(r"\s+", "", text).lower().replace("−", "-")
    compact = compact.replace("å", "a").replace("Å", "a")
    compact = compact.replace("angstrom", "a").replace("⁻¹", "^-1")
    # Require an explicit inverse-unit token. A plain column named q is unknown.
    units = set()
    if re.search(r"(?<![a-z])nm(?:\^?-1|\*\*-1)|1/nm", compact):
        units.add("nm^-1")
    if re.search(r"(?<![a-z])a(?:\^?-1|\*\*-1)|1/a(?![a-z])", compact):
        units.add("A^-1")
    if len(units) > 1:
        raise DataFormatError("Conflicting reciprocal-length units in header")
    return next(iter(units), "unknown")


def _resolve_unit(explicit: str, requested: str) -> str:
    if requested not in Q_UNITS:
        raise DataFormatError(f"q_unit must be one of {sorted(Q_UNITS)}")
    if explicit != "unknown" and requested != "unknown" and explicit != requested:
        raise DataFormatError(
            f"Header declares q in {explicit}; requested {requested}. "
            "Use the declared unit, then convert explicitly in the analysis engine."
        )
    return explicit if explicit != "unknown" else requested


def _identity(name: str) -> dict:
    sample = re.search(r"sample[\s_-]*(\d+)", name, re.I)
    frame = re.search(r"_SAXS_(\d+)(?:\.|$)", name, re.I)
    return {
        "sample": int(sample.group(1)) if sample else None,
        "frame": int(frame.group(1)) if frame else None,
    }


def _text_curve(raw: bytes, name: str, q_unit: str, metadata: dict) -> dict:
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    rows, headers = [], []
    ncols = None
    for line_number, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if not s:
            continue
        if s.startswith(("#", ";", "!")):
            headers.append(s)
            continue
        cells = next(csv.reader([s])) if "," in s else s.split()
        try:
            row = [float(v.replace("D", "E").replace("d", "e")) for v in cells]
        except (ValueError, AttributeError):
            # Descriptive preamble/column labels are permitted; corrupt numeric
            # rows and any unexpected footer after numeric data are not skipped.
            first_numeric = False
            try:
                float(cells[0])
                first_numeric = True
            except (ValueError, IndexError):
                pass
            if rows or first_numeric:
                raise DataFormatError(f"{name}: malformed numeric row at line {line_number}")
            headers.append(s)
            continue
        if len(row) not in (2, 3):
            raise DataFormatError(f"{name}: expected q,I[,sigma] at line {line_number}")
        if ncols is not None and len(row) != ncols:
            raise DataFormatError(f"{name}: inconsistent column count at line {line_number}")
        ncols = len(row)
        rows.append(row)
    if not rows:
        raise DataFormatError(f"{name}: no numeric q,I data found")
    header_text = "\n".join(headers)
    unit = _resolve_unit(_declared_unit(header_text), q_unit)
    a = np.asarray(rows, dtype=float)
    transmission = re.search(
        r"transmission\s+coefficient\s*[:=]?\s*([-+0-9.eE]+)", header_text, re.I
    )
    meta = dict(metadata, **_identity(name))
    meta.update({
        "header": header_text[:8192],
        "q_unit_origin": "header" if _declared_unit(header_text) != "unknown" else (
            "user" if unit != "unknown" else "unresolved"
        ),
        "sigma_origin": "third column; definition unconfirmed" if ncols == 3 else "not supplied",
    })
    if transmission:
        try:
            meta["transmission"] = float(transmission.group(1))
        except ValueError:
            raise DataFormatError(f"{name}: invalid transmission header")
    return {
        "label": Path(name).stem,
        "q": a[:, 0], "intensity": a[:, 1],
        "sigma": a[:, 2] if ncols == 3 else None,
        "q_unit": unit, "metadata": meta,
    }


def _xlsx_curves(path: Path, q_unit: str, source_hash: str) -> Iterator[dict]:
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        for sheet in workbook:
            rows = list(sheet.iter_rows(values_only=True))
            if not rows:
                continue
            # Discover q columns from the first 20 header rows. This supports
            # paired q/I curves and an explicit q/I/sigma three-column table.
            header_index = None
            q_columns = []
            for k, row in enumerate(rows[:20]):
                cols = [j for j, v in enumerate(row) if isinstance(v, str)
                        and re.match(r"^q(?:\s|[\[(]|$)", v.strip(), re.I)]
                if cols:
                    header_index, q_columns = k, cols
                    break
            if header_index is None:
                raise DataFormatError(f"{sheet.title}: no explicit q column header found")
            header = rows[header_index]
            for col in q_columns:
                if col + 1 >= len(header):
                    raise DataFormatError(f"{sheet.title}: q column has no intensity column")
                unit = _resolve_unit(_declared_unit(str(header[col])), q_unit)
                sigma_col = None
                if len(q_columns) == 1 and col + 2 < len(header):
                    token = str(header[col + 2] or "").lower().strip()
                    if re.match(r"^(sigma|error|err|di|uncertainty|sdi)(?:\s|[\[(]|$)", token):
                        sigma_col = col + 2
                label = f"{sheet.title}_curve_{col // 2 + 1}"
                points, started = [], False
                for k, row in enumerate(rows[header_index + 1:], header_index + 2):
                    q, intensity = row[col], row[col + 1]
                    if q is None and intensity is None:
                        continue
                    numeric = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
                    if numeric(q) and numeric(intensity):
                        started = True
                        point = [float(q), float(intensity)]
                        if sigma_col is not None:
                            if not numeric(row[sigma_col]):
                                raise DataFormatError(f"{sheet.title}: invalid sigma in row {k}")
                            point.append(float(row[sigma_col]))
                        points.append(point)
                    elif started or numeric(q) or numeric(intensity):
                        raise DataFormatError(f"{sheet.title}: incomplete numeric pair in row {k}")
                    elif q is None and isinstance(intensity, str) and intensity.strip():
                        label = intensity.strip()
                if not points:
                    raise DataFormatError(f"{sheet.title}: q column {col + 1} has no data")
                a = np.asarray(points, dtype=float)
                yield {
                    "label": label, "q": a[:, 0], "intensity": a[:, 1],
                    "sigma": a[:, 2] if sigma_col is not None else None,
                    "q_unit": unit,
                    "metadata": {
                        "source": str(path), "source_sha256": source_hash,
                        "sheet": sheet.title, "q_column": col + 1,
                        "intensity_column": col + 2, "sample": None, "frame": None,
                        "q_unit_origin": "header" if _declared_unit(str(header[col])) != "unknown"
                        else ("user" if unit != "unknown" else "unresolved"),
                        "sigma_origin": "labelled column" if sigma_col is not None else "not supplied",
                    },
                }
    finally:
        workbook.close()


def iter_curves(
    path: str | Path,
    q_unit: str = "unknown",
    *,
    max_expanded_bytes: int = 2 * 1024**3,
    max_member_bytes: int = 64 * 1024**2,
    max_members: int = 20_000,
    max_compression_ratio: float = 1000.0,
) -> Iterator[dict]:
    """Yield q/I[/sigma] curves from DAT, TXT, CSV, TSV, ZIP or XLSX.

    ``q_unit`` is a user assertion only when a file has no explicit units.
    Conflicting assertions raise rather than silently rescaling. ZIP `.out`
    aggregate files are ignored to prevent double counting their DAT sources.
    The total expanded-byte budget counts only supported text members.
    """
    _resolve_unit("unknown", q_unit)
    path = Path(path).expanduser().resolve()
    source_hash = file_hash(path)
    if path.suffix.lower() == ".xlsx":
        # XLSX is itself a ZIP: screen expansion before openpyxl reads it.
        with zipfile.ZipFile(path) as book_zip:
            infos = book_zip.infolist()
            if len(infos) > max_members or sum(i.file_size for i in infos) > max_expanded_bytes:
                raise DataFormatError("Workbook exceeds configured expansion budget")
            if any(i.file_size > max_member_bytes for i in infos):
                raise DataFormatError("Workbook member exceeds configured byte budget")
        yield from _xlsx_curves(path, q_unit, source_hash)
    elif path.suffix.lower() in TEXT_SUFFIXES:
        if path.stat().st_size > max_member_bytes:
            raise DataFormatError("Text file exceeds configured byte budget")
        yield _text_curve(path.read_bytes(), path.name, q_unit,
                          {"source": str(path), "source_sha256": source_hash})
    elif path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
            if len(infos) > max_members:
                raise DataFormatError("Archive exceeds configured member-count budget")
            selected = []
            total = 0
            for info in infos:
                if info.is_dir() or PurePosixPath(info.filename).suffix.lower() not in TEXT_SUFFIXES:
                    continue
                name = PurePosixPath(info.filename.replace("\\", "/"))
                if name.is_absolute() or ".." in name.parts or re.match(r"^[A-Za-z]:", str(name)):
                    raise DataFormatError("Archive contains an unsafe member name")
                if info.flag_bits & 1:
                    raise DataFormatError("Encrypted archive members are not supported")
                if info.file_size > max_member_bytes:
                    raise DataFormatError(f"Archive member exceeds byte budget: {info.filename}")
                if info.file_size / max(info.compress_size, 1) > max_compression_ratio:
                    raise DataFormatError(f"Archive member exceeds compression-ratio budget: {info.filename}")
                total += info.file_size
                if total > max_expanded_bytes:
                    raise DataFormatError("Archive exceeds configured expanded-byte budget")
                selected.append(info)
            if not selected:
                raise DataFormatError("Archive has no supported DAT/TXT/CSV/TSV curve members")
            if len({i.filename for i in selected}) != len(selected):
                raise DataFormatError("Archive contains duplicate curve member names")
            for info in sorted(selected, key=lambda i: _natural_key(i.filename)):
                raw = archive.read(info)
                yield _text_curve(raw, info.filename, q_unit, {
                    "source": str(path), "source_sha256": source_hash,
                    "member": info.filename, "member_crc32": f"{info.CRC:08x}",
                    "member_bytes": info.file_size,
                })
    else:
        raise DataFormatError(f"Unsupported file type: {path.suffix or '(none)'}")
