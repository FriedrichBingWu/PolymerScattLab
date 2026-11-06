"""Loader tests exercise units, provenance, masks and archive read limits."""
import hashlib
from pathlib import Path
import zipfile

import numpy as np
import pytest

from polymer_scatt.io import DataFormatError, file_hash, iter_curves


def test_text_preserves_mask_and_unresolved_units(tmp_path):
    p = tmp_path / "sample_7_SAXS_0042.dat"
    p.write_text("# Transmission coefficient 0.72\n0.1 0 0\n0.2 12 2\n0.3 -1 0.2\n")
    c, = iter_curves(p)
    np.testing.assert_array_equal(c["intensity"], [0, 12, -1])
    np.testing.assert_array_equal(c["sigma"], [0, 2, .2])
    assert c["q_unit"] == "unknown"
    assert c["metadata"]["sample"] == 7
    assert c["metadata"]["frame"] == 42
    assert c["metadata"]["transmission"] == .72
    assert c["metadata"]["source_sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()
    assert file_hash(p) == c["metadata"]["source_sha256"]


@pytest.mark.parametrize("header,unit", [("q(A-1),I", "A^-1"),
    ("q [nm^-1],I", "nm^-1"), ("q [Å⁻¹],I", "A^-1"), ("q [1/nm],I", "nm^-1")])
def test_explicit_units_and_conflicts(tmp_path, header, unit):
    p = tmp_path / "curve.csv"
    p.write_text(header + "\n0.1,1\n0.2,2\n")
    c, = iter_curves(p)
    assert c["q_unit"] == unit
    assert c["sigma"] is None
    assert c["metadata"]["q_unit_origin"] == "header"
    wrong = "nm^-1" if unit == "A^-1" else "A^-1"
    with pytest.raises(DataFormatError, match="Header declares"):
        list(iter_curves(p, q_unit=wrong))


def test_user_unit_is_assertion_without_scaling(tmp_path):
    p = tmp_path / "curve.dat"
    p.write_text("1D-1 10 1\n2D-1 8 1\n")
    c, = iter_curves(p, q_unit="nm^-1")
    assert c["q_unit"] == "nm^-1"
    assert c["q"][0] == .1
    assert c["metadata"]["q_unit_origin"] == "user"


@pytest.mark.parametrize("data", ["0.1 10 1\n0.2 X 1", "0.1 10 1\n0.2 8",
    "0.1,10,1,99", "q,I\n0.1,", "0.1 10\nfooter"])
def test_corrupt_numeric_rows_are_not_silently_dropped(tmp_path, data):
    p = tmp_path / "curve.dat"
    p.write_text(data)
    with pytest.raises(DataFormatError):
        list(iter_curves(p))


def test_archive_streams_natural_order_ignores_aggregate(tmp_path):
    p = tmp_path / "input.zip"
    with zipfile.ZipFile(p, "w") as z:
        for name in ["sample_10_SAXS_0010.dat", "sample_2_SAXS_0002.dat", "sample_2_SAXS_0010.dat"]:
            z.writestr("folder/" + name, "0.1 10 1\n0.2 8 1\n")
        z.writestr("allfiles.out", "duplicate intensity table")
    curves = list(iter_curves(p))
    assert len(curves) == 3
    assert [c["metadata"]["sample"] for c in curves] == [2, 2, 10]
    assert [c["metadata"]["frame"] for c in curves] == [2, 10, 10]
    assert not (tmp_path / "folder").exists()
    assert all("member_crc32" in c["metadata"] for c in curves)


def test_archive_budget_and_unsafe_names(tmp_path):
    p = tmp_path / "input.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("curve.dat", "0.1 10 1\n0.2 8 1\n")
    with pytest.raises(DataFormatError, match="byte budget"):
        list(iter_curves(p, max_member_bytes=10))
    with pytest.raises(DataFormatError, match="expanded-byte budget"):
        list(iter_curves(p, max_expanded_bytes=10))
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("../curve.dat", "0.1 10 1\n0.2 8 1\n")
    with pytest.raises(DataFormatError, match="unsafe"):
        list(iter_curves(p))


def test_unknown_unit_enum_rejected(tmp_path):
    p = tmp_path / "curve.dat"
    p.write_text("0.1 10 1")
    with pytest.raises(DataFormatError, match="q_unit"):
        list(iter_curves(p, "Angstrom"))


def _attached_xlsx():
    package = Path(__file__).resolve().parents[1]
    candidates = [package / "examples" / "MAXS1.xlsx", package / "sample_data" / "MAXS1.xlsx",
                  package.parent / "upload" / "MAXS1.xlsx"]
    return next((p for p in candidates if p.exists()), None)


def test_supplied_workbook_labels_units_and_values():
    p = _attached_xlsx()
    if p is None:
        pytest.skip("Original MAXS1.xlsx optional sample is not present")
    curves = list(iter_curves(p))
    assert [c["label"] for c in curves] == ["T0", "T24", "T48", "T72"]
    assert all(len(c["q"]) == 857 and c["q_unit"] == "A^-1" for c in curves)
    assert curves[0]["intensity"][0] == 215.94744
    assert curves[1]["intensity"][0] == 1383.07963
    assert all(c["sigma"] is None for c in curves)
    with pytest.raises(DataFormatError, match="Header declares"):
        list(iter_curves(p, "nm^-1"))
