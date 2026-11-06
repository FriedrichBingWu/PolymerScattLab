"""Loopback-only graphical workbench. No web framework or internet CDN required."""
from __future__ import annotations

import argparse
from collections import OrderedDict
import json
import math
from pathlib import Path
import secrets
import tempfile
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse
import webbrowser
import zipfile

import numpy as np

MAX_UPLOAD = 25 * 1024 * 1024
MAX_POINTS = 400_000
MAX_CURVES = 100
MAX_SESSIONS = 4
ALLOWED_SUFFIXES = {".dat", ".txt", ".csv", ".xlsx"}
UNITS = {"unknown", "A^-1", "nm^-1"}


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


class WorkbenchServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address):
        if address[0] not in {"127.0.0.1", "localhost"}:
            raise ValueError("The workbench may bind only to loopback.")
        super().__init__(address, WorkbenchHandler)
        self.token = secrets.token_urlsafe(32)
        self.sessions = OrderedDict()
        self.lock = threading.RLock()

    def register(self, curves, name):
        token = secrets.token_urlsafe(16)
        with self.lock:
            self.sessions[token] = {"curves": curves, "name": name, "results": {}}
            while len(self.sessions) > MAX_SESSIONS:
                self.sessions.popitem(last=False)
        return {"session": token, "filename": name, "curves": [
            {"index": i, "label": c["label"], "q_unit": c["q_unit"],
             "n_points": len(c["q"])} for i, c in enumerate(curves)]}


class WorkbenchHandler(BaseHTTPRequestHandler):
    server: WorkbenchServer

    def log_message(self, fmt, *args):
        # Do not log submitted filenames, sample identifiers or model text.
        return

    def _allowed_host(self):
        host = self.headers.get("Host", "")
        return host in {f"127.0.0.1:{self.server.server_port}",
                        f"localhost:{self.server.server_port}"}

    def _send(self, code, payload, content_type="application/json"):
        data = (json.dumps(_clean(payload), allow_nan=False).encode("utf-8")
                if content_type == "application/json" else payload)
        self.send_response(code)
        self.send_header("Content-Type", content_type + "; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "DENY")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if not self._allowed_host():
            return self._send(403, {"error": "Use the printed localhost URL."})
        if urlparse(self.path).path != "/":
            return self._send(404, {"error": "Not found"})
        page = (Path(__file__).parent / "static" / "index.html").read_text(encoding="utf-8")
        self._send(200, page.replace("__CSRF_TOKEN__", self.server.token).encode(), "text/html")

    def do_POST(self):
        if not self._allowed_host():
            return self._send(403, {"error": "Invalid host"})
        origin = self.headers.get("Origin")
        allowed = {f"http://127.0.0.1:{self.server.server_port}",
                   f"http://localhost:{self.server.server_port}"}
        if origin and origin not in allowed:
            return self._send(403, {"error": "Cross-origin requests are blocked"})
        if not secrets.compare_digest(self.headers.get("X-Workbench-Token", ""), self.server.token):
            return self._send(403, {"error": "Reload the workbench page to start a session"})
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length < 0 or length > MAX_UPLOAD:
                return self._send(413, {"error": "Upload is larger than 25 MB. Use the batch CLI."})
            data = self.rfile.read(length)
            parsed = urlparse(self.path)
            if parsed.path == "/api/upload":
                result = self._upload(data, parse_qs(parsed.query))
            elif parsed.path == "/api/demo":
                result = self._demo()
            else:
                if length > 100_000:
                    return self._send(413, {"error": "Request is too large"})
                body = json.loads(data or b"{}")
                if not isinstance(body, dict):
                    raise ValueError("Expected a JSON object")
                if parsed.path == "/api/analyze":
                    result = self._analyze(body)
                elif parsed.path == "/api/explain":
                    result = self._explain(body)
                else:
                    return self._send(404, {"error": "Not found"})
            self._send(200, result)
        except (ValueError, KeyError, IndexError, TypeError, zipfile.BadZipFile) as exc:
            self._send(400, {"error": str(exc)[:350]})
        except Exception:
            self._send(500, {"error": "The operation failed. Check the input format or local model service."})

    def _upload(self, data, query):
        from .io import iter_curves
        name = Path(query.get("name", ["curve.dat"])[0].replace("\\", "/")).name
        suffix = Path(name).suffix.lower()
        unit = query.get("unit", ["unknown"])[0]
        if suffix not in ALLOWED_SUFFIXES:
            raise ValueError("Upload DAT, TXT, CSV or XLSX. ZIP archives are processed by the batch CLI.")
        if unit not in UNITS:
            raise ValueError("Unrecognized q unit")
        with tempfile.TemporaryDirectory(prefix="polymer_scatt_") as directory:
            path = Path(directory) / ("input" + suffix)
            path.write_bytes(data)
            if suffix == ".xlsx":
                with zipfile.ZipFile(path) as archive:
                    if sum(x.file_size for x in archive.infolist()) > 100 * 1024 * 1024:
                        raise ValueError("Expanded workbook exceeds 100 MB")
            curves, points = [], 0
            for curve in iter_curves(path, q_unit=unit):
                points += len(curve["q"])
                if len(curves) >= MAX_CURVES or points > MAX_POINTS:
                    raise ValueError("This dataset exceeds the interactive limit. Use the batch CLI.")
                curves.append(curve)
            if not curves:
                raise ValueError("No numeric curves found")
        return self.server.register(curves, name)

    def _demo(self):
        q = np.linspace(.006, .24, 140)
        rg = 5.0
        sigma = np.full(len(q), .4)
        rng = np.random.default_rng(7)
        intensity = 100 * np.exp(-(q * rg)**2 / 3) + rng.normal(0, sigma)
        curve = {"q": q, "intensity": intensity, "sigma": sigma,
                 "label": "Synthetic Guinier · true Rg 5 nm", "q_unit": "nm^-1",
                 "metadata": {"synthetic": True, "true_rg_nm": 5.0}}
        return self.server.register([curve], "Synthetic validation example")

    def _session(self, body):
        with self.server.lock:
            session = self.server.sessions.get(body.get("session", ""))
            if session is None:
                raise ValueError("Session expired. Upload the file again.")
            index = int(body.get("curve", 0))
            if index < 0 or index >= len(session["curves"]):
                raise ValueError("Unknown curve")
        return session, index

    def _analyze(self, body):
        from .core import analyze_curve
        session, index = self._session(body)
        curve = session["curves"][index]
        requested = body.get("q_unit", "unknown")
        if requested not in UNITS:
            raise ValueError("Unrecognized q unit")
        unit = curve["q_unit"]
        if unit != "unknown" and requested not in {unit, "unknown"}:
            raise ValueError("Selected units conflict with the file header; convert the source explicitly.")
        unit = requested if unit == "unknown" else unit
        limits = {}
        for key in ("qmin", "qmax"):
            value = body.get(key)
            limits[key] = None if value in ("", None) else float(value)
            if limits[key] is not None and not math.isfinite(limits[key]):
                raise ValueError("Fit limits must be finite")
        if limits["qmin"] is not None and limits["qmax"] is not None and limits["qmin"] >= limits["qmax"]:
            raise ValueError("q minimum must be smaller than q maximum")
        result = analyze_curve(curve["q"], curve["intensity"], curve.get("sigma"),
                               q_unit=unit, label=curve["label"], **limits)
        result["metadata"] = curve.get("metadata", {})
        with self.server.lock:
            session["results"][index] = result
        # Only the selected curve is sent to the browser. Remote LLM calls never use this object.
        return {"result": result, "curve": {"q": curve["q"], "intensity": curve["intensity"],
                "sigma": curve.get("sigma"), "label": curve["label"], "q_unit": unit}}

    def _explain(self, body):
        from .llm import interpret
        session, index = self._session(body)
        result = session["results"].get(index)
        if result is None:
            raise ValueError("Run the numerical analysis first")
        provider = body.get("provider", "offline")
        if provider not in {"offline", "ollama"}:
            raise ValueError("The graphical workbench supports offline or local Ollama only")
        endpoint = body.get("endpoint", "http://127.0.0.1:11434")
        if provider == "ollama":
            parsed = urlparse(endpoint)
            if parsed.scheme != "http" or parsed.hostname not in {"127.0.0.1", "localhost"} or parsed.username or parsed.password:
                raise ValueError("The graphical workbench connects only to a local Ollama service")
        return interpret(result, provider=provider, endpoint=endpoint,
                         model=body.get("model") or None,
                         question=str(body.get("question", ""))[:2000] or None)


def serve(host="127.0.0.1", port=8765, open_browser=True):
    server = WorkbenchServer((host, int(port)))
    url = f"http://127.0.0.1:{server.server_port}/"
    print(f"PolymerScattLab: {url}\nRuns on this computer only. Press Ctrl+C to stop.", flush=True)
    if open_browser:
        threading.Timer(.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def main():
    parser = argparse.ArgumentParser(description="Open the local PolymerScattLab workbench")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()
    serve(port=args.port, open_browser=not args.no_browser)


if __name__ == "__main__":
    main()
