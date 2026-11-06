# Build validation — 6 November 2026

**45 automated tests passed**, Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0.

- 13 numerical tests: synthetic Guinier and Debye recovery; rejection of power-law, absent-plateau and invalid-error data; Å/nm consistency; contiguous windows around masks; WAXS peak positions and d spacing; masked gaps; no spurious peaks on a flat baseline; Debye identifiability.
- 15 input tests: text formats and units, malformed rows, masks retained, streamed archive ordering/limits, provenance and supplied workbook parsing.
- 13 LLM tests: bounded summary payloads, offline behavior, provider payloads, mocked responses, endpoint validation and sanitized errors. Aggregate counts survive prompt construction.
- 4 local HTTP integration tests: demo fitting/offline explanation, host/origin/token/file isolation, unit conflicts/model scope, and upload/analysis of all four real workbook curves.

Additional executed checks:

- Full original archive: 4,635 curves parsed, each with 2,086 rows. Counts match the independent archive inventory.
- Full supplied-data run: 4,639 curves screened individually; 0 accepted apparent Guinier fits, 4,639 rejected, 0 invalid parser/core results. No frame averaging or arbitrary q-unit conversion.
- End-to-end CLI on the included exact synthetic example: one accepted fit, Rg = 5 nm.
- Local HTTP demo with seeded noise: accepted apparent Rg = 4.98536 nm for a true 5 nm input.
- Offline CLI interpretation: full-set and group counts preserved; no model execution.
- Diagnostic figure visually inspected for measured profile, low-q curvature, candidate overlay, residuals and units.

Limits of verification:

- A browser executable was unavailable, so the GUI's visual rendering and click workflow were not tested in a real browser. Server routes and JavaScript syntax were checked.
- No live Ollama or hosted language-model call was made. Provider tests use mocks; the user must configure an available model.
- No comparison against ATSAS/SasView outputs on calibrated reference specimens was performed. This is an auditable research prototype, not validated metrology software.
- This build does not validate original experiment calibration, subtraction, uncertainty interpretation or molecular-chain contrast.

Re-run `python -m pytest -q` from the project directory after installing pytest.
