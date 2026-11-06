# PolymerScattLab 0.1.0

LLM-assisted SAXS/WAXS analysis for polymer researchers. Numerical fits, units and acceptance decisions are calculated by Python; an optional language model explains those results. This is a tested research prototype, not a model trained on the attached data or a validated commercial scattering package.

## Start here

1. Open `sample_results/report.html` in your browser for the completed analysis of all **4,639 supplied curves**.
2. See `sample_results/results.csv` for one row per curve, including rejected candidate values in explicitly labelled columns. Blank accepted-Rg cells mean that Rg was not supported by the screening.
3. Install Python 3.10 or later. In a terminal inside this folder:

```console
python -m venv .venv
```

On Windows:

```console
.venv\Scripts\activate
python -m pip install -r requirements.txt
python -m polymer_scatt serve
```

On macOS/Linux:

```console
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m polymer_scatt serve
```

The local workbench opens at `http://127.0.0.1:8765`. It runs only while the Python process is running. Upload a DAT/CSV/TXT/XLSX file, confirm q units, choose a curve and inspect the fit. ZIP files use the batch command below. The GUI has a noisy synthetic example with true Rg = 5 nm.

The interface is local; this download does not include a hosted website. It uses no external scripts or fonts. Interface HTTP behavior and calculations were tested; visual browser rendering could not be tested in the build environment because a browser executable was unavailable.

## Batch analysis

Analyze the included representative data:

```console
python -m polymer_scatt analyze examples/MAXS1.xlsx examples/representative_archive_frames.zip --output my_analysis
```

Analyze the full original archive, placing your uploaded file next to this README:

```console
python -m polymer_scatt analyze examples/MAXS1.xlsx "Sample 1-11.zip" --output full_analysis
```

The full original ZIP is not duplicated in the download. The included ZIP contains first/middle/last frames from each of the eight supplied groups (24 original frames). The report and CSV/JSON results cover the entire original archive, not just this subset.

Only set `--q-unit "nm^-1"` or `--q-unit "A^-1"` from confirmed metadata. Explicit headers are respected; conflicting overrides raise an error. To use different units for different input files, run them separately. Manual bounds `--qmin` and `--qmax` are in the native units and restrict the candidate search; they do not bypass checks.

Outputs: `report.html`, `results.csv`, `results.json`, `summary.json`, `manifest.json`, `llm_prompt.txt`, `interpretation.json`. SHA-256 input hashes, software versions, fit intervals, errors, rejection reasons and masking are retained. No original files are modified and no frame averaging is performed.

## Enable a language model

The default creates an offline deterministic summary and a grounded LLM prompt. **No separate language model is executed by default.** The sample report's numerical results were computed locally; its scientific review was prepared in ChatGPT. A local Ollama model or an OpenAI-compatible endpoint can provide new interpretations when explicitly selected.

```console
python -m polymer_scatt explain my_analysis/summary.json --provider ollama --model YOUR_INSTALLED_MODEL --endpoint http://127.0.0.1:11434 --output local_review.json
```

For a hosted provider, use `--provider openai-compatible --endpoint https://YOUR_PROVIDER/v1 --model YOUR_MODEL` and set `POLYMERSCATT_API_KEY` in your environment. The provider may charge for inference. Only compact numerical summaries and metadata are sent; raw curves are excluded. Read `README_LLM.md` for details. No API key or live model service is bundled. Mock HTTP tests verify request/response handling; a real provider/model call was not run during this build.

The LLM can discuss missing sample metadata, contrast, aggregation, background and alternative models. Numerical result files remain authoritative, even if a model's prose is wrong. This version uses a curated reference list; it does not perform live literature retrieval, fine-tuning, or autonomous model selection.

## What the supplied data show

`MAXS1.xlsx` contains T0, T24, T48 and T72, each with 857 points and q = 0.01783–0.88988 Å⁻¹. No intensity uncertainties are supplied. All four show a strong low-q upturn; a reliable Guinier Rg was not established. T labels are retained without assuming hours or temperature.

`Sample 1-11.zip` actually contains Samples 1, 2, 3, 5, 6, 7, 10 and 11. There are 4,635 DAT curves, plus seven aggregate `allfiles.out` files that are excluded to avoid double counting. Each DAT contains q, intensity and a third column treated provisionally as 1σ intensity error. The header records transmission but no q units. Zero intensity/zero error bins are treated as no-data masks, not a Guinier plateau. Physical size conversions from this archive remain unconfirmed.

The T0 workbook curve has an anomalous terminal collapse in its last four intensities (q ≥ 0.88682 Å⁻¹). Those points do not affect the low-q fit, but must be reviewed before wide-q interpretation. This anomaly is documented; raw data are not silently repaired.

The report does not equate a fitted scattering length with molecular-chain size. Polymer identity, state (solution/gel/film/bulk), concentration, solvent/contrast, preprocessing and frame meaning are still needed. A failed screening does not mean the polymer has no finite Rg; it means these data and this analysis do not establish it.

## Numerical methods and scope

- Guinier: weighted regression of ln I against q², with sigma(log I) approximated as sigma(I)/I; ordinary regression when errors are absent. Background is fixed at zero after user preprocessing. Nonpositive intensities and mask gaps break log-fit windows.
- Conservative screening: ≥10 points, qmax Rg ≤1, qmin Rg ≤0.5, qmax/qmin ≥1.5, meaningful decay, conditional relative Rg error ≤20%, endpoint sensitivity ≤15%, adequate linearity, reported-error consistency, residual curvature and log-transform signal checks. Exact criteria are in `core.py`. These are engineering heuristics, not universal physical thresholds.
- Candidate fits remain visible on rejection. Fit errors are conditional; they exclude background, calibration, q uncertainty, correlated detector bins, selection and model uncertainty.
- Kratky and residual plots are diagnostics. They do not assign polymer chemistry or conformation automatically.
- Optional Debye model: ideal monodisperse Gaussian-coil form factor, intensity-space fit, explicit range and optional constant background. Results are hypotheses with identifiability diagnostics, not substitutes for missing low-q evidence.
- Optional WAXS: descriptive peak candidates and d = 2π/q. No automatic phase identification, crystallinity fraction or Scherrer size. q units and a sensible WAXS interval must be confirmed.

To run an explicit-range Debye hypothesis (use your own justified interval):

```console
python -m polymer_scatt features YOUR_CURVE.dat --q-unit "nm^-1" --qmin 0.01 --qmax 1 --debye --output debye.json
```

To describe peaks in a confirmed WAXS curve:

```console
python -m polymer_scatt features YOUR_WAXS.dat --q-unit "A^-1" --qmin 0.8 --qmax 3 --peaks --output waxs_peaks.json
```

These commands are examples of syntax, not prescribed fit windows for the attachments. The archive's SAXS filenames and unknown q units do not establish WAXS coverage. MAXS1 does not cover a complete conventional wide-angle pattern.

The software accepts reduced 1D curves. It does not integrate 2D detector images, calibrate geometry, subtract a matched blank automatically, perform IFT/P(r), or fit semicrystalline correlation functions. These require experiment-specific inputs and can be added after metadata are established.

## Validation

```console
python -m pip install pytest
python -m pytest -q
```

Tests cover known Guinier and Debye sizes, q-unit conversion, power-law rejection, missing plateaus, masks and invalid data, archive limits and provenance, real workbook parsing, LLM payload boundaries, mocked provider responses and local server actions. See `VALIDATION.md` for this build's results.

## References

- https://www.sasview.org/docs/user/models/guinier.html
- https://www.sasview.org/docs/user/models/mono_gauss_coil.html
- https://www.sasview.org/docs/user/models/poly_gauss_coil.html
- https://bioxtas-raw.readthedocs.io/en/v2.1.4/saxs/saxs_guinier.html

Version 0.1.0 · Prepared for Bing Wu, 6 November 2026.
