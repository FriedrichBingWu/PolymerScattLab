# Independent data audit

The two attachments contain 4,639 individual reduced one-dimensional curves: 4,635 `.dat` files plus four paired-column XLSX curves. No raw detector images, wavelength/calibration file, polymer identities, solvent/concentration metadata, or experimental time mapping was supplied.

## ZIP

Eight groups are present: Samples 1, 2, 3, 5, 6, 7, 10 and 11. Samples 4, 8 and 9 are absent. Samples 1, 2, 3, 6, 7 and 10 have 480 frames each; Sample 5 has 255; Sample 11 has 1,500. Seven `allfiles.out` files aggregate the same intensity columns and should not be counted or fitted as additional experiments.

Each DAT has 2,086 rows and three numeric columns. Only a transmission-coefficient header is present. The columns appear to be q, I and uncertainty, but neither q units nor uncertainty definition is explicitly recorded. The full numeric q domain is 0.002352273–9.81133. Leading and some trailing I=error=0 rows are masked values, not physical zero scattering. First nonzero q varies from 0.1528977 to 0.2046477. There are 32–47 masked rows per curve and no negative intensities. All nonzero intensity rows have positive third-column values. Background-subtraction is suggested by folder names, not independently verified by supplied blank data. Transmission values must not be applied again without confirming reduction history.

Representative first/middle/last frames were viewed independently. Most groups rise just beyond the masked low-q edge and then form a broad, nearly flat background-like profile; Sample 7 shows a strong decreasing low-q contribution with curvature. These observations do not establish single-polymer-chain scattering or a valid Guinier range. Unknown q units independently prevent physical-length conversion. The numeric domain is compatible with a nm-inverse interpretation but this is only a hypothesis, not a justified unit assignment.

## MAXS1.xlsx

Sheet1 has paired q/I columns labelled T0, T24, T48 and T72. All four contain 857 data pairs with q=0.01783–0.88988 Å⁻¹, explicitly declared by the `q(A-1)` headers. No intensity uncertainties are provided. T labels do not establish whether these are hours, minutes, temperatures or another condition. The last four T0 intensity values (q≥0.88682 Å⁻¹) abruptly collapse to 2.9–4.7×10⁻⁸ from intensities near one; this terminal anomaly should be flagged/excluded from wide-angle peak analysis.

The four curves have a steep low-q upturn without an observed Guinier plateau. As a deliberate diagnostic, fitting ln(I) against q² over 0.02088–0.03922 Å⁻¹ gives apparent Rg values of 7.888, 8.277, 6.858 and 7.219 nm, respectively. These are **invalid Guinier estimates**, because qmaxRg=3.094, 3.246, 2.690 and 2.831. A high R² (0.91–0.97 here) is insufficient. Moving the fitted range from 0.03–0.05 to 0.04–0.06 Å⁻¹ changes these apparent dimensions substantially. They must not appear as extracted molecular Rg values.

## Interpretation boundary

No validated molecular radius of gyration is supported by this independent audit. Any automated fit should preserve model assumptions, unit provenance, selected range, qRg validity, uncertainty definitions, and rejection reasons. Polymer-chain Rg additionally requires single-chain contrast and an appropriate dilute/structure-factor-corrected regime; a contrast-domain Guinier Rg is not automatically a chain Rg. Broad maxima may support descriptive q* / d=2π/q* measurements after q units and baseline are resolved; neither crystallinity nor crystal identities should be inferred from these files alone.

`audit.json` records the full inventory and illustrative rejected fits. `overview.png` shows four XLSX curves and one middle frame per ZIP group on log–log and Guinier axes. This audit did not average time-series frames.
