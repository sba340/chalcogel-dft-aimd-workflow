# PCF vs PDF Comparison

Overlays a simulated pair correlation function (PCF, from AIMD melt-quench)
against an experimental pair distribution function (PDF) for visual
comparison, with bond-type annotations marking key coordination distances.

## Input format

Both `.dat` files are plain-text, two columns: `r` (Angstrom) and `g(r)`,
whitespace-delimited — the standard output format of common PCF/PDF tools.

## Requirements
pip install -r requirements.txt

## Usage

```bash
python pcf_pdf_comparison.py sim.dat exp.dat \
    --sim-label "Mn0.5Mo3S13 PCF Simulation" \
    --exp-label "Mn0.5Mo3S13 PDF Experiment" \
    --annotate "1.98:S-S" "2.30:Mo/Mn-S" "2.7:Mo-Mo"
```

Optional flags: `--out` (output filename), `--offset` / `--shift-right`
(align sim curve to exp curve), `--xlim` / `--ylim` (axis ranges).



## Example

`example_data/` contains a simulated PCF (`PCF_5ps.dat`) and experimental
PDF (`Mn_0_5_Mo3S13-GO_PDF.dat`) pair for a Mn0.5Mo3S13 chalcogel. Reproduce
the comparison figure with:

```bash
python pcf_pdf_comparison.py example_data/PCF_5ps.dat example_data/Mn_0_5_Mo3S13-GO_PDF.dat \
    --sim-label "Mn0.5Mo3S13 PCF Simulation" \
    --exp-label "Mn0.5Mo3S13 PDF Experiment" \
    --annotate "1.98:S-S" "2.30:Mo/Mn-S" "2.7:Mo-Mo"
```
