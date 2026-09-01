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
    --sim-label "Fe0.5Mo3S13 PCF Simulation" \
    --exp-label "Fe0.5Mo3S13 PDF Experiment" \
    --annotate "1.98:S-S" "2.30:Mo/Fe-S" "2.7:Mo-Mo"
```

Optional flags: `--out` (output filename), `--offset` / `--shift-right`
(align sim curve to exp curve), `--xlim` / `--ylim` (axis ranges).

## Example

`example_data/` contains a sample simulated PCF and experimental PDF pair
for a Fe0.5Mo3S13 chalcogel, reproducing the comparison figure used in the
manuscript.
