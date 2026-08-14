# xanes_energy_alignment

Part of the chalcogel DFT/AIMD workflow repo. Two-step energy alignment for computed core-excited-state ("XCH")
absorption edge energies from VASP total energies, generalized to any
edge (Mo L3, S K, Ni L, N K, ...) and any reference/target compound.

## The problem this solves

A raw DFT total-energy difference between a core-hole excited-state
calculation and a ground-state calculation is not directly comparable to
an experimental edge energy:

1. **VASP's total energies are referenced to an arbitrary internal zero**
   (the cell's average electrostatic potential), not to anything physical.
   This has to be corrected before the number means anything in eV
   relative to a real absorption edge.
2. **DFT has a systematic absolute-energy error** for core levels
   (self-interaction error, incomplete relaxation, pseudopotential
   approximations), so even after step 1 there's still an offset from
   experiment that has to be calibrated out empirically.

This repo does both corrections:

- **delta1** (internal re-referencing) uses the Fermi level and a
  POTCAR core-level reference energy to remove VASP's arbitrary
  potential zero.
- **delta2** (empirical calibration) compares the delta1-corrected
  energy of a **reference compound with a known experimental edge
  energy** against that literature value.
- **total_shift = delta1 + delta2** is then a single number you can add
  to the raw excitation energy of *any other compound* computed with the
  same edge, functional, and POTCAR — because delta1/delta2 correct for
  setup-dependent systematic error, not compound-specific physics.

## Usage

```bash
pip install -r requirements.txt

python src/run_alignment.py --config examples/mo_l3_edge.yaml
```

Output:
```
Reference compound: (NH4)2MoS4 -- reference, exp. value from ...
Edge: Mo L3-edge
  delta1:      70.0580 eV
  delta2:      -52.8480 eV
  total_shift: 17.2100 eV
  aligned energy: 2522.4000 eV (should match Ref_value = 2522.4)

example target complex: raw = 2498.7000 eV -> aligned = 2515.9100 eV
```

### Using it for a different edge or compound

Copy `examples/mo_l3_edge.yaml` to a new file (e.g. `examples/s_k_edge.yaml`)
and edit:

- `edge` — just a label, for your own reference
- `calibration` — your reference compound's E_fermi, E_ch, E_gs, E_core,
  E_ex_raw, and the experimental Ref_value you're calibrating against
- `targets` — any number of other compounds computed the same way, each
  needs only a raw excitation energy (no experimental value required —
  that's the point of calibrating once and reusing total_shift)

Then run:
```bash
python src/run_alignment.py --config examples/s_k_edge.yaml
```

### Using it directly in Python

```python
from src.energy_alignment import calibrate, align_target

result = calibrate(E_fermi=..., E_ch=..., E_gs=..., E_core=..., E_ex_raw=..., Ref_value=...)
aligned = align_target(new_compound_raw_energy, result.total_shift)
```

## Where each input comes from

| Variable | Source |
|---|---|
| `E_fermi` | Fermi energy (eV) from the ground-state (SCF) OUTCAR |
| `E_ch` | Total energy (eV) of the core-hole excited-state (XCH) calculation |
| `E_gs` | Total energy (eV) of the ground-state (SCF) calculation |
| `E_core` | Core-level reference energy (eV) for the excited atom, from the POTCAR |
| `E_ex_raw` | Raw computed excitation energy (eV) from your XCH workflow |
| `Ref_value` | Experimental/literature edge energy (eV) for the reference compound |

## Notebook

`notebooks/Eads_alignment.ipynb` is the original single-cell notebook this
was developed in, kept for reference. `src/` is the documented, reusable
version.
