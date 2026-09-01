# Bond & Coordination Number Analysis

Analyzes bond-length distributions and coordination numbers from a quenched
AIMD structure (CONTCAR/POSCAR) using ASE. Built for Mn-doped Mo3S13
chalcogel structures, computing Mo-Mo, Mo-S, S-S, and S-Mn bonding.

## Extracts

- All pairwise distances (Mo-Mo, Mo-S, S-S, S-Mn) with periodic boundary conditions
- Bond counts and percentages within defined cutoff distances
- Bond length statistics (average, std dev, range) within cutoffs
- Coordination numbers per atom type, with full CN distributions
- A 4-panel bond-length histogram figure

## Requirements
pip install -r requirements.txt


## Usage

```bash
python bond_coordination_analysis.py path/to/CONTCAR
python bond_coordination_analysis.py path/to/CONTCAR --out my_plot.png
```

Cutoff distances (Mo-Mo: 3.0 Å, Mo-S: 2.60 Å, S-S: 2.10 Å, S-Mn: 2.45 Å) are
defined for this Mn-doped Mo3S13 system and can be edited directly in the
script for other compositions.

## Example

Run on a quenched CONTCAR from an AIMD trajectory (`example_data/`) to
reproduce the bond distribution figure used in the chalcogel manuscript.