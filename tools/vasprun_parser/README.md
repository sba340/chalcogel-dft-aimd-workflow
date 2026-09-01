# vasprun Parser

Parses VASP `vasprun.xml` files using pymatgen and extracts the results
most commonly needed for post-processing — energy, forces, magnetic moments,
and total/element-projected density of states (DOS).

## Extracts

- Final total energy (eV)
- Final forces on each atom (eV/Angstrom)
- Final magnetic moments per atom (spin-polarized runs)
- Total DOS and element-projected DOS

## Requirements
pip install -r requirements.txt

## Usage

### Command line
```bash
python vasprun_parser.py vasprun.xml
python vasprun_parser.py vasprun.xml --dos --plot
python vasprun_parser.py vasprun.xml --csv results
```

### As a module
```python
from vasprun_parser import parse_vasprun

data = parse_vasprun("vasprun.xml")
data["final_energy"]    # float (eV)
data["forces"]          # DataFrame: atom, element, fx, fy, fz
data["magmom"]           # DataFrame: atom, element, magmom (or None)
data["total_dos"]        # DataFrame: energy, density_up, density_down
data["element_dos"]      # dict[element] -> DataFrame
```

## Example

Validated on a 144-atom MoS2 `vasprun.xml` (`example_data/`), correctly
extracting E_fermi and Mo/S-projected DOS.
