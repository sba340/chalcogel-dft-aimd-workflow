import argparse
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs, neighbor_list
from itertools import combinations

def calculate_coordination_numbers_ase(atoms, atom_indices_1, atom_indices_2, cutoff, same_type=False):
    """
    Calculate coordination numbers using ASE neighbor list
    """
    cutoffs = np.zeros(len(atoms))
    cutoffs.fill(cutoff/2)

    nl = NeighborList(cutoffs, skin=0.0, bothways=True, self_interaction=False)
    nl.update(atoms)

    coord_numbers = []

    for i, atom_idx in enumerate(atom_indices_1):
        coord_count = 0
        indices, offsets = nl.get_neighbors(atom_idx)

        for neighbor_idx, offset in zip(indices, offsets):
            if neighbor_idx in atom_indices_2:
                if same_type and atom_idx >= neighbor_idx:
                    continue

                pos1 = atoms.positions[atom_idx]
                pos2 = atoms.positions[neighbor_idx] + np.dot(offset, atoms.get_cell())
                distance = np.linalg.norm(pos1 - pos2)

                if distance <= cutoff:
                    coord_count += 1

        coord_numbers.append(coord_count)

    return np.array(coord_numbers)

def get_all_distances_ase(atoms, atom_indices_1, atom_indices_2, max_distance=5.0, same_type=False):
    """
    Calculate all distances between two atom groups using ASE
    """
    distances = []

    cutoffs = np.zeros(len(atoms))
    cutoffs.fill(max_distance/2)

    nl = NeighborList(cutoffs, skin=0.0, bothways=True, self_interaction=False)
    nl.update(atoms)

    processed_pairs = set()

    for atom_idx in atom_indices_1:
        indices, offsets = nl.get_neighbors(atom_idx)

        for neighbor_idx, offset in zip(indices, offsets):
            if neighbor_idx in atom_indices_2:
                if same_type:
                    pair = tuple(sorted([atom_idx, neighbor_idx]))
                    if pair in processed_pairs:
                        continue
                    processed_pairs.add(pair)

                pos1 = atoms.positions[atom_idx]
                pos2 = atoms.positions[neighbor_idx] + np.dot(offset, atoms.get_cell())
                distance = np.linalg.norm(pos1 - pos2)

                if distance <= max_distance:
                    distances.append(distance)

    return np.array(distances)

parser = argparse.ArgumentParser(description="Bond distribution and coordination number analysis from a quenched CONTCAR.")
parser.add_argument("contcar", help="Path to the (quenched) CONTCAR/POSCAR file")
parser.add_argument("--out", default="Bond_Distribution_ASE.png", help="Output plot filename")
args = parser.parse_args()

# Load CONTCAR file using ASE
atoms = read(args.contcar)

# Get atom indices by type
atom_symbols = atoms.get_chemical_symbols()
mo_indices = [i for i, symbol in enumerate(atom_symbols) if symbol == 'Mo']
s_indices = [i for i, symbol in enumerate(atom_symbols) if symbol == 'S']
mn_indices = [i for i, symbol in enumerate(atom_symbols) if symbol == 'Mn']

# Define cutoffs
cutoffs = {'Mo-Mo': 3.0, 'Mo-S': 2.60, 'S-S': 2.10, 'S-Mn': 2.45}

# Calculate all bond distances using ASE
bonds = {}
bonds['Mo-Mo'] = get_all_distances_ase(atoms, mo_indices, mo_indices, max_distance=5.0, same_type=True)
bonds['Mo-S'] = get_all_distances_ase(atoms, mo_indices, s_indices, max_distance=5.0)
bonds['S-S'] = get_all_distances_ase(atoms, s_indices, s_indices, max_distance=5.0, same_type=True)
bonds['S-Mn'] = get_all_distances_ase(atoms, s_indices, mn_indices, max_distance=5.0)

# Filter bonds within cutoff distances
bonds_cutoff = {key: bonds[key][bonds[key] <= cutoffs[key]] for key in bonds}

print("="*70)
print("BOND DISTRIBUTION ANALYSIS (using ASE)")
print("="*70)

total_bonds = sum(len(bonds[key]) for key in bonds)
total_cutoff = sum(len(bonds_cutoff[key]) for key in bonds_cutoff)

print(f"\nTotal atoms: {len(atoms)} (Mo: {len(mo_indices)}, S: {len(s_indices)},Mn: {len(mn_indices)})")

print("\nOverall Bond Distribution (% of all bonds):")
for key in bonds:
    percent = (len(bonds[key]) / total_bonds) * 100
    print(f"{key} Bonds: {percent:.2f}% ({len(bonds[key])} bonds)")
print(f"Total bonds: {total_bonds}")

print("\nBonds Within Cutoff Distances:")
for key in bonds:
    cutoff_count = len(bonds_cutoff[key])
    percent = (cutoff_count / len(bonds[key])) * 100 if len(bonds[key]) > 0 else 0
    print(f"{key} (<={cutoffs[key]:.2f} A): {cutoff_count} bonds ({percent:.2f}% of {key} bonds)")

print("\nBonds Within Cutoff as % of Total Bonds:")
for key in bonds_cutoff:
    percent = (len(bonds_cutoff[key]) / total_bonds) * 100
    print(f"{key}: {percent:.2f}%")
print(f"Total bonds within cutoff: {total_cutoff} ({(total_cutoff/total_bonds)*100:.2f}% of all bonds)")

# Calculate coordination numbers using ASE
cn_data = {}
cn_data['Mo-Mo'] = calculate_coordination_numbers_ase(atoms, mo_indices, mo_indices, cutoffs['Mo-Mo'], same_type=True)
cn_data['Mo-S (from Mo)'] = calculate_coordination_numbers_ase(atoms, mo_indices, s_indices, cutoffs['Mo-S'])
cn_data['S-Mo (from S)'] = calculate_coordination_numbers_ase(atoms, s_indices, mo_indices, cutoffs['Mo-S'])
cn_data['S-S'] = calculate_coordination_numbers_ase(atoms, s_indices, s_indices, cutoffs['S-S'], same_type=True)
cn_data['S-Mn (from S)'] = calculate_coordination_numbers_ase(atoms, s_indices, mn_indices, cutoffs['S-Mn'])
cn_data['Mn-S (from Mn)'] = calculate_coordination_numbers_ase(atoms, mn_indices, s_indices, cutoffs['S-Mn'])

# Bond length statistics
print("\nBond Length Statistics (within cutoff distances):")
print("="*50)
bond_averages = {}
for key in bonds_cutoff:
    if len(bonds_cutoff[key]) > 0:
        avg, std = np.mean(bonds_cutoff[key]), np.std(bonds_cutoff[key])
        bond_averages[key] = avg
        print(f"{key} bonds (<={cutoffs[key]:.2f} A):")
        print(f"  Average BL: {avg:.2f} +/- {std:.2f} A")
        print(f"  Range: {np.min(bonds_cutoff[key]):.2f} - {np.max(bonds_cutoff[key]):.2f} A")
    else:
        bond_averages[key] = None

# Coordination number statistics
print("\nCoordination Number Statistics:")
print("="*40)
atom_counts = {'Mo': len(mo_indices), 'S': len(s_indices), 'Zn': len(mn_indices)}

for key, cn in cn_data.items():
    print(f"{key} Coordination Numbers:")
    if 'from' in key:
        print(f"  Average CN = {np.mean(cn):.2f} +/- {np.std(cn):.2f}")
        print(f"    Range: {np.min(cn)} - {np.max(cn)}")
    else:
        print(f"  Average CN: {np.mean(cn):.2f} +/- {np.std(cn):.2f}")
        print(f"  Range: {np.min(cn)} - {np.max(cn)}")
        atom_type = key.split('-')[0]
        print(f"  Total {atom_type} atoms: {atom_counts[atom_type]}")

print("\nDetailed Coordination Analysis:")
print("="*35)

cn_keys = ['Mo-Mo', 'Mo-S (from Mo)', 'S-Mo (from S)', 'S-S', 'S-Mn (from S)', 'Mn-S (from Mn)']
cn_labels = ['Mo-Mo CN', 'Mo-S CN (from Mo)', 'S-Mo CN (from S)', 'S-S CN', 'S-Mn CN (from S)', 'Mn-S CN (from Mn)']

for key, label in zip(cn_keys, cn_labels):
    if key in cn_data and len(cn_data[key]) > 0:
        unique, counts = np.unique(cn_data[key], return_counts=True)
        print(f"{label} distribution: {dict(zip(unique, counts))}")

print("="*70)

# Create plots in a single horizontal row (1x4)
plt.figure(figsize=(20, 5))
bond_keys = ['Mo-S', 'S-S', 'Mo-Mo', 'S-Mn']
colors = ['#4B0082', '#FFA500', '#008080', '#b68b7c']
xlims = [(2, 5), (1.8, 5), (2.5, 5), (1.8, 5)]
xtick_ranges = [np.arange(2, 5.1, 0.5), np.arange(2, 5.1, 0.5), np.arange(3, 5.0, 0.5), np.arange(2, 5.1, 0.5)]
text_x = [3.5, 3.5, 3.8, 3.5]

bins = np.arange(1, 5, 0.05)

for i, (key, color, xlim, xticks, txt_x) in enumerate(zip(bond_keys, colors, xlims, xtick_ranges, text_x)):
    plt.subplot(1, 4, i+1)
    plt.hist(bonds[key], bins=bins, color=color, alpha=0.8, edgecolor='black')

    if bond_averages[key] is not None:
        plt.axvline(x=bond_averages[key], color='red', linestyle='--', linewidth=2,
                    label=f'Average: {bond_averages[key]:.2f} A')
    plt.axvline(x=cutoffs[key], color='orange', linestyle='--', linewidth=2,
                label=f'Cutoff: {cutoffs[key]:.2f} A')

    plt.xlabel('Atomic Distance (A)', fontsize=14, fontweight='bold')
    plt.ylabel('Count', fontsize=14, fontweight='bold')
    plt.xlim(xlim)
    plt.ylim(0, 60)
    plt.xticks(xticks, fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)

    cutoff_percent = (len(bonds_cutoff[key]) / len(bonds[key])) * 100 if len(bonds[key]) > 0 else 0
    avg_text = f'Avg: {bond_averages[key]:.2f} A' if bond_averages[key] is not None else 'No bonds'
    plt.text(txt_x, 51, f'{key}\n{cutoff_percent:.2f}% <= {cutoffs[key]:.2f} A\n{avg_text}',
             fontsize=11, fontweight='normal', color='black', ha='center')
    plt.legend(fontsize=9)

plt.suptitle(r'a-$\mathbf{Mn_{0.5}Mo_3S_{13}}$', fontsize=22, fontweight='bold', y=0.95, x=0.5)
plt.tight_layout()
plt.savefig(args.out, bbox_inches='tight', dpi=300)
plt.show()

# Additional ASE-specific analysis
print("\nASE-Specific Information:")
print("="*30)
print(f"Cell parameters: {atoms.get_cell().lengths()}")
print(f"Cell angles: {atoms.get_cell().angles()}")
print(f"Volume: {atoms.get_volume():.2f} A^3")
print(f"Density: {len(atoms)/atoms.get_volume():.4f} atoms/A^3")
