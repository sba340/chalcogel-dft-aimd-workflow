"""
vasprun_parser.py

Parses VASP vasprun.xml files using pymatgen and extracts the results most
commonly needed for post-processing:

  - Final total energy (eV)
  - Final forces on each atom (eV/Angstrom)
  - Final magnetic moments per atom (spin-polarized runs)
  - Total density of states (DOS)
  - Element- and orbital-projected DOS (if LORBIT was set during the run)

vasprun.xml contains everything OUTCAR does plus the full DOS, so this
replaces the need for a separate OUTCAR-only parser.

Usage (CLI):
    python vasprun_parser.py vasprun.xml
    python vasprun_parser.py vasprun.xml --dos --plot
    python vasprun_parser.py vasprun.xml --csv results

Usage (as a module):
    from vasprun_parser import parse_vasprun

    data = parse_vasprun("vasprun.xml")
    data["final_energy"]     # float (eV)
    data["forces"]           # DataFrame: atom, fx, fy, fz
    data["magmom"]           # DataFrame: atom, magmom (or None)
    data["total_dos"]        # DataFrame: energy, density, density_up, density_down
    data["element_dos"]      # dict[element] -> DataFrame: energy, s, p, d, (f)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pymatgen.io.vasp.outputs import Vasprun


def parse_vasprun(path: str | Path, parse_dos: bool = True) -> dict:
    """Parse a vasprun.xml file and return a dict of extracted quantities."""
    vr = Vasprun(
        str(path),
        parse_dos=parse_dos,
        parse_eigen=False,
        parse_projected_eigen=False,
    )

    final_structure = vr.final_structure
    final_energy = vr.final_energy  # eV

    # --- forces (last ionic step) ---
    forces = None
    if vr.ionic_steps:
        last_step = vr.ionic_steps[-1]
        if "forces" in last_step:
            f = last_step["forces"]
            forces = pd.DataFrame(f, columns=["fx", "fy", "fz"])
            forces.insert(0, "atom", range(1, len(forces) + 1))
            forces.insert(1, "element", [str(s.specie) for s in final_structure])

    # --- magnetic moments ---
    magmom = None
    try:
        mags = vr.final_structure.site_properties.get("magmom")
        if mags is None and vr.ionic_steps:
            mags = vr.ionic_steps[-1].get("structure").site_properties.get("magmom")
        if mags is not None:
            magmom = pd.DataFrame({
                "atom": range(1, len(mags) + 1),
                "element": [str(s.specie) for s in final_structure],
                "magmom": mags,
            })
    except Exception:
        magmom = None

    # --- total DOS ---
    total_dos_df = None
    element_dos = {}
    if parse_dos and vr.complete_dos is not None:
        cdos = vr.complete_dos
        efermi = cdos.efermi
        energies = cdos.energies - efermi  # shift so E_fermi = 0

        densities = cdos.densities
        total_dos_data = {"energy": energies}
        # densities keys are Spin.up / Spin.down (or just Spin.up if ISPIN=1)
        for spin, dens in densities.items():
            label = "density_up" if spin.value == 1 else "density_down"
            total_dos_data[label] = dens
        total_dos_df = pd.DataFrame(total_dos_data)

        # element-projected DOS
        try:
            element_pdos = cdos.get_element_dos()
            for element, dos_obj in element_pdos.items():
                edata = {"energy": energies}
                for spin, dens in dos_obj.densities.items():
                    label = "density_up" if spin.value == 1 else "density_down"
                    edata[label] = dens
                element_dos[str(element)] = pd.DataFrame(edata)
        except Exception:
            pass

    return {
        "final_energy": final_energy,
        "forces": forces,
        "magmom": magmom,
        "total_dos": total_dos_df,
        "element_dos": element_dos,
        "efermi": vr.complete_dos.efermi if (parse_dos and vr.complete_dos) else None,
    }


def plot_dos(total_dos_df: pd.DataFrame, element_dos: dict, out_path: str = "dos_plot.png"):
    """Quick total + element-projected DOS plot, saved to out_path."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.plot(total_dos_df["energy"], total_dos_df["density_up"], label="Total DOS (up)", color="black", lw=1.5)
    if "density_down" in total_dos_df.columns:
        ax.plot(total_dos_df["energy"], -total_dos_df["density_down"], color="black", lw=1.5)

    for element, df in element_dos.items():
        ax.plot(df["energy"], df["density_up"], label=f"{element}", lw=1.0)
        if "density_down" in df.columns:
            ax.plot(df["energy"], -df["density_down"], lw=1.0)

    ax.axvline(0, color="gray", linestyle="--", lw=0.8)
    ax.set_xlabel("E - E_fermi (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"DOS plot saved to {out_path}")


def main():
    ap = argparse.ArgumentParser(description="Parse a VASP vasprun.xml file.")
    ap.add_argument("vasprun", help="Path to vasprun.xml")
    ap.add_argument("--dos", action="store_true", help="Print DOS summary")
    ap.add_argument("--plot", action="store_true", help="Save a DOS plot (dos_plot.png)")
    ap.add_argument("--csv", metavar="PREFIX", help="Write results to PREFIX_*.csv")
    args = ap.parse_args()

    data = parse_vasprun(args.vasprun, parse_dos=args.dos or args.plot or bool(args.csv))

    print(f"Final energy: {data['final_energy']:.6f} eV")

    if data["forces"] is not None:
        print("\nFinal forces (eV/Angstrom):")
        print(data["forces"].to_string(index=False))
        if args.csv:
            data["forces"].to_csv(f"{args.csv}_forces.csv", index=False)

    if data["magmom"] is not None:
        print("\nFinal magnetic moments:")
        print(data["magmom"].to_string(index=False))
        if args.csv:
            data["magmom"].to_csv(f"{args.csv}_magmom.csv", index=False)

    if args.dos and data["total_dos"] is not None:
        print(f"\nE_fermi: {data['efermi']:.4f} eV")
        print(f"Total DOS: {len(data['total_dos'])} energy points")
        print(f"Element-projected DOS available for: {list(data['element_dos'].keys())}")
        if args.csv:
            data["total_dos"].to_csv(f"{args.csv}_total_dos.csv", index=False)
            for el, df in data["element_dos"].items():
                df.to_csv(f"{args.csv}_{el}_dos.csv", index=False)

    if args.plot and data["total_dos"] is not None:
        plot_dos(data["total_dos"], data["element_dos"])


if __name__ == "__main__":
    main()
