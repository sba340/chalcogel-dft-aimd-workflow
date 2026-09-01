"""
pcf_pdf_comparison.py

Overlays a simulated pair correlation function (PCF, from AIMD melt-quench,
.dat file) against an experimental pair distribution function (PDF, .dat
file) for visual comparison, with bond-type annotations.

Both input files are expected to be plain-text, two-column [r, g(r)] data
(whitespace-delimited), as produced by typical PCF/PDF analysis tools.

Usage:
    python pcf_pdf_comparison.py sim.dat exp.dat --sim-label "..." --exp-label "..."

Example:
    python pcf_pdf_comparison.py pair_correlation.dat "G14 PDF latest.dat" \\
        --sim-label "Fe0.5Mo3S13 PCF Simulation" \\
        --exp-label "Fe0.5Mo3S13 PDF Experiment" \\
        --annotate "1.98:S-S" "2.30:Mo/Fe-S" "2.7:Mo-Mo"
"""

from __future__ import annotations

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def parse_annotation(spec: str) -> tuple[float, str]:
    """Parse a 'position:label' string, e.g. '1.98:S-S'."""
    pos_str, label = spec.split(":", 1)
    return float(pos_str), label


def plot_pcf_pdf(
    sim_file: str,
    exp_file: str,
    sim_label: str,
    exp_label: str,
    out_path: str = "compared_pcf_pdf.png",
    sim_color: str = "#701705",
    offset: float = 0.25,
    shift_right: float = -0.07,
    x_limit: tuple[float, float] = (1.80, 6),
    y_limit: tuple[float, float] = (-1.2, 4),
    annotations: list[tuple[float, str]] | None = None,
):
    if annotations is None:
        annotations = [(1.98, "S-S"), (2.30, "M-S"), (2.7, "Mo-Mo")]

    sim_label_pos = (x_limit[1] - 0.5, y_limit[1] * 0.45)
    exp_label_pos = (x_limit[1] - 0.5, y_limit[1] * 0.15)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set_style("white")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 16})

    ax.grid(axis="y", linestyle="-", alpha=0.2)

    # Simulation data (optionally shifted to align with experiment)
    sim_data = np.loadtxt(sim_file)
    sim_r = sim_data[:, 0] + shift_right
    sim_g = sim_data[:, 1] + offset
    ax.plot(sim_r, sim_g, color=sim_color, linewidth=4)

    # Experimental data
    exp_data = np.loadtxt(exp_file)
    exp_r = exp_data[:, 0]
    exp_g = exp_data[:, 1]
    ax.plot(exp_r, exp_g, color="black", linewidth=4)

    ax.annotate(sim_label, xy=sim_label_pos, color=sim_color, fontsize=13, ha="right", va="center")
    ax.annotate(exp_label, xy=exp_label_pos, color="black", fontsize=13, ha="right", va="center")

    ax.set_xlim(x_limit)
    ax.set_ylim(y_limit)
    ax.set_xlabel("R(\u00c5)", fontsize=14)
    ax.set_ylabel("g(r)", fontsize=14)
    ax.set_xticks(np.arange(np.ceil(x_limit[0]), x_limit[1] + 0.1, 0.5))
    ax.set_yticks([])
    ax.tick_params(axis="x", which="major", bottom=True, top=False, length=6)

    for bond_pos, label in annotations:
        ax.axvline(x=bond_pos, color="black", linestyle="--", linewidth=1.5)
        ax.text(
            bond_pos + 0.1, y_limit[1] * 0.65, label,
            fontsize=14, fontweight="normal", rotation="vertical",
            va="bottom", ha="center", style="normal",
        )

    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1.5)

    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    print(f"Plot saved to {out_path}")
    return fig, ax


def main():
    ap = argparse.ArgumentParser(description="Overlay simulated PCF against experimental PDF.")
    ap.add_argument("sim_file", help="Path to simulated PCF .dat file (two columns: r, g(r))")
    ap.add_argument("exp_file", help="Path to experimental PDF .dat file (two columns: r, g(r))")
    ap.add_argument("--sim-label", default="Simulation", help="Legend label for simulated curve")
    ap.add_argument("--exp-label", default="Experiment", help="Legend label for experimental curve")
    ap.add_argument("--out", default="compared_pcf_pdf.png", help="Output plot filename")
    ap.add_argument("--offset", type=float, default=0.25, help="Vertical offset applied to sim curve")
    ap.add_argument("--shift-right", type=float, default=-0.07, help="Horizontal shift applied to sim curve (r-axis)")
    ap.add_argument("--xlim", type=float, nargs=2, default=[1.80, 6], metavar=("XMIN", "XMAX"))
    ap.add_argument("--ylim", type=float, nargs=2, default=[-1.2, 4], metavar=("YMIN", "YMAX"))
    ap.add_argument(
        "--annotate", nargs="*", default=["1.98:S-S", "2.30:M-S", "2.7:Mo-Mo"],
        help="Bond annotations as 'position:label' pairs, e.g. 1.98:S-S",
    )
    args = ap.parse_args()

    annotations = [parse_annotation(a) for a in args.annotate]

    plot_pcf_pdf(
        sim_file=args.sim_file,
        exp_file=args.exp_file,
        sim_label=args.sim_label,
        exp_label=args.exp_label,
        out_path=args.out,
        offset=args.offset,
        shift_right=args.shift_right,
        x_limit=tuple(args.xlim),
        y_limit=tuple(args.ylim),
        annotations=annotations,
    )


if __name__ == "__main__":
    main()
