"""
Run an energy-alignment calibration (and optionally align target compounds)
from a config file, so you don't have to hand-edit the script per edge.

Usage:
    python src/run_alignment.py --config examples/mo_l3_edge.yaml
"""
import argparse
import yaml

from energy_alignment import calibrate, align_target


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="YAML file with calibration (and optional targets) block")
    args = p.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    cal = cfg["calibration"]
    result = calibrate(
        E_fermi=cal["E_fermi"],
        E_ch=cal["E_ch"],
        E_gs=cal["E_gs"],
        E_core=cal["E_core"],
        E_ex_raw=cal["E_ex_raw"],
        Ref_value=cal["Ref_value"],
    )

    print(f"Reference compound: {cal.get('name', 'unnamed')}")
    print(f"Edge: {cfg.get('edge', 'unspecified')}")
    print(f"  delta1:      {result.delta1:.4f} eV")
    print(f"  delta2:      {result.delta2:.4f} eV")
    print(f"  total_shift: {result.total_shift:.4f} eV")
    print(f"  aligned energy: {result.aligned_energy:.4f} eV "
          f"(should match Ref_value = {cal['Ref_value']})\n")

    for target in cfg.get("targets", []):
        aligned = align_target(target["E_ex_raw"], result.total_shift)
        print(f"{target.get('name', 'unnamed target')}: "
              f"raw = {target['E_ex_raw']:.4f} eV -> aligned = {aligned:.4f} eV")


if __name__ == "__main__":
    main()
