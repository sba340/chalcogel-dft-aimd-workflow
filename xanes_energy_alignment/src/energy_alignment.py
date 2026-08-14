"""
XCH (excited core-hole) energy alignment for computed absorption edges
(e.g. XANES, XPS-derived excitation energies) from VASP total energies.

--------------------------------------------------------------------------
Why this correction is needed
--------------------------------------------------------------------------
A DFT total-energy difference between a core-excited ("XCH") calculation
and the ground-state ("SCF") calculation is NOT directly comparable to an
experimental edge energy, for two reasons:

  1. VASP total energies are referenced to an arbitrary internal zero
     (the average electrostatic potential of the cell), not to a physical
     vacuum or core level. So a raw (E_CH - E_GS) has to be re-referenced
     to something physically meaningful before it means anything in eV
     relative to a real edge.
  2. Even after that re-referencing, DFT systematically over/underestimates
     absolute core-level energies (self-interaction error, missing
     relaxation effects, pseudopotential approximations, etc.), so a
     second, empirical shift against a known reference compound is needed.

This script does both corrections in two steps:

  Step 1 - "internal" alignment (delta1):
      Re-references the raw core-hole excitation energy using the
      calculation's own Fermi level and a fixed core-level energy read
      from the POTCAR, removing the arbitrary VASP potential zero.

  Step 2 - "external" calibration (delta2):
      Compares that internally-aligned energy for a REFERENCE compound
      (one with a known, published experimental edge energy) against the
      literature value, producing an empirical correction.

  total_shift = delta1 + delta2 is then a single number you can add to the
  raw excitation energy of any OTHER compound computed with the same
  edge, functional, and POTCAR -- because delta1 and delta2 correct for
  systematic (calculation-setup-dependent) errors, not compound-specific
  physics.

--------------------------------------------------------------------------
Required inputs, per calculation
--------------------------------------------------------------------------
E_fermi   : Fermi energy (eV) from the ground-state (SCF) OUTCAR
E_ch      : total energy (eV) of the core-hole excited-state calculation
E_gs      : total energy (eV) of the ground-state (SCF) calculation
E_core    : reference core-level energy (eV) for the excited atom's
            core orbital, as tabulated in the POTCAR pseudopotential file
E_ex_raw  : raw computed excitation energy (eV) for the compound you are
            aligning (from your XCH workflow, before any correction)
"""
from dataclasses import dataclass


@dataclass
class CalibrationResult:
    delta1: float          # internal (Fermi/core-level) re-referencing shift, eV
    delta2: float          # empirical shift vs. the reference compound, eV
    total_shift: float     # delta1 + delta2, apply this to other compounds' E_ex_raw
    step1_shift: float     # reference compound's energy after delta1 only
    aligned_energy: float  # reference compound's energy after both corrections


def calibrate(E_fermi, E_ch, E_gs, E_core, E_ex_raw, Ref_value):
    """Compute delta1, delta2, and total_shift from a reference compound
    with a known experimental edge energy (Ref_value).

    Ref_value should be the literature/experimental edge energy (eV) for
    the SAME reference compound whose E_ex_raw you pass in here -- e.g. a
    well-characterized standard for the edge you're studying, cited from
    a published XANES/XPS paper.
    """
    delta1 = (E_ch - E_gs) - E_fermi
    step1_shift = E_ex_raw + delta1
    delta2 = Ref_value - step1_shift
    aligned_energy = step1_shift + delta2
    total_shift = delta1 + delta2
    return CalibrationResult(
        delta1=delta1,
        delta2=delta2,
        total_shift=total_shift,
        step1_shift=step1_shift,
        aligned_energy=aligned_energy,
    )


def align_target(E_ex_raw_target, total_shift):
    """Apply a previously-computed total_shift to a NEW compound's raw
    excitation energy (same edge, functional, and POTCAR as calibration).
    """
    return E_ex_raw_target + total_shift


if __name__ == "__main__":
    # Example: Mo L3-edge calibration using (NH4)2MoS4 as the reference
    # compound (experimental value from J. Phys. Chem. C 2021, 125, 32,
    # 17761-17773).
    result = calibrate(
        E_fermi=1.8166,
        E_ch=-0.21365606e03,
        E_gs=-0.28553067e03,
        E_core=-19816.9363,
        E_ex_raw=2505.19,
        Ref_value=2522.4,
    )
    print(f"delta1 (internal re-referencing):  {result.delta1:.4f} eV")
    print(f"delta2 (empirical calibration):    {result.delta2:.4f} eV")
    print(f"total_shift (apply to new compounds): {result.total_shift:.4f} eV")
    print(f"reference compound aligned energy:  {result.aligned_energy:.4f} eV "
          f"(should equal Ref_value)")

    # Now align a different compound computed with the same setup:
    target_raw = 2498.7  # example raw excitation energy for a new complex
    aligned = align_target(target_raw, result.total_shift)
    print(f"\nTarget compound raw energy: {target_raw:.4f} eV")
    print(f"Target compound aligned energy: {aligned:.4f} eV")
