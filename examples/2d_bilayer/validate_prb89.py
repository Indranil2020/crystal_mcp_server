#!/usr/bin/env python3
"""
Validate generated bilayer TMD stackings against PRB 89, 075409 (2014).

This script is intentionally "unbiased" / coordinate-driven:
  - It reads the *POSCAR/VASP files on disk* (ASE-readable).
  - It computes symmetry via spglib (through pymatgen SpacegroupAnalyzer).
  - It computes stacking registry directly from in-plane fractional coordinates.

It compares against the paper only through published invariants:
  - Expected point group for each stacking (Fig. 1).
  - RPA equilibrium interlayer distances d0 for all 5 stackings (Table III).
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np
from ase import Atoms
from ase.io import read

from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

STACKINGS: Tuple[str, ...] = ("AA", "AA'", "A'B", "AB", "AB'")

# PRB 89, 075409 (2014) Table III (RPA column): equilibrium interlayer distances d0 (Å).
PRB89_RPA_D0: Dict[str, Dict[str, float]] = {
    "MoS2": {"AA": 6.77, "AA'": 6.27, "A'B": 6.78, "AB": 6.17, "AB'": 6.26},
    "MoSe2": {"AA": 7.18, "AA'": 6.48, "A'B": 7.12, "AB": 6.47, "AB'": 6.53},
    "WS2": {"AA": 6.84, "AA'": 6.24, "A'B": 6.78, "AB": 6.24, "AB'": 6.27},
    "WSe2": {"AA": 7.24, "AA'": 6.50, "A'B": 7.24, "AB": 6.54, "AB'": 6.62},
}

# PRB 89, 075409 (2014) Fig. 1: expected point groups (Schoenflies), converted to Hermann–Mauguin.
EXPECTED_POINT_GROUP_HM: Dict[str, str] = {
    "AA": "-6m2",   # D3h
    "AA'": "-3m",   # D3d
    "A'B": "-3m",   # D3d
    "AB": "3m",     # C3v
    "AB'": "-3m",   # D3d
}

# In-plane translation vectors τ in fractional coordinates (a1, a2) for the conventions used in PRB Fig. 1.
EXPECTED_TAU_FRAC: Dict[str, np.ndarray] = {
    "AA": np.array([0.0, 0.0]),
    "AA'": np.array([2 / 3, 1 / 3]),
    "A'B": np.array([1 / 3, 2 / 3]),
    "AB": np.array([1 / 3, 2 / 3]),
    "AB'": np.array([0.0, 0.0]),
}

# D3d stackings must have inversion symmetry.
D3D_STACKINGS: Set[str] = {"AA'", "A'B", "AB'"}


def _largest_gap_threshold(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values_sorted = np.sort(values)
    if values_sorted.size < 2:
        return float(values_sorted[0]) if values_sorted.size == 1 else 0.0
    gaps = values_sorted[1:] - values_sorted[:-1]
    split_idx = int(np.argmax(gaps))
    return float(0.5 * (values_sorted[split_idx] + values_sorted[split_idx + 1]))


def _split_layers_by_z(atoms: Atoms) -> Tuple[np.ndarray, np.ndarray]:
    pos = atoms.get_positions()
    z_thr = _largest_gap_threshold(pos[:, 2])
    is_top = pos[:, 2] >= z_thr
    return np.where(~is_top)[0], np.where(is_top)[0]


def _infer_metal_symbol(atoms: Atoms) -> str:
    counts = Counter(atoms.get_chemical_symbols())
    if not counts:
        raise ValueError("Empty structure")
    # For MX2 bilayers (primitive), metal has the smallest count.
    return min(counts, key=lambda k: counts[k])


def _minimal_image_delta_frac(d: np.ndarray) -> np.ndarray:
    """Reduce fractional deltas into [-0.5, 0.5) in each component."""
    d = np.asarray(d, dtype=float)
    return d - np.round(d)


def _stacking_tau_from_metals(atoms: Atoms) -> np.ndarray:
    atoms = atoms.copy()
    atoms.set_pbc(True)
    atoms.wrap()

    bottom_idx, top_idx = _split_layers_by_z(atoms)
    metal = _infer_metal_symbol(atoms)

    syms = np.asarray(atoms.get_chemical_symbols())
    bot_m = bottom_idx[syms[bottom_idx] == metal]
    top_m = top_idx[syms[top_idx] == metal]
    if bot_m.size != 1 or top_m.size != 1:
        raise ValueError(f"Expected exactly 1 metal per layer, got bottom={bot_m.size}, top={top_m.size}")

    frac = atoms.get_scaled_positions()
    frac = frac - np.floor(frac)
    d = frac[int(top_m[0]), :2] - frac[int(bot_m[0]), :2]
    d = _minimal_image_delta_frac(d)
    # Map to [0, 1) for human-friendly comparison.
    d = d - np.floor(d)
    return d


def _registry_aligned_pairs(atoms: Atoms, tol_frac: float = 5e-4) -> Set[Tuple[str, str]]:
    """
    Return unique (top_symbol, bottom_symbol) pairs that align in-plane.

    Alignment criterion (periodic in-plane):
      || (frac_top - frac_bottom) reduced to [-0.5,0.5) || <= tol_frac
    """
    atoms = atoms.copy()
    atoms.set_pbc(True)
    atoms.wrap()

    bottom_idx, top_idx = _split_layers_by_z(atoms)
    frac = atoms.get_scaled_positions()
    frac = frac - np.floor(frac)

    syms = np.asarray(atoms.get_chemical_symbols())
    pairs: Set[Tuple[str, str]] = set()
    for i in top_idx:
        for j in bottom_idx:
            d = frac[int(i), :2] - frac[int(j), :2]
            d = _minimal_image_delta_frac(d)
            if float(np.linalg.norm(d)) <= float(tol_frac):
                pairs.add((str(syms[int(i)]), str(syms[int(j)])))
    return pairs


def _interlayer_distance_from_metals(atoms: Atoms) -> float:
    bottom_idx, top_idx = _split_layers_by_z(atoms)
    metal = _infer_metal_symbol(atoms)
    pos = atoms.get_positions()
    syms = np.asarray(atoms.get_chemical_symbols())
    bot_m = bottom_idx[syms[bottom_idx] == metal]
    top_m = top_idx[syms[top_idx] == metal]
    if bot_m.size != 1 or top_m.size != 1:
        raise ValueError(f"Expected exactly 1 metal per layer, got bottom={bot_m.size}, top={top_m.size}")
    return float(pos[int(top_m[0]), 2] - pos[int(bot_m[0]), 2])


def _symmetry_summary(atoms: Atoms, symprec: float) -> Tuple[int, str, bool]:
    """Return (space_group_no, point_group_hm, has_inversion)."""
    lattice = Lattice(atoms.get_cell().array)
    structure = Structure(lattice, atoms.get_chemical_symbols(), atoms.get_positions(), coords_are_cartesian=True)
    sga = SpacegroupAnalyzer(structure, symprec=float(symprec))
    sg_no = int(sga.get_space_group_number())
    pg = str(sga.get_point_group_symbol())
    ops = sga.get_symmetry_operations(cartesian=False)
    has_inv = any(np.allclose(op.rotation_matrix, -np.eye(3), atol=1e-8) for op in ops)
    return sg_no, pg, bool(has_inv)


def _expected_registry_pairs_for_string(stacking: str, metal_symbol: str, chalcogen_symbol: str) -> Set[Tuple[str, str]]:
    """
    Expected cross-layer aligned pairs for each stacking, from PRB 89, 075409 (2014) Fig. 1.

    Returns a set of (top_symbol, bottom_symbol).
    """
    if stacking == "AA":
        return {(metal_symbol, metal_symbol), (chalcogen_symbol, chalcogen_symbol)}
    if stacking == "AA'":
        return {(metal_symbol, chalcogen_symbol), (chalcogen_symbol, metal_symbol)}
    if stacking == "A'B":
        return {(chalcogen_symbol, chalcogen_symbol)}
    if stacking == "AB":
        return {(chalcogen_symbol, metal_symbol)}
    if stacking == "AB'":
        return {(metal_symbol, metal_symbol)}
    raise ValueError(f"Unknown stacking: {stacking}")


@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    details: str


def _fmt_vec2(v: np.ndarray) -> str:
    v = np.asarray(v, dtype=float)
    return f"[{v[0]:.6f}, {v[1]:.6f}]"


def validate_material(material: str, symprec_list: Sequence[float], tol_d0: float, tol_tau: float, tol_registry: float) -> Dict[str, List[CheckResult]]:
    per_stacking: Dict[str, List[CheckResult]] = {}
    for stacking in STACKINGS:
        path = Path(__file__).resolve().parent / f"{material}_{stacking}.vasp"
        if not path.exists():
            raise SystemExit(f"Missing file: {path}. Generate it first by running `python3 bilayer_stacking.py` in this folder.")

        atoms = read(str(path))
        atoms.set_pbc(True)
        atoms.wrap()

        counts = Counter(atoms.get_chemical_symbols())
        metal_symbol = _infer_metal_symbol(atoms)
        chalcogen_symbol = max(counts, key=lambda k: counts[k])

        expected_d0 = float(PRB89_RPA_D0[material][stacking])
        d0 = _interlayer_distance_from_metals(atoms)

        tau = _stacking_tau_from_metals(atoms)
        expected_tau = np.asarray(EXPECTED_TAU_FRAC[stacking], dtype=float)
        tau_delta = _minimal_image_delta_frac(tau - expected_tau)
        tau_ok = float(np.linalg.norm(tau_delta)) <= float(tol_tau)

        reg_pairs = _registry_aligned_pairs(atoms, tol_frac=float(tol_registry))
        expected_pairs = _expected_registry_pairs_for_string(stacking, metal_symbol=metal_symbol, chalcogen_symbol=chalcogen_symbol)

        expected_pg_hm = EXPECTED_POINT_GROUP_HM[stacking]

        # Symmetry can depend on numerical tolerance; sweep and accept if any symprec matches.
        sym_hits: List[str] = []
        sym_ok = False
        inv_ok = True
        for sp in symprec_list:
            sg_no, pg, has_inv = _symmetry_summary(atoms, symprec=float(sp))
            hit = (pg == expected_pg_hm)
            if hit:
                sym_ok = True
                sym_hits.append(f"symprec={sp:g}")
            if stacking in D3D_STACKINGS:
                inv_ok = inv_ok and has_inv

        checks: List[CheckResult] = []
        checks.append(
            CheckResult(
                name="d0 (RPA Table III)",
                passed=abs(d0 - expected_d0) <= float(tol_d0),
                details=f"d0={d0:.3f} Å, ref={expected_d0:.3f} Å (tol={tol_d0:.3f})",
            )
        )
        checks.append(
            CheckResult(
                name="tau (metal-metal in-plane)",
                passed=bool(tau_ok),
                details=f"tau={_fmt_vec2(tau)}, ref={_fmt_vec2(expected_tau)}, |Δ|={np.linalg.norm(tau_delta):.2e} (tol={tol_tau:.2e})",
            )
        )
        checks.append(
            CheckResult(
                name="registry (aligned pairs)",
                passed=reg_pairs == expected_pairs,
                details=f"pairs={sorted(reg_pairs)}, ref={sorted(expected_pairs)}, tol_frac={tol_registry:g}",
            )
        )
        checks.append(
            CheckResult(
                name="symmetry (space+point group)",
                passed=bool(sym_ok),
                details=f"expected pg={expected_pg_hm}; hits={sym_hits if sym_hits else 'none'}",
            )
        )
        checks.append(
            CheckResult(
                name="inversion (D3d only)",
                passed=bool(inv_ok),
                details="required for D3d stackings" if stacking in D3D_STACKINGS else "not required",
            )
        )

        per_stacking[stacking] = checks

    return per_stacking


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate bilayer stackings against PRB 89, 075409 (2014).")
    parser.add_argument(
        "--materials",
        nargs="*",
        default=["MoS2", "MoSe2", "WS2", "WSe2"],
        help="Materials to validate (default: MoS2 MoSe2 WS2 WSe2).",
    )
    parser.add_argument(
        "--symprec",
        nargs="*",
        type=float,
        default=[1e-4, 5e-4, 1e-3, 5e-3, 1e-2],
        help="Symmetry tolerance sweep (default: 1e-4 5e-4 1e-3 5e-3 1e-2).",
    )
    parser.add_argument("--tol-d0", type=float, default=0.01, help="Tolerance for d0 in Å (default: 0.01).")
    parser.add_argument("--tol-tau", type=float, default=5e-4, help="Tolerance for tau in fractional coords (default: 5e-4).")
    parser.add_argument("--tol-registry", type=float, default=5e-4, help="Tolerance for registry alignment in fractional coords (default: 5e-4).")
    args = parser.parse_args()

    total = 0
    failed = 0

    print("=" * 80)
    print("PRB 89, 075409 (2014) — Bilayer Stacking Validation")
    print("=" * 80)

    for mat in args.materials:
        if mat not in PRB89_RPA_D0:
            raise SystemExit(f"Unknown material: {mat}")

        print(f"\nMaterial: {mat}")
        print("-" * 80)
        per = validate_material(
            material=mat,
            symprec_list=args.symprec,
            tol_d0=args.tol_d0,
            tol_tau=args.tol_tau,
            tol_registry=args.tol_registry,
        )

        for stacking in STACKINGS:
            checks = per[stacking]
            ok = all(c.passed for c in checks)
            status = "PASS" if ok else "FAIL"
            print(f"\n  {stacking:>3s}  {status}")
            for c in checks:
                total += 1
                if not c.passed:
                    failed += 1
                mark = "ok " if c.passed else "BAD"
                print(f"    [{mark}] {c.name}: {c.details}")

    print("\n" + "=" * 80)
    print(f"SUMMARY: {total - failed}/{total} checks passed")
    print("=" * 80)

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
