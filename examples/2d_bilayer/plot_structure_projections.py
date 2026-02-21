#!/usr/bin/env python3
"""
Orthographic projection plots (xy/xz/yz) for structures (e.g. POSCAR/vasp files).

This is intentionally "coordinate-driven": it reads whatever coordinates are in the file
and plots them directly (no stacking-specific logic / no hardcoded registries).

Typical usage (from this folder):
  python3 plot_structure_projections.py --glob "*.vasp" --outdir plots --repeat 2 2 1
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Sequence, Tuple

import numpy as np
import warnings

warnings.filterwarnings(
    "ignore",
    message=r"Unable to import Axes3D\..*",
    category=UserWarning,
)

import matplotlib.pyplot as plt

from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.data.colors import jmol_colors
from ase.io import read


def _parse_repeat(values: Sequence[str]) -> Tuple[int, int, int]:
    if len(values) != 3:
        raise ValueError("repeat must have 3 integers: Nx Ny Nz")
    repeat = tuple(int(v) for v in values)
    if any(v < 1 for v in repeat):
        raise ValueError("repeat values must be >= 1")
    return repeat[0], repeat[1], repeat[2]


def _norm(vec: np.ndarray) -> float:
    return float(np.linalg.norm(vec))


def _unit(vec: np.ndarray) -> np.ndarray:
    length = _norm(vec)
    if length == 0.0:
        raise ValueError("Cannot normalize a zero vector")
    return np.asarray(vec, dtype=float) / length


def _lattice_aligned_basis(cell: np.ndarray) -> np.ndarray:
    """
    Return an orthonormal basis as a (3, 3) matrix with basis vectors as rows.
    x̂ is along a1, ŷ is a2 made orthogonal to a1, ẑ is x̂×ŷ.
    """
    a1 = np.asarray(cell[0], dtype=float)
    a2 = np.asarray(cell[1], dtype=float)

    x_hat = _unit(a1)
    a2_perp = a2 - np.dot(a2, x_hat) * x_hat
    y_hat = _unit(a2_perp)
    z_hat = _unit(np.cross(x_hat, y_hat))
    return np.vstack([x_hat, y_hat, z_hat])


def _largest_gap_threshold(values: np.ndarray) -> float:
    """
    Split values into two clusters by the largest gap in the sorted list.
    Intended for bilayers: the interlayer gap is typically the largest.
    """
    values = np.asarray(values, dtype=float)
    values_sorted = np.sort(values)
    if values_sorted.size < 2:
        return float(values_sorted[0]) if values_sorted.size == 1 else 0.0

    gaps = values_sorted[1:] - values_sorted[:-1]
    split_idx = int(np.argmax(gaps))
    return float(0.5 * (values_sorted[split_idx] + values_sorted[split_idx + 1]))


def _rgba_for_atoms(
    symbols: Sequence[str],
    is_top_layer: np.ndarray,
    bottom_darkening: float = 0.65,
    top_alpha: float = 0.95,
    bottom_alpha: float = 0.80,
) -> np.ndarray:
    """
    Colors are element-based (Jmol). Layer is indicated by darkening + alpha.
    """
    rgba = np.empty((len(symbols), 4), dtype=float)
    for i, sym in enumerate(symbols):
        z = atomic_numbers[sym]
        base = np.asarray(jmol_colors[z], dtype=float)
        if bool(is_top_layer[i]):
            rgb = base
            alpha = top_alpha
        else:
            rgb = np.clip(base * bottom_darkening, 0.0, 1.0)
            alpha = bottom_alpha
        rgba[i, :3] = rgb
        rgba[i, 3] = alpha
    return rgba


def _sizes_for_atoms(symbols: Sequence[str], scale: float = 40.0) -> np.ndarray:
    """
    Scatter sizes derived from covalent radii (no material/stacking hardcoding).
    Matplotlib expects marker size in points^2.
    """
    sizes = np.empty(len(symbols), dtype=float)
    for i, sym in enumerate(symbols):
        z = atomic_numbers[sym]
        radius = float(covalent_radii[z])
        sizes[i] = (radius * scale) ** 2
    return sizes


def _cell_outline_xy(cell: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Return a (5, 2) polyline for the (already repeated) cell outline in the XY plane."""
    a1 = np.asarray(cell[0], dtype=float)
    a2 = np.asarray(cell[1], dtype=float)
    a1b = a1 @ basis.T
    a2b = a2 @ basis.T
    o = np.array([0.0, 0.0])
    p1 = a1b[:2]
    p2 = (a1b + a2b)[:2]
    p3 = a2b[:2]
    return np.vstack([o, p1, p2, p3, o])


def _draw_cell_grid_xy(
    ax: plt.Axes,
    cell: np.ndarray,
    basis: np.ndarray,
    nx: int,
    ny: int,
    color: str = "0.65",
    lw: float = 1.2,
) -> None:
    """
    Draw the Nx×Ny tiling grid for a supercell defined by `cell`.

    Assumes the in-plane supercell vectors are:
      A1 = nx * a1_prim, A2 = ny * a2_prim
    so the primitive vectors are recovered by division.
    """
    if nx < 1 or ny < 1:
        return

    a1 = np.asarray(cell[0], dtype=float) / float(nx)
    a2 = np.asarray(cell[1], dtype=float) / float(ny)

    for i in range(nx + 1):
        p0 = (i * a1) @ basis.T
        p1 = (i * a1 + ny * a2) @ basis.T
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, zorder=0)

    for j in range(ny + 1):
        p0 = (j * a2) @ basis.T
        p1 = (j * a2 + nx * a1) @ basis.T
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, lw=lw, zorder=0)


def _mode_rounded(values: np.ndarray, decimals: int) -> float:
    values = np.asarray(values, dtype=float)
    rounded = np.round(values, decimals=int(decimals))
    uniq, counts = np.unique(rounded, return_counts=True)
    return float(uniq[int(np.argmax(counts))])


def _detect_cross_layer_alignments(
    atoms_prim: Atoms,
    basis: np.ndarray,
    wrap: bool,
    tol_frac: float,
) -> List[Tuple[int, int, np.ndarray]]:
    """
    Return a list of (top_index, bottom_index, shift_ij) alignments.

    `shift_ij` is the integer (a1, a2) translation that maps the top in-plane
    fractional coordinates onto the bottom (i.e. frac_top - shift ≈ frac_bottom).
    """
    atoms = atoms_prim.copy()
    atoms.set_pbc(True)
    if wrap:
        atoms.wrap()

    pos = atoms.get_positions()
    coords = pos @ basis.T
    z_threshold = _largest_gap_threshold(coords[:, 2])
    is_top = coords[:, 2] >= z_threshold

    frac = atoms.get_scaled_positions()
    frac = frac - np.floor(frac)

    top_idx = np.where(is_top)[0]
    bot_idx = np.where(~is_top)[0]

    alignments: List[Tuple[int, int, np.ndarray]] = []
    for i in top_idx:
        for j in bot_idx:
            d = frac[i, :2] - frac[j, :2]
            shift = np.round(d).astype(int)
            d_red = d - shift
            if float(np.linalg.norm(d_red)) <= float(tol_frac):
                alignments.append((int(i), int(j), shift))
    return alignments


def _detect_bonds_intralayer(
    positions: np.ndarray,
    numbers: np.ndarray,
    is_top_layer: np.ndarray,
    bond_margin: float,
    max_distance: float = 6.0,
) -> List[Tuple[int, int]]:
    """
    Detect intralayer bonds using a coordinate-derived cutoff:
      cutoff(pair_type) = d_min(pair_type) * (1 + bond_margin)
    where d_min is the minimum inter-atomic distance for that pair type inside the layer.

    This avoids radii-based cutoffs that can include too many neighbors for metals.
    """
    if bond_margin < 0.0:
        raise ValueError("bond_margin must be >= 0")

    positions = np.asarray(positions, dtype=float)
    numbers = np.asarray(numbers, dtype=int)
    is_top_layer = np.asarray(is_top_layer, dtype=bool)

    bonds: List[Tuple[int, int]] = []
    for layer in (False, True):
        idx = np.where(is_top_layer == layer)[0]
        if idx.size < 2:
            continue

        d_min_by_pair: dict[Tuple[int, int], float] = {}
        pair_dists: List[Tuple[Tuple[int, int], int, int, float]] = []
        for a in range(len(idx)):
            i = int(idx[a])
            for b in range(a + 1, len(idx)):
                j = int(idx[b])
                zi = int(numbers[i])
                zj = int(numbers[j])
                if zi == zj:
                    continue
                pair = (zi, zj) if zi < zj else (zj, zi)
                d = float(np.linalg.norm(positions[i] - positions[j]))
                if d > max_distance:
                    continue
                pair_dists.append((pair, i, j, d))
                current = d_min_by_pair.get(pair)
                if current is None or d < current:
                    d_min_by_pair[pair] = d

        if not pair_dists:
            continue

        thresholds = {pair: dmin * (1.0 + float(bond_margin)) for pair, dmin in d_min_by_pair.items()}
        for pair, i, j, d in pair_dists:
            if d <= thresholds[pair]:
                bonds.append((i, j))

    return bonds


def _closest_image_shift(
    delta: np.ndarray,
    a1_xy: np.ndarray,
    a2_xy: np.ndarray,
    search: int,
) -> np.ndarray:
    """
    Find integer (m, n) in [-search, search] that minimizes |delta + m*a1 + n*a2|.
    Returns the best shift vector m*a1 + n*a2.
    """
    best_shift = np.zeros(2, dtype=float)
    best_d2 = float("inf")
    for m in range(-search, search + 1):
        for n in range(-search, search + 1):
            shift = m * a1_xy + n * a2_xy
            d = delta + shift
            d2 = float(d[0] * d[0] + d[1] * d[1])
            if d2 < best_d2:
                best_d2 = d2
                best_shift = shift
    return best_shift


def _prepare_atoms(atoms: Atoms, repeat: Tuple[int, int, int], wrap: bool) -> Atoms:
    atoms = atoms.copy()
    atoms.set_pbc(True)
    if wrap:
        atoms.wrap()
    if repeat != (1, 1, 1):
        atoms = atoms.repeat(repeat)
    # Avoid drawing "wrap-around" bonds across the supercell boundaries.
    atoms.set_pbc(False)
    return atoms


def _tile_indices_from_scaled(
    scaled: np.ndarray,
    nx: int,
    ny: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Given scaled positions in a supercell, return the primitive-tile indices (ix, iy).

    The supercell is assumed to be built by repeating a primitive cell (nx, ny) in-plane.
    """
    scaled = np.asarray(scaled, dtype=float)
    if scaled.ndim != 2 or scaled.shape[1] < 2:
        raise ValueError("scaled positions must have shape (N, >=2)")
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be >= 1")

    scaled = scaled - np.floor(scaled)
    # Guard against floating point values like 0.49999999999999994 that should
    # represent an exact tile boundary (e.g. 0.5 for nx=2).
    eps = 1e-9
    x = np.clip(np.floor(scaled[:, 0] * float(nx) + eps).astype(int), 0, nx - 1)
    y = np.clip(np.floor(scaled[:, 1] * float(ny) + eps).astype(int), 0, ny - 1)
    return x, y


def _tile_outline_xy(
    a1_prim: np.ndarray,
    a2_prim: np.ndarray,
    ix: int,
    iy: int,
    basis: np.ndarray,
) -> np.ndarray:
    """Return a (5, 2) polyline of the primitive tile (ix, iy) projected into the XY plotting basis."""
    origin = float(ix) * np.asarray(a1_prim, dtype=float) + float(iy) * np.asarray(a2_prim, dtype=float)
    p0 = origin @ basis.T
    p1 = (origin + a1_prim) @ basis.T
    p2 = (origin + a1_prim + a2_prim) @ basis.T
    p3 = (origin + a2_prim) @ basis.T
    return np.vstack([p0[:2], p1[:2], p2[:2], p3[:2], p0[:2]])


def plot_projections(
    atoms: Atoms,
    outpath: Path,
    repeat: Tuple[int, int, int] = (2, 2, 1),
    wrap: bool = True,
    basis_mode: str = "cart",
    show_cell: bool = True,
    cell_grid: bool = True,
    grid: bool = True,
    size_scale: float = 12.0,
    bond_margin: float = 0.10,
    bond_color: str = "0.25",
    bond_lw: float = 1.0,
    bonds: Literal["all", "side", "none"] = "side",
    side_mode: Literal["cell", "slice", "all"] = "cell",
    side_cell: Tuple[int, int] = (0, 0),
    side_slice_width: float = 0.25,
    side_slice_round_decimals: int = 2,
    slice_guides: bool = True,
    align_lines: bool = True,
    align_tol_frac: float = 5e-4,
    align_color: str = "0.35",
    align_lw: float = 0.8,
    align_search: int = 1,
    dpi: int = 200,
) -> None:
    atoms_prim = atoms.copy()
    atoms = _prepare_atoms(atoms, repeat=repeat, wrap=wrap)

    cell = atoms.get_cell().array
    basis = np.eye(3) if basis_mode == "cart" else _lattice_aligned_basis(cell)

    positions = atoms.get_positions()
    coords = positions @ basis.T  # components along x̂, ŷ, ẑ

    z_threshold = _largest_gap_threshold(coords[:, 2])
    is_top = coords[:, 2] >= z_threshold

    symbols = atoms.get_chemical_symbols()
    rgba = _rgba_for_atoms(symbols, is_top_layer=is_top)
    sizes = _sizes_for_atoms(symbols, scale=size_scale)

    # Bonds (coordinate-driven): intralayer, distance-derived threshold.
    numbers = np.asarray(atoms.get_atomic_numbers(), dtype=int)
    bond_pairs = _detect_bonds_intralayer(
        positions=positions,
        numbers=numbers,
        is_top_layer=is_top,
        bond_margin=float(bond_margin),
    )

    if side_mode not in {"cell", "slice", "all"}:
        raise ValueError("side_mode must be one of: cell, slice, all")

    if side_mode == "slice":
        y0 = _mode_rounded(coords[:, 1], decimals=side_slice_round_decimals)
        x0 = _mode_rounded(coords[:, 0], decimals=side_slice_round_decimals)
        show_xz = np.abs(coords[:, 1] - y0) <= float(side_slice_width)
        show_yz = np.abs(coords[:, 0] - x0) <= float(side_slice_width)
        tile_ix = tile_iy = None
    elif side_mode == "cell":
        sx, sy = atoms.get_scaled_positions().T[:2]
        scaled = np.vstack([sx, sy]).T
        ix_all, iy_all = _tile_indices_from_scaled(scaled, nx=repeat[0], ny=repeat[1])
        tile_ix, tile_iy = int(side_cell[0]), int(side_cell[1])
        if not (0 <= tile_ix < repeat[0] and 0 <= tile_iy < repeat[1]):
            raise ValueError(f"side_cell must be within [0,{repeat[0]-1}]×[0,{repeat[1]-1}] for repeat={repeat}")
        show_xz = (ix_all == tile_ix) & (iy_all == tile_iy)
        show_yz = show_xz.copy()
        y0 = x0 = None
    else:  # all
        show_xz = np.ones(len(atoms), dtype=bool)
        show_yz = show_xz.copy()
        y0 = x0 = None
        tile_ix = tile_iy = None

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    (ax_xy, ax_xz, ax_yz) = axes

    bonds_xy = bonds == "all"
    bonds_side = bonds in {"all", "side"}

    # Draw bonds first (behind atoms) and only within the shown side-view selection.
    for i, j in bond_pairs:
        pi = coords[i]
        pj = coords[j]
        if bonds_xy:
            ax_xy.plot([pi[0], pj[0]], [pi[1], pj[1]], color=bond_color, lw=bond_lw, zorder=1)
        if bool(show_xz[i]) and bool(show_xz[j]):
            if bonds_side:
                ax_xz.plot([pi[0], pj[0]], [pi[2], pj[2]], color=bond_color, lw=bond_lw, zorder=1)
        if bool(show_yz[i]) and bool(show_yz[j]):
            if bonds_side:
                ax_yz.plot([pi[1], pj[1]], [pi[2], pj[2]], color=bond_color, lw=bond_lw, zorder=1)

    # XY projection
    # Plot bottom layer first for readability
    order = np.argsort(is_top.astype(int))
    ax_xy.scatter(coords[order, 0], coords[order, 1], s=sizes[order], c=rgba[order], linewidths=0.5, edgecolors="k", zorder=2)
    if slice_guides and side_mode == "slice":
        ax_xy.axhline(float(y0), color="0.55", lw=0.8, ls=(0, (3, 3)), zorder=0)
        ax_xy.axvline(float(x0), color="0.55", lw=0.8, ls=(0, (3, 3)), zorder=0)
    if show_cell:
        if cell_grid:
            _draw_cell_grid_xy(ax_xy, cell=cell, basis=basis, nx=repeat[0], ny=repeat[1])
        outline = _cell_outline_xy(cell, basis=basis)
        ax_xy.plot(outline[:, 0], outline[:, 1], color="0.15", lw=2.2, zorder=0)
        if side_mode == "cell" and tile_ix is not None and tile_iy is not None:
            a1_prim = cell[0] / float(repeat[0])
            a2_prim = cell[1] / float(repeat[1])
            tile_outline = _tile_outline_xy(a1_prim=a1_prim, a2_prim=a2_prim, ix=tile_ix, iy=tile_iy, basis=basis)
            ax_xy.plot(tile_outline[:, 0], tile_outline[:, 1], color="0.35", lw=2.0, ls=(0, (4, 3)), zorder=1)
    ax_xy.set_title("XY")
    ax_xy.set_xlabel("x (Å)")
    ax_xy.set_ylabel("y (Å)")
    ax_xy.set_aspect("equal", adjustable="box")

    # XZ projection
    order_xz = order[show_xz[order]]
    ax_xz.scatter(coords[order_xz, 0], coords[order_xz, 2], s=sizes[order_xz], c=rgba[order_xz], linewidths=0.5, edgecolors="k", zorder=2)
    ax_xz.set_title("XZ")
    ax_xz.set_xlabel("x (Å)")
    ax_xz.set_ylabel("z (Å)")
    ax_xz.set_aspect("equal", adjustable="box")

    # YZ projection
    order_yz = order[show_yz[order]]
    ax_yz.scatter(coords[order_yz, 1], coords[order_yz, 2], s=sizes[order_yz], c=rgba[order_yz], linewidths=0.5, edgecolors="k", zorder=2)
    ax_yz.set_title("YZ")
    ax_yz.set_xlabel("y (Å)")
    ax_yz.set_ylabel("z (Å)")
    ax_yz.set_aspect("equal", adjustable="box")

    # Dotted alignment guides (cross-layer in-plane registry)
    if align_lines:
        basis_prim = np.eye(3) if basis_mode == "cart" else _lattice_aligned_basis(atoms_prim.get_cell().array)
        alignments = _detect_cross_layer_alignments(
            atoms_prim=atoms_prim, basis=basis_prim, wrap=wrap, tol_frac=float(align_tol_frac)
        )
        pos_prim = atoms_prim.get_positions()
        a1_prim = atoms_prim.get_cell().array[0]
        a2_prim = atoms_prim.get_cell().array[1]

        a1_xy = (a1_prim @ basis.T)[:2]
        a2_xy = (a2_prim @ basis.T)[:2]

        xz_lines: Dict[float, Tuple[float, float]] = {}
        yz_lines: Dict[float, Tuple[float, float]] = {}

        def _accum_line(lines: Dict[float, Tuple[float, float]], key: float, z0: float, z1: float) -> None:
            zmin, zmax = lines.get(key, (float("inf"), float("-inf")))
            zmin = min(zmin, float(z0), float(z1))
            zmax = max(zmax, float(z0), float(z1))
            lines[key] = (zmin, zmax)

        for top_i, bot_j, shift_ij in alignments:
            for m in range(repeat[0]):
                for n in range(repeat[1]):
                    if side_mode == "cell" and tile_ix is not None and tile_iy is not None:
                        if m != tile_ix or n != tile_iy:
                            continue

                    shift = m * a1_prim + n * a2_prim
                    bot_cart = pos_prim[bot_j] + shift
                    top_cart = pos_prim[top_i] + shift - shift_ij[0] * a1_prim - shift_ij[1] * a2_prim

                    bot_b = bot_cart @ basis.T
                    top_b = top_cart @ basis.T

                    delta_xy = top_b[:2] - bot_b[:2]
                    top_b[:2] = top_b[:2] + _closest_image_shift(delta_xy, a1_xy=a1_xy, a2_xy=a2_xy, search=int(align_search))

                    if side_mode == "slice":
                        if abs(bot_b[1] - float(y0)) <= float(side_slice_width):
                            key = float(np.round(bot_b[0], 3))
                            _accum_line(xz_lines, key=key, z0=float(bot_b[2]), z1=float(top_b[2]))
                        if abs(bot_b[0] - float(x0)) <= float(side_slice_width):
                            key = float(np.round(bot_b[1], 3))
                            _accum_line(yz_lines, key=key, z0=float(bot_b[2]), z1=float(top_b[2]))
                    elif side_mode == "cell":
                        key_x = float(np.round(bot_b[0], 3))
                        _accum_line(xz_lines, key=key_x, z0=float(bot_b[2]), z1=float(top_b[2]))
                        key_y = float(np.round(bot_b[1], 3))
                        _accum_line(yz_lines, key=key_y, z0=float(bot_b[2]), z1=float(top_b[2]))
                    else:  # all
                        key_x = float(np.round(bot_b[0], 3))
                        _accum_line(xz_lines, key=key_x, z0=float(bot_b[2]), z1=float(top_b[2]))
                        key_y = float(np.round(bot_b[1], 3))
                        _accum_line(yz_lines, key=key_y, z0=float(bot_b[2]), z1=float(top_b[2]))

        for x_key, (zmin, zmax) in xz_lines.items():
            ax_xz.plot(
                [x_key, x_key],
                [zmin, zmax],
                color=align_color,
                lw=align_lw,
                ls=(0, (1.2, 2.4)),
                zorder=0,
            )

        for y_key, (zmin, zmax) in yz_lines.items():
            ax_yz.plot(
                [y_key, y_key],
                [zmin, zmax],
                color=align_color,
                lw=align_lw,
                ls=(0, (1.2, 2.4)),
                zorder=0,
            )

    # Tighten limits to the drawn content
    pad = 0.6
    if show_cell:
        outline = _cell_outline_xy(cell, basis=basis)
        ax_xy.set_xlim(outline[:, 0].min() - pad, outline[:, 0].max() + pad)
        ax_xy.set_ylim(outline[:, 1].min() - pad, outline[:, 1].max() + pad)
    else:
        ax_xy.set_xlim(coords[:, 0].min() - pad, coords[:, 0].max() + pad)
        ax_xy.set_ylim(coords[:, 1].min() - pad, coords[:, 1].max() + pad)
    if order_xz.size:
        ax_xz.set_xlim(coords[order_xz, 0].min() - pad, coords[order_xz, 0].max() + pad)
        ax_xz.set_ylim(coords[order_xz, 2].min() - pad, coords[order_xz, 2].max() + pad)
    if order_yz.size:
        ax_yz.set_xlim(coords[order_yz, 1].min() - pad, coords[order_yz, 1].max() + pad)
        ax_yz.set_ylim(coords[order_yz, 2].min() - pad, coords[order_yz, 2].max() + pad)

    if grid:
        ax_xz.grid(True, color="0.92", lw=0.6)
        ax_yz.grid(True, color="0.92", lw=0.6)
    ax_xy.grid(False)

    fig.suptitle(outpath.stem)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def _iter_input_files(inputs: Sequence[str], glob_pattern: str) -> Iterable[Path]:
    if not inputs:
        inputs = ["."]

    for raw in inputs:
        p = Path(raw)
        if p.is_dir():
            yield from sorted(p.glob(glob_pattern))
        else:
            yield p


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot XY/XZ/YZ projections for VASP/ASE-readable structures."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="Input files or directories (default: current directory).",
    )
    parser.add_argument(
        "--glob",
        default="*.vasp",
        help="When an input is a directory, include files matching this glob (default: *.vasp).",
    )
    parser.add_argument(
        "--outdir",
        default="plots",
        help="Output directory (default: plots).",
    )
    parser.add_argument(
        "--repeat",
        nargs=3,
        default=("2", "2", "1"),
        metavar=("NX", "NY", "NZ"),
        help="Repeat supercell before plotting (default: 2 2 1).",
    )
    parser.add_argument(
        "--no-wrap",
        action="store_true",
        help="Do not wrap atoms into the cell before plotting.",
    )
    parser.add_argument(
        "--basis",
        choices=("lattice", "cart"),
        default="cart",
        help="Coordinate basis for plotting (default: Cartesian).",
    )
    parser.add_argument(
        "--no-cell",
        action="store_true",
        help="Do not draw the in-plane unit-cell grid in the XY panel.",
    )
    parser.add_argument(
        "--no-cell-grid",
        action="store_true",
        help="Draw only the outer cell outline (no internal grid).",
    )
    parser.add_argument(
        "--no-grid",
        action="store_true",
        help="Disable the background plot grid.",
    )
    parser.add_argument(
        "--bonds",
        choices=("all", "side", "none"),
        default="side",
        help="Where to draw bonds (default: side).",
    )
    parser.add_argument(
        "--side-mode",
        choices=("cell", "slice", "all"),
        default="cell",
        help="How to populate side views: one primitive cell, a thin slice, or full projection (default: cell).",
    )
    parser.add_argument(
        "--side-cell",
        nargs=2,
        default=("0", "0"),
        metavar=("IX", "IY"),
        help="Primitive tile index (ix iy) to show when --side-mode=cell (default: 0 0).",
    )
    parser.add_argument(
        "--size-scale",
        type=float,
        default=12.0,
        help="Marker size scaling factor (default: 12.0).",
    )
    parser.add_argument(
        "--bond-margin",
        type=float,
        default=0.10,
        help="Bond threshold margin above nearest-neighbor distance (default: 0.10).",
    )
    parser.add_argument(
        "--bond-lw",
        type=float,
        default=1.0,
        help="Bond line width (default: 1.0).",
    )
    parser.add_argument(
        "--bond-color",
        default="0.25",
        help="Bond line color (default: 0.25).",
    )
    parser.add_argument(
        "--side-slice-width",
        type=float,
        default=0.25,
        help="Side-view slice half-width in Å (only for --side-mode=slice; default: 0.25).",
    )
    parser.add_argument(
        "--side-slice-round",
        type=int,
        default=2,
        help="Decimals used to pick the densest side-view slice (only for --side-mode=slice; default: 2).",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="Disable dotted cross-layer alignment guides in side views.",
    )
    parser.add_argument(
        "--align-tol",
        type=float,
        default=5e-4,
        help="Alignment tolerance in fractional in-plane coords (default: 5e-4).",
    )
    parser.add_argument(
        "--align-search",
        type=int,
        default=1,
        help="Closest-image search radius for alignment guides (default: 1).",
    )
    parser.add_argument(
        "--no-slice-guides",
        action="store_true",
        help="Do not draw slice guides (x=x0, y=y0) in the XY panel.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Output DPI (default: 200).",
    )

    args = parser.parse_args()

    repeat = _parse_repeat(args.repeat)
    outdir = Path(args.outdir)
    wrap = not args.no_wrap
    show_cell = not args.no_cell
    cell_grid = not args.no_cell_grid
    grid = not args.no_grid
    side_cell = (int(args.side_cell[0]), int(args.side_cell[1]))

    files = list(_iter_input_files(args.inputs, glob_pattern=args.glob))
    if not files:
        raise SystemExit(f"No input files found for glob {args.glob!r}")

    for path in files:
        atoms = read(str(path))
        outpath = outdir / f"{path.stem}.projections.png"
        plot_projections(
            atoms,
            outpath=outpath,
            repeat=repeat,
            wrap=wrap,
            basis_mode=args.basis,
            show_cell=show_cell,
            cell_grid=cell_grid,
            grid=grid,
            size_scale=args.size_scale,
            bonds=args.bonds,
            side_mode=args.side_mode,
            side_cell=side_cell,
            bond_margin=args.bond_margin,
            bond_color=args.bond_color,
            bond_lw=args.bond_lw,
            side_slice_width=args.side_slice_width,
            side_slice_round_decimals=args.side_slice_round,
            slice_guides=not args.no_slice_guides,
            align_lines=not args.no_align,
            align_tol_frac=args.align_tol,
            align_search=args.align_search,
            dpi=args.dpi,
        )


if __name__ == "__main__":
    main()
