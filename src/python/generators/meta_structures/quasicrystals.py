"""
meta_structures/quasicrystals.py - Quasicrystal Structure Generation

Provides generation of aperiodic crystalline structures:
- Penrose tilings (2D)
- Fibonacci approximants
- Icosahedral quasicrystals (3D)
- Decagonal quasicrystals
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


def generate_quasicrystal(
    quasicrystal_type: str = "penrose",
    approximant_order: int = 5,
    element: str = "Al",
    vacuum: float = 15.0
) -> Dict[str, Any]:
    """
    Generate quasicrystal structure or approximant.

    Quasicrystals have long-range order but no translational periodicity.
    This generates periodic approximants for DFT calculations.

    Args:
        quasicrystal_type: "penrose", "fibonacci", "icosahedral", "decagonal"
        approximant_order: Order of Fibonacci approximant (higher = closer to true QC)
        element: Element for atomic positions
        vacuum: Vacuum for 2D structures (Angstroms)

    Returns:
        Quasicrystal approximant structure
    """
    if quasicrystal_type == "penrose":
        return generate_penrose_tiling(order=approximant_order, element=element, vacuum=vacuum)
    elif quasicrystal_type == "fibonacci":
        return generate_fibonacci_approximant(order=approximant_order, element=element)
    elif quasicrystal_type == "icosahedral":
        return generate_icosahedral_quasicrystal(approximant_order=approximant_order, elements=[element])
    elif quasicrystal_type == "decagonal":
        return generate_decagonal_quasicrystal(order=approximant_order, element=element, vacuum=vacuum)
    else:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_TYPE",
                "message": f"Unknown quasicrystal type: {quasicrystal_type}",
                "available": ["penrose", "fibonacci", "icosahedral", "decagonal"]
            }
        }


def generate_penrose_tiling(
    order: int = 5,
    element: str = "C",
    bond_length: float = 1.42,
    vacuum: float = 15.0,
    vertex_decoration: bool = True
) -> Dict[str, Any]:
    """
    Generate Penrose tiling (P3 rhombus tiling).

    Creates a 2D quasicrystalline structure based on Penrose tiling
    using the inflation/deflation method.

    Args:
        order: Number of inflation steps (determines size)
        element: Element for vertex atoms
        bond_length: Distance between vertices (Angstroms)
        vacuum: Vacuum above/below 2D layer
        vertex_decoration: Place atoms at tile vertices

    Returns:
        Penrose tiling structure
    """
    from pymatgen.core import Structure, Lattice

    # Golden ratio
    phi = (1 + np.sqrt(5)) / 2

    # Penrose tiling angles
    angles = [i * 2 * np.pi / 5 for i in range(5)]  # 72 degree rotations

    # Unit vectors for Penrose tiling
    e = [np.array([np.cos(a), np.sin(a)]) for a in angles]

    # Generate vertices using de Bruijn's method (multigrid)
    # Use 5 families of parallel lines
    def get_vertices(n_lines: int = 10):
        vertices = set()
        offsets = [0.2, 0.4, 0.6, 0.8, 1.0]  # Offsets for 5 grids

        for i in range(5):
            for j in range(i + 1, 5):
                for ki in range(-n_lines, n_lines + 1):
                    for kj in range(-n_lines, n_lines + 1):
                        # Intersection of two grid lines
                        # Line i: r · e[i] = ki + offset[i]
                        # Line j: r · e[j] = kj + offset[j]
                        det = e[i][0] * e[j][1] - e[i][1] * e[j][0]
                        if abs(det) > 1e-10:
                            x = ((ki + offsets[i]) * e[j][1] - (kj + offsets[j]) * e[i][1]) / det
                            y = ((kj + offsets[j]) * e[i][0] - (ki + offsets[i]) * e[j][0]) / det
                            # Round to avoid duplicates
                            vertices.add((round(x, 6), round(y, 6)))
        return list(vertices)

    n_lines = order + 3
    vertices = get_vertices(n_lines)

    # Filter vertices within a reasonable range
    max_radius = order * 3 * bond_length
    filtered = [(x, y) for x, y in vertices if x**2 + y**2 < max_radius**2]

    if len(filtered) == 0:
        filtered = vertices[:100]  # Fallback

    # Scale by bond length
    coords_2d = [(x * bond_length, y * bond_length) for x, y in filtered]

    # Convert to 3D with vacuum
    coords_3d = [[x, y, vacuum/2] for x, y in coords_2d]

    # Create bounding box lattice
    if coords_3d:
        xs = [c[0] for c in coords_3d]
        ys = [c[1] for c in coords_3d]
        a = max(xs) - min(xs) + 5
        b = max(ys) - min(ys) + 5

        # Shift to positive coordinates
        min_x, min_y = min(xs), min(ys)
        coords_3d = [[c[0] - min_x + 2.5, c[1] - min_y + 2.5, c[2]] for c in coords_3d]
    else:
        a, b = 20, 20

    lattice = Lattice.orthorhombic(a, b, vacuum)

    struct = Structure(
        lattice=lattice,
        species=[element] * len(coords_3d),
        coords=coords_3d,
        coords_are_cartesian=True
    )

    return {
        "success": True,
        "type": "penrose_tiling",
        "order": order,
        "symmetry": "5-fold (quasiperiodic)",
        "n_atoms": len(struct),
        "golden_ratio": phi,
        "lattice_type": "approximant_orthorhombic",
        "notes": [
            "This is a periodic approximant to true Penrose tiling",
            "Higher order = better approximation to quasicrystal",
            "True Penrose tiling is aperiodic"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_fibonacci_approximant(
    order: int = 6,
    element: str = "Al",
    a_short: float = 2.5,
    a_long: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate 1D Fibonacci chain approximant extended to 3D.

    The Fibonacci sequence generates quasiperiodic ordering.
    Approximants use F(n)/F(n-1) ratio approaching golden ratio.

    Args:
        order: Fibonacci order (F_n atoms in chain)
        element: Element symbol
        a_short: Short bond length (S)
        a_long: Long bond length (L), defaults to phi * a_short

    Returns:
        Fibonacci approximant structure
    """
    from pymatgen.core import Structure, Lattice

    phi = (1 + np.sqrt(5)) / 2

    if a_long is None:
        a_long = phi * a_short

    # Generate Fibonacci sequence
    def fib_sequence(n):
        """Generate Fibonacci word: L -> LS, S -> L"""
        if n == 0:
            return "L"
        elif n == 1:
            return "LS"
        else:
            prev2 = "L"
            prev1 = "LS"
            for _ in range(2, n + 1):
                current = prev1 + prev2
                prev2 = prev1
                prev1 = current
            return current

    seq = fib_sequence(order)

    # Calculate positions along chain
    positions = [0.0]
    for char in seq:
        if char == "L":
            positions.append(positions[-1] + a_long)
        else:
            positions.append(positions[-1] + a_short)

    # Remove last (it's beyond the cell)
    positions = positions[:-1]

    # Cell length
    chain_length = sum(a_long if c == "L" else a_short for c in seq)

    # Create 3D structure (chain along x, vacuum in y and z)
    b = 10.0  # Vacuum
    c = 10.0

    lattice = Lattice.orthorhombic(chain_length, b, c)

    coords = [[x, b/2, c/2] for x in positions]

    struct = Structure(
        lattice=lattice,
        species=[element] * len(coords),
        coords=coords,
        coords_are_cartesian=True
    )

    # Fibonacci numbers
    fib_nums = [1, 1]
    for _ in range(order):
        fib_nums.append(fib_nums[-1] + fib_nums[-2])

    return {
        "success": True,
        "type": "fibonacci_approximant",
        "order": order,
        "fibonacci_numbers": fib_nums[:order+2],
        "sequence": seq[:50] + ("..." if len(seq) > 50 else ""),
        "sequence_length": len(seq),
        "n_L": seq.count("L"),
        "n_S": seq.count("S"),
        "ratio_L_S": seq.count("L") / max(seq.count("S"), 1),
        "golden_ratio": phi,
        "a_short": a_short,
        "a_long": a_long,
        "chain_length_A": chain_length,
        "n_atoms": len(struct),
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_icosahedral_quasicrystal(
    approximant_order: int = 1,
    elements: List[str] = ["Al", "Mn"],
    a: float = 12.0
) -> Dict[str, Any]:
    """
    Generate icosahedral quasicrystal approximant.

    Creates periodic approximants to icosahedral quasicrystals
    like Al-Mn, Al-Cu-Fe, etc.

    Args:
        approximant_order: 1/1, 2/1, 3/2 approximant (Fibonacci ratios)
        elements: Elements for decoration
        a: Lattice parameter (Angstroms)

    Returns:
        Icosahedral approximant structure
    """
    from pymatgen.core import Structure, Lattice

    phi = (1 + np.sqrt(5)) / 2

    # Icosahedral approximants have BCC-related structures
    # 1/1 approximant is related to W structure (Im-3)

    if approximant_order == 1:
        # 1/1 approximant - alpha-AlMnSi type
        # Body-centered structure with icosahedral clusters

        lattice = Lattice.cubic(a)

        # Icosahedron vertices (normalized)
        ico_verts = []
        for i in [-1, 1]:
            for j in [-1, 1]:
                ico_verts.append([0, i/phi, j])
                ico_verts.append([i/phi, j, 0])
                ico_verts.append([i, 0, j/phi])

        # Normalize and scale
        ico_verts = np.array(ico_verts)
        ico_verts = ico_verts / np.linalg.norm(ico_verts[0]) * 0.2  # Scale to fit in cell

        # Place icosahedra at BCC positions
        centers = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

        coords = []
        species = []

        for center in centers:
            # Add shell atoms
            for vert in ico_verts:
                coord = np.array(center) + vert
                coords.append(coord.tolist())
                species.append(elements[0])

            # Add center atom
            coords.append(list(center))
            species.append(elements[-1] if len(elements) > 1 else elements[0])

    elif approximant_order == 2:
        # 2/1 approximant - larger cell
        a_21 = a * phi

        lattice = Lattice.cubic(a_21)

        # More complex icosahedral cluster packing
        coords = []
        species = []

        # Multiple cluster centers
        cluster_centers = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.5, 0.5, 0.5]
        ]

        # Simplified icosahedron
        ico_verts = []
        for i in [-1, 1]:
            ico_verts.append([0, 0.1*i, 0.16*i])
            ico_verts.append([0.1*i, 0.16*i, 0])
            ico_verts.append([0.16*i, 0, 0.1*i])

        for center in cluster_centers:
            for vert in ico_verts:
                coord = np.array(center) + np.array(vert)
                coords.append(coord.tolist())
                species.append(elements[0])
            coords.append(list(center))
            species.append(elements[-1] if len(elements) > 1 else elements[0])

    else:
        # Higher order - use general formula
        lattice = Lattice.cubic(a * phi ** (approximant_order - 1))
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        species = [elements[0], elements[-1] if len(elements) > 1 else elements[0]]

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "type": "icosahedral_approximant",
        "approximant_order": f"{approximant_order}/1",
        "symmetry": "Im-3 (approximant)",
        "true_symmetry": "icosahedral (m-3-5)",
        "n_atoms": len(struct),
        "golden_ratio": phi,
        "elements": elements,
        "notes": [
            "Periodic approximant to true icosahedral quasicrystal",
            "Real quasicrystals have sharp diffraction spots at irrational positions",
            f"This {approximant_order}/1 approximant captures local icosahedral order"
        ],
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_decagonal_quasicrystal(
    order: int = 3,
    element: str = "Al",
    c: float = 4.0,
    vacuum: float = 0.0
) -> Dict[str, Any]:
    """
    Generate decagonal quasicrystal approximant.

    Decagonal quasicrystals are periodic along one axis (c)
    but quasiperiodic in the ab plane.

    Args:
        order: Approximant order
        element: Element symbol
        c: Periodic axis length
        vacuum: Vacuum (0 for bulk, >0 for slab)

    Returns:
        Decagonal approximant structure
    """
    from pymatgen.core import Structure, Lattice

    phi = (1 + np.sqrt(5)) / 2

    # Decagonal symmetry - 10-fold in plane
    angles = [i * np.pi / 5 for i in range(10)]

    # Generate vertices on concentric rings
    coords = []
    radii = [2.5, 2.5 * phi, 2.5 * phi ** 2]

    for r in radii[:order]:
        for angle in angles:
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            coords.append([x, y])

    # Add center
    coords.insert(0, [0, 0])

    # Create 3D coords
    coords_3d = [[x, y, c/2] for x, y in coords]

    # Bounding box
    max_r = max(radii[:order]) + 2
    a = 2 * max_r
    b = 2 * max_r

    # Shift to positive
    coords_3d = [[x + max_r, y + max_r, z] for x, y, z in coords_3d]

    if vacuum > 0:
        lattice = Lattice.orthorhombic(a, b, c + vacuum)
    else:
        lattice = Lattice.orthorhombic(a, b, c)

    struct = Structure(
        lattice=lattice,
        species=[element] * len(coords_3d),
        coords=coords_3d,
        coords_are_cartesian=True
    )

    return {
        "success": True,
        "type": "decagonal_approximant",
        "order": order,
        "symmetry": "10-fold in plane, periodic along c",
        "point_group": "10/mmm",
        "n_atoms": len(struct),
        "c_axis_A": c,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }
