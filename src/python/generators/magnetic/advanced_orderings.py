"""
magnetic/advanced_orderings.py - Advanced Magnetic Orderings

Multi-q magnetic structures and DMI-induced spin textures:
- Single-q, 2-q, 3-q, 4-q magnetic orderings
- DMI cycloids and spin spirals
- Skyrmion lattices
- Conical phases
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """Convert pymatgen Structure to dictionary."""
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_multi_q_structure(
    material: str = "Fe",
    q_type: str = "single_q",
    q_vectors: List[List[float]] = None,
    moment_magnitude: float = 2.2,
    supercell: List[int] = None,
    a: float = 2.87
) -> Dict[str, Any]:
    """
    Generate multi-q magnetic structure.

    Multi-q magnetic structures arise from superposition of spin density waves
    with different propagation vectors. Common in frustrated magnets and
    heavy-fermion systems.

    Args:
        material: Element or compound name
        q_type: 'single_q', '2q', '3q', '4q'
        q_vectors: Custom q-vectors (fractional coordinates in reciprocal space)
        moment_magnitude: Magnetic moment per atom (Bohr magnetons)
        supercell: Supercell dimensions [nx, ny, nz]
        a: Lattice parameter (Angstrom)

    Returns:
        Structure with multi-q magnetic ordering

    Examples:
        >>> result = generate_multi_q_structure('Fe', q_type='3q')
        >>> result['n_q_vectors']
        3
    """
    if supercell is None:
        supercell = [4, 4, 4]

    # Default q-vectors for different orderings
    Q_VECTORS = {
        "single_q": [[0.5, 0, 0]],  # Simple AFM
        "2q": [[0.5, 0, 0], [0, 0.5, 0]],  # 2-q ordering
        "3q": [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],  # Triple-q (tetrahedral)
        "4q": [[0.25, 0.25, 0], [0.25, -0.25, 0],
               [-0.25, 0.25, 0], [-0.25, -0.25, 0]],  # 4-q (square skyrmion lattice)
    }

    if q_type not in Q_VECTORS and q_vectors is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_Q_TYPE",
                "message": f"Unknown q_type: {q_type}",
                "available": list(Q_VECTORS.keys())
            }
        }

    q_vecs = q_vectors if q_vectors else Q_VECTORS[q_type]

    # Create BCC structure
    lattice = Lattice.cubic(a)
    base_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]

    nx, ny, nz = supercell
    species = []
    coords = []
    spins = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for base in base_coords:
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz

                    species.append(material)
                    coords.append([x % 1, y % 1, z % 1])

                    # Position in real space (in units of lattice parameter)
                    r = np.array([i + base[0], j + base[1], k + base[2]])

                    # Calculate multi-q spin texture
                    spin = np.zeros(3)
                    for q_idx, q in enumerate(q_vecs):
                        q_arr = np.array(q)
                        phase = 2 * np.pi * np.dot(q_arr, r)

                        if len(q_vecs) == 1:
                            # Single-q: simple spiral
                            spin[0] += moment_magnitude * np.cos(phase)
                            spin[1] += moment_magnitude * np.sin(phase)
                        elif len(q_vecs) == 2:
                            # 2-q: orthogonal modulation
                            if q_idx == 0:
                                spin[0] += moment_magnitude * np.cos(phase) / np.sqrt(2)
                            else:
                                spin[1] += moment_magnitude * np.cos(phase) / np.sqrt(2)
                        elif len(q_vecs) == 3:
                            # 3-q: tetrahedral (hedgehog) structure
                            axis = q_arr / np.linalg.norm(q_arr) if np.linalg.norm(q_arr) > 0 else np.array([1, 0, 0])
                            spin += axis * moment_magnitude * np.cos(phase) / np.sqrt(3)
                        else:
                            # 4-q: complex texture
                            angle = 2 * np.pi * q_idx / 4
                            spin[0] += moment_magnitude * np.cos(phase) * np.cos(angle) / 2
                            spin[1] += moment_magnitude * np.cos(phase) * np.sin(angle) / 2
                            spin[2] += moment_magnitude * np.sin(phase) / 4

                    # Normalize to moment magnitude
                    spin_mag = np.linalg.norm(spin)
                    if spin_mag > 0:
                        spin = spin * moment_magnitude / spin_mag

                    spins.append(spin.tolist())

    # Create structure
    supercell_lattice = Lattice(lattice.matrix * np.diag(supercell))
    structure = Structure(supercell_lattice, species, coords)

    return {
        "success": True,
        "material": material,
        "q_type": q_type,
        "q_vectors": q_vecs,
        "n_q_vectors": len(q_vecs),
        "moment_magnitude_muB": moment_magnitude,
        "supercell": supercell,
        "n_atoms": len(structure),
        "magnetic_moments": spins,
        "structure": structure_to_dict(structure),
        "description": _get_multi_q_description(q_type)
    }


def _get_multi_q_description(q_type: str) -> str:
    """Get description for multi-q ordering type."""
    descriptions = {
        "single_q": "Single-q spin density wave or spiral",
        "2q": "Double-q ordering with orthogonal modulation vectors",
        "3q": "Triple-q (tetrahedral/hedgehog) magnetic structure",
        "4q": "Quadruple-q ordering (square skyrmion lattice precursor)"
    }
    return descriptions.get(q_type, "Multi-q magnetic ordering")


def generate_dmi_cycloid(
    material: str = "BiFeO3",
    propagation_direction: List[float] = None,
    cycloid_period_nm: float = 62.0,
    moment_magnitude: float = 4.0,
    dmi_strength: float = 0.1,
    supercell_periods: int = 2,
    easy_axis: str = "z"
) -> Dict[str, Any]:
    """
    Generate DMI-induced spin cycloid structure.

    Dzyaloshinskii-Moriya Interaction (DMI) arises from spin-orbit coupling
    in non-centrosymmetric systems. It causes spin canting and can stabilize
    cycloidal spin spirals, as in multiferroic BiFeO3.

    Args:
        material: Material name ('BiFeO3', 'FeGe', 'MnSi', etc.)
        propagation_direction: Cycloid propagation direction [h, k, l]
        cycloid_period_nm: Period of the cycloid in nanometers
        moment_magnitude: Magnetic moment per atom (Bohr magnetons)
        dmi_strength: Relative DMI strength (affects cone angle)
        supercell_periods: Number of cycloid periods in supercell
        easy_axis: Magnetic easy axis ('x', 'y', 'z')

    Returns:
        Structure with DMI cycloid spin texture

    Examples:
        >>> result = generate_dmi_cycloid('BiFeO3', cycloid_period_nm=62)
        >>> result['cycloid_type']
        'spin_cycloid'
    """
    if propagation_direction is None:
        propagation_direction = [1, 1, 0]  # BiFeO3 default

    # Material database for DMI materials
    DMI_MATERIALS = {
        "BiFeO3": {"a": 5.58, "c": 13.87, "structure": "rhombohedral",
                   "moment": 4.0, "period_nm": 62, "multiferroic": True},
        "FeGe": {"a": 4.70, "structure": "B20", "moment": 1.0,
                 "period_nm": 70, "skyrmion": True},
        "MnSi": {"a": 4.56, "structure": "B20", "moment": 0.4,
                 "period_nm": 18, "skyrmion": True},
        "Cu2OSeO3": {"a": 8.92, "structure": "cubic", "moment": 1.0,
                     "period_nm": 60, "insulator": True},
        "Cr1/3NbS2": {"a": 5.74, "c": 12.10, "structure": "hexagonal",
                      "moment": 3.0, "period_nm": 48, "chiral": True},
    }

    if material in DMI_MATERIALS:
        mat_info = DMI_MATERIALS[material]
        a = mat_info["a"]
        moment_magnitude = mat_info.get("moment", moment_magnitude)
        cycloid_period_nm = mat_info.get("period_nm", cycloid_period_nm)
    else:
        a = 4.0
        mat_info = {}

    # Calculate supercell needed for cycloid
    # period in lattice units
    period_lattice = cycloid_period_nm * 10 / a  # Convert nm to Angstrom
    n_cells = max(int(period_lattice * supercell_periods), 8)

    # Propagation vector normalized
    prop_dir = np.array(propagation_direction, dtype=float)
    prop_dir = prop_dir / np.linalg.norm(prop_dir)

    # Create base structure
    lattice = Lattice.cubic(a)
    base_coords = [[0, 0, 0]]  # Simple cubic for now

    species = []
    coords = []
    spins = []

    # Determine supercell shape based on propagation direction
    if abs(prop_dir[0]) > 0.5:
        nx, ny, nz = n_cells, 2, 2
    elif abs(prop_dir[1]) > 0.5:
        nx, ny, nz = 2, n_cells, 2
    else:
        nx, ny, nz = 2, 2, n_cells

    # Easy axis vector
    easy_axis_vec = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}[easy_axis]
    easy = np.array(easy_axis_vec, dtype=float)

    # Calculate rotation plane (perpendicular to propagation)
    perp1 = np.cross(prop_dir, [0, 0, 1])
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(prop_dir, [0, 1, 0])
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(prop_dir, perp1)

    elem = material.replace("BiFeO3", "Fe").replace("FeGe", "Fe").replace("MnSi", "Mn")
    if len(elem) > 2:
        elem = elem[:2]

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for base in base_coords:
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz

                    species.append(elem)
                    coords.append([x % 1, y % 1, z % 1])

                    # Position along propagation direction
                    r = np.array([i + base[0], j + base[1], k + base[2]]) * a
                    r_prop = np.dot(r, prop_dir)

                    # Cycloid phase (convert period to Angstrom)
                    phase = 2 * np.pi * r_prop / (cycloid_period_nm * 10)

                    # Cycloid: spins rotate in plane containing easy axis and propagation
                    # m = m0 * (cos(qr) * e_easy + sin(qr) * e_prop)
                    spin = moment_magnitude * (np.cos(phase) * easy +
                                               np.sin(phase) * prop_dir)

                    # Add small canting from DMI
                    canting = dmi_strength * moment_magnitude * np.sin(phase) * perp1
                    spin = spin + canting

                    # Normalize
                    spin = spin * moment_magnitude / np.linalg.norm(spin)
                    spins.append(spin.tolist())

    # Create structure
    supercell_lattice = Lattice(lattice.matrix * np.diag([nx, ny, nz]))
    structure = Structure(supercell_lattice, species, coords)

    return {
        "success": True,
        "material": material,
        "cycloid_type": "spin_cycloid",
        "propagation_direction": propagation_direction,
        "cycloid_period_nm": cycloid_period_nm,
        "moment_magnitude_muB": moment_magnitude,
        "dmi_strength": dmi_strength,
        "easy_axis": easy_axis,
        "supercell": [nx, ny, nz],
        "n_atoms": len(structure),
        "magnetic_moments": spins,
        "is_multiferroic": mat_info.get("multiferroic", False),
        "can_host_skyrmions": mat_info.get("skyrmion", False),
        "structure": structure_to_dict(structure),
        "description": f"DMI-induced spin cycloid with {cycloid_period_nm:.1f} nm period"
    }


def generate_conical_phase(
    material: str = "MnSi",
    cone_angle_deg: float = 45.0,
    helix_period_nm: float = 18.0,
    field_direction: List[float] = None,
    moment_magnitude: float = 0.4,
    supercell_periods: int = 2
) -> Dict[str, Any]:
    """
    Generate conical magnetic phase.

    Conical phase arises when a magnetic field is applied along the helix
    axis in chiral magnets. Spins precess on a cone surface combining
    helical rotation with field-aligned component.

    Args:
        material: Material name
        cone_angle_deg: Half-angle of the spin cone (0=FM, 90=flat helix)
        helix_period_nm: Helical modulation period
        field_direction: Applied field direction [h, k, l]
        moment_magnitude: Magnetic moment per atom
        supercell_periods: Number of helix periods

    Returns:
        Structure with conical spin texture
    """
    if field_direction is None:
        field_direction = [0, 0, 1]

    # Material parameters
    CONICAL_MATERIALS = {
        "MnSi": {"a": 4.56, "moment": 0.4, "period_nm": 18},
        "FeGe": {"a": 4.70, "moment": 1.0, "period_nm": 70},
        "Fe0.5Co0.5Si": {"a": 4.50, "moment": 0.6, "period_nm": 30},
    }

    if material in CONICAL_MATERIALS:
        mat_info = CONICAL_MATERIALS[material]
        a = mat_info["a"]
        moment_magnitude = mat_info.get("moment", moment_magnitude)
        helix_period_nm = mat_info.get("period_nm", helix_period_nm)
    else:
        a = 4.5
        mat_info = {}

    # Field direction normalized
    H_dir = np.array(field_direction, dtype=float)
    H_dir = H_dir / np.linalg.norm(H_dir)

    # Perpendicular basis
    perp1 = np.cross(H_dir, [1, 0, 0])
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(H_dir, [0, 1, 0])
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(H_dir, perp1)

    cone_angle = np.radians(cone_angle_deg)

    # Supercell for helix
    period_lattice = helix_period_nm * 10 / a
    n_cells = max(int(period_lattice * supercell_periods), 8)

    # Determine supercell shape
    if abs(H_dir[2]) > 0.5:
        nx, ny, nz = 2, 2, n_cells
    elif abs(H_dir[1]) > 0.5:
        nx, ny, nz = 2, n_cells, 2
    else:
        nx, ny, nz = n_cells, 2, 2

    lattice = Lattice.cubic(a)
    base_coords = [[0, 0, 0]]

    species = []
    coords = []
    spins = []

    elem = material[:2] if len(material) > 2 else material

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for base in base_coords:
                    x = (i + base[0]) / nx
                    y = (j + base[1]) / ny
                    z = (k + base[2]) / nz

                    species.append(elem)
                    coords.append([x % 1, y % 1, z % 1])

                    # Position along helix axis
                    r = np.array([i + base[0], j + base[1], k + base[2]]) * a
                    r_H = np.dot(r, H_dir)

                    # Helix phase
                    phase = 2 * np.pi * r_H / (helix_period_nm * 10)

                    # Conical spin: component along H + rotating perpendicular
                    spin_H = moment_magnitude * np.cos(cone_angle) * H_dir
                    spin_perp = moment_magnitude * np.sin(cone_angle) * (
                        np.cos(phase) * perp1 + np.sin(phase) * perp2)
                    spin = spin_H + spin_perp

                    spins.append(spin.tolist())

    supercell_lattice = Lattice(lattice.matrix * np.diag([nx, ny, nz]))
    structure = Structure(supercell_lattice, species, coords)

    return {
        "success": True,
        "material": material,
        "phase": "conical",
        "cone_angle_deg": cone_angle_deg,
        "helix_period_nm": helix_period_nm,
        "field_direction": field_direction,
        "moment_magnitude_muB": moment_magnitude,
        "supercell": [nx, ny, nz],
        "n_atoms": len(structure),
        "magnetic_moments": spins,
        "structure": structure_to_dict(structure),
        "description": f"Conical phase with {cone_angle_deg}° cone angle"
    }


def generate_skyrmion_lattice(
    material: str = "MnSi",
    skyrmion_type: str = "bloch",
    lattice_type: str = "hexagonal",
    skyrmion_radius_nm: float = 10.0,
    moment_magnitude: float = 0.4,
    supercell_size: int = 20
) -> Dict[str, Any]:
    """
    Generate magnetic skyrmion lattice.

    Skyrmions are topologically protected spin textures. They form
    lattices in chiral magnets under applied magnetic field.

    Args:
        material: Material name
        skyrmion_type: 'bloch' (chiral), 'neel' (polar), 'antiskyrmion'
        lattice_type: 'hexagonal' (triangular), 'square'
        skyrmion_radius_nm: Characteristic radius
        moment_magnitude: Magnetic moment per atom
        supercell_size: Supercell size (atoms per side)

    Returns:
        Structure with skyrmion lattice
    """
    # Material parameters
    SKYRMION_MATERIALS = {
        "MnSi": {"a": 4.56, "moment": 0.4, "radius_nm": 9},
        "FeGe": {"a": 4.70, "moment": 1.0, "radius_nm": 35},
        "Cu2OSeO3": {"a": 8.92, "moment": 1.0, "radius_nm": 30},
        "PdFe_bilayer": {"a": 2.75, "moment": 2.5, "radius_nm": 1},
    }

    if material in SKYRMION_MATERIALS:
        mat_info = SKYRMION_MATERIALS[material]
        a = mat_info["a"]
        moment_magnitude = mat_info.get("moment", moment_magnitude)
        skyrmion_radius_nm = mat_info.get("radius_nm", skyrmion_radius_nm)
    else:
        a = 4.0
        mat_info = {}

    n = supercell_size
    lattice = Lattice.cubic(a)

    # Skyrmion lattice vectors
    if lattice_type == "hexagonal":
        # Triangular lattice of skyrmions
        sk_period = skyrmion_radius_nm * 10 * 2.5  # Angstrom
        sk_a1 = np.array([sk_period, 0, 0])
        sk_a2 = np.array([sk_period * 0.5, sk_period * np.sqrt(3) / 2, 0])
    else:  # square
        sk_period = skyrmion_radius_nm * 10 * 2
        sk_a1 = np.array([sk_period, 0, 0])
        sk_a2 = np.array([0, sk_period, 0])

    species = []
    coords = []
    spins = []

    elem = material[:2] if len(material) > 2 else material
    radius_A = skyrmion_radius_nm * 10  # Convert to Angstrom

    for i in range(n):
        for j in range(n):
            x = float(i) / n
            y = float(j) / n

            species.append(elem)
            coords.append([x, y, 0.5])

            # Real space position
            r_x = i * a
            r_y = j * a

            # Find distance to nearest skyrmion center
            min_dist = float('inf')
            for ni in range(-1, 2):
                for nj in range(-1, 2):
                    center = ni * sk_a1 + nj * sk_a2
                    dist = np.sqrt((r_x - center[0])**2 + (r_y - center[1])**2)
                    if dist < min_dist:
                        min_dist = dist
                        dx = r_x - center[0]
                        dy = r_y - center[1]

            # Skyrmion profile
            rho = min_dist / radius_A
            theta = np.pi * np.exp(-rho)  # Spin angle from z-axis
            phi = np.arctan2(dy, dx)  # Azimuthal angle

            if skyrmion_type == "bloch":
                # Bloch skyrmion: phi_spin = phi + pi/2
                phi_spin = phi + np.pi / 2
            elif skyrmion_type == "neel":
                # Neel skyrmion: phi_spin = phi
                phi_spin = phi
            else:  # antiskyrmion
                # Antiskyrmion: phi_spin = -phi
                phi_spin = -phi

            spin = moment_magnitude * np.array([
                np.sin(theta) * np.cos(phi_spin),
                np.sin(theta) * np.sin(phi_spin),
                np.cos(theta)
            ])

            spins.append(spin.tolist())

    supercell_lattice = Lattice(lattice.matrix * np.diag([n, n, 1]))
    structure = Structure(supercell_lattice, species, coords)

    return {
        "success": True,
        "material": material,
        "phase": "skyrmion_lattice",
        "skyrmion_type": skyrmion_type,
        "lattice_type": lattice_type,
        "skyrmion_radius_nm": skyrmion_radius_nm,
        "moment_magnitude_muB": moment_magnitude,
        "supercell": [n, n, 1],
        "n_atoms": len(structure),
        "magnetic_moments": spins,
        "topological_charge_per_skyrmion": -1 if skyrmion_type != "antiskyrmion" else 1,
        "structure": structure_to_dict(structure),
        "description": f"{skyrmion_type.capitalize()} skyrmion lattice ({lattice_type})"
    }


def get_available_advanced_orderings() -> Dict[str, Any]:
    """Get list of available advanced magnetic orderings."""
    return {
        "success": True,
        "orderings": {
            "multi_q": {
                "types": ["single_q", "2q", "3q", "4q"],
                "description": "Multi-q spin density wave structures"
            },
            "dmi_cycloid": {
                "materials": ["BiFeO3", "FeGe", "MnSi", "Cu2OSeO3", "Cr1/3NbS2"],
                "description": "DMI-induced spin cycloids"
            },
            "conical": {
                "materials": ["MnSi", "FeGe", "Fe0.5Co0.5Si"],
                "description": "Conical phase under applied field"
            },
            "skyrmion_lattice": {
                "types": ["bloch", "neel", "antiskyrmion"],
                "lattices": ["hexagonal", "square"],
                "materials": ["MnSi", "FeGe", "Cu2OSeO3", "PdFe_bilayer"],
                "description": "Topological skyrmion lattices"
            }
        }
    }
