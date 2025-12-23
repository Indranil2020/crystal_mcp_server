"""
output_formats/visualization.py - Visualization Export Formats

Exports structures to visualization-specific formats:
- XSF (XCrySDen Structure Format)
- CUBE (Gaussian Cube for volumetric data)
- VESTA (Enhanced format for VESTA viewer)
- JSON (for web-based 3D viewers like Three.js)

Scientific References:
- XSF format: http://www.xcrysden.org/doc/XSF.html
- Gaussian CUBE: https://gaussian.com/cubegen/
- VESTA: K. Momma and F. Izumi, J. Appl. Crystallogr. 44, 1272 (2011)
"""

from typing import Dict, Any, List, Optional, Union, Tuple
import numpy as np
import json


# Element data for visualization
ELEMENT_COLORS = {
    "H": "#FFFFFF", "He": "#D9FFFF", "Li": "#CC80FF", "Be": "#C2FF00",
    "B": "#FFB5B5", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "Ne": "#B3E3F5", "Na": "#AB5CF2", "Mg": "#8AFF00",
    "Al": "#BFA6A6", "Si": "#F0C8A0", "P": "#FF8000", "S": "#FFFF30",
    "Cl": "#1FF01F", "Ar": "#80D1E3", "K": "#8F40D4", "Ca": "#3DFF00",
    "Sc": "#E6E6E6", "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7",
    "Mn": "#9C7AC7", "Fe": "#E06633", "Co": "#F090A0", "Ni": "#50D050",
    "Cu": "#C88033", "Zn": "#7D80B0", "Ga": "#C28F8F", "Ge": "#668F8F",
    "As": "#BD80E3", "Se": "#FFA100", "Br": "#A62929", "Kr": "#5CB8D1",
    "Rb": "#702EB0", "Sr": "#00FF00", "Y": "#94FFFF", "Zr": "#94E0E0",
    "Nb": "#73C2C9", "Mo": "#54B5B5", "Tc": "#3B9E9E", "Ru": "#248F8F",
    "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F",
    "In": "#A67573", "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00",
    "I": "#940094", "Xe": "#429EB0", "Cs": "#57178F", "Ba": "#00C900",
    "La": "#70D4FF", "Ce": "#FFFFC7", "Pr": "#D9FFC7", "Nd": "#C7FFC7",
    "Pm": "#A3FFC7", "Sm": "#8FFFC7", "Eu": "#61FFC7", "Gd": "#45FFC7",
    "Tb": "#30FFC7", "Dy": "#1FFFC7", "Ho": "#00FF9C", "Er": "#00E675",
    "Tm": "#00D452", "Yb": "#00BF38", "Lu": "#00AB24", "Hf": "#4DC2FF",
    "Ta": "#4DA6FF", "W": "#2194D6", "Re": "#267DAB", "Os": "#266696",
    "Ir": "#175487", "Pt": "#D0D0E0", "Au": "#FFD123", "Hg": "#B8B8D0",
    "Tl": "#A6544D", "Pb": "#575961", "Bi": "#9E4FB5", "Po": "#AB5C00",
    "At": "#754F45", "Rn": "#428296", "Fr": "#420066", "Ra": "#007D00",
    "Ac": "#70ABFA", "Th": "#00BAFF", "Pa": "#00A1FF", "U": "#008FFF",
    "Np": "#0080FF", "Pu": "#006BFF", "Am": "#545CF2", "Cm": "#785CE3",
}

ELEMENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40, "Cs": 2.44,
    "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99,
    "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92,
    "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75, "Ta": 1.70,
    "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40, "At": 1.50,
    "Rn": 1.50, "Fr": 2.60, "Ra": 2.21, "Ac": 2.15, "Th": 2.06, "Pa": 2.00,
    "U": 1.96, "Np": 1.90, "Pu": 1.87, "Am": 1.80, "Cm": 1.69,
}


def _get_structure_data(structure: Union[Dict, Any]) -> Tuple[List, Dict]:
    """Extract atoms and lattice from structure dict."""
    if hasattr(structure, 'as_dict'):
        structure = structure.as_dict()

    # Handle nested structure
    if "structure" in structure:
        struct = structure["structure"]
    else:
        struct = structure

    atoms = struct.get("atoms", [])
    lattice = struct.get("lattice", {})

    return atoms, lattice


def _get_lattice_matrix(lattice: Dict) -> np.ndarray:
    """Get 3x3 lattice matrix from lattice dict."""
    if "matrix" in lattice:
        return np.array(lattice["matrix"])

    a = lattice.get("a", 5.0)
    b = lattice.get("b", a)
    c = lattice.get("c", a)
    alpha = np.radians(lattice.get("alpha", 90))
    beta = np.radians(lattice.get("beta", 90))
    gamma = np.radians(lattice.get("gamma", 90))

    # Crystallographic to Cartesian conversion
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)

    val = (cos_alpha - cos_beta * cos_gamma) / sin_gamma

    matrix = np.array([
        [a, 0, 0],
        [b * cos_gamma, b * sin_gamma, 0],
        [c * cos_beta, c * val, c * np.sqrt(1 - cos_beta**2 - val**2)]
    ])

    return matrix


# =============================================================================
# XSF FORMAT (XCrySDen)
# =============================================================================

def export_xsf(
    structure: Dict[str, Any],
    include_forces: bool = False,
    forces: Optional[List[List[float]]] = None,
    animated: bool = False,
    frames: Optional[List[Dict]] = None
) -> Dict[str, Any]:
    """
    Export structure to XCrySDen XSF format.

    XSF format is widely used for visualization in XCrySDen, VESTA,
    and other crystallography software.

    Args:
        structure: Structure dictionary with atoms and lattice
        include_forces: Include atomic forces if available
        forces: List of force vectors [fx, fy, fz] per atom
        animated: Export as animated XSF (AXSF)
        frames: List of structure frames for animation

    Returns:
        Dictionary with success status and XSF content

    Reference: http://www.xcrysden.org/doc/XSF.html
    """
    atoms, lattice = _get_structure_data(structure)

    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}

    matrix = _get_lattice_matrix(lattice)
    lines = []

    # Header for animated XSF
    if animated and frames:
        lines.append("ANIMSTEPS {}".format(len(frames)))
        lines.append("")

    # Crystal structure
    lines.append("CRYSTAL")
    lines.append("")

    # Primitive lattice vectors (in Angstrom)
    lines.append("PRIMVEC")
    for i in range(3):
        lines.append("  {:15.10f} {:15.10f} {:15.10f}".format(
            matrix[i][0], matrix[i][1], matrix[i][2]
        ))
    lines.append("")

    # Conventional lattice vectors (same as primitive for now)
    lines.append("CONVVEC")
    for i in range(3):
        lines.append("  {:15.10f} {:15.10f} {:15.10f}".format(
            matrix[i][0], matrix[i][1], matrix[i][2]
        ))
    lines.append("")

    # Atomic coordinates
    lines.append("PRIMCOORD")
    lines.append("{} 1".format(len(atoms)))

    for i, atom in enumerate(atoms):
        elem = atom.get("element", "X")

        # Get Cartesian coordinates
        if "cartesian" in atom:
            x, y, z = atom["cartesian"]
        elif "coords" in atom:
            frac = np.array(atom["coords"])
            cart = frac @ matrix
            x, y, z = cart
        else:
            x, y, z = 0, 0, 0

        if include_forces and forces and i < len(forces):
            fx, fy, fz = forces[i]
            lines.append("{:3s} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}".format(
                elem, x, y, z, fx, fy, fz
            ))
        else:
            lines.append("{:3s} {:15.10f} {:15.10f} {:15.10f}".format(
                elem, x, y, z
            ))

    content = "\n".join(lines)

    return {
        "success": True,
        "format": "xsf",
        "content": content,
        "extension": ".xsf",
        "n_atoms": len(atoms),
        "description": "XCrySDen Structure Format"
    }


# =============================================================================
# CUBE FORMAT (Gaussian Volumetric)
# =============================================================================

def export_cube(
    structure: Dict[str, Any],
    volumetric_data: Optional[np.ndarray] = None,
    grid_shape: Tuple[int, int, int] = (50, 50, 50),
    data_type: str = "electron_density",
    origin: Optional[List[float]] = None
) -> Dict[str, Any]:
    """
    Export structure to Gaussian CUBE format.

    CUBE format is standard for volumetric data visualization
    (electron density, electrostatic potential, orbitals).

    Args:
        structure: Structure dictionary with atoms and lattice
        volumetric_data: 3D numpy array with volumetric data (optional)
        grid_shape: Grid dimensions (nx, ny, nz) if no data provided
        data_type: Type of volumetric data (for header comment)
        origin: Origin point [x, y, z] in Angstrom

    Returns:
        Dictionary with success status and CUBE content

    Reference: https://gaussian.com/cubegen/
    """
    atoms, lattice = _get_structure_data(structure)

    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}

    matrix = _get_lattice_matrix(lattice)

    # Convert to Bohr (CUBE uses atomic units)
    ANGSTROM_TO_BOHR = 1.8897259886
    matrix_bohr = matrix * ANGSTROM_TO_BOHR

    # Grid setup
    nx, ny, nz = grid_shape

    # Calculate grid vectors
    voxel_x = matrix_bohr[0] / nx
    voxel_y = matrix_bohr[1] / ny
    voxel_z = matrix_bohr[2] / nz

    # Origin
    if origin is None:
        origin_bohr = [0.0, 0.0, 0.0]
    else:
        origin_bohr = [o * ANGSTROM_TO_BOHR for o in origin]

    lines = []

    # Header comments
    lines.append("Crystal MCP Server - {} data".format(data_type))
    lines.append("Generated volumetric data in Gaussian CUBE format")

    # Number of atoms and origin
    lines.append("{:5d} {:12.6f} {:12.6f} {:12.6f}".format(
        len(atoms), origin_bohr[0], origin_bohr[1], origin_bohr[2]
    ))

    # Grid vectors
    lines.append("{:5d} {:12.6f} {:12.6f} {:12.6f}".format(
        nx, voxel_x[0], voxel_x[1], voxel_x[2]
    ))
    lines.append("{:5d} {:12.6f} {:12.6f} {:12.6f}".format(
        ny, voxel_y[0], voxel_y[1], voxel_y[2]
    ))
    lines.append("{:5d} {:12.6f} {:12.6f} {:12.6f}".format(
        nz, voxel_z[0], voxel_z[1], voxel_z[2]
    ))

    # Atomic numbers and positions
    ELEMENT_Z = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
        "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
        "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
        "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
        "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
        "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
        "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
        "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
        "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
        "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "U": 92,
    }

    for atom in atoms:
        elem = atom.get("element", "X")
        z = ELEMENT_Z.get(elem, 0)

        # Get Cartesian coordinates in Bohr
        if "cartesian" in atom:
            cart = np.array(atom["cartesian"]) * ANGSTROM_TO_BOHR
        elif "coords" in atom:
            frac = np.array(atom["coords"])
            cart = (frac @ matrix) * ANGSTROM_TO_BOHR
        else:
            cart = np.array([0, 0, 0])

        lines.append("{:5d} {:12.6f} {:12.6f} {:12.6f} {:12.6f}".format(
            z, float(z), cart[0], cart[1], cart[2]
        ))

    # Volumetric data
    if volumetric_data is not None:
        data = volumetric_data.flatten()
    else:
        # Generate placeholder data (zeros)
        data = np.zeros(nx * ny * nz)

    # Write data (6 values per line)
    for i in range(0, len(data), 6):
        chunk = data[i:i+6]
        line = " ".join("{:13.5E}".format(v) for v in chunk)
        lines.append(line)

    content = "\n".join(lines)

    return {
        "success": True,
        "format": "cube",
        "content": content,
        "extension": ".cube",
        "n_atoms": len(atoms),
        "grid_shape": grid_shape,
        "description": "Gaussian CUBE volumetric format"
    }


# =============================================================================
# VESTA FORMAT (Enhanced CIF with VESTA-specific settings)
# =============================================================================

def export_vesta(
    structure: Dict[str, Any],
    bond_search: bool = True,
    bond_radius: float = 0.2,
    polyhedral: bool = False,
    style: str = "ball_and_stick"
) -> Dict[str, Any]:
    """
    Export structure to VESTA-compatible format with visualization settings.

    VESTA can read CIF files, so we output enhanced CIF with comments
    for VESTA-specific visualization settings.

    Args:
        structure: Structure dictionary with atoms and lattice
        bond_search: Enable automatic bond search
        bond_radius: Bond cylinder radius
        polyhedral: Show coordination polyhedra
        style: Visualization style (ball_and_stick, space_filling, wireframe)

    Returns:
        Dictionary with success status and VESTA-ready CIF content

    Reference: K. Momma and F. Izumi, J. Appl. Crystallogr. 44, 1272 (2011)
    """
    atoms, lattice = _get_structure_data(structure)

    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}

    lines = []

    # CIF header with VESTA optimization notes
    lines.append("# VESTA-optimized CIF file")
    lines.append("# Generated by Crystal MCP Server")
    lines.append("# Recommended VESTA settings:")
    lines.append("#   Bond search: {}".format("enabled" if bond_search else "disabled"))
    lines.append("#   Style: {}".format(style))
    lines.append("")
    lines.append("data_crystal_structure")
    lines.append("")

    # Symmetry (P1 for general)
    lines.append("_symmetry_space_group_name_H-M   'P 1'")
    lines.append("_symmetry_Int_Tables_number      1")
    lines.append("")

    # Cell parameters
    a = lattice.get("a", 5.0)
    b = lattice.get("b", a)
    c = lattice.get("c", a)
    alpha = lattice.get("alpha", 90.0)
    beta = lattice.get("beta", 90.0)
    gamma = lattice.get("gamma", 90.0)

    lines.append("_cell_length_a     {:10.6f}".format(a))
    lines.append("_cell_length_b     {:10.6f}".format(b))
    lines.append("_cell_length_c     {:10.6f}".format(c))
    lines.append("_cell_angle_alpha  {:10.6f}".format(alpha))
    lines.append("_cell_angle_beta   {:10.6f}".format(beta))
    lines.append("_cell_angle_gamma  {:10.6f}".format(gamma))
    lines.append("")

    # Atom sites
    lines.append("loop_")
    lines.append("_atom_site_label")
    lines.append("_atom_site_type_symbol")
    lines.append("_atom_site_fract_x")
    lines.append("_atom_site_fract_y")
    lines.append("_atom_site_fract_z")
    lines.append("_atom_site_occupancy")

    # Track element counts for labels
    elem_counts = {}

    for atom in atoms:
        elem = atom.get("element", "X")

        # Generate unique label
        if elem not in elem_counts:
            elem_counts[elem] = 0
        elem_counts[elem] += 1
        label = "{}{}".format(elem, elem_counts[elem])

        # Get fractional coordinates
        if "coords" in atom:
            x, y, z = atom["coords"]
        else:
            x, y, z = 0, 0, 0

        occupancy = atom.get("occupancy", 1.0)

        lines.append("{:6s} {:3s} {:10.6f} {:10.6f} {:10.6f} {:6.4f}".format(
            label, elem, x, y, z, occupancy
        ))

    content = "\n".join(lines)

    # Also create VESTA settings file content
    vesta_settings = {
        "bond_search": bond_search,
        "bond_radius": bond_radius,
        "polyhedral": polyhedral,
        "style": style,
        "element_colors": {elem: ELEMENT_COLORS.get(elem, "#808080")
                         for elem in elem_counts.keys()},
        "element_radii": {elem: ELEMENT_RADII.get(elem, 1.0)
                         for elem in elem_counts.keys()}
    }

    return {
        "success": True,
        "format": "vesta",
        "content": content,
        "extension": ".cif",
        "n_atoms": len(atoms),
        "vesta_settings": vesta_settings,
        "description": "VESTA-compatible CIF with visualization hints"
    }


# =============================================================================
# JSON FORMAT (for Web 3D Viewers like Three.js)
# =============================================================================

def export_threejs_json(
    structure: Dict[str, Any],
    include_bonds: bool = True,
    bond_cutoff: float = 2.5,
    supercell: Tuple[int, int, int] = (1, 1, 1)
) -> Dict[str, Any]:
    """
    Export structure to JSON format optimized for Three.js visualization.

    This format is designed for web-based 3D crystal viewers using
    Three.js or similar WebGL libraries.

    Args:
        structure: Structure dictionary with atoms and lattice
        include_bonds: Calculate and include bond information
        bond_cutoff: Maximum bond length in Angstrom
        supercell: Supercell repetitions (nx, ny, nz)

    Returns:
        Dictionary with Three.js-ready JSON data
    """
    atoms, lattice = _get_structure_data(structure)

    if not atoms:
        return {"success": False, "error": {"code": "NO_ATOMS", "message": "No atoms in structure"}}

    matrix = _get_lattice_matrix(lattice)

    # Build atom list with Cartesian coordinates
    atom_data = []
    nx, ny, nz = supercell

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                offset = np.array([ix, iy, iz])

                for atom in atoms:
                    elem = atom.get("element", "X")

                    # Get fractional coordinates
                    if "coords" in atom:
                        frac = np.array(atom["coords"]) + offset
                    else:
                        frac = offset

                    # Convert to Cartesian
                    cart = frac @ matrix

                    atom_data.append({
                        "element": elem,
                        "position": cart.tolist(),
                        "color": ELEMENT_COLORS.get(elem, "#808080"),
                        "radius": ELEMENT_RADII.get(elem, 1.0)
                    })

    # Calculate bonds
    bonds = []
    if include_bonds:
        positions = np.array([a["position"] for a in atom_data])
        n_atoms = len(positions)

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist < bond_cutoff:
                    bonds.append({
                        "atom1": i,
                        "atom2": j,
                        "length": float(dist)
                    })

    # Unit cell edges for visualization
    cell_edges = []
    corners = [
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]
    ]
    edges = [
        (0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4),
        (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)
    ]

    # Scale to supercell
    full_matrix = matrix * np.array([[nx], [ny], [nz]])

    corner_positions = [np.array(c) @ full_matrix for c in corners]

    for i, j in edges:
        cell_edges.append({
            "start": corner_positions[i].tolist(),
            "end": corner_positions[j].tolist()
        })

    return {
        "success": True,
        "format": "threejs_json",
        "data": {
            "atoms": atom_data,
            "bonds": bonds,
            "unitCell": {
                "matrix": matrix.tolist(),
                "edges": cell_edges,
                "a": float(lattice.get("a", 5.0)),
                "b": float(lattice.get("b", 5.0)),
                "c": float(lattice.get("c", 5.0)),
                "alpha": float(lattice.get("alpha", 90)),
                "beta": float(lattice.get("beta", 90)),
                "gamma": float(lattice.get("gamma", 90))
            },
            "supercell": list(supercell),
            "n_atoms": len(atom_data),
            "n_bonds": len(bonds)
        },
        "description": "Three.js-optimized JSON for web visualization"
    }


# =============================================================================
# HTML VISUALIZATION (Self-contained Three.js viewer)
# =============================================================================

def export_html_viewer(
    structure: Dict[str, Any],
    title: str = "Crystal Structure Viewer",
    background_color: str = "#0a0a0a",
    auto_rotate: bool = True
) -> Dict[str, Any]:
    """
    Export structure as self-contained HTML file with Three.js viewer.

    Creates a standalone HTML file that can be opened in any browser
    to view the crystal structure in 3D.

    Args:
        structure: Structure dictionary with atoms and lattice
        title: Page title
        background_color: Background color (hex)
        auto_rotate: Enable auto-rotation

    Returns:
        Dictionary with success status and HTML content
    """
    # Get JSON data first
    json_result = export_threejs_json(structure)
    if not json_result["success"]:
        return json_result

    structure_data = json.dumps(json_result["data"], indent=2)

    html_template = '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: Arial, sans-serif; }}
        #info {{ position: absolute; top: 10px; left: 10px; color: white;
                 background: rgba(0,0,0,0.7); padding: 10px; border-radius: 5px; }}
        #controls {{ position: absolute; bottom: 10px; left: 10px; color: white;
                    background: rgba(0,0,0,0.7); padding: 10px; border-radius: 5px; font-size: 12px; }}
    </style>
</head>
<body>
    <div id="info">
        <h3 style="margin:0 0 10px 0">{title}</h3>
        <div id="atom-count"></div>
    </div>
    <div id="controls">
        Drag to rotate | Scroll to zoom | Right-click to pan
    </div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <script>
        const structureData = {structure_data};

        // Scene setup
        const scene = new THREE.Scene();
        scene.background = new THREE.Color("{background_color}");

        const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.1, 1000);
        camera.position.set(20, 20, 20);

        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.appendChild(renderer.domElement);

        // Controls
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.autoRotate = {auto_rotate};
        controls.autoRotateSpeed = 1.0;

        // Lighting
        scene.add(new THREE.AmbientLight(0xffffff, 0.5));
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight.position.set(10, 10, 10);
        scene.add(dirLight);

        // Center calculation
        let center = new THREE.Vector3();
        structureData.atoms.forEach(atom => {{
            center.add(new THREE.Vector3(...atom.position));
        }});
        center.divideScalar(structureData.atoms.length);

        // Add atoms
        structureData.atoms.forEach(atom => {{
            const geometry = new THREE.SphereGeometry(atom.radius * 0.3, 32, 32);
            const material = new THREE.MeshPhongMaterial({{
                color: atom.color,
                shininess: 80
            }});
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(
                atom.position[0] - center.x,
                atom.position[1] - center.y,
                atom.position[2] - center.z
            );
            scene.add(sphere);
        }});

        // Add bonds
        structureData.bonds.forEach(bond => {{
            const atom1 = structureData.atoms[bond.atom1];
            const atom2 = structureData.atoms[bond.atom2];
            const start = new THREE.Vector3(...atom1.position).sub(center);
            const end = new THREE.Vector3(...atom2.position).sub(center);

            const dir = new THREE.Vector3().subVectors(end, start);
            const len = dir.length();

            const geometry = new THREE.CylinderGeometry(0.08, 0.08, len, 8);
            const material = new THREE.MeshPhongMaterial({{ color: 0x666666 }});
            const cylinder = new THREE.Mesh(geometry, material);

            cylinder.position.copy(start).add(dir.multiplyScalar(0.5));
            cylinder.quaternion.setFromUnitVectors(
                new THREE.Vector3(0, 1, 0),
                dir.clone().normalize()
            );
            scene.add(cylinder);
        }});

        // Add unit cell
        const cellMaterial = new THREE.LineBasicMaterial({{ color: 0x00ffff }});
        structureData.unitCell.edges.forEach(edge => {{
            const points = [
                new THREE.Vector3(...edge.start).sub(center),
                new THREE.Vector3(...edge.end).sub(center)
            ];
            const geometry = new THREE.BufferGeometry().setFromPoints(points);
            scene.add(new THREE.Line(geometry, cellMaterial));
        }});

        // Update info
        document.getElementById('atom-count').textContent =
            `Atoms: ${{structureData.n_atoms}} | Bonds: ${{structureData.n_bonds}}`;

        // Animation
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }}
        animate();

        // Resize handler
        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});
    </script>
</body>
</html>'''

    html_content = html_template.format(
        title=title,
        structure_data=structure_data,
        background_color=background_color,
        auto_rotate="true" if auto_rotate else "false"
    )

    return {
        "success": True,
        "format": "html",
        "content": html_content,
        "extension": ".html",
        "n_atoms": json_result["data"]["n_atoms"],
        "description": "Self-contained HTML viewer with Three.js"
    }


# =============================================================================
# CONVENIENCE FUNCTION
# =============================================================================

def export_visualization(
    structure: Dict[str, Any],
    format: str = "html",
    **kwargs
) -> Dict[str, Any]:
    """
    Export structure to any visualization format.

    Args:
        structure: Structure dictionary
        format: Output format (xsf, cube, vesta, threejs, html)
        **kwargs: Format-specific options

    Returns:
        Dictionary with export result
    """
    exporters = {
        "xsf": export_xsf,
        "cube": export_cube,
        "vesta": export_vesta,
        "threejs": export_threejs_json,
        "html": export_html_viewer
    }

    if format not in exporters:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_FORMAT",
                "message": f"Unknown format: {format}. Supported: {list(exporters.keys())}"
            }
        }

    return exporters[format](structure, **kwargs)
