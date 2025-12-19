"""
photonic/metamaterials.py - Metamaterials and Metasurfaces

Comprehensive metamaterial generation per structure_catalogue.md Category 13:
(iii) Photonic & phononic crystals with deliberately broken symmetries
(iv) Metasurface unit cells (V-antenna, split-ring, nano-fin, catenary)
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice


# Metamaterial unit cell database
METAMATERIAL_DATABASE = {
    # Split-ring resonators
    "SRR": {
        "type": "magnetic", "shape": "ring_with_gap",
        "resonance": "magnetic_dipole", "frequency_range": "GHz-THz"
    },
    "double_SRR": {
        "type": "magnetic", "shape": "concentric_rings",
        "resonance": "coupled_magnetic", "frequency_range": "GHz"
    },
    "complementary_SRR": {
        "type": "babinet", "shape": "SRR_negative",
        "resonance": "electric_dipole", "frequency_range": "THz"
    },
    
    # Antenna elements
    "V_antenna": {
        "type": "phase_gradient", "shape": "V",
        "resonance": "localized_plasmon", "frequency_range": "IR-visible"
    },
    "dipole_antenna": {
        "type": "resonator", "shape": "rod",
        "resonance": "electric_dipole", "frequency_range": "IR"
    },
    "bow_tie": {
        "type": "field_enhancement", "shape": "double_triangle",
        "resonance": "gap_plasmon", "frequency_range": "visible"
    },
    "cross_antenna": {
        "type": "polarization_independent", "shape": "cross",
        "resonance": "dual_polarization", "frequency_range": "THz-IR"
    },
    
    # Geometric phase elements
    "nano_fin": {
        "type": "Pancharatnam-Berry", "shape": "rectangular_pillar",
        "resonance": "waveguide", "frequency_range": "visible"
    },
    "nano_brick": {
        "type": "birefringent", "shape": "rectangular_block",
        "resonance": "shape_resonance", "frequency_range": "visible-IR"
    },
    "catenary": {
        "type": "catenary_curve", "shape": "hanging_chain",
        "resonance": "continuous_phase", "frequency_range": "THz"
    },
    
    # Plasmonic structures
    "nano_hole_array": {
        "type": "EOT", "shape": "periodic_holes",
        "resonance": "extraordinary_transmission", "frequency_range": "visible"
    },
    "nano_slit_array": {
        "type": "SPP", "shape": "periodic_slits",
        "resonance": "surface_plasmon", "frequency_range": "visible-IR"
    },
    "fishnet": {
        "type": "negative_index", "shape": "double_mesh",
        "resonance": "magnetic_electric", "frequency_range": "IR"
    },
    
    # Huygens metasurfaces
    "huygens_element": {
        "type": "Huygens", "shape": "disk_pair",
        "resonance": "Kerker_condition", "frequency_range": "visible"
    },
}


# Photonic crystal types
PHOTONIC_CRYSTAL_DATABASE = {
    "woodpile": {
        "structure": "3D", "symmetry": "FCC", "bandgap": True,
        "fabrication": "layer_by_layer", "materials": ["Si", "GaAs"]
    },
    "inverse_opal": {
        "structure": "3D", "symmetry": "FCC", "bandgap": True,
        "fabrication": "self_assembly", "materials": ["Si", "TiO2"]
    },
    "yablonovite": {
        "structure": "3D", "symmetry": "FCC", "bandgap": True,
        "fabrication": "drilling", "materials": ["Al2O3", "GaAs"]
    },
    "triangular_hole": {
        "structure": "2D", "symmetry": "hexagonal", "bandgap": "TE only",
        "fabrication": "lithography", "materials": ["Si", "GaAs", "InP"]
    },
    "square_rod": {
        "structure": "2D", "symmetry": "square", "bandgap": "TM only",
        "fabrication": "lithography", "materials": ["Si", "polymer"]
    },
    "hexagonal_rod": {
        "structure": "2D", "symmetry": "hexagonal", "bandgap": "TM only",
        "fabrication": "self_assembly", "materials": ["SiO2", "PS"]
    },
    "kagome": {
        "structure": "2D", "symmetry": "kagome", "bandgap": "flat_band",
        "fabrication": "lithography", "materials": ["Si", "graphene"]
    },
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_metasurface_unit_cell(
    element_type: str = "V_antenna",
    material: str = "Au",
    period_nm: float = 300.0,
    element_size_nm: float = 100.0,
    rotation_deg: float = 0.0,
    substrate: str = "SiO2"
) -> Dict[str, Any]:
    """
    Generate metasurface unit cell.
    
    Args:
        element_type: Metamaterial element type
        material: Metal material (Au, Ag, Al)
        period_nm: Unit cell period
        element_size_nm: Element characteristic size
        rotation_deg: Element rotation for phase control
        substrate: Substrate material
    
    Returns:
        Metasurface unit cell structure
    """
    if element_type not in METAMATERIAL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_ELEMENT", "message": f"Unknown element '{element_type}'",
                      "available": list(METAMATERIAL_DATABASE.keys())}
        }
    
    info = METAMATERIAL_DATABASE[element_type]
    
    # Convert to Angstrom
    period = period_nm * 10
    size = element_size_nm * 10
    
    lattice = Lattice.orthorhombic(period, period, 500)  # 50nm height
    
    species = []
    coords = []
    
    # Add substrate layer
    n_sub = 5
    for i in range(n_sub):
        for j in range(n_sub):
            x = (i + 0.5) / n_sub
            y = (j + 0.5) / n_sub
            z = 0.1  # Substrate region
            
            if "SiO2" in substrate:
                species.extend(["Si", "O", "O"])
                coords.extend([[x, y, z], [x + 0.02, y, z], [x, y + 0.02, z]])
            else:
                species.append(substrate[0] if substrate else "Si")
                coords.append([x, y, z])
    
    # Add metamaterial element
    rotation_rad = np.radians(rotation_deg)
    
    if element_type == "V_antenna":
        # V-shaped antenna
        arm_length = size / period
        arm_angle = np.pi / 6  # 30 degrees
        
        n_atoms = 20
        for arm in [-1, 1]:
            for i in range(n_atoms // 2):
                t = i / (n_atoms // 2 - 1)
                x = 0.5 + arm * t * arm_length * np.cos(arm_angle + rotation_rad)
                y = 0.5 + t * arm_length * np.sin(arm_angle * arm)
                z = 0.5
                species.append(material)
                coords.append([x, y, z])
    
    elif element_type == "SRR":
        # Split-ring resonator
        n_atoms = 30
        gap_angle = np.pi / 6
        ring_radius = size / (2 * period)
        
        for i in range(n_atoms):
            angle = 2 * np.pi * i / n_atoms
            # Skip gap region
            if abs(angle - np.pi) < gap_angle / 2:
                continue
            
            x = 0.5 + ring_radius * np.cos(angle + rotation_rad)
            y = 0.5 + ring_radius * np.sin(angle + rotation_rad)
            z = 0.5
            species.append(material)
            coords.append([x, y, z])
    
    elif element_type == "nano_fin":
        # Rectangular nano-fin (pillar)
        fin_width = size / (3 * period)
        fin_length = size / period
        
        n_w = 3
        n_l = 8
        n_h = 5
        
        for i in range(n_w):
            for j in range(n_l):
                for k in range(n_h):
                    # Rotate coordinates
                    x0 = (i - n_w/2) * fin_width / n_w
                    y0 = (j - n_l/2) * fin_length / n_l
                    
                    x = 0.5 + x0 * np.cos(rotation_rad) - y0 * np.sin(rotation_rad)
                    y = 0.5 + x0 * np.sin(rotation_rad) + y0 * np.cos(rotation_rad)
                    z = 0.4 + k * 0.1 / n_h
                    
                    species.append(material)
                    coords.append([x, y, z])
    
    elif element_type == "bow_tie":
        # Bow-tie antenna
        n_atoms = 20
        triangle_size = size / period
        gap = 0.01
        
        for triangle in [-1, 1]:
            for i in range(n_atoms // 2):
                t = i / (n_atoms // 2 - 1)
                x = 0.5 + triangle * (gap / 2 + t * triangle_size / 2)
                # Triangle shape
                width = (1 - t) * triangle_size / 4
                for w in [-1, 1]:
                    y = 0.5 + w * width / 2
                    z = 0.5
                    species.append(material)
                    coords.append([x, y, z])
    
    else:
        # Generic circular element
        n_atoms = 20
        radius = size / (2 * period)
        for i in range(n_atoms):
            angle = 2 * np.pi * i / n_atoms
            x = 0.5 + radius * np.cos(angle)
            y = 0.5 + radius * np.sin(angle)
            z = 0.5
            species.append(material)
            coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "element_type": element_type,
        "resonance_type": info["resonance"],
        "frequency_range": info["frequency_range"],
        "material": material,
        "period_nm": period_nm,
        "element_size_nm": element_size_nm,
        "rotation_deg": rotation_deg,
        "substrate": substrate,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_photonic_crystal(
    crystal_type: str = "woodpile",
    material: str = "Si",
    lattice_constant_nm: float = 500.0,
    fill_fraction: float = 0.3,
    n_periods: List[int] = [2, 2, 2]
) -> Dict[str, Any]:
    """
    Generate photonic crystal structure.
    
    Args:
        crystal_type: Photonic crystal type
        material: Dielectric material
        lattice_constant_nm: Lattice constant in nm
        fill_fraction: Volume fraction of high-index material
        n_periods: Number of periods in each direction
    
    Returns:
        Photonic crystal structure
    """
    if crystal_type not in PHOTONIC_CRYSTAL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown crystal type",
                      "available": list(PHOTONIC_CRYSTAL_DATABASE.keys())}
        }
    
    info = PHOTONIC_CRYSTAL_DATABASE[crystal_type]
    
    a = lattice_constant_nm * 10  # Convert to Angstrom
    nx, ny, nz = n_periods
    
    if info["structure"] == "3D":
        lattice = Lattice.cubic(a * nx)
    else:
        lattice = Lattice.orthorhombic(a * nx, a * ny, 500)  # 2D with vacuum
    
    species = []
    coords = []
    
    np.random.seed(42)
    
    if crystal_type == "woodpile":
        # Log-pile structure
        rod_width = fill_fraction * 0.5
        n_rods = 4
        
        for layer in range(nz * 4):
            z = (layer + 0.5) / (nz * 4)
            direction = layer % 2  # Alternating x and y
            
            for rod in range(n_rods * (nx if direction == 0 else ny)):
                if direction == 0:
                    y = (rod + 0.5) / (n_rods * ny)
                    for xi in range(10):
                        x = xi / 10
                        species.append(material)
                        coords.append([x, y, z])
                else:
                    x = (rod + 0.5) / (n_rods * nx)
                    for yi in range(10):
                        y = yi / 10
                        species.append(material)
                        coords.append([x, y, z])
    
    elif crystal_type in ["triangular_hole", "hexagonal_rod"]:
        # 2D hexagonal photonic crystal
        for i in range(nx * 3):
            for j in range(ny * 3):
                x = (i + (j % 2) * 0.5) / (nx * 3)
                y = j * np.sqrt(3) / 2 / (ny * 3)
                z = 0.5
                
                species.append(material)
                coords.append([x % 1, y % 1, z])
    
    elif crystal_type == "inverse_opal":
        # FCC arrangement of spherical voids
        fcc_sites = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for site in fcc_sites:
                        x = (i + site[0]) / nx
                        y = (j + site[1]) / ny
                        z = (k + site[2]) / nz
                        
                        # Add shell atoms around each site
                        for _ in range(8):
                            dx = np.random.uniform(-0.05, 0.05)
                            dy = np.random.uniform(-0.05, 0.05)
                            dz = np.random.uniform(-0.05, 0.05)
                            r = np.sqrt(dx**2 + dy**2 + dz**2)
                            if r > 0.02 and r < 0.08:
                                species.append(material)
                                coords.append([(x + dx) % 1, (y + dy) % 1, (z + dz) % 1])
    
    else:
        # Generic square lattice
        for i in range(nx * 3):
            for j in range(ny * 3):
                x = i / (nx * 3)
                y = j / (ny * 3)
                z = 0.5
                species.append(material)
                coords.append([x, y, z])
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "crystal_type": crystal_type,
        "dimension": info["structure"],
        "symmetry": info["symmetry"],
        "has_bandgap": info["bandgap"],
        "fabrication_method": info["fabrication"],
        "material": material,
        "lattice_constant_nm": lattice_constant_nm,
        "fill_fraction": fill_fraction,
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_broken_symmetry_photonic(
    base_type: str = "triangular_hole",
    symmetry_breaking: str = "missing_hole",
    defect_position: str = "center",
    material: str = "Si"
) -> Dict[str, Any]:
    """
    Generate photonic crystal with broken symmetry (defect modes).
    
    Args:
        base_type: Base photonic crystal type
        symmetry_breaking: Type of symmetry breaking
        defect_position: Position of defect
        material: Dielectric material
    
    Returns:
        Photonic crystal with defect
    """
    symmetry_breaking_types = {
        "missing_hole": "point_defect_acceptor",
        "missing_rod": "point_defect_donor",
        "line_defect": "waveguide",
        "L3_cavity": "high_Q_cavity",
        "H1_cavity": "single_defect",
        "heterostructure": "interface_states",
    }
    
    if symmetry_breaking not in symmetry_breaking_types:
        return {"success": False, "error": {"code": "INVALID_BREAKING", "message": "Unknown symmetry breaking type"}}
    
    # Generate base crystal
    base = generate_photonic_crystal(base_type, material)
    
    if not base["success"]:
        return base
    
    atoms = list(base["structure"]["atoms"])
    
    # Apply symmetry breaking
    if defect_position == "center":
        center = [0.5, 0.5, 0.5]
    else:
        center = [0.25, 0.25, 0.5]
    
    if symmetry_breaking in ["missing_hole", "missing_rod", "H1_cavity"]:
        # Remove atoms near center
        filtered_atoms = []
        for atom in atoms:
            dist = np.sqrt(sum((atom["coords"][i] - center[i])**2 for i in range(3)))
            if dist > 0.1:
                filtered_atoms.append(atom)
        atoms = filtered_atoms
    
    elif symmetry_breaking == "L3_cavity":
        # Remove 3 atoms in a line
        filtered_atoms = []
        for atom in atoms:
            in_defect = (abs(atom["coords"][0] - 0.5) < 0.15 and 
                        abs(atom["coords"][1] - 0.5) < 0.05)
            if not in_defect:
                filtered_atoms.append(atom)
        atoms = filtered_atoms
    
    elif symmetry_breaking == "line_defect":
        # Remove a row to create waveguide
        filtered_atoms = []
        for atom in atoms:
            if abs(atom["coords"][1] - 0.5) > 0.05:
                filtered_atoms.append(atom)
        atoms = filtered_atoms
    
    return {
        "success": True,
        "base_type": base_type,
        "symmetry_breaking": symmetry_breaking,
        "defect_type": symmetry_breaking_types[symmetry_breaking],
        "defect_position": defect_position,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms, "lattice": base["structure"]["lattice"]}
    }
