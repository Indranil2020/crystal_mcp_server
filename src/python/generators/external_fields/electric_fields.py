"""
external_fields/electric_fields.py - Electric Field Effects on Crystal Structures

Provides structure modifications for:
- Electric field polarized cells
- Ferroelectric domain structures
- Polar surface configurations
- Field-induced phase transitions
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np


def generate_polarized_cell(
    structure,
    field_direction: List[float] = [0, 0, 1],
    field_magnitude_V_per_A: float = 0.1,
    polarization_mode: str = "ionic",
    max_displacement_A: float = 0.5
) -> Dict[str, Any]:
    """
    Generate structure with electric field-induced polarization.

    Simulates the effect of an external electric field by displacing
    ions according to their Born effective charges.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        field_direction: Electric field direction [x, y, z] (will be normalized)
        field_magnitude_V_per_A: Field strength in V/Angstrom
        polarization_mode: "ionic" (displace ions), "electronic" (add dipole info)
        max_displacement_A: Maximum allowed ionic displacement

    Returns:
        Polarized structure with field information
    """
    from pymatgen.core import Structure, Lattice

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure.copy()

    # Normalize field direction
    field_dir = np.array(field_direction)
    field_dir = field_dir / np.linalg.norm(field_dir)

    # Approximate Born effective charges (Z*) for common elements
    # In reality, these should be calculated from DFT
    born_charges = {
        # Cations (positive)
        "Li": 1.0, "Na": 1.0, "K": 1.0,
        "Mg": 2.0, "Ca": 2.0, "Sr": 2.0, "Ba": 2.0,
        "Ti": 4.0, "Zr": 4.0, "Hf": 4.0,
        "Pb": 2.5, "Bi": 3.0,
        "Fe": 2.5, "Co": 2.5, "Ni": 2.5,
        "Zn": 2.0, "Cd": 2.0,
        # Anions (negative)
        "O": -2.0, "S": -2.0, "Se": -2.0, "Te": -2.0,
        "F": -1.0, "Cl": -1.0, "Br": -1.0, "I": -1.0,
        "N": -3.0, "P": -3.0, "As": -3.0,
        # Special cases for ferroelectrics
        "C": 0.0, "Si": 0.0, "Ge": 0.0
    }

    # Calculate displacements
    new_coords = []
    displacements = []
    total_dipole = np.zeros(3)

    for site in struct.sites:
        element = str(site.specie)
        z_star = born_charges.get(element, 0.0)

        # Displacement = Z* * E (simplified model)
        # Real calculation: u = Z* * E / (m * omega^2) but we use max_displacement constraint
        displacement_magnitude = abs(z_star) * field_magnitude_V_per_A * 0.1  # Scale factor
        displacement_magnitude = min(displacement_magnitude, max_displacement_A)

        # Direction based on charge sign
        displacement = displacement_magnitude * field_dir * np.sign(z_star)

        new_cart = site.coords + displacement
        new_coords.append(new_cart)
        displacements.append(displacement.tolist())

        # Accumulate dipole moment
        total_dipole += z_star * displacement

    # Create polarized structure
    polarized_struct = Structure(
        lattice=struct.lattice,
        species=[str(s.specie) for s in struct.sites],
        coords=new_coords,
        coords_are_cartesian=True
    )

    # Calculate polarization (dipole per unit volume)
    polarization = total_dipole / struct.lattice.volume  # e*A/A^3 = e/A^2
    polarization_uC_cm2 = polarization * 1.602e-19 * 1e16  # Convert to μC/cm²

    return {
        "success": True,
        "mode": polarization_mode,
        "field_direction": field_dir.tolist(),
        "field_magnitude_V_per_A": field_magnitude_V_per_A,
        "n_atoms": len(struct),
        "displacements_A": displacements,
        "max_displacement_A": max(np.linalg.norm(d) for d in displacements),
        "total_dipole_moment_eA": total_dipole.tolist(),
        "polarization_e_per_A2": polarization.tolist(),
        "polarization_uC_per_cm2": polarization_uC_cm2.tolist(),
        "structure": {
            "lattice": {"matrix": polarized_struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in polarized_struct.sites
            ]
        },
        "pymatgen_structure": polarized_struct.as_dict()
    }


def generate_ferroelectric_domain(
    material: str = "BaTiO3",
    domain_type: str = "180_degree",
    n_domains: int = 2,
    domain_width_unit_cells: int = 5,
    wall_type: str = "sharp"
) -> Dict[str, Any]:
    """
    Generate ferroelectric domain structure.

    Creates multi-domain ferroelectric structures with specified
    domain wall configurations.

    Args:
        material: Ferroelectric material (BaTiO3, PbTiO3, PZT, BiFeO3, etc.)
        domain_type: "180_degree", "90_degree", "71_degree", "109_degree"
        n_domains: Number of domains to create
        domain_width_unit_cells: Width of each domain in unit cells
        wall_type: "sharp" (abrupt) or "gradual" (Bloch/Néel-like)

    Returns:
        Multi-domain ferroelectric structure
    """
    from pymatgen.core import Structure, Lattice

    # Material parameters for common ferroelectrics
    ferroelectrics = {
        "BaTiO3": {
            "a": 3.99, "c": 4.04,  # Tetragonal phase
            "species": ["Ba", "Ti", "O", "O", "O"],
            "frac_coords": [
                [0.0, 0.0, 0.0],      # Ba
                [0.5, 0.5, 0.52],     # Ti (displaced along c)
                [0.5, 0.5, 0.02],     # O (apical, displaced)
                [0.5, 0.0, 0.49],     # O (equatorial)
                [0.0, 0.5, 0.49],     # O (equatorial)
            ],
            "Ti_displacement": 0.02,  # Fractional displacement
            "polarization_uC_cm2": 26
        },
        "PbTiO3": {
            "a": 3.90, "c": 4.15,
            "species": ["Pb", "Ti", "O", "O", "O"],
            "frac_coords": [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.54],
                [0.5, 0.5, 0.11],
                [0.5, 0.0, 0.62],
                [0.0, 0.5, 0.62],
            ],
            "Ti_displacement": 0.04,
            "polarization_uC_cm2": 75
        },
        "BiFeO3": {
            "a": 5.58, "c": 13.87,  # Hexagonal setting
            "species": ["Bi", "Bi", "Fe", "Fe", "O", "O", "O", "O", "O", "O"],
            "frac_coords": [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.5],
                [0.0, 0.0, 0.221],
                [0.0, 0.0, 0.721],
                [0.443, 0.0, 0.083],
                [0.0, 0.443, 0.083],
                [-0.443, -0.443, 0.083],
                [0.443, 0.0, 0.583],
                [0.0, 0.443, 0.583],
                [-0.443, -0.443, 0.583],
            ],
            "Ti_displacement": 0.02,  # Fe displacement
            "polarization_uC_cm2": 100
        },
        "PZT": {  # Pb(Zr0.5Ti0.5)O3
            "a": 4.05, "c": 4.10,
            "species": ["Pb", "Ti", "O", "O", "O"],  # Using Ti site
            "frac_coords": [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.53],
                [0.5, 0.5, 0.05],
                [0.5, 0.0, 0.55],
                [0.0, 0.5, 0.55],
            ],
            "Ti_displacement": 0.03,
            "polarization_uC_cm2": 50
        }
    }

    if material not in ferroelectrics:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_MATERIAL",
                "message": f"Unknown ferroelectric material '{material}'",
                "available": list(ferroelectrics.keys())
            }
        }

    params = ferroelectrics[material]

    # Create base unit cell
    lattice = Lattice.tetragonal(params["a"], params["c"])
    base_struct = Structure(
        lattice=lattice,
        species=params["species"],
        coords=params["frac_coords"],
        coords_are_cartesian=False
    )

    # Create supercell for domains
    supercell_size = [domain_width_unit_cells * n_domains, 1, 1]
    supercell = base_struct * supercell_size

    # Apply domain structure
    domain_info = []
    new_coords = []

    displacement = params["Ti_displacement"] * params["c"]  # Actual displacement in Angstroms

    for i, site in enumerate(supercell.sites):
        element = str(site.specie)
        coords = site.coords.copy()

        # Determine which domain this atom belongs to
        domain_idx = int(site.frac_coords[0] * n_domains / 1.0) % n_domains

        # Determine polarization direction for this domain
        if domain_type == "180_degree":
            polarization_sign = 1 if domain_idx % 2 == 0 else -1

            # Displace ferroactive cation (Ti, Fe, etc.)
            if element in ["Ti", "Fe"]:
                if wall_type == "sharp":
                    coords[2] += polarization_sign * displacement
                else:
                    # Gradual wall - smooth transition
                    wall_position = (domain_idx + 0.5) / n_domains
                    distance_from_wall = abs(site.frac_coords[0] - wall_position)
                    smoothing = np.tanh(distance_from_wall * 10)
                    coords[2] += polarization_sign * displacement * smoothing

        elif domain_type == "90_degree":
            # 90° domains alternate between [001] and [100] polarization
            if domain_idx % 2 == 0:
                if element in ["Ti", "Fe"]:
                    coords[2] += displacement
            else:
                if element in ["Ti", "Fe"]:
                    coords[0] += displacement

        new_coords.append(coords)
        domain_info.append({
            "site_index": i,
            "element": element,
            "domain": domain_idx
        })

    # Create domain structure
    domain_struct = Structure(
        lattice=supercell.lattice,
        species=[str(s.specie) for s in supercell.sites],
        coords=new_coords,
        coords_are_cartesian=True
    )

    return {
        "success": True,
        "material": material,
        "domain_type": domain_type,
        "n_domains": n_domains,
        "domain_width_unit_cells": domain_width_unit_cells,
        "wall_type": wall_type,
        "n_atoms": len(domain_struct),
        "supercell_size": supercell_size,
        "bulk_polarization_uC_cm2": params["polarization_uC_cm2"],
        "domains": [
            {"domain_id": i, "polarization_direction": [0, 0, 1 if i % 2 == 0 else -1]}
            for i in range(n_domains)
        ],
        "structure": {
            "lattice": {"matrix": domain_struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in domain_struct.sites
            ]
        },
        "pymatgen_structure": domain_struct.as_dict()
    }


def generate_polar_surface(
    structure,
    miller_index: List[int] = [0, 0, 1],
    n_layers: int = 6,
    vacuum: float = 15.0,
    compensation: str = "stoichiometric",
    surface_charge: float = 0.0
) -> Dict[str, Any]:
    """
    Generate polar surface with charge compensation.

    Creates polar surfaces with various compensation mechanisms:
    - Stoichiometric reconstruction
    - Hydroxylation
    - Surface charge
    - Vacancy formation

    Args:
        structure: Bulk crystal structure
        miller_index: Surface orientation [h, k, l]
        n_layers: Number of atomic layers
        vacuum: Vacuum thickness in Angstroms
        compensation: "stoichiometric", "hydroxylated", "charged", "vacancy"
        surface_charge: Surface charge density in e/A² (for charged surfaces)

    Returns:
        Polar surface structure
    """
    from pymatgen.core import Structure, Lattice
    from pymatgen.core.surface import SlabGenerator

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Generate slab
    slabgen = SlabGenerator(
        struct,
        miller_index=tuple(miller_index),
        min_slab_size=n_layers * 2,  # Approximate
        min_vacuum_size=vacuum,
        center_slab=True
    )

    slabs = slabgen.get_slabs()

    if not slabs:
        return {
            "success": False,
            "error": {
                "code": "SLAB_GENERATION_FAILED",
                "message": f"Could not generate slab for {miller_index} surface"
            }
        }

    slab = slabs[0]  # Take first slab

    # Analyze polarity
    # Check if top and bottom terminations have different compositions
    top_species = []
    bottom_species = []

    for site in slab.sites:
        z = site.frac_coords[2]
        if z > 0.7:  # Top region
            top_species.append(str(site.specie))
        elif z < 0.3:  # Bottom region
            bottom_species.append(str(site.specie))

    is_polar = set(top_species) != set(bottom_species)

    # Calculate surface dipole (simplified)
    dipole_moment = 0.0
    if is_polar:
        # Estimate based on ionic charges
        for site in slab.sites:
            z = site.coords[2]
            element = str(site.specie)
            # Use oxidation state approximation
            charge = {"O": -2, "N": -3, "Ti": 4, "Ba": 2, "Pb": 2}.get(element, 0)
            dipole_moment += charge * z

    result = {
        "success": True,
        "miller_index": list(miller_index),
        "n_layers": n_layers,
        "vacuum_A": vacuum,
        "compensation": compensation,
        "is_polar": is_polar,
        "dipole_moment_eA": dipole_moment,
        "top_termination": list(set(top_species)),
        "bottom_termination": list(set(bottom_species)),
        "n_atoms": len(slab),
        "slab_thickness_A": slab.lattice.c - vacuum,
        "structure": {
            "lattice": {"matrix": slab.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in slab.sites
            ]
        },
        "pymatgen_structure": slab.as_dict()
    }

    # Apply compensation if polar
    if is_polar and compensation == "hydroxylated":
        result["compensation_note"] = "Hydroxyl groups should be added to polar surface"
        result["suggested_adsorbate"] = "OH"
    elif is_polar and compensation == "charged":
        result["surface_charge_e_per_A2"] = surface_charge

    return result


def generate_field_induced_phase(
    material: str = "SrTiO3",
    electric_field_kV_cm: float = 100.0,
    temperature_K: float = 300.0
) -> Dict[str, Any]:
    """
    Generate field-induced phase structure.

    Simulates electric field-induced phase transitions like:
    - Paraelectric to ferroelectric
    - AFE to FE
    - Quantum paraelectric softening

    Args:
        material: Material undergoing transition
        electric_field_kV_cm: Electric field strength
        temperature_K: Temperature in Kelvin

    Returns:
        Field-induced phase structure
    """
    from pymatgen.core import Structure, Lattice

    # Material-specific phase transition parameters
    materials = {
        "SrTiO3": {
            # Quantum paraelectric - field can induce ferroelectricity
            "cubic_a": 3.905,
            "tc": 0,  # No bulk Tc (quantum paraelectric)
            "field_threshold_kV_cm": 50,
            "soft_mode_freq_cm-1": 15,
            "species": ["Sr", "Ti", "O", "O", "O"],
            "coords_para": [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
            ],
            "displacement_scaling": 0.0002  # Displacement per kV/cm
        },
        "PbZrO3": {
            # Antiferroelectric - field can induce FE
            "ortho_a": 5.87, "ortho_b": 11.74, "ortho_c": 8.19,
            "tc": 503,
            "field_threshold_kV_cm": 200,
            "species": ["Pb", "Pb", "Pb", "Pb", "Zr", "Zr", "Zr", "Zr",
                       "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O", "O"],
            "phase": "antiferroelectric"
        },
        "BaTiO3": {
            "tetra_a": 3.99, "tetra_c": 4.04,
            "tc": 393,
            "species": ["Ba", "Ti", "O", "O", "O"],
            "coords_ferro": [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.52],
                [0.5, 0.5, 0.02],
                [0.5, 0.0, 0.49],
                [0.0, 0.5, 0.49],
            ],
            "polarization_uC_cm2": 26
        }
    }

    if material not in materials:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_MATERIAL",
                "message": f"Unknown material '{material}'",
                "available": list(materials.keys())
            }
        }

    params = materials[material]

    if material == "SrTiO3":
        # Field-induced ferroelectricity in STO
        a = params["cubic_a"]

        # Calculate induced displacement
        if electric_field_kV_cm > params["field_threshold_kV_cm"]:
            displacement = params["displacement_scaling"] * electric_field_kV_cm
            phase = "field_induced_ferroelectric"
        else:
            displacement = 0.0
            phase = "paraelectric"

        # Create structure with displacement
        coords = [list(c) for c in params["coords_para"]]
        coords[1][2] += displacement  # Displace Ti
        coords[2][2] -= displacement * 0.5  # Oxygen counter-displacement

        lattice = Lattice.cubic(a)
        struct = Structure(
            lattice=lattice,
            species=params["species"],
            coords=coords,
            coords_are_cartesian=False
        )

        return {
            "success": True,
            "material": material,
            "electric_field_kV_cm": electric_field_kV_cm,
            "temperature_K": temperature_K,
            "initial_phase": "quantum_paraelectric",
            "final_phase": phase,
            "Ti_displacement_A": displacement * a,
            "induced_polarization": displacement > 0,
            "structure": {
                "lattice": {"matrix": struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in struct.sites
                ]
            },
            "pymatgen_structure": struct.as_dict()
        }

    elif material == "BaTiO3":
        # Enhanced polarization under field
        a = params["tetra_a"]
        c = params["tetra_c"]

        # Field enhancement of polarization
        enhancement_factor = 1 + electric_field_kV_cm / 1000.0

        coords = [list(c) for c in params["coords_ferro"]]
        # Enhance Ti displacement
        coords[1][2] = 0.5 + (coords[1][2] - 0.5) * enhancement_factor

        lattice = Lattice.tetragonal(a, c)
        struct = Structure(
            lattice=lattice,
            species=params["species"],
            coords=coords,
            coords_are_cartesian=False
        )

        return {
            "success": True,
            "material": material,
            "electric_field_kV_cm": electric_field_kV_cm,
            "temperature_K": temperature_K,
            "initial_phase": "ferroelectric_tetragonal",
            "final_phase": "enhanced_ferroelectric",
            "polarization_enhancement": enhancement_factor,
            "structure": {
                "lattice": {"matrix": struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in struct.sites
                ]
            },
            "pymatgen_structure": struct.as_dict()
        }

    return {"success": False, "error": {"message": "Material handling not implemented"}}


def apply_electric_field_metadata(
    structure,
    field_vector_V_per_A: List[float] = [0, 0, 0.1],
    field_type: str = "sawtooth",
    reference_potential: float = 0.0
) -> Dict[str, Any]:
    """
    Add electric field metadata for DFT calculations.

    Prepares structure with electric field information for codes that
    support field calculations (VASP, QE with sawtooth potential).

    Args:
        structure: Crystal structure
        field_vector_V_per_A: Electric field vector [Ex, Ey, Ez] in V/Angstrom
        field_type: "sawtooth" (dipole correction), "finite" (finite field)
        reference_potential: Reference potential for sawtooth

    Returns:
        Structure with field metadata and calculation settings
    """
    from pymatgen.core import Structure, Lattice

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    field_vec = np.array(field_vector_V_per_A)
    field_magnitude = np.linalg.norm(field_vec)

    # Convert to other units
    field_Ry_per_Bohr = field_magnitude * 0.019446  # V/A to Ry/Bohr
    field_Ha_per_Bohr = field_Ry_per_Bohr / 2  # Ry to Ha

    # VASP settings for electric field
    vasp_settings = {
        "LDIPOL": True,
        "IDIPOL": 3 if abs(field_vec[2]) > abs(field_vec[0]) and abs(field_vec[2]) > abs(field_vec[1]) else (
            2 if abs(field_vec[1]) > abs(field_vec[0]) else 1
        ),
        "DIPOL": [0.5, 0.5, 0.5],  # Dipole center
        "EFIELD": field_magnitude,  # V/Angstrom
    }

    # QE settings for sawtooth potential
    qe_settings = {
        "tefield": True,
        "dipfield": True,
        "edir": 3 if abs(field_vec[2]) > abs(field_vec[0]) and abs(field_vec[2]) > abs(field_vec[1]) else (
            2 if abs(field_vec[1]) > abs(field_vec[0]) else 1
        ),
        "eamp": field_Ry_per_Bohr,
        "emaxpos": 0.5,
        "eopreg": 0.1
    }

    return {
        "success": True,
        "field_type": field_type,
        "field_vector_V_per_A": field_vec.tolist(),
        "field_magnitude_V_per_A": float(field_magnitude),
        "field_magnitude_kV_cm": float(field_magnitude * 1e8),  # V/A = 10^8 V/cm
        "field_magnitude_Ry_per_Bohr": float(field_Ry_per_Bohr),
        "field_magnitude_Ha_per_Bohr": float(field_Ha_per_Bohr),
        "vasp_settings": vasp_settings,
        "qe_settings": qe_settings,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "notes": [
            "For VASP: Set LDIPOL=.TRUE., IDIPOL, EFIELD in INCAR",
            "For QE: Set tefield=.true., dipfield=.true., edir, eamp in &CONTROL",
            "Ensure sufficient vacuum for slab calculations"
        ]
    }
