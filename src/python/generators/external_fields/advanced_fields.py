"""
external_fields/advanced_fields.py - Advanced External Field Effects

Provides structure generation for:
- Floquet-engineered systems
- Optical phonon pumped structures
- Warm dense matter conditions
- Dynamic strain fields
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


def generate_floquet_cell(
    structure,
    drive_frequency_THz: float = 10.0,
    drive_amplitude: float = 0.1,
    drive_type: str = "circularly_polarized",
    n_floquet_zones: int = 3
) -> Dict[str, Any]:
    """
    Generate structure with Floquet drive parameters.

    Floquet engineering uses periodic driving to modify band topology
    and create effective time-independent Hamiltonians.

    Args:
        structure: Base crystal structure
        drive_frequency_THz: Driving frequency in THz
        drive_amplitude: Electric field amplitude (dimensionless A·ω/c)
        drive_type: "circularly_polarized", "linearly_polarized", "elliptical"
        n_floquet_zones: Number of Floquet zones to consider

    Returns:
        Structure with Floquet parameters for calculations
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

    # Convert frequency to energy
    hbar_eV_s = 6.582e-16  # eV·s
    omega_eV = drive_frequency_THz * 1e12 * 2 * np.pi * hbar_eV_s

    # Calculate Floquet zone boundaries
    floquet_zones = []
    for n in range(-n_floquet_zones, n_floquet_zones + 1):
        floquet_zones.append({
            "zone": n,
            "energy_shift_eV": n * omega_eV
        })

    # Drive polarization vector
    if drive_type == "circularly_polarized":
        # Right-handed circular polarization
        polarization = {
            "type": "circular",
            "handedness": "right",
            "components": {
                "x": {"amplitude": drive_amplitude, "phase": 0},
                "y": {"amplitude": drive_amplitude, "phase": np.pi/2}
            }
        }
        effective_field = drive_amplitude * np.sqrt(2)
    elif drive_type == "linearly_polarized":
        polarization = {
            "type": "linear",
            "direction": [1, 0, 0],
            "components": {
                "x": {"amplitude": drive_amplitude, "phase": 0}
            }
        }
        effective_field = drive_amplitude
    else:  # elliptical
        polarization = {
            "type": "elliptical",
            "components": {
                "x": {"amplitude": drive_amplitude, "phase": 0},
                "y": {"amplitude": drive_amplitude * 0.5, "phase": np.pi/4}
            }
        }
        effective_field = drive_amplitude * np.sqrt(1.25)

    # Calculate dimensionless driving parameter
    # For graphene-like systems: A = eE/(ℏω)
    driving_parameter = effective_field  # Already dimensionless in input

    # Expected effects based on driving regime
    if driving_parameter < 0.1:
        regime = "perturbative"
        expected_effects = ["Weak gap opening", "Small band renormalization"]
    elif driving_parameter < 1.0:
        regime = "intermediate"
        expected_effects = ["Floquet-Bloch bands", "Dynamic localization possible"]
    else:
        regime = "strong"
        expected_effects = ["Strong hybridization", "Floquet topological transitions"]

    return {
        "success": True,
        "drive_type": drive_type,
        "frequency_THz": drive_frequency_THz,
        "frequency_eV": omega_eV,
        "amplitude": drive_amplitude,
        "polarization": polarization,
        "driving_parameter": driving_parameter,
        "regime": regime,
        "expected_effects": expected_effects,
        "n_floquet_zones": n_floquet_zones,
        "floquet_zones": floquet_zones,
        "calculation_notes": [
            "Use Floquet-TDDFT or Floquet Green's function methods",
            "Time-average over one period for effective Hamiltonian",
            f"Floquet quasi-energy range: ±{n_floquet_zones * omega_eV:.3f} eV"
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


def generate_phonon_pumped_cell(
    structure,
    phonon_mode: str = "IR_active",
    pump_amplitude_A: float = 0.1,
    pump_frequency_THz: float = 5.0,
    target_atom: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate structure with optical phonon pumping.

    Simulates coherent phonon excitation that can induce
    phase transitions or modify electronic structure.

    Args:
        structure: Base crystal structure
        phonon_mode: "IR_active", "Raman_active", "soft_mode"
        pump_amplitude_A: Atomic displacement amplitude in Angstroms
        pump_frequency_THz: Pump frequency
        target_atom: Specific atom to displace (None = use mode pattern)

    Returns:
        Structure at phonon-displaced geometry
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

    # Generate mode pattern based on type
    n_atoms = len(struct)

    if phonon_mode == "IR_active":
        # Dipolar mode - opposite displacements for +/- charges
        # Estimate charges from electronegativity difference
        mode_pattern = []
        for site in struct.sites:
            # Simple heuristic: heavier atoms move less
            mass = site.specie.atomic_mass
            weight = 1.0 / np.sqrt(mass)
            # Alternate direction for different species
            species_list = list(set(str(s.specie) for s in struct.sites))
            idx = species_list.index(str(site.specie))
            direction = np.array([0, 0, 1]) * (1 if idx % 2 == 0 else -1)
            mode_pattern.append(weight * direction)

    elif phonon_mode == "Raman_active":
        # Symmetric breathing mode
        mode_pattern = []
        center = np.mean([s.coords for s in struct.sites], axis=0)
        for site in struct.sites:
            direction = site.coords - center
            norm = np.linalg.norm(direction)
            if norm > 0.1:
                direction = direction / norm
            else:
                direction = np.array([0, 0, 1])
            mode_pattern.append(direction)

    elif phonon_mode == "soft_mode":
        # Ferroelectric soft mode - polar displacement
        mode_pattern = []
        for site in struct.sites:
            element = str(site.specie)
            # Typical soft mode: cation vs anion
            if element in ["Ti", "Zr", "Nb", "Ta", "Fe", "Pb", "Bi"]:
                mode_pattern.append(np.array([0, 0, 1]))
            elif element in ["O", "S", "Se", "F", "Cl"]:
                mode_pattern.append(np.array([0, 0, -0.5]))
            else:
                mode_pattern.append(np.array([0, 0, 0]))

    else:
        mode_pattern = [np.array([0, 0, 1])] * n_atoms

    # Apply displacement
    new_coords = []
    displacements = []

    for i, site in enumerate(struct.sites):
        if target_atom is not None and i != target_atom:
            displacement = np.zeros(3)
        else:
            displacement = pump_amplitude_A * np.array(mode_pattern[i])

        new_coords.append(site.coords + displacement)
        displacements.append(displacement.tolist())

    # Create displaced structure
    displaced_struct = Structure(
        lattice=struct.lattice,
        species=[str(s.specie) for s in struct.sites],
        coords=new_coords,
        coords_are_cartesian=True
    )

    # Calculate energy scale
    hbar = 6.582e-16  # eV·s
    phonon_energy_eV = pump_frequency_THz * 1e12 * 2 * np.pi * hbar

    return {
        "success": True,
        "phonon_mode": phonon_mode,
        "pump_amplitude_A": pump_amplitude_A,
        "pump_frequency_THz": pump_frequency_THz,
        "phonon_energy_eV": phonon_energy_eV,
        "phonon_energy_meV": phonon_energy_eV * 1000,
        "displacements_A": displacements,
        "max_displacement_A": max(np.linalg.norm(d) for d in displacements),
        "applications": [
            "Light-induced phase transitions",
            "Nonlinear phononics",
            "Hidden phase access"
        ],
        "structure": {
            "lattice": {"matrix": displaced_struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in displaced_struct.sites
            ]
        },
        "pymatgen_structure": displaced_struct.as_dict()
    }


def generate_warm_dense_matter(
    element: str = "Al",
    temperature_eV: float = 10.0,
    density_compression: float = 1.0,
    n_atoms: int = 64,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate warm dense matter configuration.

    Creates disordered high-temperature, high-density structures
    for WDM simulations (plasmas, ICF, planetary interiors).

    Args:
        element: Element symbol
        temperature_eV: Electronic temperature in eV
        density_compression: Density relative to ambient (ρ/ρ₀)
        n_atoms: Number of atoms
        seed: Random seed

    Returns:
        WDM configuration
    """
    from pymatgen.core import Structure, Lattice, Element

    if seed is not None:
        np.random.seed(seed)

    elem = Element(element)
    atomic_mass = elem.atomic_mass
    ambient_density = elem.density_of_solid  # g/cm³

    # Calculate compressed density
    compressed_density = ambient_density * density_compression

    # Calculate box size from density
    # ρ = N * m / V => V = N * m / ρ
    mass_per_atom_g = atomic_mass / 6.022e23  # g
    volume_cm3 = n_atoms * mass_per_atom_g / compressed_density
    volume_A3 = volume_cm3 * 1e24  # Convert to Å³

    # Cubic box
    box_length = volume_A3 ** (1/3)

    lattice = Lattice.cubic(box_length)

    # Generate random positions (uniform, then apply thermal jitter)
    # Start with simple cubic as reference, then randomize
    n_per_side = int(np.ceil(n_atoms ** (1/3)))

    coords = []
    for i in range(n_atoms):
        # Random uniform position
        coord = np.random.uniform(0, 1, 3)
        coords.append(coord.tolist())

    # Check for overlaps and adjust (minimum distance criterion)
    min_distance = 1.5  # Angstroms - typical atomic radius
    for iteration in range(100):  # Max iterations
        adjusted = False
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                diff = np.array(coords[i]) - np.array(coords[j])
                # Apply PBC
                diff = diff - np.round(diff)
                dist = np.linalg.norm(diff * box_length)
                if dist < min_distance:
                    # Push apart
                    push = (min_distance - dist) / 2 / box_length
                    direction = diff / max(np.linalg.norm(diff), 1e-10)
                    coords[i] = (np.array(coords[i]) + push * direction).tolist()
                    coords[j] = (np.array(coords[j]) - push * direction).tolist()
                    adjusted = True
        if not adjusted:
            break

    # Wrap coordinates
    coords = [[c % 1 for c in coord] for coord in coords]

    struct = Structure(
        lattice=lattice,
        species=[element] * n_atoms,
        coords=coords,
        coords_are_cartesian=False
    )

    # Calculate WDM parameters
    # Wigner-Seitz radius
    r_s_bohr = (3 * volume_A3 / (4 * np.pi * n_atoms)) ** (1/3) / 0.529177  # in Bohr

    # Coupling parameter Γ = e²/(r_s * k_B * T)
    # For T in eV: Γ = 27.2 / (r_s * T)
    gamma = 27.2 / (r_s_bohr * temperature_eV)

    # Degeneracy parameter θ = T / T_F
    # T_F ≈ 0.5 * (9π/4)^(2/3) / r_s² in eV
    T_fermi = 0.5 * (9 * np.pi / 4) ** (2/3) / (r_s_bohr ** 2) * 27.2
    theta = temperature_eV / T_fermi

    # Classify regime
    if gamma > 1 and theta < 1:
        regime = "strongly_coupled_degenerate"
    elif gamma > 1 and theta > 1:
        regime = "strongly_coupled_classical"
    elif gamma < 1 and theta < 1:
        regime = "weakly_coupled_degenerate"
    else:
        regime = "weakly_coupled_classical"

    return {
        "success": True,
        "element": element,
        "n_atoms": n_atoms,
        "temperature_eV": temperature_eV,
        "temperature_K": temperature_eV * 11604.5,
        "density_compression": density_compression,
        "density_g_cm3": compressed_density,
        "box_length_A": box_length,
        "wigner_seitz_radius_bohr": r_s_bohr,
        "coupling_parameter_Gamma": gamma,
        "degeneracy_parameter_theta": theta,
        "fermi_temperature_eV": T_fermi,
        "regime": regime,
        "is_wdm": 0.1 < theta < 10 and 0.1 < gamma < 10,
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }


def generate_laser_shocked_cell(
    element: str = "Fe",
    shock_pressure_GPa: float = 100.0,
    electron_temperature_eV: float = 1.0,
    ion_temperature_eV: float = 0.1,
    compression_ratio: float = 1.5,
    n_atoms: int = 64,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate laser-shocked two-temperature structure.

    Creates structures for simulating laser-matter interaction
    with separate electron and ion temperatures.

    Args:
        element: Element symbol
        shock_pressure_GPa: Shock pressure
        electron_temperature_eV: Hot electron temperature
        ion_temperature_eV: Cooler ion temperature
        compression_ratio: Density compression (1.0 = ambient)
        n_atoms: Number of atoms
        seed: Random seed

    Returns:
        Two-temperature shocked structure
    """
    from pymatgen.core import Structure, Lattice, Element

    if seed is not None:
        np.random.seed(seed)

    elem = Element(element)
    ambient_density = elem.density_of_solid

    # Hugoniot compression
    compressed_density = ambient_density * compression_ratio

    # Calculate box
    atomic_mass = elem.atomic_mass
    mass_per_atom = atomic_mass / 6.022e23
    volume_cm3 = n_atoms * mass_per_atom / compressed_density
    volume_A3 = volume_cm3 * 1e24
    box_length = volume_A3 ** (1/3)

    lattice = Lattice.cubic(box_length)

    # Generate thermally displaced positions
    # Ion thermal displacement ~ sqrt(k_B * T_ion / (m * ω²))
    # Approximate with Debye model
    debye_T_K = 400  # Approximate Debye temperature
    debye_T_eV = debye_T_K / 11604.5

    # RMS displacement
    thermal_amp = 0.1 * np.sqrt(ion_temperature_eV / debye_T_eV)  # Å

    # Start from BCC if Fe-like, FCC otherwise
    if element in ["Fe", "Cr", "W", "Mo", "V", "Nb", "Ta"]:
        # BCC arrangement
        base_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    else:
        # FCC arrangement
        base_coords = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    # Replicate and add thermal displacements
    n_per_side = int(np.ceil((n_atoms / len(base_coords)) ** (1/3)))
    coords = []

    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                for base in base_coords:
                    if len(coords) >= n_atoms:
                        break
                    coord = np.array([
                        (base[0] + i) / n_per_side,
                        (base[1] + j) / n_per_side,
                        (base[2] + k) / n_per_side
                    ])
                    # Add thermal displacement
                    thermal_disp = np.random.normal(0, thermal_amp / box_length, 3)
                    coord += thermal_disp
                    coord = coord % 1  # Wrap
                    coords.append(coord.tolist())

    coords = coords[:n_atoms]

    struct = Structure(
        lattice=lattice,
        species=[element] * len(coords),
        coords=coords,
        coords_are_cartesian=False
    )

    # Calculate particle velocity from Hugoniot
    # Up ≈ √(P * V₀ / ρ₀) for simple estimate
    particle_velocity = np.sqrt(shock_pressure_GPa * 1e9 / (ambient_density * 1e3))  # m/s

    return {
        "success": True,
        "element": element,
        "n_atoms": len(coords),
        "shock_pressure_GPa": shock_pressure_GPa,
        "compression_ratio": compression_ratio,
        "electron_temperature_eV": electron_temperature_eV,
        "electron_temperature_K": electron_temperature_eV * 11604.5,
        "ion_temperature_eV": ion_temperature_eV,
        "ion_temperature_K": ion_temperature_eV * 11604.5,
        "temperature_ratio_Te_Ti": electron_temperature_eV / ion_temperature_eV,
        "is_two_temperature": electron_temperature_eV / ion_temperature_eV > 2,
        "particle_velocity_km_s": particle_velocity / 1000,
        "thermal_displacement_rms_A": thermal_amp,
        "applications": [
            "Laser ablation",
            "Inertial confinement fusion",
            "High-energy-density physics"
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


def generate_bimeron(
    material: str = "CrI3",
    size_nm: float = 10.0,
    supercell_size: int = 20
) -> Dict[str, Any]:
    """
    Generate bimeron magnetic texture.

    Bimerons are in-plane analogs of skyrmions, relevant for
    easy-plane magnetic materials.

    Args:
        material: Magnetic material with easy-plane anisotropy
        size_nm: Bimeron size in nm
        supercell_size: Supercell dimension

    Returns:
        Structure with bimeron spin texture
    """
    from pymatgen.core import Structure, Lattice

    # Material database
    materials = {
        "CrI3": {"a": 6.87, "c": 20.0, "magnetic_element": "Cr"},
        "Fe3GeTe2": {"a": 3.99, "c": 16.4, "magnetic_element": "Fe"},
        "MnSi": {"a": 4.56, "c": 4.56, "magnetic_element": "Mn"}
    }

    if material not in materials:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_MATERIAL",
                "message": f"Unknown material: {material}",
                "available": list(materials.keys())
            }
        }

    params = materials[material]
    a = params["a"]
    mag_elem = params["magnetic_element"]

    # Create supercell
    lattice = Lattice.tetragonal(a * supercell_size, params["c"])

    # Place magnetic atoms on grid
    coords = []
    spins = []
    species = []

    center = supercell_size / 2
    r_bimeron = size_nm * 10 / a  # Convert to lattice units

    for i in range(supercell_size):
        for j in range(supercell_size):
            # Position
            x = (i + 0.5) / supercell_size
            y = (j + 0.5) / supercell_size
            coords.append([x, y, 0.5])
            species.append(mag_elem)

            # Calculate spin for bimeron
            dx = i - center + 0.5
            dy = j - center + 0.5
            r = np.sqrt(dx**2 + dy**2)
            phi = np.arctan2(dy, dx)

            # Bimeron profile (in-plane vortex with out-of-plane core)
            if r < r_bimeron:
                theta = np.pi * r / r_bimeron
                mx = np.cos(phi) * np.sin(theta)
                my = np.sin(phi) * np.sin(theta)
                mz = np.cos(theta)
            else:
                # Uniform in-plane background
                mx = 1.0
                my = 0.0
                mz = 0.0

            spins.append([mx, my, mz])

    struct = Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "texture_type": "bimeron",
        "material": material,
        "size_nm": size_nm,
        "n_atoms": len(struct),
        "magnetic_element": mag_elem,
        "spin_texture": spins,
        "topological_charge": 1,  # Bimeron has Q=1
        "structure": {
            "lattice": {"matrix": struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in struct.sites
            ]
        },
        "pymatgen_structure": struct.as_dict()
    }
