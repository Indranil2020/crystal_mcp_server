"""
Parametric Scan System for Structure Generation

Provides capabilities for generating series of structures with systematically
varying parameters, supporting:
- Z-direction scans for molecule-surface systems
- Distance scans between structures
- Angular rotation scans
- General parameter sweeps
- Trajectory/animation output

Example usage:
    # Z-scan of PTCDA molecule on NaCl(100) surface from 3 to 20 Å
    result = generate_z_scan(
        surface_material="NaCl",
        surface_miller=(1, 0, 0),
        molecule="PTCDA",
        z_start=3.0,
        z_end=20.0,
        z_step=0.1
    )
"""

from typing import Dict, Any, List, Optional, Tuple, Union
from copy import deepcopy
import importlib.util
import inspect
import numpy as np


def generate_z_scan(
    surface_structure: Optional[Dict[str, Any]] = None,
    surface_material: Optional[str] = None,
    surface_miller: Tuple[int, int, int] = (1, 0, 0),
    molecule: Optional[Union[str, Dict[str, Any]]] = None,
    z_start: float = 3.0,
    z_end: float = 20.0,
    z_step: float = 0.1,
    site: str = "ontop",
    site_index: Optional[int] = None,
    surface_layers: int = 4,
    surface_vacuum: float = 30.0,
    supercell: Tuple[int, int, int] = (3, 3, 1),
    output_format: str = "trajectory",
    include_energy_template: bool = True
) -> Dict[str, Any]:
    """
    Generate a z-direction scan of a molecule on a surface.

    Creates a series of structures with the molecule placed at different
    heights above the surface, suitable for:
    - Potential energy surface (PES) scanning
    - Adsorption energy curves
    - Binding distance optimization

    Parameters
    ----------
    surface_structure : dict, optional
        Pre-generated surface structure dict. If None, generates from surface_material.
    surface_material : str, optional
        Material for surface generation (e.g., "NaCl", "Au", "TiO2").
        Used if surface_structure is None.
    surface_miller : tuple
        Miller indices for surface (default: (1, 0, 0))
    molecule : str or dict
        Molecule name (e.g., "PTCDA", "H2O", "CO") or pre-generated molecule dict
    z_start : float
        Starting z-distance in Angstroms (default: 3.0)
    z_end : float
        Ending z-distance in Angstroms (default: 20.0)
    z_step : float
        Step size in Angstroms (default: 0.1)
    site : str
        Adsorption site type ("ontop", "bridge", "hollow", "center")
    site_index : int, optional
        Specific surface atom index for placement
    surface_layers : int
        Number of surface layers (default: 4)
    surface_vacuum : float
        Vacuum above surface in Angstroms (default: 30.0)
    supercell : tuple
        Surface supercell size (default: (3, 3, 1))
    output_format : str
        Output format: "trajectory", "list", "animated_xyz", "xsf_animation"
    include_energy_template : bool
        Include template for energy calculations (default: True)

    Returns
    -------
    dict
        Result containing:
        - structures: List of structure dicts at each z-height
        - z_values: Array of z-values sampled
        - n_frames: Number of frames in trajectory
        - trajectory_xyz: XYZ trajectory string (if output_format includes trajectory)
        - metadata: Scan parameters and info
    """
    from ase import Atoms

    # Generate z-values
    z_values = np.arange(z_start, z_end + z_step/2, z_step)
    n_frames = len(z_values)

    # Get or generate surface
    if surface_structure is not None:
        surface_atoms = _dict_to_atoms(surface_structure)
    elif surface_material is not None:
        surface_atoms = _generate_surface(
            material=surface_material,
            miller=surface_miller,
            layers=surface_layers,
            vacuum=surface_vacuum,
            supercell=supercell
        )
    else:
        return {
            "success": False,
            "error": "Must provide either surface_structure or surface_material"
        }

    # Get or generate molecule
    if isinstance(molecule, dict):
        mol_atoms = _dict_to_atoms(molecule)
    elif isinstance(molecule, str):
        mol_atoms = _generate_molecule(molecule)
    else:
        return {
            "success": False,
            "error": "Must provide molecule name (str) or structure (dict)"
        }

    # Find adsorption site position (xy-coordinates)
    if site_index is not None:
        site_pos = surface_atoms.positions[site_index, :2]
    else:
        site_pos = _find_adsorption_site(surface_atoms, site)

    # Generate trajectory of structures
    structures = []
    trajectory_xyz_lines = [f"{len(surface_atoms) + len(mol_atoms)}\n"]

    for i, z in enumerate(z_values):
        # Create combined structure
        combined = surface_atoms.copy()
        mol_copy = mol_atoms.copy()

        # Center molecule and position above surface
        mol_center = mol_copy.get_center_of_mass()
        surface_top_z = surface_atoms.positions[:, 2].max()

        # Translate molecule to site position and z-height
        mol_copy.translate([
            site_pos[0] - mol_center[0],
            site_pos[1] - mol_center[1],
            surface_top_z + z - mol_center[2]
        ])

        # Combine surface and molecule
        combined = combined + mol_copy

        # Convert to dict and store
        struct_dict = _atoms_to_dict(combined)
        struct_dict["z_height"] = float(z)
        struct_dict["frame_index"] = i
        structures.append(struct_dict)

        # Build XYZ trajectory frame
        comment = f"z_height={z:.3f} frame={i}"
        trajectory_xyz_lines.append(f"{comment}\n")
        for atom in combined:
            trajectory_xyz_lines.append(
                f"{atom.symbol:2s} {atom.position[0]:12.6f} {atom.position[1]:12.6f} {atom.position[2]:12.6f}\n"
            )
        if i < n_frames - 1:
            trajectory_xyz_lines.append(f"{len(combined)}\n")

    # Build result
    result = {
        "success": True,
        "operation": "z_scan",
        "structures": structures,
        "z_values": z_values.tolist(),
        "n_frames": n_frames,
        "z_start": z_start,
        "z_end": z_end,
        "z_step": z_step,
        "molecule": molecule if isinstance(molecule, str) else "custom",
        "surface_material": surface_material,
        "surface_miller": list(surface_miller),
        "site": site,
        "trajectory_xyz": "".join(trajectory_xyz_lines),
    }

    # Add XSF animation if requested
    if output_format in ["xsf_animation", "all"]:
        result["trajectory_xsf"] = _generate_xsf_animation(structures)

    # Add energy calculation template
    if include_energy_template:
        result["energy_template"] = {
            "description": "Template for computing binding energy curve",
            "calculation_steps": [
                "1. Compute E_surface (isolated surface)",
                "2. Compute E_molecule (isolated molecule in vacuum)",
                "3. For each z: E_binding(z) = E_combined(z) - E_surface - E_molecule",
                "4. Find z_min where E_binding is minimum"
            ],
            "expected_output": {
                "binding_energies": "Array of binding energies at each z",
                "equilibrium_distance": "z value at minimum binding energy",
                "binding_energy_min": "Minimum binding energy"
            }
        }

    return result


def generate_distance_scan(
    structure1: Dict[str, Any],
    structure2: Dict[str, Any],
    distance_start: float = 2.0,
    distance_end: float = 10.0,
    distance_step: float = 0.1,
    approach_axis: str = "z",
    align_centers: bool = True,
    output_format: str = "trajectory"
) -> Dict[str, Any]:
    """
    Generate a distance scan between two structures.

    Creates a series of structures with varying separation distance,
    useful for:
    - Dimer interaction curves
    - Layer separation scans
    - Interface formation studies

    Parameters
    ----------
    structure1 : dict
        First structure (kept fixed)
    structure2 : dict
        Second structure (translated)
    distance_start : float
        Starting distance in Angstroms (default: 2.0)
    distance_end : float
        Ending distance in Angstroms (default: 10.0)
    distance_step : float
        Step size in Angstroms (default: 0.1)
    approach_axis : str
        Axis for approach ("x", "y", "z") (default: "z")
    align_centers : bool
        Align structure centers perpendicular to approach axis (default: True)
    output_format : str
        Output format: "trajectory", "list"

    Returns
    -------
    dict
        Result with trajectory of structures at each distance
    """
    from ase import Atoms

    # Parse axis
    axis_map = {"x": 0, "y": 1, "z": 2}
    axis_idx = axis_map.get(approach_axis.lower(), 2)

    # Convert to Atoms
    atoms1 = _dict_to_atoms(structure1)
    atoms2 = _dict_to_atoms(structure2)

    # Generate distances
    distances = np.arange(distance_start, distance_end + distance_step/2, distance_step)
    n_frames = len(distances)

    structures = []
    trajectory_xyz_lines = []

    for i, d in enumerate(distances):
        # Copy structures
        a1 = atoms1.copy()
        a2 = atoms2.copy()

        if align_centers:
            # Align centers in perpendicular directions
            c1 = a1.get_center_of_mass()
            c2 = a2.get_center_of_mass()

            shift = np.zeros(3)
            for j in range(3):
                if j != axis_idx:
                    shift[j] = c1[j] - c2[j]
            a2.translate(shift)

        # Get reference points for distance
        max1 = a1.positions[:, axis_idx].max()
        min2 = a2.positions[:, axis_idx].min()

        # Translate structure2 to target distance
        translation = np.zeros(3)
        translation[axis_idx] = max1 + d - min2
        a2.translate(translation)

        # Combine
        combined = a1 + a2

        # Store
        struct_dict = _atoms_to_dict(combined)
        struct_dict["distance"] = float(d)
        struct_dict["frame_index"] = i
        structures.append(struct_dict)

        # XYZ frame
        n_atoms = len(combined)
        trajectory_xyz_lines.append(f"{n_atoms}\n")
        trajectory_xyz_lines.append(f"distance={d:.3f} frame={i}\n")
        for atom in combined:
            trajectory_xyz_lines.append(
                f"{atom.symbol:2s} {atom.position[0]:12.6f} {atom.position[1]:12.6f} {atom.position[2]:12.6f}\n"
            )

    return {
        "success": True,
        "operation": "distance_scan",
        "structures": structures,
        "distances": distances.tolist(),
        "n_frames": n_frames,
        "approach_axis": approach_axis,
        "trajectory_xyz": "".join(trajectory_xyz_lines)
    }


def generate_rotation_scan(
    structure: Dict[str, Any],
    rotation_axis: str = "z",
    angle_start: float = 0.0,
    angle_end: float = 360.0,
    angle_step: float = 5.0,
    rotation_center: Optional[List[float]] = None,
    rotate_cell: bool = False,
    output_format: str = "trajectory"
) -> Dict[str, Any]:
    """
    Generate a rotation scan of a structure.

    Creates a series of structures with incremental rotations,
    useful for:
    - Molecular orientation studies
    - Rotational barrier calculations
    - Symmetry analysis

    Parameters
    ----------
    structure : dict
        Structure to rotate
    rotation_axis : str
        Axis of rotation ("x", "y", "z") (default: "z")
    angle_start : float
        Starting angle in degrees (default: 0.0)
    angle_end : float
        Ending angle in degrees (default: 360.0)
    angle_step : float
        Step size in degrees (default: 5.0)
    rotation_center : list, optional
        Center point for rotation. If None, uses center of mass.
    rotate_cell : bool
        Also rotate the cell vectors (default: False)
    output_format : str
        Output format: "trajectory", "list"

    Returns
    -------
    dict
        Result with trajectory of rotated structures
    """
    from ase import Atoms
    from scipy.spatial.transform import Rotation as R

    atoms = _dict_to_atoms(structure)

    # Generate angles
    angles = np.arange(angle_start, angle_end + angle_step/2, angle_step)
    n_frames = len(angles)

    # Get rotation center
    if rotation_center is None:
        center = atoms.get_center_of_mass()
    else:
        center = np.array(rotation_center)

    # Axis vector
    axis_map = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
    axis_vec = axis_map.get(rotation_axis.lower(), [0, 0, 1])

    structures = []
    trajectory_xyz_lines = []

    for i, angle in enumerate(angles):
        rotated = atoms.copy()

        # Create rotation
        rot = R.from_rotvec(np.radians(angle) * np.array(axis_vec))

        # Apply rotation around center
        positions = rotated.positions - center
        rotated_positions = rot.apply(positions) + center
        rotated.set_positions(rotated_positions)

        if rotate_cell:
            cell = rotated.get_cell()
            rotated_cell = rot.apply(cell)
            rotated.set_cell(rotated_cell)

        # Store
        struct_dict = _atoms_to_dict(rotated)
        struct_dict["rotation_angle"] = float(angle)
        struct_dict["frame_index"] = i
        structures.append(struct_dict)

        # XYZ frame
        n_atoms = len(rotated)
        trajectory_xyz_lines.append(f"{n_atoms}\n")
        trajectory_xyz_lines.append(f"angle={angle:.1f} frame={i}\n")
        for atom in rotated:
            trajectory_xyz_lines.append(
                f"{atom.symbol:2s} {atom.position[0]:12.6f} {atom.position[1]:12.6f} {atom.position[2]:12.6f}\n"
            )

    return {
        "success": True,
        "operation": "rotation_scan",
        "structures": structures,
        "angles": angles.tolist(),
        "n_frames": n_frames,
        "rotation_axis": rotation_axis,
        "trajectory_xyz": "".join(trajectory_xyz_lines)
    }


def generate_parameter_sweep(
    generator_name: str,
    base_params: Dict[str, Any],
    sweep_param: str,
    sweep_values: List[Any],
    output_format: str = "list"
) -> Dict[str, Any]:
    """
    Generate a parameter sweep using any generator function.

    Creates structures by varying a single parameter across specified values,
    enabling systematic exploration of parameter space.

    Parameters
    ----------
    generator_name : str
        Name of the generator function (e.g., "generate_perovskite")
    base_params : dict
        Base parameters for the generator
    sweep_param : str
        Parameter name to sweep
    sweep_values : list
        Values to sweep over
    output_format : str
        Output format: "list", "trajectory"

    Returns
    -------
    dict
        Result with structures for each parameter value

    Example
    -------
    >>> # Sweep lattice parameter for perovskite
    >>> result = generate_parameter_sweep(
    ...     generator_name="generate_perovskite",
    ...     base_params={"a_site": "Ba", "b_site": "Ti", "x_site": "O"},
    ...     sweep_param="a",
    ...     sweep_values=[3.9, 3.95, 4.0, 4.05, 4.1]
    ... )
    """
    from generators import GENERATOR_REGISTRY
    import importlib

    # Find the generator
    generator_info = None
    for category, data in GENERATOR_REGISTRY.items():
        if generator_name in data["operations"]:
            generator_info = data["operations"][generator_name]
            break

    if generator_info is None:
        return {
            "success": False,
            "error": f"Generator '{generator_name}' not found in registry"
        }

    # Import the function
    module_path = generator_info["module"]
    func_name = generator_info["function"]

    if importlib.util.find_spec(module_path) is None:
        return {
            "success": False,
            "error": f"Module '{module_path}' not found for generator '{generator_name}'"
        }

    module = importlib.import_module(module_path)
    generator_func = getattr(module, func_name, None)
    if generator_func is None or not callable(generator_func):
        return {
            "success": False,
            "error": f"Generator function '{func_name}' not found in '{module_path}'"
        }

    required, allowed, accepts_kwargs = _get_generator_spec(generator_func)

    # Generate structures
    structures = []
    errors = []

    for i, value in enumerate(sweep_values):
        params = base_params.copy()
        params[sweep_param] = value

        valid, error = _validate_generator_params(required, allowed, accepts_kwargs, params)
        if not valid:
            errors.append({"index": i, "value": value, "error": error})
            continue

        result = generator_func(**params)
        if isinstance(result, dict):
            result["sweep_value"] = value
            result["sweep_index"] = i
            structures.append(result)
        else:
            errors.append({"index": i, "value": value, "error": "Unexpected result type"})

    return {
        "success": True,
        "operation": "parameter_sweep",
        "generator": generator_name,
        "sweep_param": sweep_param,
        "sweep_values": sweep_values,
        "n_structures": len(structures),
        "structures": structures,
        "errors": errors if errors else None
    }


def generate_multi_parameter_scan(
    generator_name: str,
    base_params: Dict[str, Any],
    scan_params: Dict[str, List[Any]],
    scan_type: str = "grid",
    output_format: str = "list"
) -> Dict[str, Any]:
    """
    Generate a multi-parameter scan using any generator function.

    Creates structures by varying multiple parameters, either as a grid
    or along a path in parameter space.

    Parameters
    ----------
    generator_name : str
        Name of the generator function
    base_params : dict
        Base parameters for the generator
    scan_params : dict
        Dictionary mapping parameter names to lists of values
    scan_type : str
        "grid" for full grid, "path" for sequential variation
    output_format : str
        Output format: "list", "trajectory"

    Returns
    -------
    dict
        Result with structures for each parameter combination

    Example
    -------
    >>> # 2D grid scan of strain
    >>> result = generate_multi_parameter_scan(
    ...     generator_name="apply_strain",
    ...     base_params={"structure_dict": si_structure},
    ...     scan_params={
    ...         "strain": [-0.02, -0.01, 0, 0.01, 0.02],
    ...         "strain_type": ["biaxial", "uniaxial"]
    ...     },
    ...     scan_type="grid"
    ... )
    """
    from itertools import product
    from generators import GENERATOR_REGISTRY
    import importlib

    # Find and import generator (same as above)
    generator_info = None
    for category, data in GENERATOR_REGISTRY.items():
        if generator_name in data["operations"]:
            generator_info = data["operations"][generator_name]
            break

    if generator_info is None:
        return {
            "success": False,
            "error": f"Generator '{generator_name}' not found"
        }

    module_path = generator_info["module"]
    func_name = generator_info["function"]
    if importlib.util.find_spec(module_path) is None:
        return {
            "success": False,
            "error": f"Module '{module_path}' not found for generator '{generator_name}'"
        }

    module = importlib.import_module(module_path)
    generator_func = getattr(module, func_name, None)
    if generator_func is None or not callable(generator_func):
        return {
            "success": False,
            "error": f"Generator function '{func_name}' not found in '{module_path}'"
        }

    required, allowed, accepts_kwargs = _get_generator_spec(generator_func)

    # Generate parameter combinations
    param_names = list(scan_params.keys())
    param_values = list(scan_params.values())

    if scan_type == "grid":
        combinations = list(product(*param_values))
    elif scan_type == "path":
        # Sequential variation - all lists must be same length
        n_points = len(param_values[0])
        combinations = [tuple(pv[i] for pv in param_values) for i in range(n_points)]
    else:
        return {"success": False, "error": f"Unknown scan_type: {scan_type}"}

    structures = []
    errors = []

    for i, combo in enumerate(combinations):
        params = base_params.copy()
        param_dict = dict(zip(param_names, combo))
        params.update(param_dict)

        valid, error = _validate_generator_params(required, allowed, accepts_kwargs, params)
        if not valid:
            errors.append({"index": i, "params": param_dict, "error": error})
            continue

        result = generator_func(**params)
        if isinstance(result, dict):
            result["scan_params"] = param_dict
            result["scan_index"] = i
            structures.append(result)

    return {
        "success": True,
        "operation": "multi_parameter_scan",
        "generator": generator_name,
        "scan_params": scan_params,
        "scan_type": scan_type,
        "n_combinations": len(combinations),
        "n_structures": len(structures),
        "structures": structures,
        "errors": errors if errors else None
    }


def generate_trajectory_animation(
    structures: List[Dict[str, Any]],
    output_format: str = "xyz",
    frame_labels: Optional[List[str]] = None,
    include_cell: bool = True
) -> Dict[str, Any]:
    """
    Generate animation files from a list of structures.

    Parameters
    ----------
    structures : list
        List of structure dictionaries
    output_format : str
        Output format: "xyz", "xsf", "extxyz", "html"
    frame_labels : list, optional
        Labels for each frame
    include_cell : bool
        Include cell information (default: True)

    Returns
    -------
    dict
        Result with animation content in requested format
    """
    n_frames = len(structures)

    if output_format == "xyz":
        content = _generate_xyz_trajectory(structures, frame_labels)
    elif output_format == "xsf":
        content = _generate_xsf_animation(structures)
    elif output_format == "extxyz":
        content = _generate_extxyz_trajectory(structures, frame_labels)
    elif output_format == "html":
        content = _generate_html_animation(structures)
    else:
        return {"success": False, "error": f"Unknown format: {output_format}"}

    return {
        "success": True,
        "operation": "trajectory_animation",
        "format": output_format,
        "n_frames": n_frames,
        "content": content
    }


# =============================================================================
# Helper Functions
# =============================================================================

def _get_generator_spec(generator_func):
    signature = inspect.signature(generator_func)
    required = []
    allowed = set()
    accepts_kwargs = False

    for name, param in signature.parameters.items():
        if param.kind == inspect.Parameter.VAR_KEYWORD:
            accepts_kwargs = True
        if param.kind == inspect.Parameter.VAR_POSITIONAL:
            continue
        allowed.add(name)
        if param.default is inspect.Parameter.empty and param.kind in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        ):
            required.append(name)

    return required, allowed, accepts_kwargs


def _validate_generator_params(required, allowed, accepts_kwargs, params):
    missing = [name for name in required if name not in params]
    if missing:
        return False, f"Missing parameters: {', '.join(missing)}"

    if not accepts_kwargs:
        unknown = [name for name in params if name not in allowed]
        if unknown:
            return False, f"Unknown parameters: {', '.join(unknown)}"

    return True, None


def _ase_molecule_or_none(name: str):
    from ase.collections import g2
    import ase.build.molecule as ase_molecule_module

    known = set(getattr(g2, "names", []))
    extra = getattr(ase_molecule_module, "extra", {})
    if isinstance(extra, dict):
        known.update(extra.keys())

    if name in known:
        return ase_molecule_module.molecule(name)

    return None


def _dict_to_atoms(structure: Dict[str, Any]):
    """Convert structure dict to ASE Atoms object."""
    from ase import Atoms

    if "atoms" in structure:
        # Already ASE-compatible format
        atoms_data = structure["atoms"]
    else:
        atoms_data = structure

    species = atoms_data.get("species", atoms_data.get("symbols", []))
    positions = atoms_data.get("positions", atoms_data.get("coords", []))
    cell = atoms_data.get("cell", atoms_data.get("lattice", None))
    pbc = atoms_data.get("pbc", True if cell else False)

    atoms = Atoms(
        symbols=species,
        positions=positions,
        cell=cell,
        pbc=pbc
    )

    return atoms


def _atoms_to_dict(atoms) -> Dict[str, Any]:
    """Convert ASE Atoms to structure dict."""
    return {
        "species": list(atoms.get_chemical_symbols()),
        "positions": atoms.positions.tolist(),
        "cell": atoms.get_cell().tolist(),
        "pbc": list(atoms.pbc),
        "n_atoms": len(atoms),
        "formula": atoms.get_chemical_formula()
    }


def _generate_surface(
    material: str,
    miller: Tuple[int, int, int],
    layers: int,
    vacuum: float,
    supercell: Tuple[int, int, int]
):
    """Generate a surface structure."""
    from ase.build import surface, bulk, fcc111, fcc100, bcc110, bcc100
    from ase.build import add_vacuum

    # Common materials lattice parameters (Å)
    lattice_params = {
        "NaCl": {"a": 5.64, "structure": "rocksalt"},
        "Au": {"a": 4.08, "structure": "fcc"},
        "Ag": {"a": 4.09, "structure": "fcc"},
        "Cu": {"a": 3.61, "structure": "fcc"},
        "Pt": {"a": 3.92, "structure": "fcc"},
        "Pd": {"a": 3.89, "structure": "fcc"},
        "Fe": {"a": 2.87, "structure": "bcc"},
        "W": {"a": 3.16, "structure": "bcc"},
        "TiO2": {"a": 4.59, "c": 2.96, "structure": "rutile"},
        "Al2O3": {"a": 4.76, "c": 12.99, "structure": "corundum"},
        "Si": {"a": 5.43, "structure": "diamond"},
        "Ge": {"a": 5.66, "structure": "diamond"},
    }

    params = lattice_params.get(material, {"a": 4.0, "structure": "fcc"})

    if params["structure"] == "rocksalt":
        # NaCl-type structure
        from ase.spacegroup import crystal
        a = params["a"]
        atoms = crystal(
            ['Na', 'Cl'],
            [(0, 0, 0), (0.5, 0.5, 0.5)],
            spacegroup=225,
            cellpar=[a, a, a, 90, 90, 90]
        )
        slab = surface(atoms, miller, layers, vacuum=vacuum)
    elif params["structure"] in ["fcc", "bcc", "diamond"]:
        from ase.build import bulk as ase_bulk
        bulk_atoms = ase_bulk(material, params["structure"], a=params["a"])
        slab = surface(bulk_atoms, miller, layers, vacuum=vacuum)
    else:
        # Fallback
        from ase.build import bulk as ase_bulk
        bulk_atoms = ase_bulk(material, "fcc", a=params.get("a", 4.0))
        slab = surface(bulk_atoms, miller, layers, vacuum=vacuum)

    # Apply supercell
    slab = slab * supercell

    return slab


def _generate_molecule(name: str):
    """Generate a molecule structure."""
    mol = _ase_molecule_or_none(name)
    if mol is not None:
        return mol

    # Extended molecule database
    molecules = {
        "PTCDA": {
            "formula": "C24H8O6",
            "description": "Perylene-3,4,9,10-tetracarboxylic dianhydride",
            # Simplified planar structure
            "atoms": _build_ptcda()
        },
        "pentacene": {
            "atoms": _build_pentacene()
        },
        "C60": {
            "atoms": _build_c60()
        }
    }

    if name in molecules:
        return molecules[name]["atoms"]

    # Try building from formula
    return _build_molecule_from_formula(name)


def _build_ptcda():
    """Build PTCDA molecule (simplified planar model)."""
    from ase import Atoms
    import numpy as np

    # PTCDA: C24H8O6 - perylene core with anhydride groups
    # Using approximate coordinates for a planar molecule
    # Perylene backbone (16 carbons) + 4 C for anhydride rings + 4 C bridging
    # + 8 H + 6 O (4 C=O + 2 bridge O)

    a = 1.42  # C-C bond length
    h = 2.46  # Perylene width parameter

    # Simplified planar model - perylene core
    positions = []
    symbols = []

    # Central perylene unit (simplified as connected hexagons)
    # Row of benzene rings along x-axis
    for ix in range(-2, 3):
        x_off = ix * 2.46
        for iy in [-1, 0, 1]:
            if abs(ix) == 2 and abs(iy) == 1:
                continue  # Skip corners
            positions.append([x_off, iy * 1.42, 0.0])
            symbols.append('C')

    # Add anhydride groups at ends
    for sign in [-1, 1]:
        x = sign * 5.5
        for y_off, sym in [(0, 'C'), (1.2, 'O'), (-1.2, 'O'), (0, 'O')]:
            positions.append([x, y_off, 0.0])
            symbols.append(sym)

    # Add hydrogens around perimeter (simplified)
    for i in range(8):
        angle = i * np.pi / 4
        r = 5.0
        positions.append([r * np.cos(angle), r * np.sin(angle), 0.0])
        symbols.append('H')

    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.center()

    return atoms


def _build_pentacene():
    """Build pentacene molecule."""
    from ase import Atoms
    import numpy as np

    a = 1.42  # C-C bond
    positions = []
    symbols = []

    # 5 fused benzene rings
    for ring in range(5):
        x_off = ring * 2.46
        for i in range(6):
            angle = i * np.pi / 3
            x = x_off + 0.71 * np.cos(angle)
            y = 0.71 * np.sin(angle)
            positions.append([x, y, 0.0])
            symbols.append('C')

    # Add H atoms
    for i, (x, y, z) in enumerate(positions[:]):
        if len([p for p in positions if np.sqrt((p[0]-x)**2 + (p[1]-y)**2) < 1.6]) < 3:
            dx = x - np.mean([p[0] for p in positions])
            dy = y - np.mean([p[1] for p in positions])
            norm = np.sqrt(dx**2 + dy**2)
            positions.append([x + dx/norm * 1.08, y + dy/norm * 1.08, 0.0])
            symbols.append('H')

    atoms = Atoms(symbols=symbols[:30], positions=positions[:30])  # Truncate to reasonable size
    atoms.center()

    return atoms


def _build_c60():
    """Build C60 fullerene."""
    mol = _ase_molecule_or_none("C60")
    if mol is not None:
        return mol

    # Fallback: icosahedral approximation
    from ase import Atoms
    import numpy as np

    # Golden ratio
    phi = (1 + np.sqrt(5)) / 2

    # Icosahedron vertices
    vertices = []
    for s1 in [-1, 1]:
        for s2 in [-1, 1]:
            vertices.append([0, s1 * 1, s2 * phi])
            vertices.append([s1 * 1, s2 * phi, 0])
            vertices.append([s1 * phi, 0, s2 * 1])

    # Scale to C60 radius (~3.5 Å)
    vertices = np.array(vertices) * 3.5 / np.linalg.norm(vertices[0])

    atoms = Atoms('C' * len(vertices), positions=vertices)
    return atoms


def _build_molecule_from_formula(formula: str):
    """Attempt to build molecule from chemical formula."""
    from ase import Atoms

    # Very basic fallback - single atom
    return Atoms(formula)


def _find_adsorption_site(surface_atoms, site_type: str) -> np.ndarray:
    """Find adsorption site on surface."""
    import numpy as np

    # Get top layer atoms
    z_coords = surface_atoms.positions[:, 2]
    z_max = z_coords.max()
    top_mask = z_coords > z_max - 0.5
    top_positions = surface_atoms.positions[top_mask]

    if site_type == "ontop" or site_type == "top":
        # On top of first atom in top layer
        return top_positions[0, :2]
    elif site_type == "center":
        # Center of cell
        cell = surface_atoms.get_cell()
        return np.array([cell[0, 0]/2, cell[1, 1]/2])
    elif site_type == "bridge":
        # Midpoint between first two atoms
        if len(top_positions) >= 2:
            return (top_positions[0, :2] + top_positions[1, :2]) / 2
    elif site_type == "hollow":
        # Centroid of first three atoms
        if len(top_positions) >= 3:
            return np.mean(top_positions[:3, :2], axis=0)

    # Default to first atom
    return top_positions[0, :2]


def _generate_xyz_trajectory(structures: List[Dict[str, Any]], labels: Optional[List[str]] = None) -> str:
    """Generate XYZ trajectory string."""
    lines = []

    for i, struct in enumerate(structures):
        n_atoms = struct.get("n_atoms", len(struct.get("species", [])))
        label = labels[i] if labels else f"frame_{i}"

        lines.append(f"{n_atoms}\n")
        lines.append(f"{label}\n")

        species = struct.get("species", [])
        positions = struct.get("positions", [])

        for sym, pos in zip(species, positions):
            lines.append(f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")

    return "".join(lines)


def _generate_xsf_animation(structures: List[Dict[str, Any]]) -> str:
    """Generate XSF animation format."""
    n_frames = len(structures)
    lines = [f"ANIMSTEPS {n_frames}\n", "CRYSTAL\n"]

    for i, struct in enumerate(structures):
        lines.append(f"PRIMVEC {i+1}\n")
        cell = struct.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])
        for vec in cell:
            lines.append(f"  {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")

        lines.append(f"PRIMCOORD {i+1}\n")
        species = struct.get("species", [])
        positions = struct.get("positions", [])
        n_atoms = len(species)
        lines.append(f"  {n_atoms} 1\n")

        # Element number mapping
        from ase.data import atomic_numbers
        for sym, pos in zip(species, positions):
            z = atomic_numbers.get(sym, 1)
            lines.append(f"  {z:3d} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")

    return "".join(lines)


def _generate_extxyz_trajectory(structures: List[Dict[str, Any]], labels: Optional[List[str]] = None) -> str:
    """Generate extended XYZ trajectory."""
    lines = []

    for i, struct in enumerate(structures):
        n_atoms = struct.get("n_atoms", len(struct.get("species", [])))
        cell = struct.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])
        label = labels[i] if labels else f"frame_{i}"

        # Format cell as Lattice="..."
        cell_str = " ".join([f"{v:.6f}" for row in cell for v in row])

        lines.append(f"{n_atoms}\n")
        lines.append(f'Lattice="{cell_str}" Properties=species:S:1:pos:R:3 comment="{label}"\n')

        species = struct.get("species", [])
        positions = struct.get("positions", [])

        for sym, pos in zip(species, positions):
            lines.append(f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")

    return "".join(lines)


def _generate_html_animation(structures: List[Dict[str, Any]]) -> str:
    """Generate HTML animation viewer using Three.js."""
    n_frames = len(structures)

    # Serialize structure data for JavaScript
    import json
    frames_json = json.dumps([{
        "species": s.get("species", []),
        "positions": s.get("positions", []),
        "cell": s.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    } for s in structures])

    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>Structure Animation - {n_frames} Frames</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: sans-serif; }}
        #controls {{
            position: absolute;
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(255,255,255,0.9);
            padding: 10px 20px;
            border-radius: 5px;
            display: flex;
            gap: 10px;
            align-items: center;
        }}
        button {{ padding: 8px 16px; cursor: pointer; }}
        #slider {{ width: 300px; }}
        #frame-info {{ min-width: 80px; }}
    </style>
</head>
<body>
    <div id="container"></div>
    <div id="controls">
        <button id="play-pause">Play</button>
        <input type="range" id="slider" min="0" max="{n_frames-1}" value="0">
        <span id="frame-info">Frame: 0/{n_frames-1}</span>
    </div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <script>
        const frames = {frames_json};
        let currentFrame = 0;
        let playing = false;
        let atomMeshes = [];

        // Element colors
        const colors = {{
            'H': 0xFFFFFF, 'C': 0x404040, 'N': 0x0000FF, 'O': 0xFF0000,
            'S': 0xFFFF00, 'P': 0xFFA500, 'Na': 0xAB82FF, 'Cl': 0x00FF00,
            'Fe': 0xE25822, 'Au': 0xFFD700, 'Ag': 0xC0C0C0, 'Cu': 0xB87333,
            'Pt': 0xD4D4D4, 'Ti': 0x878787, 'Si': 0xF0C8A0, 'default': 0x808080
        }};

        const radii = {{
            'H': 0.31, 'C': 0.77, 'N': 0.75, 'O': 0.73, 'S': 1.02, 'P': 1.06,
            'Na': 1.66, 'Cl': 0.99, 'Fe': 1.32, 'Au': 1.44, 'Ag': 1.45, 'Cu': 1.28,
            'Pt': 1.39, 'Ti': 1.47, 'Si': 1.17, 'default': 1.0
        }};

        // Scene setup
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x1a1a2e);

        const camera = new THREE.PerspectiveCamera(75, window.innerWidth/window.innerHeight, 0.1, 1000);
        camera.position.set(20, 20, 20);

        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.getElementById('container').appendChild(renderer.domElement);

        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;

        // Lighting
        scene.add(new THREE.AmbientLight(0x404040, 0.5));
        const light = new THREE.DirectionalLight(0xffffff, 0.8);
        light.position.set(10, 10, 10);
        scene.add(light);

        function createAtoms(frame) {{
            // Clear existing
            atomMeshes.forEach(m => scene.remove(m));
            atomMeshes = [];

            const data = frames[frame];
            for (let i = 0; i < data.species.length; i++) {{
                const elem = data.species[i];
                const pos = data.positions[i];
                const color = colors[elem] || colors.default;
                const radius = (radii[elem] || radii.default) * 0.4;

                const geo = new THREE.SphereGeometry(radius, 16, 16);
                const mat = new THREE.MeshPhongMaterial({{ color: color }});
                const mesh = new THREE.Mesh(geo, mat);
                mesh.position.set(pos[0], pos[1], pos[2]);
                scene.add(mesh);
                atomMeshes.push(mesh);
            }}
        }}

        function updateFrame(frame) {{
            currentFrame = frame;
            createAtoms(frame);
            document.getElementById('slider').value = frame;
            document.getElementById('frame-info').textContent = `Frame: ${{frame}}/{n_frames-1}`;
        }}

        document.getElementById('play-pause').onclick = () => {{
            playing = !playing;
            document.getElementById('play-pause').textContent = playing ? 'Pause' : 'Play';
        }};

        document.getElementById('slider').oninput = (e) => {{
            updateFrame(parseInt(e.target.value));
        }};

        let lastTime = 0;
        function animate(time) {{
            requestAnimationFrame(animate);

            if (playing && time - lastTime > 100) {{
                updateFrame((currentFrame + 1) % {n_frames});
                lastTime = time;
            }}

            controls.update();
            renderer.render(scene, camera);
        }}

        updateFrame(0);
        animate(0);

        window.onresize = () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }};
    </script>
</body>
</html>'''

    return html


# =============================================================================
# Convenience functions for common scans
# =============================================================================

def generate_ptcda_nacl_zscan(
    z_start: float = 3.0,
    z_end: float = 20.0,
    z_step: float = 0.1,
    surface_layers: int = 4,
    supercell: Tuple[int, int, int] = (4, 4, 1)
) -> Dict[str, Any]:
    """
    Convenience function for PTCDA on NaCl(100) z-scan.

    This is a commonly requested structure scan for studying
    molecule-surface interactions.

    Parameters
    ----------
    z_start : float
        Starting z-distance (default: 3.0 Å)
    z_end : float
        Ending z-distance (default: 20.0 Å)
    z_step : float
        Step size (default: 0.1 Å)
    surface_layers : int
        Number of NaCl layers (default: 4)
    supercell : tuple
        Surface supercell size (default: (4, 4, 1))

    Returns
    -------
    dict
        Complete z-scan result with trajectory
    """
    return generate_z_scan(
        surface_material="NaCl",
        surface_miller=(1, 0, 0),
        molecule="PTCDA",
        z_start=z_start,
        z_end=z_end,
        z_step=z_step,
        site="center",
        surface_layers=surface_layers,
        supercell=supercell,
        output_format="all"
    )


def generate_directional_scan(
    structure: Dict[str, Any],
    direction: List[float],
    distance_start: float = 0.0,
    distance_end: float = 10.0,
    distance_step: float = 0.1,
    move_indices: Optional[List[int]] = None,
    normalize_direction: bool = True,
    output_format: str = "trajectory"
) -> Dict[str, Any]:
    """
    Scan structure along an arbitrary direction vector.

    Translates selected atoms (or all atoms) along any specified direction,
    useful for:
    - Reaction coordinate scans
    - Diffusion path studies
    - Custom approach trajectories

    Parameters
    ----------
    structure : dict
        Input structure
    direction : list
        Direction vector [dx, dy, dz] for translation
    distance_start : float
        Starting distance along direction (default: 0.0)
    distance_end : float
        Ending distance along direction (default: 10.0)
    distance_step : float
        Step size in Angstroms (default: 0.1)
    move_indices : list, optional
        Indices of atoms to move. If None, moves all atoms.
    normalize_direction : bool
        Normalize direction to unit vector (default: True)
    output_format : str
        Output format: "trajectory", "list"

    Returns
    -------
    dict
        Result with structures at each position along direction
    """
    direction = np.array(direction, dtype=float)
    if normalize_direction:
        direction = direction / np.linalg.norm(direction)

    distances = np.arange(distance_start, distance_end + distance_step/2, distance_step)
    n_frames = len(distances)

    structures = []
    trajectory_xyz_lines = []

    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))
    cell = structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    for i, d in enumerate(distances):
        positions = base_positions.copy()
        translation = d * direction

        if move_indices is None:
            positions += translation
        else:
            for idx in move_indices:
                positions[idx] += translation

        struct_dict = {
            "species": species,
            "positions": positions.tolist(),
            "cell": cell,
            "n_atoms": len(species),
            "distance_along_direction": float(d),
            "frame_index": i
        }
        structures.append(struct_dict)

        # XYZ frame
        n_atoms = len(species)
        trajectory_xyz_lines.append(f"{n_atoms}\n")
        trajectory_xyz_lines.append(f"distance={d:.3f} frame={i}\n")
        for sym, pos in zip(species, positions):
            trajectory_xyz_lines.append(
                f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n"
            )

    return {
        "success": True,
        "operation": "directional_scan",
        "direction": direction.tolist(),
        "structures": structures,
        "distances": distances.tolist(),
        "n_frames": n_frames,
        "trajectory_xyz": "".join(trajectory_xyz_lines)
    }


def generate_planar_scan(
    structure: Dict[str, Any],
    plane_normal: List[float] = [0, 0, 1],
    x_range: Tuple[float, float, float] = (-5.0, 5.0, 0.5),
    y_range: Tuple[float, float, float] = (-5.0, 5.0, 0.5),
    move_indices: Optional[List[int]] = None,
    output_format: str = "grid"
) -> Dict[str, Any]:
    """
    2D planar scan of structure positions.

    Generates a grid of structures with positions varied across a plane,
    useful for:
    - Potential energy surface (PES) mapping
    - Lateral position optimization
    - Diffusion barrier calculations

    Parameters
    ----------
    structure : dict
        Input structure
    plane_normal : list
        Normal vector to the scan plane (default: [0, 0, 1] for XY plane)
    x_range : tuple
        (start, end, step) for first in-plane direction
    y_range : tuple
        (start, end, step) for second in-plane direction
    move_indices : list, optional
        Indices of atoms to move. If None, moves all atoms.
    output_format : str
        Output format: "grid", "list"

    Returns
    -------
    dict
        Result with 2D grid of structures
    """
    # Build orthonormal basis for the plane
    normal = np.array(plane_normal, dtype=float)
    normal = normal / np.linalg.norm(normal)

    # Find two perpendicular vectors in the plane
    if abs(normal[2]) < 0.9:
        u = np.cross(normal, [0, 0, 1])
    else:
        u = np.cross(normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)

    x_values = np.arange(x_range[0], x_range[1] + x_range[2]/2, x_range[2])
    y_values = np.arange(y_range[0], y_range[1] + y_range[2]/2, y_range[2])

    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))
    cell = structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    structures = []
    grid_shape = (len(x_values), len(y_values))
    frame_idx = 0

    for ix, x in enumerate(x_values):
        for iy, y in enumerate(y_values):
            positions = base_positions.copy()
            translation = x * u + y * v

            if move_indices is None:
                positions += translation
            else:
                for idx in move_indices:
                    positions[idx] += translation

            struct_dict = {
                "species": species,
                "positions": positions.tolist(),
                "cell": cell,
                "n_atoms": len(species),
                "x_position": float(x),
                "y_position": float(y),
                "grid_index": (ix, iy),
                "frame_index": frame_idx
            }
            structures.append(struct_dict)
            frame_idx += 1

    return {
        "success": True,
        "operation": "planar_scan",
        "plane_normal": normal.tolist(),
        "x_range": list(x_range),
        "y_range": list(y_range),
        "grid_shape": grid_shape,
        "n_structures": len(structures),
        "structures": structures
    }


def generate_radial_scan(
    structure: Dict[str, Any],
    center: List[float],
    radius_start: float = 0.0,
    radius_end: float = 10.0,
    radius_step: float = 0.2,
    n_angles: int = 12,
    plane_normal: List[float] = [0, 0, 1],
    move_indices: Optional[List[int]] = None
) -> Dict[str, Any]:
    """
    Radial scan around a center point.

    Generates structures at different radii and angles from a center,
    useful for:
    - Adsorption site mapping
    - Radial distribution studies
    - Polar coordinate PES

    Parameters
    ----------
    structure : dict
        Input structure
    center : list
        Center point [x, y, z] for radial scan
    radius_start : float
        Starting radius (default: 0.0)
    radius_end : float
        Ending radius (default: 10.0)
    radius_step : float
        Radial step size (default: 0.2)
    n_angles : int
        Number of angular divisions (default: 12 for 30° steps)
    plane_normal : list
        Normal to the plane of rotation (default: [0, 0, 1])
    move_indices : list, optional
        Indices of atoms to move

    Returns
    -------
    dict
        Result with radial grid of structures
    """
    from scipy.spatial.transform import Rotation as R

    center = np.array(center)
    normal = np.array(plane_normal, dtype=float)
    normal = normal / np.linalg.norm(normal)

    radii = np.arange(radius_start, radius_end + radius_step/2, radius_step)
    angles = np.linspace(0, 360, n_angles, endpoint=False)

    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))
    cell = structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    # Find initial radial direction
    if abs(normal[2]) < 0.9:
        radial_dir = np.cross(normal, [0, 0, 1])
    else:
        radial_dir = np.cross(normal, [1, 0, 0])
    radial_dir = radial_dir / np.linalg.norm(radial_dir)

    structures = []
    frame_idx = 0

    for radius in radii:
        for angle in angles:
            # Rotate radial direction around normal
            rot = R.from_rotvec(np.radians(angle) * normal)
            direction = rot.apply(radial_dir)

            positions = base_positions.copy()
            translation = center + radius * direction - base_positions.mean(axis=0)

            if move_indices is None:
                positions += translation
            else:
                for idx in move_indices:
                    positions[idx] += translation

            struct_dict = {
                "species": species,
                "positions": positions.tolist(),
                "cell": cell,
                "n_atoms": len(species),
                "radius": float(radius),
                "angle_deg": float(angle),
                "frame_index": frame_idx
            }
            structures.append(struct_dict)
            frame_idx += 1

    return {
        "success": True,
        "operation": "radial_scan",
        "center": center.tolist(),
        "radii": radii.tolist(),
        "angles": angles.tolist(),
        "n_structures": len(structures),
        "structures": structures
    }


def generate_lattice_scan(
    structure: Dict[str, Any],
    a_range: Optional[Tuple[float, float, float]] = None,
    b_range: Optional[Tuple[float, float, float]] = None,
    c_range: Optional[Tuple[float, float, float]] = None,
    alpha_range: Optional[Tuple[float, float, float]] = None,
    beta_range: Optional[Tuple[float, float, float]] = None,
    gamma_range: Optional[Tuple[float, float, float]] = None,
    volume_range: Optional[Tuple[float, float, float]] = None,
    scale_positions: bool = True
) -> Dict[str, Any]:
    """
    Scan over lattice parameters.

    Generates structures with varying cell parameters, useful for:
    - Equation of state calculations
    - Phase transition studies
    - Lattice relaxation

    Parameters
    ----------
    structure : dict
        Input structure
    a_range, b_range, c_range : tuple, optional
        (start, end, step) for lattice vectors a, b, c in Angstroms
    alpha_range, beta_range, gamma_range : tuple, optional
        (start, end, step) for angles in degrees
    volume_range : tuple, optional
        (start, end, step) for volume scaling factor
    scale_positions : bool
        Scale atomic positions with cell (default: True)

    Returns
    -------
    dict
        Result with structures at each lattice parameter set
    """
    base_cell = np.array(structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))

    # Convert to cell parameters
    a = np.linalg.norm(base_cell[0])
    b = np.linalg.norm(base_cell[1])
    c = np.linalg.norm(base_cell[2])

    structures = []
    frame_idx = 0

    # Simple case: volume scaling
    if volume_range is not None:
        scales = np.arange(volume_range[0], volume_range[1] + volume_range[2]/2, volume_range[2])
        for scale in scales:
            # Volume scales as cube of linear scaling
            linear_scale = scale ** (1/3)
            new_cell = base_cell * linear_scale

            if scale_positions:
                new_positions = base_positions * linear_scale
            else:
                new_positions = base_positions.copy()

            struct_dict = {
                "species": species,
                "positions": new_positions.tolist(),
                "cell": new_cell.tolist(),
                "n_atoms": len(species),
                "volume_scale": float(scale),
                "frame_index": frame_idx
            }
            structures.append(struct_dict)
            frame_idx += 1

    # Lattice parameter a scan
    elif a_range is not None:
        a_values = np.arange(a_range[0], a_range[1] + a_range[2]/2, a_range[2])
        for new_a in a_values:
            scale = new_a / a
            new_cell = base_cell.copy()
            new_cell[0] = base_cell[0] * scale

            if scale_positions:
                # Scale positions in a direction only
                inv_cell = np.linalg.inv(base_cell)
                frac = base_positions @ inv_cell
                new_positions = frac @ new_cell
            else:
                new_positions = base_positions.copy()

            struct_dict = {
                "species": species,
                "positions": new_positions.tolist(),
                "cell": new_cell.tolist(),
                "n_atoms": len(species),
                "a": float(new_a),
                "frame_index": frame_idx
            }
            structures.append(struct_dict)
            frame_idx += 1

    return {
        "success": True,
        "operation": "lattice_scan",
        "n_structures": len(structures),
        "structures": structures
    }


def generate_strain_tensor_scan(
    structure: Dict[str, Any],
    strain_component: str = "xx",
    strain_range: Tuple[float, float, float] = (-0.05, 0.05, 0.01),
    strain_type: str = "engineering"
) -> Dict[str, Any]:
    """
    Scan over strain tensor components.

    Generates structures with systematically varied strain, useful for:
    - Elastic constant calculations
    - Piezoelectric response
    - Stress-strain curves

    Parameters
    ----------
    structure : dict
        Input structure
    strain_component : str
        Strain component to vary: "xx", "yy", "zz", "xy", "xz", "yz",
        "biaxial", "uniaxial_x", "uniaxial_y", "uniaxial_z", "hydrostatic"
    strain_range : tuple
        (start, end, step) for strain values (e.g., -0.05 to 0.05 for ±5%)
    strain_type : str
        "engineering" or "true" strain

    Returns
    -------
    dict
        Result with strained structures
    """
    strains = np.arange(strain_range[0], strain_range[1] + strain_range[2]/2, strain_range[2])

    base_cell = np.array(structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))
    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))

    structures = []
    frame_idx = 0

    for eps in strains:
        # Build strain tensor
        if strain_component == "xx":
            F = np.array([[1+eps, 0, 0], [0, 1, 0], [0, 0, 1]])
        elif strain_component == "yy":
            F = np.array([[1, 0, 0], [0, 1+eps, 0], [0, 0, 1]])
        elif strain_component == "zz":
            F = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1+eps]])
        elif strain_component == "xy" or strain_component == "yx":
            F = np.array([[1, eps/2, 0], [eps/2, 1, 0], [0, 0, 1]])
        elif strain_component == "xz" or strain_component == "zx":
            F = np.array([[1, 0, eps/2], [0, 1, 0], [eps/2, 0, 1]])
        elif strain_component == "yz" or strain_component == "zy":
            F = np.array([[1, 0, 0], [0, 1, eps/2], [0, eps/2, 1]])
        elif strain_component == "biaxial":
            F = np.array([[1+eps, 0, 0], [0, 1+eps, 0], [0, 0, 1]])
        elif strain_component == "uniaxial_x":
            F = np.array([[1+eps, 0, 0], [0, 1, 0], [0, 0, 1]])
        elif strain_component == "uniaxial_y":
            F = np.array([[1, 0, 0], [0, 1+eps, 0], [0, 0, 1]])
        elif strain_component == "uniaxial_z":
            F = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1+eps]])
        elif strain_component == "hydrostatic":
            scale = 1 + eps
            F = np.array([[scale, 0, 0], [0, scale, 0], [0, 0, scale]])
        else:
            F = np.eye(3)

        # Apply deformation
        new_cell = (F @ base_cell.T).T
        new_positions = (F @ base_positions.T).T

        struct_dict = {
            "species": species,
            "positions": new_positions.tolist(),
            "cell": new_cell.tolist(),
            "n_atoms": len(species),
            "strain": float(eps),
            "strain_component": strain_component,
            "frame_index": frame_idx
        }
        structures.append(struct_dict)
        frame_idx += 1

    return {
        "success": True,
        "operation": "strain_tensor_scan",
        "strain_component": strain_component,
        "strain_values": strains.tolist(),
        "n_structures": len(structures),
        "structures": structures
    }


def generate_coverage_scan(
    surface_structure: Dict[str, Any],
    adsorbate: Union[str, Dict[str, Any]],
    coverage_range: Tuple[float, float, float] = (0.1, 1.0, 0.1),
    site_type: str = "ontop",
    height: float = 2.0,
    random_placement: bool = False,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Scan over adsorbate coverage.

    Generates structures with varying numbers of adsorbates, useful for:
    - Coverage-dependent adsorption studies
    - Catalytic activity vs coverage
    - Saturation behavior

    Parameters
    ----------
    surface_structure : dict
        Surface structure
    adsorbate : str or dict
        Adsorbate molecule
    coverage_range : tuple
        (start, end, step) for fractional coverage (0 to 1)
    site_type : str
        Adsorption site type
    height : float
        Height above surface
    random_placement : bool
        Randomly place adsorbates vs ordered
    seed : int, optional
        Random seed

    Returns
    -------
    dict
        Result with structures at each coverage
    """
    if seed is not None:
        np.random.seed(seed)

    coverages = np.arange(coverage_range[0], coverage_range[1] + coverage_range[2]/2, coverage_range[2])

    surface_species = surface_structure.get("species", [])
    surface_positions = np.array(surface_structure.get("positions", []))
    cell = np.array(surface_structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]]))

    # Find surface sites
    z_max = surface_positions[:, 2].max()
    top_mask = surface_positions[:, 2] > z_max - 0.5
    top_positions = surface_positions[top_mask]
    n_sites = len(top_positions)

    # Get adsorbate
    if isinstance(adsorbate, str):
        if not adsorbate.strip():
            return {
                "success": False,
                "error": "Adsorbate name cannot be empty"
            }
        mol = _ase_molecule_or_none(adsorbate)
        if mol is None:
            from ase import Atoms
            mol = Atoms(adsorbate)
        mol_species = mol.get_chemical_symbols()
        mol_positions = mol.positions - mol.positions.mean(axis=0)
    else:
        mol_species = adsorbate.get("species", [])
        mol_positions = np.array(adsorbate.get("positions", [])) - np.mean(adsorbate.get("positions", []), axis=0)

    structures = []
    frame_idx = 0

    for coverage in coverages:
        n_adsorbates = max(1, int(n_sites * coverage))

        if random_placement:
            site_indices = np.random.choice(n_sites, n_adsorbates, replace=False)
        else:
            # Evenly spaced
            site_indices = np.linspace(0, n_sites - 1, n_adsorbates, dtype=int)

        all_species = list(surface_species)
        all_positions = surface_positions.tolist()

        for site_idx in site_indices:
            site_pos = top_positions[site_idx]
            for ms, mp in zip(mol_species, mol_positions):
                all_species.append(ms)
                all_positions.append([
                    site_pos[0] + mp[0],
                    site_pos[1] + mp[1],
                    z_max + height + mp[2]
                ])

        struct_dict = {
            "species": all_species,
            "positions": all_positions,
            "cell": cell.tolist(),
            "n_atoms": len(all_species),
            "coverage": float(coverage),
            "n_adsorbates": n_adsorbates,
            "frame_index": frame_idx
        }
        structures.append(struct_dict)
        frame_idx += 1

    return {
        "success": True,
        "operation": "coverage_scan",
        "coverages": coverages.tolist(),
        "n_structures": len(structures),
        "structures": structures
    }


def generate_neb_path(
    initial_structure: Dict[str, Any],
    final_structure: Dict[str, Any],
    n_images: int = 9,
    interpolation: str = "linear",
    climbing_image: bool = False
) -> Dict[str, Any]:
    """
    Generate initial path for Nudged Elastic Band (NEB) calculation.

    Creates interpolated structures between initial and final states,
    useful for:
    - Transition state searches
    - Reaction path studies
    - Diffusion barrier calculations

    Parameters
    ----------
    initial_structure : dict
        Initial (reactant) structure
    final_structure : dict
        Final (product) structure
    n_images : int
        Number of intermediate images (default: 9)
    interpolation : str
        Interpolation method: "linear", "idpp" (requires ASE)
    climbing_image : bool
        Mark middle image as climbing image

    Returns
    -------
    dict
        Result with NEB path structures
    """
    pos_initial = np.array(initial_structure.get("positions", []))
    pos_final = np.array(final_structure.get("positions", []))
    species = initial_structure.get("species", [])
    cell = initial_structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    if len(pos_initial) != len(pos_final):
        return {
            "success": False,
            "error": "Initial and final structures must have same number of atoms"
        }

    structures = []

    # Linear interpolation
    for i in range(n_images + 2):  # +2 for endpoints
        t = i / (n_images + 1)
        positions = (1 - t) * pos_initial + t * pos_final

        is_climbing = climbing_image and i == (n_images + 2) // 2

        struct_dict = {
            "species": species,
            "positions": positions.tolist(),
            "cell": cell,
            "n_atoms": len(species),
            "image_index": i,
            "interpolation_parameter": float(t),
            "is_endpoint": i == 0 or i == n_images + 1,
            "is_climbing_image": is_climbing,
            "frame_index": i
        }
        structures.append(struct_dict)

    return {
        "success": True,
        "operation": "neb_path",
        "n_images": n_images,
        "interpolation": interpolation,
        "climbing_image": climbing_image,
        "n_structures": len(structures),
        "structures": structures,
        "neb_info": {
            "description": "Initial NEB path ready for optimization",
            "next_steps": [
                "1. Run NEB optimization with DFT/ML force field",
                "2. Analyze transition state (maximum energy image)",
                "3. Compute activation barrier"
            ]
        }
    }


def generate_bond_scan(
    structure: Dict[str, Any],
    atom_indices: Tuple[int, int],
    bond_range: Tuple[float, float, float] = (0.8, 3.0, 0.05),
    fix_center_of_mass: bool = True
) -> Dict[str, Any]:
    """
    Scan bond distance between two atoms.

    Generates structures with varying bond length, useful for:
    - Bond dissociation curves
    - Potential energy curves
    - Force constant extraction

    Parameters
    ----------
    structure : dict
        Input structure
    atom_indices : tuple
        Indices of two atoms (i, j) for bond scan
    bond_range : tuple
        (start, end, step) for bond distance in Angstroms
    fix_center_of_mass : bool
        Keep center of mass fixed

    Returns
    -------
    dict
        Result with structures at each bond distance
    """
    i, j = atom_indices
    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))
    cell = structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    # Current bond vector
    bond_vec = base_positions[j] - base_positions[i]
    current_length = np.linalg.norm(bond_vec)
    bond_dir = bond_vec / current_length

    distances = np.arange(bond_range[0], bond_range[1] + bond_range[2]/2, bond_range[2])

    structures = []
    frame_idx = 0

    for d in distances:
        positions = base_positions.copy()

        # Scale factor
        scale = d / current_length
        new_bond = bond_dir * d

        if fix_center_of_mass:
            # Move both atoms symmetrically
            midpoint = (positions[i] + positions[j]) / 2
            positions[i] = midpoint - new_bond / 2
            positions[j] = midpoint + new_bond / 2
        else:
            # Move only atom j
            positions[j] = positions[i] + new_bond

        struct_dict = {
            "species": species,
            "positions": positions.tolist(),
            "cell": cell,
            "n_atoms": len(species),
            "bond_distance": float(d),
            "atom_pair": list(atom_indices),
            "frame_index": frame_idx
        }
        structures.append(struct_dict)
        frame_idx += 1

    return {
        "success": True,
        "operation": "bond_scan",
        "atom_indices": list(atom_indices),
        "bond_distances": distances.tolist(),
        "n_structures": len(structures),
        "structures": structures
    }


def generate_dihedral_scan(
    structure: Dict[str, Any],
    atom_indices: Tuple[int, int, int, int],
    angle_range: Tuple[float, float, float] = (0.0, 360.0, 10.0)
) -> Dict[str, Any]:
    """
    Scan dihedral angle between four atoms.

    Generates structures with varying dihedral angle, useful for:
    - Torsional barrier calculations
    - Conformational analysis
    - Rotational isomerism studies

    Parameters
    ----------
    structure : dict
        Input structure
    atom_indices : tuple
        Indices of four atoms (i, j, k, l) defining dihedral
    angle_range : tuple
        (start, end, step) for dihedral angle in degrees

    Returns
    -------
    dict
        Result with structures at each dihedral angle
    """
    from scipy.spatial.transform import Rotation as R

    i, j, k, l = atom_indices
    species = structure.get("species", [])
    base_positions = np.array(structure.get("positions", []))
    cell = structure.get("cell", [[10, 0, 0], [0, 10, 0], [0, 0, 10]])

    # Define rotation axis (bond j-k)
    axis = base_positions[k] - base_positions[j]
    axis = axis / np.linalg.norm(axis)

    # Get current dihedral angle
    def get_dihedral(p):
        b1 = p[j] - p[i]
        b2 = p[k] - p[j]
        b3 = p[l] - p[k]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        return np.degrees(np.arctan2(y, x))

    current_angle = get_dihedral(base_positions)
    angles = np.arange(angle_range[0], angle_range[1] + angle_range[2]/2, angle_range[2])

    structures = []
    frame_idx = 0

    for target_angle in angles:
        positions = base_positions.copy()

        # Rotate atoms after k around j-k axis
        rotation_angle = target_angle - current_angle
        rot = R.from_rotvec(np.radians(rotation_angle) * axis)

        # Find atoms connected to l side (simplified: just rotate l and beyond)
        pivot = positions[k]
        for idx in range(l, len(positions)):  # Simplified: rotate all after k
            positions[idx] = rot.apply(positions[idx] - pivot) + pivot

        struct_dict = {
            "species": species,
            "positions": positions.tolist(),
            "cell": cell,
            "n_atoms": len(species),
            "dihedral_angle": float(target_angle),
            "atom_indices": list(atom_indices),
            "frame_index": frame_idx
        }
        structures.append(struct_dict)
        frame_idx += 1

    return {
        "success": True,
        "operation": "dihedral_scan",
        "atom_indices": list(atom_indices),
        "angles": angles.tolist(),
        "n_structures": len(structures),
        "structures": structures
    }


def get_available_scans() -> Dict[str, Any]:
    """Get list of available parametric scan operations."""
    return {
        "z_scan": {
            "function": "generate_z_scan",
            "description": "Scan molecule height above surface",
            "typical_use": "Binding energy curves, adsorption studies"
        },
        "distance_scan": {
            "function": "generate_distance_scan",
            "description": "Scan separation between two structures",
            "typical_use": "Dimer interactions, layer separation"
        },
        "rotation_scan": {
            "function": "generate_rotation_scan",
            "description": "Rotate structure around axis",
            "typical_use": "Rotational barriers, orientation studies"
        },
        "directional_scan": {
            "function": "generate_directional_scan",
            "description": "Scan along arbitrary direction vector",
            "typical_use": "Reaction coordinates, custom trajectories"
        },
        "planar_scan": {
            "function": "generate_planar_scan",
            "description": "2D grid scan on a plane",
            "typical_use": "PES mapping, lateral position optimization"
        },
        "radial_scan": {
            "function": "generate_radial_scan",
            "description": "Radial scan around a center point",
            "typical_use": "Adsorption site mapping, radial PES"
        },
        "lattice_scan": {
            "function": "generate_lattice_scan",
            "description": "Scan lattice parameters (a, b, c, angles, volume)",
            "typical_use": "Equation of state, phase transitions"
        },
        "strain_tensor_scan": {
            "function": "generate_strain_tensor_scan",
            "description": "Scan strain tensor components",
            "typical_use": "Elastic constants, stress-strain curves"
        },
        "coverage_scan": {
            "function": "generate_coverage_scan",
            "description": "Scan adsorbate coverage on surface",
            "typical_use": "Coverage-dependent properties"
        },
        "neb_path": {
            "function": "generate_neb_path",
            "description": "Generate NEB interpolation path",
            "typical_use": "Transition state searches, diffusion barriers"
        },
        "bond_scan": {
            "function": "generate_bond_scan",
            "description": "Scan bond distance between atoms",
            "typical_use": "Bond dissociation, potential curves"
        },
        "dihedral_scan": {
            "function": "generate_dihedral_scan",
            "description": "Scan dihedral angle",
            "typical_use": "Torsional barriers, conformational analysis"
        },
        "parameter_sweep": {
            "function": "generate_parameter_sweep",
            "description": "Sweep single parameter for any generator",
            "typical_use": "Systematic parameter exploration"
        },
        "multi_parameter_scan": {
            "function": "generate_multi_parameter_scan",
            "description": "Grid or path scan over multiple parameters",
            "typical_use": "Multi-dimensional phase space exploration"
        },
        "trajectory_animation": {
            "function": "generate_trajectory_animation",
            "description": "Create animation from structure list",
            "typical_use": "Visualization of scan results"
        }
    }
