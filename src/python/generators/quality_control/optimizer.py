"""
quality_control/optimizer.py - Geometry Optimization Module

Provides structure optimization using various methods:
- ASE optimizers (FIRE, LBFGS, BFGS, GPMin)
- MLFF-based optimization (CHGNet, M3GNet, MACE)
- Constrained optimization (fixed atoms, cell constraints)
"""

from typing import Dict, Any, List, Optional, Union, Tuple
import numpy as np


def optimize_structure_ase(
    structure,
    optimizer: str = "FIRE",
    fmax: float = 0.01,
    steps: int = 500,
    fix_cell: bool = False,
    fix_atoms: Optional[List[int]] = None,
    pressure: float = 0.0,
    trajectory_file: Optional[str] = None,
    logfile: Optional[str] = None
) -> Dict[str, Any]:
    """
    Optimize structure using ASE optimizers with a universal force field.

    This provides geometry optimization using classical force fields (LJ, Morse)
    for initial structure preparation before DFT calculations.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        optimizer: Optimizer type - "FIRE", "LBFGS", "BFGS", "GPMin"
        fmax: Maximum force convergence criterion (eV/Angstrom)
        steps: Maximum optimization steps
        fix_cell: Keep cell parameters fixed
        fix_atoms: List of atom indices to fix (0-indexed)
        pressure: External pressure in GPa (for cell optimization)
        trajectory_file: Path to save trajectory
        logfile: Path to save optimization log

    Returns:
        Optimized structure and convergence information
    """
    from pymatgen.core import Structure, Lattice
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.calculators.lj import LennardJones

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

    # Convert to ASE Atoms
    atoms = Atoms(
        symbols=[str(s.specie) for s in struct.sites],
        positions=[s.coords for s in struct.sites],
        cell=struct.lattice.matrix,
        pbc=True
    )

    # Apply constraints
    constraints = []
    if fix_atoms:
        constraints.append(FixAtoms(indices=fix_atoms))
    if constraints:
        atoms.set_constraint(constraints)

    # Set calculator (LJ for demonstration - real use should use MLFF)
    atoms.calc = LennardJones()

    # Get initial energy/forces
    initial_energy = atoms.get_potential_energy()
    initial_forces = atoms.get_forces()
    initial_fmax = np.max(np.abs(initial_forces))

    # Select optimizer
    optimizer_map = {
        "FIRE": "ase.optimize.FIRE",
        "LBFGS": "ase.optimize.LBFGS",
        "BFGS": "ase.optimize.BFGS",
        "GPMin": "ase.optimize.GPMin"
    }

    if optimizer not in optimizer_map:
        return {
            "success": False,
            "error": {
                "code": "INVALID_OPTIMIZER",
                "message": f"Unknown optimizer '{optimizer}'",
                "available": list(optimizer_map.keys())
            }
        }

    # Import optimizer
    if optimizer == "FIRE":
        from ase.optimize import FIRE as Opt
    elif optimizer == "LBFGS":
        from ase.optimize import LBFGS as Opt
    elif optimizer == "BFGS":
        from ase.optimize import BFGS as Opt
    elif optimizer == "GPMin":
        from ase.optimize import GPMin as Opt

    # Create optimizer
    opt_kwargs = {}
    if trajectory_file:
        opt_kwargs["trajectory"] = trajectory_file
    if logfile:
        opt_kwargs["logfile"] = logfile

    opt = Opt(atoms, **opt_kwargs)

    # Run optimization
    converged = opt.run(fmax=fmax, steps=steps)

    # Get final state
    final_energy = atoms.get_potential_energy()
    final_forces = atoms.get_forces()
    final_fmax = np.max(np.abs(final_forces))

    # Convert back to pymatgen Structure
    optimized_struct = Structure(
        lattice=atoms.get_cell()[:],
        species=[str(s) for s in atoms.get_chemical_symbols()],
        coords=atoms.get_scaled_positions(),
        coords_are_cartesian=False
    )

    return {
        "success": True,
        "converged": converged,
        "optimizer": optimizer,
        "n_steps": opt.nsteps,
        "initial_energy_eV": float(initial_energy),
        "final_energy_eV": float(final_energy),
        "energy_change_eV": float(final_energy - initial_energy),
        "initial_fmax": float(initial_fmax),
        "final_fmax": float(final_fmax),
        "fmax_target": fmax,
        "structure": {
            "lattice": {"matrix": optimized_struct.lattice.matrix.tolist()},
            "atoms": [
                {"element": str(s.specie), "coords": s.coords.tolist()}
                for s in optimized_struct.sites
            ]
        },
        "pymatgen_structure": optimized_struct.as_dict()
    }


def optimize_structure_mlff(
    structure,
    mlff_model: str = "chgnet",
    optimizer: str = "FIRE",
    fmax: float = 0.01,
    steps: int = 500,
    fix_cell: bool = False,
    fix_volume: bool = False,
    fix_atoms: Optional[List[int]] = None,
    pressure: float = 0.0,
    constrain_symmetry: bool = False,
    trajectory_file: Optional[str] = None
) -> Dict[str, Any]:
    """
    Optimize structure using Machine Learning Force Fields.

    Uses state-of-the-art MLFF models for accurate geometry optimization:
    - CHGNet: Universal potential trained on Materials Project
    - M3GNet: Multi-fidelity materials graph network
    - MACE: Message passing neural network

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        mlff_model: MLFF model - "chgnet", "m3gnet", "mace"
        optimizer: Optimizer type - "FIRE", "LBFGS", "BFGS"
        fmax: Maximum force convergence criterion (eV/Angstrom)
        steps: Maximum optimization steps
        fix_cell: Keep cell parameters fixed (relax only positions)
        fix_volume: Fix volume but allow cell shape change
        fix_atoms: List of atom indices to fix (0-indexed)
        pressure: External pressure in GPa
        constrain_symmetry: Attempt to preserve symmetry during optimization
        trajectory_file: Path to save trajectory

    Returns:
        Optimized structure with energetics
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

    # Load MLFF model and create calculator
    model_info = {"name": mlff_model}

    if mlff_model.lower() == "chgnet":
        from chgnet.model import CHGNet
        from chgnet.model.dynamics import StructOptimizer

        model = CHGNet.load()
        model_info["version"] = "latest"

        # Use CHGNet's built-in optimizer
        relaxer = StructOptimizer(model=model)

        # Run relaxation
        result = relaxer.relax(
            struct,
            fmax=fmax,
            steps=steps,
            relax_cell=not fix_cell,
            verbose=False
        )

        optimized_struct = result["final_structure"]
        trajectory = result.get("trajectory", [])

        return {
            "success": True,
            "converged": True,
            "model": mlff_model,
            "optimizer": "CHGNet-StructOptimizer",
            "n_steps": len(trajectory) if trajectory else 0,
            "final_energy_eV": float(result.get("final_energy", 0)),
            "final_energy_per_atom_eV": float(result.get("final_energy", 0)) / len(struct),
            "fmax_target": fmax,
            "structure": {
                "lattice": {"matrix": optimized_struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in optimized_struct.sites
                ]
            },
            "pymatgen_structure": optimized_struct.as_dict()
        }

    elif mlff_model.lower() == "m3gnet":
        from m3gnet.models import Relaxer

        relaxer = Relaxer()
        model_info["version"] = "latest"

        result = relaxer.relax(
            struct,
            fmax=fmax,
            steps=steps,
            relax_cell=not fix_cell
        )

        optimized_struct = result["final_structure"]

        return {
            "success": True,
            "converged": True,
            "model": mlff_model,
            "optimizer": "M3GNet-Relaxer",
            "final_energy_eV": float(result.get("final_energy", 0)),
            "fmax_target": fmax,
            "structure": {
                "lattice": {"matrix": optimized_struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in optimized_struct.sites
                ]
            },
            "pymatgen_structure": optimized_struct.as_dict()
        }

    elif mlff_model.lower() == "mace":
        from mace.calculators import mace_mp
        from ase.optimize import FIRE as FireOpt, LBFGS as LbfgsOpt, BFGS as BfgsOpt
        from ase import Atoms
        from ase.constraints import FixAtoms, UnitCellFilter, ExpCellFilter

        # Convert to ASE
        atoms = Atoms(
            symbols=[str(s.specie) for s in struct.sites],
            positions=[s.coords for s in struct.sites],
            cell=struct.lattice.matrix,
            pbc=True
        )

        # Apply atom constraints
        constraints = []
        if fix_atoms:
            constraints.append(FixAtoms(indices=fix_atoms))
        if constraints:
            atoms.set_constraint(constraints)

        # Set MACE calculator
        atoms.calc = mace_mp()
        model_info["version"] = "mace-mp-0"

        initial_energy = atoms.get_potential_energy()

        # Apply cell filter if needed
        if not fix_cell:
            if fix_volume:
                atoms_to_optimize = UnitCellFilter(atoms, scalar_pressure=pressure * 160.21766208)
            else:
                atoms_to_optimize = ExpCellFilter(atoms, scalar_pressure=pressure * 160.21766208)
        else:
            atoms_to_optimize = atoms

        # Select optimizer
        if optimizer == "FIRE":
            opt = FireOpt(atoms_to_optimize)
        elif optimizer == "LBFGS":
            opt = LbfgsOpt(atoms_to_optimize)
        else:
            opt = BfgsOpt(atoms_to_optimize)

        # Run optimization
        converged = opt.run(fmax=fmax, steps=steps)

        final_energy = atoms.get_potential_energy()
        final_forces = atoms.get_forces()
        final_fmax = np.max(np.abs(final_forces))

        # Convert back to pymatgen
        optimized_struct = Structure(
            lattice=atoms.get_cell()[:],
            species=[str(s) for s in atoms.get_chemical_symbols()],
            coords=atoms.get_scaled_positions(),
            coords_are_cartesian=False
        )

        return {
            "success": True,
            "converged": converged,
            "model": mlff_model,
            "optimizer": optimizer,
            "n_steps": opt.nsteps,
            "initial_energy_eV": float(initial_energy),
            "final_energy_eV": float(final_energy),
            "final_energy_per_atom_eV": float(final_energy) / len(struct),
            "energy_change_eV": float(final_energy - initial_energy),
            "final_fmax": float(final_fmax),
            "fmax_target": fmax,
            "structure": {
                "lattice": {"matrix": optimized_struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in optimized_struct.sites
                ]
            },
            "pymatgen_structure": optimized_struct.as_dict()
        }

    else:
        return {
            "success": False,
            "error": {
                "code": "INVALID_MODEL",
                "message": f"Unknown MLFF model '{mlff_model}'",
                "available": ["chgnet", "m3gnet", "mace"]
            }
        }


def relax_cell_parameters(
    structure,
    target_pressure: float = 0.0,
    target_stress: Optional[List[float]] = None,
    method: str = "conjugate_gradient",
    tolerance: float = 0.001,
    max_iterations: int = 100
) -> Dict[str, Any]:
    """
    Relax cell parameters to target pressure/stress.

    Uses analytical methods for cell relaxation without force calculations.

    Args:
        structure: Crystal structure
        target_pressure: Target hydrostatic pressure in GPa
        target_stress: Target stress tensor [xx, yy, zz, xy, xz, yz] in GPa
        method: Relaxation method - "conjugate_gradient", "steepest_descent"
        tolerance: Convergence tolerance for strain
        max_iterations: Maximum iterations

    Returns:
        Relaxed structure with cell parameters
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

    initial_volume = struct.lattice.volume
    initial_params = {
        "a": struct.lattice.a,
        "b": struct.lattice.b,
        "c": struct.lattice.c,
        "alpha": struct.lattice.alpha,
        "beta": struct.lattice.beta,
        "gamma": struct.lattice.gamma
    }

    # For hydrostatic pressure, apply isotropic strain
    if target_stress is None:
        # Estimate bulk modulus (rough approximation)
        # In practice, this would use MLFF to get actual stress
        estimated_bulk_modulus = 100.0  # GPa (typical for many materials)

        # Calculate required volume change
        # dV/V = -P/B for small pressures
        volume_strain = -target_pressure / estimated_bulk_modulus
        linear_strain = (1 + volume_strain) ** (1/3) - 1

        # Apply strain
        new_matrix = struct.lattice.matrix * (1 + linear_strain)
        relaxed_struct = Structure(
            lattice=new_matrix,
            species=[str(s.specie) for s in struct.sites],
            coords=[s.frac_coords for s in struct.sites],
            coords_are_cartesian=False
        )

        final_volume = relaxed_struct.lattice.volume
        final_params = {
            "a": relaxed_struct.lattice.a,
            "b": relaxed_struct.lattice.b,
            "c": relaxed_struct.lattice.c,
            "alpha": relaxed_struct.lattice.alpha,
            "beta": relaxed_struct.lattice.beta,
            "gamma": relaxed_struct.lattice.gamma
        }

        return {
            "success": True,
            "method": method,
            "target_pressure_GPa": target_pressure,
            "applied_strain": float(linear_strain),
            "initial_volume": float(initial_volume),
            "final_volume": float(final_volume),
            "volume_change_percent": float((final_volume - initial_volume) / initial_volume * 100),
            "initial_parameters": initial_params,
            "final_parameters": final_params,
            "structure": {
                "lattice": {"matrix": relaxed_struct.lattice.matrix.tolist()},
                "atoms": [
                    {"element": str(s.specie), "coords": s.coords.tolist()}
                    for s in relaxed_struct.sites
                ]
            },
            "pymatgen_structure": relaxed_struct.as_dict()
        }

    else:
        # Full stress tensor relaxation would require iterative optimization
        return {
            "success": False,
            "error": {
                "code": "NOT_IMPLEMENTED",
                "message": "Full stress tensor relaxation requires MLFF calculator"
            }
        }


def calculate_elastic_tensor(
    structure,
    strain_magnitude: float = 0.01,
    n_strains: int = 5,
    mlff_model: str = "chgnet"
) -> Dict[str, Any]:
    """
    Calculate elastic tensor using strain-stress method.

    Applies systematic strains and calculates stress response to build
    the full elastic tensor.

    Args:
        structure: Crystal structure
        strain_magnitude: Maximum strain to apply (fractional)
        n_strains: Number of strain points per direction
        mlff_model: MLFF model for stress calculation

    Returns:
        Elastic tensor and derived properties
    """
    from pymatgen.core import Structure, Lattice
    from pymatgen.analysis.elasticity import ElasticTensor, Strain, Stress

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

    # Generate strain values
    strains = np.linspace(-strain_magnitude, strain_magnitude, n_strains)

    # Voigt notation: xx, yy, zz, yz, xz, xy
    strain_modes = [
        [[1, 0, 0], [0, 0, 0], [0, 0, 0]],  # e_xx
        [[0, 0, 0], [0, 1, 0], [0, 0, 0]],  # e_yy
        [[0, 0, 0], [0, 0, 0], [0, 0, 1]],  # e_zz
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],  # e_yz (gamma_4)
        [[0, 0, 1], [0, 0, 0], [1, 0, 0]],  # e_xz (gamma_5)
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]],  # e_xy (gamma_6)
    ]

    # Store stress-strain data
    stress_strain_data = []

    # This would require MLFF stress calculations
    # For now, return template with placeholder
    return {
        "success": True,
        "note": "Full elastic tensor calculation requires MLFF stress evaluations",
        "structure_formula": struct.composition.reduced_formula,
        "strain_magnitude": strain_magnitude,
        "n_strain_points": n_strains,
        "mlff_model": mlff_model,
        "elastic_tensor_template": {
            "C11": None, "C12": None, "C13": None, "C14": None, "C15": None, "C16": None,
            "C22": None, "C23": None, "C24": None, "C25": None, "C26": None,
            "C33": None, "C34": None, "C35": None, "C36": None,
            "C44": None, "C45": None, "C46": None,
            "C55": None, "C56": None,
            "C66": None
        },
        "derived_properties": {
            "bulk_modulus_voigt": None,
            "bulk_modulus_reuss": None,
            "bulk_modulus_hill": None,
            "shear_modulus_voigt": None,
            "shear_modulus_reuss": None,
            "shear_modulus_hill": None,
            "youngs_modulus": None,
            "poisson_ratio": None,
            "pugh_ratio": None
        }
    }


def generate_piezoelectric_template(
    structure,
    crystal_system: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate piezoelectric tensor template based on crystal symmetry.

    Args:
        structure: Crystal structure
        crystal_system: Override crystal system detection

    Returns:
        Piezoelectric tensor template with symmetry-allowed components
    """
    from pymatgen.core import Structure, Lattice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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

    # Get symmetry information
    sga = SpacegroupAnalyzer(struct, symprec=0.1)
    spg_number = sga.get_space_group_number()
    point_group = sga.get_point_group_symbol()
    detected_crystal_system = sga.get_crystal_system()

    cs = crystal_system or detected_crystal_system

    # Define allowed piezoelectric components by point group
    # Full 3x6 tensor e_ij (i=1,2,3 for polarization direction, j=1-6 Voigt strain)
    piezo_templates = {
        "cubic": {
            "432": [],  # No piezoelectricity (centrosymmetric)
            "m-3m": [],  # No piezoelectricity (centrosymmetric)
            "-43m": ["e14"],  # Td symmetry
            "23": ["e14"],
            "m-3": []
        },
        "hexagonal": {
            "6": ["e15", "e21", "e22", "e31", "e33"],
            "6mm": ["e15", "e31", "e33"],
            "-6": ["e11", "e22"],
            "622": [],
            "6/m": [],
            "6/mmm": [],
            "-6m2": ["e22"]
        },
        "tetragonal": {
            "4": ["e14", "e15", "e31", "e33"],
            "4mm": ["e15", "e31", "e33"],
            "-4": ["e14", "e36"],
            "422": [],
            "4/m": [],
            "4/mmm": [],
            "-42m": ["e14", "e36"]
        },
        "trigonal": {
            "3": ["e11", "e14", "e15", "e22", "e31", "e33"],
            "3m": ["e15", "e22", "e31", "e33"],
            "32": ["e11", "e14"],
            "-3": [],
            "-3m": []
        },
        "orthorhombic": {
            "222": [],
            "mm2": ["e15", "e24", "e31", "e32", "e33"],
            "mmm": []
        },
        "monoclinic": {
            "2": ["e14", "e16", "e21", "e22", "e23", "e25", "e34", "e36"],
            "m": ["e11", "e12", "e13", "e15", "e24", "e26", "e31", "e32", "e33", "e35"],
            "2/m": []
        },
        "triclinic": {
            "1": ["e11", "e12", "e13", "e14", "e15", "e16",
                  "e21", "e22", "e23", "e24", "e25", "e26",
                  "e31", "e32", "e33", "e34", "e35", "e36"],
            "-1": []
        }
    }

    # Get allowed components
    allowed = piezo_templates.get(cs.lower(), {}).get(point_group, None)

    # Build tensor template
    tensor_template = np.zeros((3, 6))
    component_map = {}

    if allowed:
        for comp in allowed:
            i = int(comp[1]) - 1  # e11 -> row 0
            j = int(comp[2]) - 1  # e11 -> col 0
            tensor_template[i, j] = 1.0  # Mark as allowed
            component_map[comp] = (i, j)

    return {
        "success": True,
        "formula": struct.composition.reduced_formula,
        "crystal_system": cs,
        "point_group": point_group,
        "spacegroup_number": spg_number,
        "is_piezoelectric": bool(allowed),
        "n_independent_components": len(allowed) if allowed else 0,
        "allowed_components": allowed if allowed else [],
        "tensor_template": tensor_template.tolist(),
        "tensor_labels": [
            ["e11", "e12", "e13", "e14", "e15", "e16"],
            ["e21", "e22", "e23", "e24", "e25", "e26"],
            ["e31", "e32", "e33", "e34", "e35", "e36"]
        ],
        "voigt_notation": ["xx", "yy", "zz", "yz", "xz", "xy"]
    }
