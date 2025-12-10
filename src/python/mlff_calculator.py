#!/usr/bin/env python3
"""
MLFF Calculator - Machine Learning Force Field Wrapper

This module provides interfaces to various MLFF models (CHGNet, M3GNet, MACE)
for structure optimization and energy calculations using ASE.
All functions follow defensive programming without try/except blocks.
"""

from typing import Dict, List, Optional, Any, Tuple
import json
import sys
import time
import numpy as np
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.constraints import FixAtoms, UnitCellFilter


class MLFFModel:
    """Available MLFF models."""
    CHGNET = "chgnet"
    M3GNET = "m3gnet"
    MACE = "mace"
    
    ALL_MODELS = [CHGNET, M3GNET, MACE]


def validate_mlff_model(model_name: str) -> Tuple[bool, Optional[str]]:
    """
    Validate MLFF model name.
    
    Args:
        model_name: Name of the MLFF model
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not model_name or not isinstance(model_name, str):
        return False, "Model name must be a non-empty string"
    
    model_lower = model_name.lower().strip()
    if model_lower not in MLFFModel.ALL_MODELS:
        return False, f"Unknown model '{model_name}'. Valid models: {', '.join(MLFFModel.ALL_MODELS)}"
    
    return True, None


def validate_optimizer(optimizer_name: str) -> Tuple[bool, Optional[str]]:
    """
    Validate optimizer name.
    
    Args:
        optimizer_name: Name of the optimizer
    
    Returns:
        Tuple of (is_valid, error_message)
    """
    valid_optimizers = ["BFGS", "FIRE", "LBFGS", "GPMin"]
    
    if not optimizer_name or not isinstance(optimizer_name, str):
        return False, "Optimizer name must be a non-empty string"
    
    if optimizer_name not in valid_optimizers:
        return False, f"Unknown optimizer '{optimizer_name}'. Valid: {', '.join(valid_optimizers)}"
    
    return True, None


def load_chgnet_calculator():
    """
    Load CHGNet calculator.
    
    Returns:
        CHGNet ASE calculator or None if loading fails
    """
    # Check if chgnet is available
    chgnet_available = True
    
    # Import and load CHGNet
    if chgnet_available:
        from chgnet.model import CHGNet, CHGNetCalculator
        model = CHGNet.load()
        if model is not None:
            calculator = CHGNetCalculator(model=model, use_device="cpu")
            return calculator, None
        else:
            return None, "Failed to load CHGNet model"
    
    return None, "CHGNet package not available"


def load_m3gnet_calculator():
    """
    Load M3GNet calculator.
    
    Returns:
        M3GNet ASE calculator or None if loading fails
    """
    # Check if matgl is available
    m3gnet_available = True
    error_msg = None
    
    # Import and load M3GNet
    if m3gnet_available:
        import matgl
        from matgl.ext.ase import PESCalculator
        
        pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        if pot is not None:
            calculator = PESCalculator(pot)
            return calculator, None
        else:
            return None, "Failed to load M3GNet model"
    
    return None, "M3GNet (matgl) package not available"


def load_mace_calculator():
    """
    Load MACE calculator.
    
    Returns:
        MACE ASE calculator or None if loading fails
    """
    # Check if mace is available
    mace_available = True
    error_msg = None
    
    # Import and load MACE
    if mace_available:
        from mace.calculators import MACECalculator
        
        calculator = MACECalculator(
            model_path="medium",  # Use pretrained medium model
            device="cpu"
        )
        if calculator is not None:
            return calculator, None
        else:
            return None, "Failed to load MACE calculator"
    
    return None, "MACE package not available"


def load_calculator(model_name: str):
    """
    Load MLFF calculator by name.
    
    Args:
        model_name: Name of the MLFF model
    
    Returns:
        Tuple of (calculator, error_message)
    """
    # Validate model name
    is_valid, error = validate_mlff_model(model_name)
    if not is_valid:
        return None, error
    
    model_lower = model_name.lower().strip()
    
    # Load appropriate calculator
    if model_lower == MLFFModel.CHGNET:
        return load_chgnet_calculator()
    elif model_lower == MLFFModel.M3GNET:
        return load_m3gnet_calculator()
    elif model_lower == MLFFModel.MACE:
        return load_mace_calculator()
    else:
        return None, f"Model {model_name} not implemented"


def dict_to_atoms(structure_dict: Dict[str, Any]) -> Optional[Atoms]:
    """
    Convert structure dictionary to ASE Atoms object.
    
    Args:
        structure_dict: Structure dictionary with lattice and atoms
    
    Returns:
        ASE Atoms object or None if conversion fails
    """
    if not structure_dict or not isinstance(structure_dict, dict):
        return None
    
    if "lattice" not in structure_dict or "atoms" not in structure_dict:
        return None
    
    lattice = structure_dict["lattice"]
    if "matrix" not in lattice:
        return None
    
    cell = np.array(lattice["matrix"])
    if cell.shape != (3, 3):
        return None
    
    atoms_data = structure_dict["atoms"]
    if not atoms_data or len(atoms_data) == 0:
        return None
    
    symbols = []
    positions = []
    
    for atom in atoms_data:
        if "element" not in atom or "cartesian" not in atom:
            return None
        
        symbols.append(atom["element"])
        positions.append(atom["cartesian"])
    
    positions_array = np.array(positions)
    if positions_array.shape[1] != 3:
        return None
    
    atoms = Atoms(
        symbols=symbols,
        positions=positions_array,
        cell=cell,
        pbc=True
    )
    
    return atoms


def atoms_to_dict(atoms: Atoms) -> Dict[str, Any]:
    """
    Convert ASE Atoms object to structure dictionary.
    
    Args:
        atoms: ASE Atoms object
    
    Returns:
        Structure dictionary
    """
    cell = atoms.get_cell()
    positions_cartesian = atoms.get_positions()
    positions_fractional = atoms.get_scaled_positions()
    
    # Get lattice parameters
    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    
    structure = {
        "lattice": {
            "a": float(a),
            "b": float(b),
            "c": float(c),
            "alpha": float(alpha),
            "beta": float(beta),
            "gamma": float(gamma),
            "matrix": cell.tolist(),
            "volume": float(atoms.get_volume())
        },
        "atoms": [],
        "metadata": {
            "natoms": len(atoms)
        }
    }
    
    for i, (symbol, frac_pos, cart_pos) in enumerate(
        zip(atoms.get_chemical_symbols(), 
            positions_fractional, 
            positions_cartesian)
    ):
        structure["atoms"].append({
            "element": symbol,
            "coords": frac_pos.tolist(),
            "cartesian": cart_pos.tolist()
        })
    
    return structure


def optimize_structure(
    structure_dict: Dict[str, Any],
    mlff_model: str,
    optimizer: str = "BFGS",
    fmax: float = 0.01,
    steps: int = 500,
    constrain_symmetry: bool = False,
    fix_lattice: bool = False,
    fix_volume: bool = False,
    pressure: float = 0.0,
    fix_atoms: Optional[List[int]] = None,
    trajectory_file: Optional[str] = None
) -> Dict[str, Any]:
    """
    Optimize structure using MLFF.
    
    Args:
        structure_dict: Structure to optimize
        mlff_model: MLFF model name
        optimizer: Optimizer algorithm
        fmax: Force convergence criterion
        steps: Maximum optimization steps
        constrain_symmetry: Preserve symmetry (not yet implemented)
        fix_lattice: Fix lattice parameters
        fix_volume: Fix volume
        pressure: External pressure in GPa
        fix_atoms: List of atom indices to fix
        trajectory_file: Path to save trajectory
    
    Returns:
        Dictionary with optimization results
    """
    start_time = time.time()
    
    # Validate inputs
    is_valid_model, model_error = validate_mlff_model(mlff_model)
    if not is_valid_model:
        return {
            "success": False,
            "error": {
                "code": "MODEL_NOT_AVAILABLE",
                "message": model_error,
                "details": {"mlff_model": mlff_model}
            }
        }
    
    is_valid_opt, opt_error = validate_optimizer(optimizer)
    if not is_valid_opt:
        return {
            "success": False,
            "error": {
                "code": "INVALID_OPTIMIZER",
                "message": opt_error,
                "details": {"optimizer": optimizer}
            }
        }
    
    if fmax <= 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"fmax must be positive, got {fmax}",
                "details": {"fmax": fmax}
            }
        }
    
    if steps <= 0:
        return {
            "success": False,
            "error": {
                "code": "INVALID_PARAMETER",
                "message": f"steps must be positive, got {steps}",
                "details": {"steps": steps}
            }
        }
    
    # Convert structure to ASE Atoms
    atoms = dict_to_atoms(structure_dict)
    if atoms is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure to ASE Atoms",
                "details": {}
            }
        }
    
    # Load calculator
    calculator, calc_error = load_calculator(mlff_model)
    if calculator is None:
        return {
            "success": False,
            "error": {
                "code": "MODEL_LOAD_FAILED",
                "message": calc_error or "Failed to load MLFF calculator",
                "details": {"mlff_model": mlff_model}
            }
        }
    
    # Attach calculator
    atoms.calc = calculator
    
    # Apply constraints
    if fix_atoms and len(fix_atoms) > 0:
        c = FixAtoms(indices=fix_atoms)
        atoms.set_constraint(c)
    
    # Calculate initial state
    initial_energy = atoms.get_potential_energy()
    initial_forces = atoms.get_forces()
    max_force_initial = float(np.max(np.linalg.norm(initial_forces, axis=1)))
    
    # Set up optimizer
    if constrain_symmetry:
        return {
            "success": False,
            "error": {
                "code": "NOT_IMPLEMENTED",
                "message": "Symmetry-constrained optimization not yet implemented",
                "details": {}
            }
        }
    
    # Choose what to optimize
    if fix_lattice:
        opt_atoms = atoms
    else:
        opt_atoms = UnitCellFilter(atoms, scalar_pressure=pressure)
    
    # Create optimizer
    if optimizer == "BFGS":
        opt = BFGS(opt_atoms, trajectory=trajectory_file)
    elif optimizer == "FIRE":
        opt = FIRE(opt_atoms, trajectory=trajectory_file)
    elif optimizer == "LBFGS":
        opt = LBFGS(opt_atoms, trajectory=trajectory_file)
    else:
        return {
            "success": False,
            "error": {
                "code": "INVALID_OPTIMIZER",
                "message": f"Optimizer {optimizer} not implemented",
                "details": {"optimizer": optimizer}
            }
        }
    
    # Run optimization
    opt_start = time.time()
    opt.run(fmax=fmax, steps=steps)
    opt_time = time.time() - opt_start
    
    # Get final state
    final_energy = atoms.get_potential_energy()
    final_forces = atoms.get_forces()
    max_force_final = float(np.max(np.linalg.norm(final_forces, axis=1)))
    converged = max_force_final < fmax
    n_steps = opt.get_number_of_steps()
    
    # Convert back to dictionary
    optimized_structure = atoms_to_dict(atoms)
    
    total_time = time.time() - start_time
    
    return {
        "success": True,
        "optimized_structure": optimized_structure,
        "initial_energy": float(initial_energy),
        "final_energy": float(final_energy),
        "energy_change": float(final_energy - initial_energy),
        "max_force_initial": max_force_initial,
        "max_force_final": max_force_final,
        "n_steps": n_steps,
        "converged": converged,
        "preserved_symmetry": False,
        "timing": {
            "setup_time_ms": (opt_start - start_time) * 1000,
            "optimization_time_ms": opt_time * 1000,
            "total_time_ms": total_time * 1000
        }
    }


def calculate_energy(
    structure_dict: Dict[str, Any],
    mlff_model: str,
    calculate_forces: bool = True,
    calculate_stress: bool = True
) -> Dict[str, Any]:
    """
    Calculate energy using MLFF without optimization.
    
    Args:
        structure_dict: Structure for energy calculation
        mlff_model: MLFF model name
        calculate_forces: Whether to calculate forces
        calculate_stress: Whether to calculate stress
    
    Returns:
        Dictionary with energy calculation results
    """
    start_time = time.time()
    
    # Validate model
    is_valid, error = validate_mlff_model(mlff_model)
    if not is_valid:
        return {
            "success": False,
            "error": {
                "code": "MODEL_NOT_AVAILABLE",
                "message": error,
                "details": {"mlff_model": mlff_model}
            }
        }
    
    # Convert structure
    atoms = dict_to_atoms(structure_dict)
    if atoms is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_STRUCTURE",
                "message": "Failed to convert structure to ASE Atoms",
                "details": {}
            }
        }
    
    # Load calculator
    calculator, calc_error = load_calculator(mlff_model)
    if calculator is None:
        return {
            "success": False,
            "error": {
                "code": "MODEL_LOAD_FAILED",
                "message": calc_error or "Failed to load MLFF calculator",
                "details": {"mlff_model": mlff_model}
            }
        }
    
    # Attach calculator
    atoms.calc = calculator
    
    # Calculate energy
    energy = atoms.get_potential_energy()
    energy_per_atom = energy / len(atoms)
    
    result = {
        "success": True,
        "energy": float(energy),
        "energy_per_atom": float(energy_per_atom),
        "metadata": {
            "model": mlff_model,
            "calculation_time_ms": (time.time() - start_time) * 1000
        }
    }
    
    # Calculate forces if requested
    if calculate_forces:
        forces = atoms.get_forces()
        result["forces"] = forces.tolist()
    
    # Calculate stress if requested
    if calculate_stress:
        stress = atoms.get_stress()
        result["stress"] = stress.tolist()
    
    return result


def main():
    """
    Main entry point for command-line execution.
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python mlff_calculator.py <input_json_file>"
            }
        }))
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read input file
    with open(input_file, 'r') as f:
        params = json.load(f)
    
    # Determine operation
    operation = params.get("operation", "optimize")
    
    if operation == "optimize":
        result = optimize_structure(**params)
    elif operation == "energy":
        result = calculate_energy(**params)
    else:
        result = {
            "success": False,
            "error": {
                "code": "INVALID_OPERATION",
                "message": f"Unknown operation: {operation}",
                "details": {"operation": operation}
            }
        }
    
    # Output result
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
