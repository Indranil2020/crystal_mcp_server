"""
Ground State Searcher

Primary backend for ground state search via MCP tools.
Generates structures across multiple space groups and finds the lowest energy state using MLFF.

Related modular generators (for advanced use):
  - generators.bulk.spacegroups.generate_space_group_scan
  - generators.quality_control.optimizer.optimize_structure_mlff
"""
import sys
import json
import numpy as np
from typing import Dict, Any, List, Optional

# Import local modules
# We assume these are in the same directory or python path
from crystal_generator import generate_crystal
from mlff_calculator import optimize_structure, validate_mlff_model

def search_ground_state(
    composition: List[str],
    space_groups: List[int],
    num_structures_per_sg: int = 1,
    mlff_model: str = "chgnet",
    max_attempts: int = 10,
    fmax: float = 0.05
) -> Dict[str, Any]:
    """
    Search for ground state structure.
    
    Args:
        composition: List of elements (e.g. ["C"])
        space_groups: List of space group numbers to search
        num_structures_per_sg: Number of random structures to generate per space group
        mlff_model: MLFF model to use for optimization
        max_attempts: Max attempts for generation
        fmax: Optimization convergence criteria
        
    Returns:
        Dictionary with search results
    """
    # Validate inputs
    if not composition:
        return {"success": False, "error": "Composition cannot be empty"}
    
    if not space_groups:
        return {"success": False, "error": "Space groups list cannot be empty"}
        
    valid_model, error = validate_mlff_model(mlff_model)
    if not valid_model:
        return {"success": False, "error": error}

    results = []
    
    for sg in space_groups:
        for i in range(num_structures_per_sg):
            # 1. Generate random structure
            # Note: num_atoms logic needs to be handled carefully. 
            # crystal_generator usually takes total atoms if provided, or infers.
            # Here we might need to assume a stoichiometry or just pass the species list 
            # and let pyxtal decide multiplicity based on Wyckoff positions, 
            # but generate_crystal expects explicit instructions.
            # For simplicity, we'll try to generate with default settings.
            
            gen_params = {
                "composition": composition, 
                "spacegroup": sg,
                "numIons": [4] * len(set(composition)) if not isinstance(composition[0], str) else None # Placeholder logic
                # Better: Let the user specify exactly what they want or infer from composition length if it's explicit atoms
            }
            
            # If composition is ["Si", "Si"], we likely mean 2 atoms of Si.
            # generate_crystal handles this.
            
            gen_result = generate_crystal(composition, sg)
            
            if not gen_result["success"]:
                results.append({
                    "success": False,
                    "space_group": sg,
                    "error": gen_result.get("error", "Generation failed"),
                    "stage": "generation"
                })
                continue
                
            structure = gen_result["structure"]
            
            # 2. Optimize structure
            opt_result = optimize_structure(
                structure_dict=structure,
                mlff_model=mlff_model,
                fmax=fmax,
                steps=100
            )
            
            if not opt_result["success"]:
                results.append({
                    "success": False,
                    "space_group": sg,
                    "error": opt_result.get("error", {}).get("message", "Optimization failed"),
                    "stage": "optimization"
                })
                continue
                
            # 3. Store result
            results.append({
                "success": True,
                "space_group": sg,
                "initial_structure": structure,
                "optimized_structure": opt_result["optimized_structure"],
                "energy": opt_result["final_energy"],
                "energy_per_atom": opt_result["final_energy"] / structure["metadata"]["natoms"],
                "stage": "complete"
            })

    # Sort results by energy
    valid_results = [r for r in results if r["success"]]
    if not valid_results:
        return {
            "success": False, 
            "error": "No valid structures could be generated or optimized",
            "details": results
        }
        
    valid_results.sort(key=lambda x: x["energy_per_atom"])
    
    return {
        "success": True,
        "best_structure": valid_results[0]["optimized_structure"],
        "min_energy_per_atom": valid_results[0]["energy_per_atom"],
        "all_results": results
    }

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "No input file provided"}))
        sys.exit(1)
        
    input_file = sys.argv[1]
    with open(input_file, 'r') as f:
        data = json.load(f)
        
    result = search_ground_state(
        composition=data.get("composition", []),
        space_groups=data.get("space_groups", []),
        num_structures_per_sg=data.get("num_structures_per_sg", 1),
        mlff_model=data.get("mlff_model", "chgnet"),
        fmax=data.get("fmax", 0.05)
    )
    
    print(json.dumps(result))

if __name__ == "__main__":
    main()
