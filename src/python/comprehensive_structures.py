"""
comprehensive_structures.py - Unified Entry Point for All Generator Modules

This file provides a single entry point for TypeScript to access all Python generators.
Routes requests through the central registry in generators/__init__.py.
"""

import sys
import json
import importlib
from typing import Dict, Any


def handle_request(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    Route requests to appropriate generator module.
    
    Args:
        args: Dictionary with 'operation' key and parameter keys
    
    Returns:
        Result from the generator function
    """
    from generators import GENERATOR_REGISTRY, get_all_operations, get_operation_info, list_category_operations
    
    operation = args.get("operation", "")
    
    # Special operations
    if operation == "list_all":
        return {
            "success": True,
            "categories": list(GENERATOR_REGISTRY.keys()),
            "operations": get_all_operations(),
            "total_operations": sum(len(ops) for ops in get_all_operations().values())
        }
    
    if operation == "list_category":
        category = args.get("category", "")
        return list_category_operations(category)
    
    if operation == "operation_info":
        op_name = args.get("operation_name", "")
        return get_operation_info(op_name)
    
    # Find operation in registry
    # Prioritize requested category if provided
    search_categories = list(GENERATOR_REGISTRY.keys())
    preferred_cat = args.get("category")
    if preferred_cat and preferred_cat in GENERATOR_REGISTRY:
        search_categories.remove(preferred_cat)
        search_categories.insert(0, preferred_cat)
        
    for category in search_categories:
        data = GENERATOR_REGISTRY[category]
        if operation in data["operations"]:
            op_info = data["operations"][operation]
            
            # Import module dynamically
            module_path = op_info["module"]
            function_name = op_info["function"]
            
            module = importlib.import_module(module_path)
            func = getattr(module, function_name)
            
            # Extract parameters (remove 'operation' and 'category' keys)
            params = {k: v for k, v in args.items() if k not in ["operation", "category"]}
            
            # AUTO-CONVERSION: Convert dicts to Structure objects if needed
            # This enables iterative workflows where output of one tool is input to another
            from pymatgen.core import Structure, Lattice
            
            for key, value in params.items():
                if isinstance(value, dict):
                    # Check if it looks like a structure
                    if "lattice" in value and ("sites" in value or "atoms" in value):
                        try:
                            # Handle our custom dict format
                            if "matrix" in value["lattice"] and "atoms" in value:
                                lat = Lattice(value["lattice"]["matrix"])
                                species = [s["element"] for s in value["atoms"]]
                                coords = [s["coords"] for s in value["atoms"]]
                                params[key] = Structure(lat, species, coords)
                            # Handle pymatgen as_dict format
                            elif "sites" in value:
                                params[key] = Structure.from_dict(value)
                        except Exception as e:
                            # If conversion fails, pass as dict (some functions might expect dict)
                            pass
            
            # Call the function with error handling
            try:
                result = func(**params)
            except TypeError as e:
                # Handle missing required parameters
                error_msg = str(e)
                if "missing" in error_msg and "required" in error_msg:
                    return {
                        "success": False,
                        "error": {
                            "code": "MISSING_REQUIRED_PARAMETER",
                            "message": error_msg,
                            "operation": operation
                        }
                    }
                raise
            except Exception as e:
                return {
                    "success": False,
                    "error": {
                        "code": "GENERATOR_ERROR",
                        "message": str(e),
                        "operation": operation
                    }
                }

            # Ensure success key exists
            if isinstance(result, dict) and "success" not in result:
                result["success"] = True

            return result
    
    # Operation not found - return helpful error
    return {
        "success": False,
        "error": {
            "code": "UNKNOWN_OPERATION",
            "message": f"Unknown operation: '{operation}'",
            "hint": "Use operation='list_all' to see all available operations",
            "available_categories": list(GENERATOR_REGISTRY.keys())
        }
    }


# =============================================================================
# DIRECT OPERATION HANDLERS (for backward compatibility)
# =============================================================================

def _handle_battery(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle battery operations."""
    op = args.get("operation", "")
    
    if op == "generate_cathode":
        from generators.battery.cathodes import generate_cathode
        return generate_cathode(
            args.get("cathode_type", "NCM811"),
            lithiation=args.get("lithiation", 1.0),
            supercell=args.get("supercell", [2, 2, 1])
        )
    
    if op == "generate_electrolyte":
        from generators.battery.electrolytes import generate_electrolyte
        return generate_electrolyte(
            args.get("electrolyte_type", "LLZO"),
            supercell=args.get("supercell", [1, 1, 1])
        )
    
    return {"success": False, "error": {"message": f"Unknown battery operation: {op}"}}


def _handle_catalyst(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle catalyst operations."""
    op = args.get("operation", "")
    
    if op == "generate_saa":
        from generators.catalyst.saa import generate_saa_surface
        return generate_saa_surface(
            args.get("saa_type", "Pt/Cu"),
            supercell=args.get("supercell", [4, 4, 4]),
            n_dopant=args.get("n_dopant", 1)
        )
    
    if op == "generate_supported_cluster":
        from generators.catalyst.clusters import generate_supported_cluster
        return generate_supported_cluster(
            args.get("cluster", "Pt4"),
            args.get("support", "CeO2"),
            supercell=args.get("supercell", [3, 3, 4])
        )
    
    return {"success": False, "error": {"message": f"Unknown catalyst operation: {op}"}}


def _handle_electronic(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle electronic operations."""
    op = args.get("operation", "")
    
    if op == "generate_topological_insulator":
        from generators.electronic.topological import generate_topological_insulator
        return generate_topological_insulator(
            args.get("material", "Bi2Se3"),
            thickness_QL=args.get("thickness_QL", 6)
        )
    
    if op == "generate_quantum_well":
        from generators.electronic.semiconductors import generate_quantum_well
        return generate_quantum_well(
            args.get("well_material", "GaAs"),
            args.get("barrier_material", "AlAs"),
            well_thickness_nm=args.get("well_thickness_nm", 10.0),
            barrier_thickness_nm=args.get("barrier_thickness_nm", 20.0),
            n_wells=args.get("n_wells", 1)
        )
    
    if op == "generate_superconductor":
        from generators.electronic.superconductors import generate_superconductor
        return generate_superconductor(
            args.get("material", "YBCO"),
            supercell=args.get("supercell", [1, 1, 1])
        )
    
    return {"success": False, "error": {"message": f"Unknown electronic operation: {op}"}}


def _handle_bulk(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle bulk structure operations."""
    op = args.get("operation", "")
    
    if op == "generate_from_spacegroup":
        from generators.bulk.spacegroups import generate_from_spacegroup
        return generate_from_spacegroup(
            args.get("spacegroup", 225),
            args.get("elements", ["Na", "Cl"]),
            args.get("composition", [4, 4]),
            a=args.get("a", 5.0),
            factor=args.get("factor", 1.0)
        )
    
    if op == "generate_prototype":
        from generators.bulk.spacegroups import generate_prototype
        return generate_prototype(
            args.get("prototype", "rocksalt"),
            args.get("elements", {"cation": "Na", "anion": "Cl"}),
            a=args.get("a", 5.0)
        )
    
    if op == "generate_zeolite":
        from generators.bulk.zeolites import generate_zeolite
        return generate_zeolite(
            args.get("framework", "MFI"),
            cation=args.get("cation", "Na"),
            Si_Al_ratio=args.get("Si_Al_ratio", 30.0)
        )
    
    if op == "generate_clathrate":
        from generators.bulk.clathrates import generate_clathrate
        return generate_clathrate(
            args.get("clathrate", "Ba8Ga16Ge30"),
            empty_cages=args.get("empty_cages", False)
        )
    
    return {"success": False, "error": {"message": f"Unknown bulk operation: {op}"}}


def _handle_defect(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle defect operations."""
    op = args.get("operation", "")
    
    if op == "generate_dislocation":
        from generators.defect.extended_defects import generate_dislocation
        return generate_dislocation(
            args.get("element", "Cu"),
            args.get("dislocation_type", "FCC_edge_110"),
            supercell=args.get("supercell", [10, 10, 3])
        )
    
    if op == "generate_grain_boundary":
        from generators.defect.extended_defects import generate_grain_boundary
        return generate_grain_boundary(
            args.get("element", "Cu"),
            args.get("gb_type", "Sigma3_twin_FCC"),
            supercell=args.get("supercell", [4, 4, 8])
        )
    
    if op == "generate_amorphous":
        from generators.defect.amorphous import generate_amorphous_structure
        return generate_amorphous_structure(
            args.get("material", "SiO2"),
            n_atoms=args.get("n_atoms", 100)
        )
    
    return {"success": False, "error": {"message": f"Unknown defect operation: {op}"}}


def _handle_molecule(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle molecule operations."""
    op = args.get("operation", "")
    
    if op == "generate_mof":
        from generators.molecule.frameworks import generate_mof
        return generate_mof(
            args.get("mof_name", "MOF-5"),
            supercell=args.get("supercell", [1, 1, 1])
        )
    
    if op == "generate_cof":
        from generators.molecule.frameworks import generate_cof
        return generate_cof(
            args.get("cof_name", "COF-5"),
            n_layers=args.get("n_layers", 3)
        )
    
    if op == "generate_fullerene":
        from generators.molecule.cages import generate_fullerene
        return generate_fullerene(args.get("fullerene", "C60"))
    
    if op == "generate_metal_cluster":
        from generators.molecule.cages import generate_metal_cluster
        return generate_metal_cluster(args.get("cluster", "Au13"))
    
    return {"success": False, "error": {"message": f"Unknown molecule operation: {op}"}}


def _handle_surface(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle surface operations."""
    op = args.get("operation", "")
    
    if op == "generate_nanoparticle":
        from generators.surface.nanoparticles import generate_nanoparticle
        return generate_nanoparticle(
            args.get("element", "Au"),
            args.get("shape", "cuboctahedron"),
            size_nm=args.get("size_nm", 3.0)
        )
    
    if op == "generate_reconstruction":
        from generators.surface.reconstructions import generate_reconstruction
        return generate_reconstruction(args.get("reconstruction_type", "Si(100)-2x1"))
    
    return {"success": False, "error": {"message": f"Unknown surface operation: {op}"}}


def _handle_twist(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle twist operations."""
    op = args.get("operation", "")
    
    if op == "generate_multilayer_stack":
        from generators.twist.multilayers import generate_multilayer_stack
        return generate_multilayer_stack(
            args.get("system", "twisted_trilayer_ALT"),
            supercell=args.get("supercell", [5, 5])
        )
    
    return {"success": False, "error": {"message": f"Unknown twist operation: {op}"}}


def _handle_two_d(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle 2D material operations."""
    op = args.get("operation", "")
    
    if op == "generate_vertical_heterostructure":
        from generators.two_d.he_2d import generate_vertical_heterostructure
        return generate_vertical_heterostructure(
            args.get("materials", ["graphene", "hBN"]),
            supercell=args.get("supercell", [5, 5])
        )
    
    return {"success": False, "error": {"message": f"Unknown 2D operation: {op}"}}


def _handle_magnetic(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle magnetic operations."""
    op = args.get("operation", "")
    
    if op == "generate_magnetic_structure":
        from generators.magnetic.materials import generate_magnetic_structure
        return generate_magnetic_structure(
            args.get("material", "Fe_bcc"),
            ordering=args.get("ordering", "FM"),
            supercell=args.get("supercell", [2, 2, 2])
        )
    
    return {"success": False, "error": {"message": f"Unknown magnetic operation: {op}"}}


def _handle_quantum(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle quantum operations."""
    op = args.get("operation", "")
    
    if op == "generate_qubit_structure":
        from generators.quantum.qubits import generate_qubit_structure
        return generate_qubit_structure(
            args.get("platform", "NV_diamond"),
            supercell=args.get("supercell", [3, 3, 3])
        )
    
    return {"success": False, "error": {"message": f"Unknown quantum operation: {op}"}}


def _handle_photonic(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle photonic operations."""
    op = args.get("operation", "")
    
    if op == "generate_photonic_crystal":
        from generators.photonic.crystals import generate_2d_photonic_crystal
        return generate_2d_photonic_crystal(
            args.get("crystal_type", "triangular_holes"),
            material=args.get("material", "Si"),
            a=args.get("a", 500.0)
        )
    
    return {"success": False, "error": {"message": f"Unknown photonic operation: {op}"}}


def _handle_nanotube(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle nanotube operations."""
    op = args.get("operation", "")
    
    if op == "generate_cnt":
        from generators.nanotube.cnt import generate_cnt
        return generate_cnt(
            args.get("n", 10),
            args.get("m", 10),
            length_nm=args.get("length_nm", 5.0)
        )
    
    if op == "generate_nanowire":
        from generators.nanotube.nanowires import generate_nanowire
        return generate_nanowire(
            args.get("material", "Si"),
            diameter_nm=args.get("diameter_nm", 3.0),
            length_nm=args.get("length_nm", 10.0)
        )
    
    return {"success": False, "error": {"message": f"Unknown nanotube operation: {op}"}}


def _handle_high_pressure(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle high pressure operations."""
    op = args.get("operation", "")
    
    if op == "generate_superhydride":
        from generators.high_pressure.phases import generate_superhydride
        return generate_superhydride(
            args.get("hydride", "LaH10"),
            pressure_GPa=args.get("pressure_GPa", 200)
        )
    
    return {"success": False, "error": {"message": f"Unknown high pressure operation: {op}"}}


def _handle_thermoelectric(args: Dict[str, Any]) -> Dict[str, Any]:
    """Handle thermoelectric operations."""
    op = args.get("operation", "")
    
    if op == "generate_thermoelectric":
        from generators.thermoelectric.materials import generate_thermoelectric
        return generate_thermoelectric(
            args.get("material", "Bi2Te3"),
            supercell=args.get("supercell", [1, 1, 1])
        )
    
    return {"success": False, "error": {"message": f"Unknown thermoelectric operation: {op}"}}


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # Check if input file provided as argument
    if len(sys.argv) > 1 and sys.argv[1]:
        input_file = sys.argv[1]
        with open(input_file, 'r') as f:
            input_data = json.load(f)
    else:
        # Fallback to stdin
        input_data = json.load(sys.stdin)
    
    # Process request
    result = handle_request(input_data)
    
    # Output result
    print(json.dumps(result))
