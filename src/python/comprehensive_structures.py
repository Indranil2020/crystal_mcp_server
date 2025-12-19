"""
comprehensive_structures.py - Unified Entry Point for All Generator Modules

This file provides a single entry point for TypeScript to access all Python generators.
"""

import sys
import json
from typing import Dict, Any


def handle_request(args: Dict[str, Any]) -> Dict[str, Any]:
    """Route requests to appropriate generator module."""
    
    operation = args.get("operation", "")
    
    # Battery operations
    if operation == "generate_cathode":
        from generators.battery import generate_cathode, CATHODE_DATABASE
        result = generate_cathode(
            args.get("cathode_type", "NCM811"),
            lithiation=args.get("lithiation", 1.0),
            supercell=args.get("supercell", [2, 2, 1])
        )
        if args.get("list_available", False):
            result["available_cathodes"] = list(CATHODE_DATABASE.keys())
        return result
    
    elif operation == "generate_electrolyte":
        from generators.battery import generate_electrolyte, ELECTROLYTE_DATABASE
        result = generate_electrolyte(
            args.get("electrolyte_type", "LLZO"),
            supercell=args.get("supercell", [1, 1, 1])
        )
        if args.get("list_available", False):
            result["available_electrolytes"] = list(ELECTROLYTE_DATABASE.keys())
        return result
    
    # Catalyst operations
    elif operation == "generate_saa":
        from generators.catalyst import generate_saa_surface, SAA_DATABASE
        result = generate_saa_surface(
            args.get("saa_type", "Pt/Cu"),
            supercell=args.get("supercell", [4, 4, 4]),
            n_dopant=args.get("n_dopant", 1),
            dopant_layer=args.get("dopant_layer", "surface")
        )
        if args.get("list_available", False):
            result["available_saa"] = list(SAA_DATABASE.keys())
        return result
    
    elif operation == "generate_supported_cluster":
        from generators.catalyst import generate_supported_cluster, CLUSTER_DATABASE, OXIDE_SUPPORTS
        result = generate_supported_cluster(
            args.get("cluster", "Pt4"),
            args.get("support", "CeO2"),
            supercell=args.get("supercell", [3, 3, 4])
        )
        if args.get("list_available", False):
            result["available_clusters"] = list(CLUSTER_DATABASE.keys())
            result["available_supports"] = list(OXIDE_SUPPORTS.keys())
        return result
    
    elif operation == "generate_zeolite_cluster":
        from generators.catalyst import generate_zeolite_cluster, ZEOLITE_FRAMEWORKS
        result = generate_zeolite_cluster(
            args.get("metal", "Cu"),
            args.get("zeolite", "CHA"),
            args.get("cluster_size", 2)
        )
        if args.get("list_available", False):
            result["available_zeolites"] = list(ZEOLITE_FRAMEWORKS.keys())
        return result
    
    # Photonic operations
    elif operation == "generate_metasurface":
        from generators.photonic import generate_metasurface_unit_cell, METAMATERIAL_DATABASE
        result = generate_metasurface_unit_cell(
            args.get("element_type", "V_antenna"),
            material=args.get("material", "Au"),
            period_nm=args.get("period_nm", 300.0),
            rotation_deg=args.get("rotation_deg", 0.0)
        )
        if args.get("list_available", False):
            result["available_elements"] = list(METAMATERIAL_DATABASE.keys())
        return result
    
    elif operation == "generate_photonic_crystal":
        from generators.photonic import generate_photonic_crystal, PHOTONIC_CRYSTAL_DATABASE
        result = generate_photonic_crystal(
            args.get("crystal_type", "woodpile"),
            material=args.get("material", "Si"),
            lattice_constant_nm=args.get("lattice_constant_nm", 500.0)
        )
        if args.get("list_available", False):
            result["available_types"] = list(PHOTONIC_CRYSTAL_DATABASE.keys())
        return result
    
    # Electronic operations
    elif operation == "generate_topological_insulator":
        from generators.electronic.topological import generate_topological_insulator, TOPOLOGICAL_INSULATOR_DATABASE
        result = generate_topological_insulator(
            args.get("material", "Bi2Se3"),
            thickness_QL=args.get("thickness_QL", 6)
        )
        if args.get("list_available", False):
            result["available_materials"] = list(TOPOLOGICAL_INSULATOR_DATABASE.keys())
        return result
    
    elif operation == "generate_majorana_nanowire":
        from generators.electronic.topological import generate_majorana_nanowire, MAJORANA_DATABASE
        result = generate_majorana_nanowire(
            args.get("platform", "InSb_Al"),
            wire_length_nm=args.get("wire_length_nm", 1000),
            wire_diameter_nm=args.get("wire_diameter_nm", 100)
        )
        if args.get("list_available", False):
            result["available_platforms"] = list(MAJORANA_DATABASE.keys())
        return result
    
    # Nanotube operations
    elif operation == "generate_cnt":
        from generators.nanotube import generate_cnt, CNT_DATABASE
        result = generate_cnt(
            args.get("n", 10),
            args.get("m", 10),
            length_nm=args.get("length_nm", 5.0),
            capped=args.get("capped", False)
        )
        if args.get("list_available", False):
            result["common_chiralities"] = list(CNT_DATABASE.keys())
        return result
    
    elif operation == "generate_nanowire":
        from generators.nanotube import generate_nanowire, NANOWIRE_DATABASE
        result = generate_nanowire(
            args.get("material", "Si"),
            diameter_nm=args.get("diameter_nm", 3.0),
            length_nm=args.get("length_nm", 10.0)
        )
        if args.get("list_available", False):
            result["available_materials"] = list(NANOWIRE_DATABASE.keys())
        return result
    
    # High pressure operations
    elif operation == "generate_superhydride":
        from generators.high_pressure import generate_superhydride, SUPERHYDRIDE_DATABASE
        result = generate_superhydride(
            args.get("hydride", "LaH10"),
            pressure_GPa=args.get("pressure_GPa", 200)
        )
        if args.get("list_available", False):
            result["available_hydrides"] = list(SUPERHYDRIDE_DATABASE.keys())
        return result
    
    # List all available generators
    elif operation == "list_generators":
        return {
            "success": True,
            "generators": {
                "battery": ["generate_cathode", "generate_electrolyte"],
                "catalyst": ["generate_saa", "generate_supported_cluster", "generate_zeolite_cluster"],
                "photonic": ["generate_metasurface", "generate_photonic_crystal"],
                "electronic": ["generate_topological_insulator", "generate_majorana_nanowire"],
                "nanotube": ["generate_cnt", "generate_nanowire"],
                "high_pressure": ["generate_superhydride"],
            }
        }
    
    else:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_OPERATION",
                "message": f"Unknown operation: {operation}",
                "available_operations": [
                    "generate_cathode", "generate_electrolyte",  # Battery
                    "generate_saa", "generate_supported_cluster", "generate_zeolite_cluster",  # Catalyst
                    "generate_metasurface", "generate_photonic_crystal",  # Photonic
                    "generate_topological_insulator", "generate_majorana_nanowire",  # Electronic
                    "generate_cnt", "generate_nanowire",  # Nanotube
                    "generate_superhydride",  # High pressure
                    "list_generators"
                ]
            }
        }


if __name__ == "__main__":
    try:
        input_data = json.load(sys.stdin)
        result = handle_request(input_data)
        print(json.dumps(result))
    except Exception as e:
        print(json.dumps({
            "success": False,
            "error": {"code": "PYTHON_ERROR", "message": str(e)}
        }))
