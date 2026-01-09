#!/usr/bin/env python3
"""
Dependency Check Module

Provides graceful degradation for optional dependencies.
Checks availability of all optional Python packages used by Crystal MCP Server
and returns capability flags with installation instructions.

Usage:
    from dependency_check import check_dependencies, get_capabilities
    
    caps = get_capabilities()
    if caps["rdkit"]["available"]:
        from rdkit import Chem
        # use RDKit
    else:
        print(caps["rdkit"]["install_hint"])
"""

import sys
import importlib.util
import time
from typing import Dict, Any

# Debug logging
def _debug_log(action: str, details: str = "") -> None:
    """Print debug message with timestamp."""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[DEPENDENCY_CHECK] {timestamp} {action}: {details}", file=sys.stderr)


# List of optional dependencies with metadata
OPTIONAL_DEPENDENCIES = {
    "rdkit": {
        "module": "rdkit",
        "description": "Chemistry toolkit for SMILES parsing and 3D generation",
        "install_hint": "pip install rdkit",
        "features": ["SMILES parsing", "3D coordinate generation", "Molecular optimization"],
    },
    "pubchempy": {
        "module": "pubchempy",
        "description": "PubChem API client for molecule lookups",
        "install_hint": "pip install pubchempy",
        "features": ["PubChem compound search", "CID lookup", "IUPAC name resolution"],
    },
    "ase": {
        "module": "ase",
        "description": "Atomic Simulation Environment",
        "install_hint": "pip install ase",
        "features": ["G2 molecule database", "Atomic structure manipulation"],
    },
    "pymatgen": {
        "module": "pymatgen",
        "description": "Python Materials Genomics",
        "install_hint": "pip install pymatgen",
        "features": ["Crystal structure handling", "Space group analysis", "CIF export"],
    },
    "pyxtal": {
        "module": "pyxtal",
        "description": "Crystal structure generation",
        "install_hint": "pip install pyxtal",
        "features": ["Random crystal generation", "Molecular crystal packing"],
    },
    "chgnet": {
        "module": "chgnet",
        "description": "Machine learning force field for crystals",
        "install_hint": "pip install chgnet",
        "features": ["Structure relaxation", "Energy prediction"],
    },
    "mace": {
        "module": "mace",
        "description": "MACE machine learning potentials",
        "install_hint": "pip install mace-torch",
        "features": ["Structure relaxation", "Energy prediction"],
    },
    "opsin": {
        "module": "opsin",
        "description": "IUPAC name to structure converter",
        "install_hint": "pip install opsin",
        "features": ["IUPAC name parsing", "SMILES generation from names"],
    },
}


def check_dependency(name: str) -> Dict[str, Any]:
    """
    Check if a single dependency is available.
    
    Args:
        name: Name of the dependency (key in OPTIONAL_DEPENDENCIES)
        
    Returns:
        Dictionary with availability status and metadata
    """
    if name not in OPTIONAL_DEPENDENCIES:
        return {
            "available": False,
            "error": f"Unknown dependency: {name}",
        }
    
    dep_info = OPTIONAL_DEPENDENCIES[name]
    module_name = dep_info["module"]
    
    is_available = importlib.util.find_spec(module_name) is not None
    
    result = {
        "available": is_available,
        "module": module_name,
        "description": dep_info["description"],
        "install_hint": dep_info["install_hint"],
        "features": dep_info["features"],
    }
    
    if is_available:
        # Try to get version
        imported = importlib.import_module(module_name)
        version = getattr(imported, "__version__", None)
        if version:
            result["version"] = version
    
    return result


def check_all_dependencies() -> Dict[str, Dict[str, Any]]:
    """
    Check all optional dependencies.
    
    Returns:
        Dictionary mapping dependency names to their availability status
    """
    _debug_log("CHECK", "Checking all optional dependencies...")
    
    results = {}
    for name in OPTIONAL_DEPENDENCIES:
        results[name] = check_dependency(name)
        status = "✓" if results[name]["available"] else "✗"
        _debug_log("STATUS", f"{name}: {status}")
    
    return results


def get_capabilities() -> Dict[str, Any]:
    """
    Get system capabilities based on available dependencies.
    
    Returns:
        Dictionary with capability flags and dependency status
    """
    deps = check_all_dependencies()
    
    # Calculate high-level capabilities
    capabilities = {
        "molecule_generation": {
            "smiles_to_3d": deps["rdkit"]["available"],
            "pubchem_lookup": deps["pubchempy"]["available"],
            "iupac_names": deps["opsin"]["available"] or deps["rdkit"]["available"],
            "g2_database": deps["ase"]["available"],
        },
        "crystal_generation": {
            "random_crystals": deps["pyxtal"]["available"],
            "molecular_crystals": deps["pyxtal"]["available"],
            "structure_analysis": deps["pymatgen"]["available"],
        },
        "optimization": {
            "ml_relaxation": deps["chgnet"]["available"] or deps["mace"]["available"],
            "force_field": deps["rdkit"]["available"],
        },
        "dependencies": deps,
    }
    
    # Count available vs total
    available_count = sum(1 for d in deps.values() if d["available"])
    total_count = len(deps)
    capabilities["summary"] = {
        "available": available_count,
        "total": total_count,
        "percentage": round(100 * available_count / total_count, 1),
    }
    
    _debug_log("SUMMARY", f"{available_count}/{total_count} dependencies available")
    
    return capabilities


def print_dependency_report():
    """Print a human-readable dependency report to stderr."""
    caps = get_capabilities()
    
    print("\n" + "=" * 60, file=sys.stderr)
    print("  Crystal MCP Server - Dependency Report", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    print(f"\n  Status: {caps['summary']['available']}/{caps['summary']['total']} dependencies available ({caps['summary']['percentage']}%)\n", file=sys.stderr)
    
    deps = caps["dependencies"]
    
    print("  Package          Available   Version", file=sys.stderr)
    print("  " + "-" * 45, file=sys.stderr)
    
    for name, info in deps.items():
        status = "✓ Yes" if info["available"] else "✗ No"
        version = info.get("version", "-")
        print(f"  {name:16} {status:11} {version}", file=sys.stderr)
    
    # Print missing dependencies with install hints
    missing = [name for name, info in deps.items() if not info["available"]]
    if missing:
        print("\n  Missing dependencies:", file=sys.stderr)
        for name in missing:
            print(f"    {name}: {deps[name]['install_hint']}", file=sys.stderr)
    
    print("\n" + "=" * 60 + "\n", file=sys.stderr)


if __name__ == "__main__":
    import json
    
    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        # JSON output for programmatic use
        print(json.dumps(get_capabilities(), indent=2))
    else:
        # Human-readable report
        print_dependency_report()
