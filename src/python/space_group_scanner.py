#!/usr/bin/env python3
"""
Space Group Scanner - Scan multiple space groups to find valid structures

This module provides functionality to generate structures across multiple
space groups for a given composition, useful for polymorph searching and
phase stability analysis.
"""

import json
import sys
import os
import importlib.util
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import generators (no try/except; pre-check module path)
if importlib.util.find_spec("generators.bulk.spacegroups") is None:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from generators.bulk.spacegroups import generate_from_spacegroup, CRYSTAL_SYSTEMS


def scan_space_groups(
    composition: List[str],
    space_groups: Optional[List[int]] = None,
    space_group_range: Optional[List[int]] = None,
    crystal_systems: Optional[List[str]] = None,
    volume_factor: float = 1.0,
    num_atoms: Optional[int] = None,
    parallel: bool = False,
    output_directory: Optional[str] = None,
    naming_scheme: Optional[str] = None,
    max_attempts: int = 50
) -> Dict[str, Any]:
    """
    Scan multiple space groups to generate structures for a given composition.

    Args:
        composition: List of element symbols (e.g., ['Si', 'O', 'O'] for SiO2)
        space_groups: Specific space group numbers to try
        space_group_range: Range [start, end] of space groups (inclusive)
        crystal_systems: Filter by crystal system names
        volume_factor: Volume scaling factor
        num_atoms: Target number of atoms
        parallel: Enable parallel generation
        output_directory: Directory to save structures
        naming_scheme: Filename template (e.g., '{formula}_{spg}.cif')
        max_attempts: Max attempts per space group

    Returns:
        Dictionary with scan results including successful and failed structures
    """
    # Determine which space groups to scan
    sgs_to_scan = []

    if space_groups:
        sgs_to_scan = space_groups
    elif space_group_range:
        start, end = space_group_range
        sgs_to_scan = list(range(max(1, start), min(231, end + 1)))
    elif crystal_systems:
        for system in crystal_systems:
            if system.lower() in CRYSTAL_SYSTEMS:
                info = CRYSTAL_SYSTEMS[system.lower()]
                sgs_to_scan.extend(range(info["min_sg"], info["max_sg"] + 1))
    else:
        # Default: all 230 space groups
        sgs_to_scan = list(range(1, 231))

    # Remove duplicates and sort
    sgs_to_scan = sorted(set(sgs_to_scan))

    # Count elements for composition
    element_counts = {}
    for elem in composition:
        element_counts[elem] = element_counts.get(elem, 0) + 1
    elements = list(element_counts.keys())
    comp = list(element_counts.values())

    results = {
        "success": True,
        "generated_structures": [],
        "failed_space_groups": [],
        "statistics": {
            "total_attempted": len(sgs_to_scan),
            "successful": 0,
            "failed": 0
        }
    }

    def try_space_group(sg: int) -> Dict[str, Any]:
        """Try to generate structure in a single space group."""
        if not isinstance(sg, int):
            return {
                "space_group": sg,
                "success": False,
                "error": "Space group must be an integer"
            }
        if sg < 1 or sg > 230:
            return {
                "space_group": sg,
                "success": False,
                "error": f"Space group out of range: {sg}"
            }

        result = generate_from_spacegroup(
            spacegroup=sg,
            elements=elements,
            composition=comp,
            factor=volume_factor  # Map volume_factor to factor parameter
        )

        if result.get("success"):
            return {
                "space_group": sg,
                "success": True,
                "structure": result.get("structure"),
                "spacegroup_symbol": result.get("spacegroup_symbol", ""),
                "crystal_system": result.get("crystal_system", "")
            }

        return {
            "space_group": sg,
            "success": False,
            "error": result.get("error", {}).get("message", "Generation failed")
        }

    if parallel:
        # Parallel execution
        with ThreadPoolExecutor(max_workers=min(8, len(sgs_to_scan))) as executor:
            futures = {executor.submit(try_space_group, sg): sg for sg in sgs_to_scan}
            for future in as_completed(futures):
                result = future.result()
                if result["success"]:
                    results["generated_structures"].append(result)
                    results["statistics"]["successful"] += 1
                else:
                    results["failed_space_groups"].append(result)
                    results["statistics"]["failed"] += 1
    else:
        # Sequential execution
        for sg in sgs_to_scan:
            result = try_space_group(sg)
            if result["success"]:
                results["generated_structures"].append(result)
                results["statistics"]["successful"] += 1
            else:
                results["failed_space_groups"].append(result)
                results["statistics"]["failed"] += 1

    # Sort results by space group number
    results["generated_structures"].sort(key=lambda x: x["space_group"])
    results["failed_space_groups"].sort(key=lambda x: x["space_group"])

    # Save structures if output directory specified
    if output_directory and results["generated_structures"]:
        os.makedirs(output_directory, exist_ok=True)
        pymatgen_available = importlib.util.find_spec("pymatgen.core") is not None

        if pymatgen_available:
            from pymatgen.core import Structure, Lattice

            for struct_info in results["generated_structures"]:
                sg = struct_info["space_group"]
                structure = struct_info["structure"]

                if not structure or "lattice" not in structure:
                    continue

                lat = structure["lattice"]
                required = ["a", "b", "c", "alpha", "beta", "gamma"]
                if any(key not in lat for key in required):
                    continue

                sites = structure.get("sites", structure.get("atoms", []))
                if not sites:
                    continue

                # Generate filename
                if naming_scheme:
                    formula = "".join(f"{e}{c}" for e, c in element_counts.items())
                    filename = naming_scheme.format(formula=formula, spg=sg)
                else:
                    formula = "".join(f"{e}{c}" for e, c in element_counts.items())
                    filename = f"{formula}_sg{sg}.cif"

                filepath = os.path.join(output_directory, filename)

                lattice = Lattice.from_parameters(
                    lat["a"], lat["b"], lat["c"],
                    lat["alpha"], lat["beta"], lat["gamma"]
                )
                species = [s["element"] for s in sites]
                coords = [s["coords"] for s in sites]
                pmg_structure = Structure(lattice, species, coords)
                pmg_structure.to(filename=filepath)

    return results


def main():
    """Main entry point for CLI usage."""
    if len(sys.argv) < 2:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "INVALID_USAGE",
                "message": "Usage: python space_group_scanner.py <input.json>"
            }
        }))
        sys.exit(1)

    input_file = sys.argv[1]

    with open(input_file, 'r') as f:
        params = json.load(f)

    # Execute scan
    result = scan_space_groups(
        composition=params.get("composition", []),
        space_groups=params.get("space_groups"),
        space_group_range=params.get("space_group_range"),
        crystal_systems=params.get("crystal_systems"),
        volume_factor=params.get("volume_factor", 1.0),
        num_atoms=params.get("num_atoms"),
        parallel=params.get("parallel", False),
        output_directory=params.get("output_directory"),
        naming_scheme=params.get("naming_scheme"),
        max_attempts=params.get("max_attempts", 50)
    )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
