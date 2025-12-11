#!/usr/bin/env python3
"""
Symmetry Transformations Backend

This module provides tools for exploring group-subgroup relationships
and symmetry transformations using PyXtal and Pymatgen.
"""

import sys
import json
from typing import Dict, Any, List, Optional
import numpy as np
from pyxtal.symmetry import Group

def get_subgroups(space_group: int, strict: bool = False) -> Dict[str, Any]:
    """
    Get maximal subgroups of a space group.
    
    Args:
        space_group: Space group number (1-230)
        strict: If True, only return direct maximal subgroups
        
    Returns:
        Dictionary with subgroup information
    """
    if space_group < 1 or space_group > 230:
         return {
            "success": False,
            "error": {
                "code": "INVALID_SPACE_GROUP",
                "message": f"Space group must be 1-230, got {space_group}"
            }
        }
    
    g = Group(space_group)
    # PyXtal functionality for subgroups
    # get_max_subgroup_numbers() returns list of list? or list of ints?
    # Based on inspection: [141, 141, ...]
    subgroup_numbers = g.get_max_subgroup_numbers()
    
    # We can also get detailed relation?
    # g.get_max_t_subgroup() ? g.get_max_k_subgroup() ?
    # Let's collect unique subgroups and count types if possible
    
    unique_subgroups = sorted(list(set(subgroup_numbers)))
    subgroups_data = []
    
    for num in unique_subgroups:
        sub_g = Group(num)
        subgroups_data.append({
            "number": int(num),
            "symbol": sub_g.symbol,
            "crystal_system": sub_g.lattice_type,
            "point_group": sub_g.point_group
        })
        
    return {
        "success": True,
        "space_group": {
            "number": space_group,
            "symbol": g.symbol,
            "hall_number": g.hall_number
        },
        "subgroups": subgroups_data,
        "raw_numbers": subgroup_numbers
    }

def get_symmetry_path(start_spg: int, end_spg: int, max_depth: int = 5) -> Dict[str, Any]:
    """
    Find a path from supergroup to subgroup.
    Using simple BFS on PyXtal group relations.
    """
    if start_spg == end_spg:
         return {"success": True, "path": [start_spg]}
         
    queue = [(start_spg, [start_spg])]
    visited = {start_spg}
    
    while queue:
        current, path = queue.pop(0)
        
        if len(path) > max_depth:
            continue
            
        if current == end_spg:
             # Found path
             formatted_path = []
             for spg_num in path:
                 g = Group(spg_num)
                 formatted_path.append({
                     "number": spg_num,
                     "symbol": g.symbol
                 })
             return {
                 "success": True, 
                 "path": formatted_path,
                 "length": len(path) - 1
             }
        
        # Get subgroups
        # Get subgroups
        g = Group(current)
        subs = g.get_max_subgroup_numbers()
        for sub in subs:
            if sub not in visited:
                visited.add(sub)
                new_path = list(path)
                new_path.append(sub)
                queue.append((sub, new_path))
            
    return {
        "success": False, 
        "error": {
            "code": "PATH_NOT_FOUND", 
            "message": f"No subgroup path found from {start_spg} to {end_spg} within depth {max_depth}"
        }
    }

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "Usage: python symmetry_transformations.py <input.json>"}))
        sys.exit(1)
        
    input_file = sys.argv[1]
    with open(input_file, 'r') as f:
        params = json.load(f)
        
    operation = params.get("operation")
    
    if operation == "get_subgroups":
        result = get_subgroups(params.get("space_group"))
    elif operation == "get_path":
        result = get_symmetry_path(params.get("start_group"), params.get("end_group"), params.get("max_depth", 5))
    else:
        result = {"success": False, "error": f"Unknown operation {operation}"}
        
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
