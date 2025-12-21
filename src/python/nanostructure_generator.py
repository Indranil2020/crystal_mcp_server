"""
Nanostructure Generator using ASE

DEPRECATED: This module is maintained for backward compatibility.
New code should use the modular generators:
  - generators.nanotube.cnt.generate_cnt
  - generators.two_d.graphene.generate_graphene
  - generators.molecule.cages.generate_fullerene

Supports generation of nanotubes, nanoribbons, graphene sheets, and fullerenes.
"""
import sys
import json
import json
import numpy as np
from ase.build import nanotube, graphene_nanoribbon, molecule, mx2, bulk
from ase.build.molecule import extra
from ase.lattice.hexagonal import Graphene
from ase.collections import g2
from typing import Dict, Any, List

def generate_nanostructure(
    structure_type: str,
    params: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Generate nanostructures using ASE.
    
    Args:
        structure_type: Type of nanostructure ('nanotube', 'graphene', 'nanoribbon', 'fullerene', 'mos2')
        params: Parameters specific to the structure type
    
    Returns:
        Dictionary containing structure data
    """
    atoms = None
    
    # 1. Generate core structure
    if structure_type == "nanotube":
        n = int(params.get("n", 5))
        m = int(params.get("m", 5))
        length = int(params.get("length", 1))
        # Note: ASE nanotube might return 1D or 3D pbc depending on version, 
        # usually 1D along z for 'nanotube' function
        atoms = nanotube(n, m, length=length, bond=params.get("bond", 1.42))
        
        # Nanotubes need vacuum only in x/y, z is periodic
        # But for 3D compatibility we might want to ensure a big cell
        atoms.center(vacuum=params.get("vacuum", 10.0), axis=(0, 1))

    elif structure_type == "graphene":
        # ASE Graphene expects an element symbol (e.g., 'C'), not a formula like 'C2'
        formula = str(params.get("formula", "C"))
        if formula == "C2": formula = "C" # Handle common user input for unit cell formula
        
        a = float(params.get("a", 2.46))
        c = float(params.get("c", 10.0)) # Vacuum separation
        size = [int(x) for x in params.get("size", [1, 1, 1])]
        
        # Graphene class creates a 3D cell with 'c' as separation
        atoms = Graphene(symbol=formula, latticeconstant={'a': a, 'c': c}, size=size)

    elif structure_type == "nanoribbon":
        w = int(params.get("width", 5))
        l = int(params.get("length", 5))
        ribbon_type = str(params.get("type", "armchair"))
        saturated = bool(params.get("saturated", True))
        vacuum = float(params.get("vacuum", 10.0))
        
        # graphene_nanoribbon creates structure periodic in one direction
        atoms = graphene_nanoribbon(w, l, type=ribbon_type, saturated=saturated, vacuum=vacuum)
        
        # Ensure vacuum in non-periodic directions (usually y, z for ribbon where x is periodic)
        # However, check pbc to be sure. ASE docs say periodic along Z usually?
        # Actually standard graphene_nanoribbon puts it in XY plane, periodic along X (armchair) or Y?
        # Let's trust .center(vacuum) to adding padding where PBC is false
        if not all(atoms.pbc):
             atoms.center(vacuum=vacuum)

    elif structure_type == "fullerene":
        name = str(params.get("name", "C60"))
        
        # Check if molecule exists in database before calling molecule()
        if name not in g2.names and name not in extra:
            return {
                "success": False, 
                "error": f"Molecule '{name}' not found in ASE database. Try 'C60'."
            }
            
        atoms = molecule(name)
        # 0D object, needs full 3D vacuum box
        atoms.center(vacuum=5.0)

    elif structure_type == "mos2":
        formula = str(params.get("formula", "MoS2"))
        kind = str(params.get("kind", "2H"))
        a = float(params.get("a", 3.16))
        thickness = float(params.get("thickness", 3.19))
        vacuum = float(params.get("vacuum", 10.0))
        size = [int(x) for x in params.get("size", [1, 1, 1])]
        
        atoms = mx2(formula=formula, kind=kind, a=a, thickness=thickness, size=size, vacuum=vacuum)

    elif structure_type == "nanowire":
        # Carve a nanowire from bulk
        formula = str(params.get("formula", "Au"))
        radius = float(params.get("radius", 5.0))
        length_cells = int(params.get("length", 4))
        vacuum = float(params.get("vacuum", 10.0))
        
        # 1. Create bulk
        # Check if formula is supported by 'bulk' (simple elements/compounds)
        # We assume valid input or let it raise error if formula is unknown.
        atoms_bulk = bulk(formula, cubic=True) 
            
        # 2. Supercell
        # Estimate size needed for radius
        # Volume = (2*r)^2 * L
        # We need a box of at least 2*radius in x/y
        # Get cell dimensions
        cell = atoms_bulk.get_cell_lengths_and_angles()[:3]
        nx = int(np.ceil(2 * radius / cell[0])) + 2
        ny = int(np.ceil(2 * radius / cell[1])) + 2
        nz = length_cells
        
        atoms = atoms_bulk * (nx, ny, nz)
        atoms.center()
        
        # 3. Carve cylinder
        # Center in XY plane is at cell center? 
        # Actually .center() puts it in middle of cell.
        # We want to keep atoms within radius of center axis (Z)
        
        com = atoms.get_center_of_mass()
        center_x, center_y = com[0], com[1]
        
        # Filter atoms
        # positions = atoms.get_positions()
        # mask = (x-cx)^2 + (y-cy)^2 <= r^2
        keep_indices = []
        for i, pos in enumerate(atoms.get_positions()):
            r_sq = (pos[0] - center_x)**2 + (pos[1] - center_y)**2
            if r_sq <= radius**2:
                keep_indices.append(i)
        
        if not keep_indices:
             return {
                "success": False,
                "error": f"Radius {radius} is too small, no atoms remained."
            }
            
        atoms = atoms[keep_indices]
        
        # 4. Set vacuum
        # Make it 1D along Z? 
        # ASE standard for wires is periodic in Z, vacuum in X/Y
        atoms.set_pbc([False, False, True])
        atoms.center(vacuum=vacuum, axis=(0, 1))

    else:
        return {
            "success": False, 
            "error": f"Unknown nanostructure type: {structure_type}"
        }

    if atoms is None:
        return {
            "success": False,
            "error": "Failed to generate structure (atoms object is None)"
        }

    # 2. Ensure fully 3D cell for compatibility
    # If cell is missing 3 vectors (e.g. 1D), make it 3D
    # .center() usually upgrades the cell to 3D if vacuum is provided
    # but we must be explicit
    
    if len(atoms.get_cell()) < 3:
         atoms.center(vacuum=10.0) # Force 3D box if somehow still lower dim
         
    # Force full PBC for standard visualizers if it's a "simulated crystal"
    # Or keep it as is?
    # For compatibility with tool "extract_structure_data", we usually expect a defined volume.
    # atoms.get_volume() raises error if not fully 3D.
    
    # Check dimensions
    # Explicitly check if cell volume is non-zero/definable
    cell = atoms.get_cell()
    if cell.shape != (3, 3) or np.abs(np.linalg.det(cell)) < 1e-6:
        # If we still don't have a volume, force a generic cubic box fitting the atoms
        atoms.center(vacuum=10.0)
    
    return {
        "success": True, 
        "structure": atoms_to_dict(atoms)
    }

def atoms_to_dict(atoms) -> Dict[str, Any]:
    """Convert ASE atoms to dictionary format."""
    cell = atoms.get_cell()
    
    # Ensure 3x3 matrix
    if cell.shape != (3, 3):
        # Build a 3x3 cell from what we have or identity
        # This shouldn't happen if we ran center(vacuum) correctly
        new_cell = np.eye(3) * (np.max(atoms.positions) - np.min(atoms.positions) + 10.0)
        atoms.set_cell(new_cell)
        atoms.center()
        cell = atoms.get_cell()

    # Calculate lattice params manually to be safe
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    # Angles
    # Very basic angle calc if needed, or assume orthogonal if simple box
    # For general case:
    def angle(v1, v2):
        return np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
        
    alpha = angle(cell[1], cell[2])
    beta = angle(cell[0], cell[2])
    gamma = angle(cell[0], cell[1])

    return {
        "lattice": {
            "matrix": cell.tolist(),
            "volume": float(atoms.get_volume()),
            "a": float(a),
            "b": float(b),
            "c": float(c),
            "alpha": float(alpha),
            "beta": float(beta),
            "gamma": float(gamma)
        },
        "atoms": [
            {
                "element": atom.symbol,
                "coords": atoms.get_scaled_positions()[i].tolist(),
                "cartesian": atom.position.tolist()
            }
            for i, atom in enumerate(atoms)
        ],
        "metadata": {
            "formula": atoms.get_chemical_formula(),
            "natoms": len(atoms),
            "pbc": atoms.pbc.tolist()
        }
    }

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "No input file provided"}))
        sys.exit(1)
        
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)
            
    # We deliberately let Python errors crash the script if they occur outside our logic,
    # or handle them if they are expected.
    # But as per request "i hate try: except", we avoid wrapping the whole block.
    # However, for the JSON-RPC contract, we must output JSON. 
    # If the logic above crashes, the MCP wrapper will catch the stderr/exit code.
    
    result = generate_nanostructure(data.get("type"), data.get("params", {}))
    print(json.dumps(result))

if __name__ == "__main__":
    main()
