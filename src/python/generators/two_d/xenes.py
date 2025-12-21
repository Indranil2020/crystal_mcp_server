"""
two_d/xenes.py - X-ene 2D Materials (Graphene Analogs)

Generates X-ene structures:
- Graphene (flat honeycomb)
- Silicene (buckled honeycomb, 0.44 Å)
- Germanene (buckled honeycomb, 0.66 Å)
- Stanene (buckled honeycomb, 0.85 Å)
- Plumbene (buckled honeycomb, 0.93 Å)
- Borophene (various polymorphs: α, β12, χ3)
- Phosphorene (puckered orthorhombic)
- Arsenene, antimonene, bismuthene

Scientific basis:
- Experimental lattice constants
- DFT-optimized buckling heights
- ASE graphene_nanoribbon for ribbon generation
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from pymatgen.core import Structure, Lattice
from ase.build import graphene

from .base import structure_to_dict, atoms_to_dict


# X-ene parameters from DFT and experiments
XENE_PARAMS = {
    "graphene": {
        "element": "C",
        "a": 2.46,  # Å
        "buckling": 0.0,
        "bond_length": 1.42,
        "description": "Flat sp2 honeycomb"
    },
    "silicene": {
        "element": "Si",
        "a": 3.86,
        "buckling": 0.44,
        "bond_length": 2.28,
        "description": "Buckled sp3/sp2 hybrid"
    },
    "germanene": {
        "element": "Ge",
        "a": 4.02,
        "buckling": 0.66,
        "bond_length": 2.38,
        "description": "Buckled honeycomb"
    },
    "stanene": {
        "element": "Sn",
        "a": 4.67,
        "buckling": 0.85,
        "bond_length": 2.83,
        "description": "Buckled honeycomb, topological insulator"
    },
    "plumbene": {
        "element": "Pb",
        "a": 4.93,
        "buckling": 0.93,
        "bond_length": 3.00,
        "description": "Buckled honeycomb"
    },
    "borophene_beta12": {
        "element": "B",
        "a": 2.87,
        "b": 5.07,
        "buckling": 0.0,
        "puckered": False,
        "description": "β12 phase with vacancy pattern"
    },
    "borophene_chi3": {
        "element": "B",
        "a": 4.45,
        "b": 4.45,
        "buckling": 0.0,
        "description": "χ3 triangular phase"
    },
    "phosphorene": {
        "element": "P",
        "a": 3.31,
        "b": 4.38,
        "buckling": 2.10,
        "puckered": True,
        "description": "Puckered orthorhombic (black phosphorus monolayer)"
    },
    "arsenene": {
        "element": "As",
        "a": 3.61,
        "buckling": 1.40,
        "description": "Buckled honeycomb"
    },
    "antimonene": {
        "element": "Sb",
        "a": 4.12,
        "buckling": 1.65,
        "description": "Buckled honeycomb"
    },
    "bismuthene": {
        "element": "Bi",
        "a": 4.33,
        "buckling": 1.73,
        "description": "Buckled honeycomb, strong SOC"
    }
}


def generate_xene(
    xene_type: str,
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0,
    passivate: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate X-ene monolayer structure.
    
    Args:
        xene_type: Type of X-ene (graphene, silicene, germanene, etc.)
        size: Supercell size [nx, ny, 1]
        vacuum: Vacuum padding in Angstroms
        passivate: Element for surface passivation (H, F, OH)
    
    Returns:
        X-ene structure dictionary
    
    Examples:
        >>> result = generate_xene('silicene')
        >>> result["buckling"]
        0.44
    """
    xene_lower = xene_type.lower().replace("-", "_").replace(" ", "_")
    
    if xene_lower not in XENE_PARAMS:
        return {
            "success": False,
            "error": {
                "code": "INVALID_XENE",
                "message": f"Unknown X-ene '{xene_type}'",
                "available": list(XENE_PARAMS.keys())
            }
        }
    
    params = XENE_PARAMS[xene_lower]
    element = params["element"]
    a = params["a"]
    buckling = params.get("buckling", 0.0)
    is_puckered = params.get("puckered", False)
    
    if is_puckered:
        # Phosphorene-type puckered structure
        b = params.get("b", a * 1.3)
        c = vacuum
        lattice = Lattice.orthorhombic(a, b, c)
        
        # 4 atoms in unit cell (ridge-trough pattern)
        species = [element] * 4
        dz = buckling / c / 2
        coords = [
            [0.0, 0.08, 0.5 - dz],
            [0.0, 0.42, 0.5 + dz],
            [0.5, 0.58, 0.5 + dz],
            [0.5, 0.92, 0.5 - dz],
        ]
    elif "borophene" in xene_lower:
        # Borophene structures
        b = params.get("b", a)
        c = vacuum
        if "beta12" in xene_lower:
            lattice = Lattice.orthorhombic(a, b, c)
            # β12 has 5 atoms with hexagonal hole pattern
            species = [element] * 5
            coords = [
                [0.0, 0.0, 0.5], [0.33, 0.17, 0.5], [0.67, 0.17, 0.5],
                [0.33, 0.5, 0.5], [0.67, 0.5, 0.5]
            ]
        else:
            lattice = Lattice.hexagonal(a, c)
            species = [element] * 3
            coords = [[0, 0, 0.5], [1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]
    else:
        # Standard honeycomb (buckled or flat)
        c = vacuum
        lattice = Lattice.hexagonal(a, c)
        dz = buckling / c / 2 if c > 0 else 0
        species = [element] * 2
        coords = [
            [1/3, 2/3, 0.5 - dz],
            [2/3, 1/3, 0.5 + dz],
        ]
    
    structure = Structure(lattice, species, coords)
    
    # Apply supercell
    if size != [1, 1, 1]:
        structure.make_supercell(size)
    
    # Add passivation
    if passivate:
        passivation_z = buckling / vacuum + 0.02 if vacuum > 0 else 0.02
        for site in list(structure):
            z = site.frac_coords[2]
            # Add passivant above and below
            structure.append(passivate, [site.frac_coords[0], site.frac_coords[1], z + passivation_z])
            structure.append(passivate, [site.frac_coords[0], site.frac_coords[1], z - passivation_z])
    
    return {
        "success": True,
        "xene_type": xene_type,
        "element": element,
        "buckling_angstrom": buckling,
        "bond_length_angstrom": params.get("bond_length", a / np.sqrt(3)),
        "description": params["description"],
        "passivated": passivate is not None,
        "structure": structure_to_dict(structure, vacuum)
    }


def generate_graphene(
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0,
    stacking: str = "AA",
    n_layers: int = 1
) -> Dict[str, Any]:
    """
    Generate graphene using ASE for maximum accuracy.
    
    Args:
        size: Supercell size
        vacuum: Vacuum padding
        stacking: Stacking for multilayer (AA, AB, ABC)
        n_layers: Number of layers
    
    Returns:
        Graphene structure
    """
    # Use ASE's graphene builder
    atoms = graphene(vacuum=vacuum)
    
    if size[0] > 1 or size[1] > 1:
        atoms = atoms.repeat((size[0], size[1], 1))
    
    if n_layers > 1:
        # Add additional layers
        layer_sep = 3.35  # Å
        original = atoms.copy()
        for i in range(1, n_layers):
            layer = original.copy()
            positions = layer.get_positions()
            if stacking == "AB" and i % 2 == 1:
                # Shift by a/sqrt(3) in x
                positions[:, 0] += 1.42
            elif stacking == "ABC":
                shift = (i % 3) * 1.42
                positions[:, 0] += shift
            positions[:, 2] += i * layer_sep
            layer.set_positions(positions)
            atoms += layer
        
        # Update cell
        cell = atoms.get_cell()
        cell[2, 2] = n_layers * layer_sep + vacuum
        atoms.set_cell(cell)
    
    return {
        "success": True,
        "material": "graphene",
        "n_layers": n_layers,
        "stacking": stacking,
        "structure": atoms_to_dict(atoms, vacuum)
    }


def generate_silicene(
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0,
    passivate: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate silicene monolayer structure.

    Silicene is the silicon analogue of graphene with buckled honeycomb structure.

    Args:
        size: Supercell size [nx, ny, 1]
        vacuum: Vacuum padding in Angstroms
        passivate: Element for surface passivation (H, F)

    Returns:
        Silicene structure dictionary
    """
    result = generate_xene("silicene", size, vacuum, passivate)
    if result.get("success"):
        result["material"] = "silicene"
        result["is_buckled"] = True
        result["buckling_angstrom"] = 0.44
    return result


def generate_germanene(
    size: List[int] = [1, 1, 1],
    vacuum: float = 15.0,
    passivate: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate germanene monolayer structure.

    Germanene is the germanium analogue of graphene with buckled honeycomb structure.

    Args:
        size: Supercell size [nx, ny, 1]
        vacuum: Vacuum padding in Angstroms
        passivate: Element for surface passivation (H, F)

    Returns:
        Germanene structure dictionary
    """
    result = generate_xene("germanene", size, vacuum, passivate)
    if result.get("success"):
        result["material"] = "germanene"
        result["is_buckled"] = True
        result["buckling_angstrom"] = 0.66
    return result
