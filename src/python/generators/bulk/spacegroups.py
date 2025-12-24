"""
bulk/spacegroups.py - Space Group Based Structure Generation using PyXtal

Reliable generation of structures in all 230 space groups using PyXtal.
Per project requirements: Create any material in any of 230 space groups.
"""

from typing import Dict, Any, List, Optional, Union
import importlib.util
import numpy as np
from pymatgen.core import Structure, Lattice
from .base import validate_space_group

# Import PyXtal for space group generation
PYXTAL_AVAILABLE = importlib.util.find_spec("pyxtal") is not None
if PYXTAL_AVAILABLE:
    from pyxtal import pyxtal
    from pyxtal.symmetry import Group
else:
    pyxtal = None
    Group = None


# Crystal system to space group range mapping
CRYSTAL_SYSTEMS = {
    "triclinic": {"range": (1, 2), "min_sg": 1, "max_sg": 2},
    "monoclinic": {"range": (3, 15), "min_sg": 3, "max_sg": 15},
    "orthorhombic": {"range": (16, 74), "min_sg": 16, "max_sg": 74},
    "tetragonal": {"range": (75, 142), "min_sg": 75, "max_sg": 142},
    "trigonal": {"range": (143, 167), "min_sg": 143, "max_sg": 167},
    "hexagonal": {"range": (168, 194), "min_sg": 168, "max_sg": 194},
    "cubic": {"range": (195, 230), "min_sg": 195, "max_sg": 230},
}


# Common space groups with examples
COMMON_SPACEGROUPS = {
    # Cubic
    225: {"symbol": "Fm-3m", "system": "cubic", "examples": ["NaCl", "Cu", "Au", "Ag", "MgO"]},
    227: {"symbol": "Fd-3m", "system": "cubic", "examples": ["diamond", "Si", "Ge", "spinel"]},
    221: {"symbol": "Pm-3m", "system": "cubic", "examples": ["CsCl", "SrTiO3", "BaTiO3"]},
    229: {"symbol": "Im-3m", "system": "cubic", "examples": ["bcc Fe", "W", "Mo", "Cr", "V"]},
    216: {"symbol": "F-43m", "system": "cubic", "examples": ["zincblende GaAs", "ZnS", "CdTe"]},
    223: {"symbol": "Pm-3n", "system": "cubic", "examples": ["A15 Nb3Sn", "Cr3Si"]},
    
    # Hexagonal
    194: {"symbol": "P6_3/mmc", "system": "hexagonal", "examples": ["hcp Mg", "graphite", "MoS2"]},
    191: {"symbol": "P6/mmm", "system": "hexagonal", "examples": ["AlB2", "MgB2"]},
    186: {"symbol": "P6_3mc", "system": "hexagonal", "examples": ["wurtzite ZnO", "GaN", "AlN"]},
    164: {"symbol": "P-3m1", "system": "trigonal", "examples": ["CdI2", "1T-MoS2"]},
    
    # Trigonal
    166: {"symbol": "R-3m", "system": "trigonal", "examples": ["Bi2Te3", "Bi2Se3", "calcite"]},
    167: {"symbol": "R-3c", "system": "trigonal", "examples": ["corundum Al2O3", "Fe2O3"]},
    
    # Tetragonal
    139: {"symbol": "I4/mmm", "system": "tetragonal", "examples": ["ThCr2Si2", "CaC2"]},
    140: {"symbol": "I4/mcm", "system": "tetragonal", "examples": ["BaFe2As2"]},
    123: {"symbol": "P4/mmm", "system": "tetragonal", "examples": ["CuO2 planes", "cuprates"]},
    136: {"symbol": "P4_2/mnm", "system": "tetragonal", "examples": ["rutile TiO2", "SnO2"]},
    
    # Orthorhombic
    62: {"symbol": "Pnma", "system": "orthorhombic", "examples": ["perovskite GdFeO3", "SnSe"]},
    63: {"symbol": "Cmcm", "system": "orthorhombic", "examples": ["CrB", "TiB"]},
    
    # Monoclinic
    14: {"symbol": "P2_1/c", "system": "monoclinic", "examples": ["molecular crystals"]},
    15: {"symbol": "C2/c", "system": "monoclinic", "examples": ["pyroxene"]},
    12: {"symbol": "C2/m", "system": "monoclinic", "examples": ["layered LiCoO2"]},
    
    # Triclinic
    2: {"symbol": "P-1", "system": "triclinic", "examples": ["low symmetry phases"]},
    1: {"symbol": "P1", "system": "triclinic", "examples": ["no symmetry"]},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    """Convert pymatgen Structure to dict with both 'sites' and 'atoms' keys."""
    lattice = structure.lattice
    sites_data = []
    for s in structure:
        site_info = {
            "element": str(s.specie),
            "coords": list(s.frac_coords),
            "species": [{"element": str(s.specie), "occupation": 1.0}]  # For test compatibility
        }
        sites_data.append(site_info)
    return {
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "matrix": lattice.matrix.tolist()
        },
        "sites": sites_data,  # Standard key expected by tests
        "atoms": sites_data,  # Alias for user convenience
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def _coerce_spacegroup(spacegroup: Any) -> Optional[int]:
    """Convert a space group value to int if it is integer-like."""
    if isinstance(spacegroup, bool):
        return None
    if isinstance(spacegroup, (int, np.integer)):
        return int(spacegroup)
    if isinstance(spacegroup, str):
        cleaned = spacegroup.strip()
        if not cleaned:
            return None
        sign = ""
        if cleaned[0] in "+-":
            sign = cleaned[0]
            cleaned = cleaned[1:]
        if cleaned.isdigit():
            return int(f"{sign}{cleaned}")
    return None


def generate_from_spacegroup(
    spacegroup: int,
    elements: List[str],
    composition: Optional[List[int]] = None,
    a: float = 5.0,
    b: Optional[float] = None,
    c: Optional[float] = None,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
    factor: float = 1.0,
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate structure in any of 230 space groups using PyXtal.
    
    Args:
        spacegroup: Space group number (1-230)
        elements: List of element symbols, e.g., ["Na", "Cl"]
        composition: Number of each element, e.g., [4, 4] for 4 Na and 4 Cl.
                    If None, defaults to [1] * len(elements)
        a, b, c: Lattice parameters in Angstrom
        alpha, beta, gamma: Lattice angles in degrees
        factor: Volume scaling factor
    
    Returns:
        Structure in the specified space group
    """
    if not PYXTAL_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "PYXTAL_NOT_INSTALLED", 
                      "message": "PyXtal is required for space group generation. Install with: pip install pyxtal"}
        }
    
    sg_number = _coerce_spacegroup(spacegroup)
    if sg_number is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SPACEGROUP_TYPE",
                      "message": f"Space group must be a number, got {type(spacegroup)}"}
        }
    spacegroup = sg_number

    if not 1 <= spacegroup <= 230:
        return {
            "success": False,
            "error": {"code": "INVALID_SPACEGROUP", 
                      "message": f"Space group must be 1-230, got {spacegroup}"}
        }
    
    if composition is None:
        composition = [1] * len(elements)
    
    if len(elements) != len(composition):
        return {
            "success": False,
            "error": {"code": "LENGTH_MISMATCH", 
                      "message": "elements and composition must have same length"}
        }
    
    # Validate lattice parameters are positive
    if a <= 0 or (b is not None and b <= 0) or (c is not None and c <= 0):
        return {
            "success": False,
            "error": {"code": "INVALID_LATTICE_PARAMETER",
                      "message": f"Lattice parameters must be positive. Got a={a}, b={b}, c={c}"}
        }
    
    # Validate element symbols using pymatgen's valid symbol list
    from pymatgen.core.periodic_table import Element
    valid_symbols = {e.symbol for e in Element}
    for elem in elements:
        if elem not in valid_symbols:
            return {
                "success": False,
                "error": {"code": "INVALID_ELEMENT",
                          "message": f"Invalid element symbol: '{elem}'. Must be a valid periodic table element."}
            }
    
    # Determine crystal system from space group
    crystal_system = None
    for system, info in CRYSTAL_SYSTEMS.items():
        if info["min_sg"] <= spacegroup <= info["max_sg"]:
            crystal_system = system
            break
    
    # Generate structure using PyXtal
    crystal = pyxtal()
    
    # Set lattice parameters based on crystal system
    if b is None:
        b = a
    if c is None:
        c = a
    
    # Pre-validate composition compatibility with space group
    group = Group(spacegroup)

    # Check if composition is compatible with available Wyckoff positions
    available_multiplicities = sorted(set([wp.multiplicity for wp in group.Wyckoff_positions]))
    min_multiplicity = min(available_multiplicities)

    # Adjust composition if needed to fit Wyckoff constraints
    adjusted_composition = list(composition)
    composition_adjusted = False

    for i, (elem, count) in enumerate(zip(elements, composition)):
        if count < min_multiplicity and count != 0:
            # Find the smallest valid multiplicity >= count
            valid_count = min(m for m in available_multiplicities if m >= count)
            adjusted_composition[i] = valid_count
            composition_adjusted = True
        elif count > 0:
            # Check if count is achievable with available multiplicities
            # Simple check: can we sum to this count using available multiplicities?
            achievable = False
            for mult in available_multiplicities:
                if count % mult == 0 or count in available_multiplicities:
                    achievable = True
                    break
                # Check if sum of multiplicities can reach count
                for m1 in available_multiplicities:
                    for m2 in available_multiplicities:
                        if m1 + m2 == count or m1 * 2 == count or m2 * 2 == count:
                            achievable = True
                            break
                    if achievable:
                        break
            if not achievable:
                # Find nearest achievable count
                valid_counts = set()
                for m in available_multiplicities:
                    valid_counts.add(m)
                    valid_counts.add(m * 2)
                    for m2 in available_multiplicities:
                        valid_counts.add(m + m2)
                nearest = min(valid_counts, key=lambda x: abs(x - count))
                adjusted_composition[i] = nearest
                composition_adjusted = True

    composition = adjusted_composition
    
    # Generate the crystal - PyXtal will find valid Wyckoff assignments
    # Set seed for reproducibility if provided
    if seed is not None:
        np.random.seed(seed)
    
    crystal.from_random(
        dim=3,
        group=spacegroup,
        species=elements,
        numIons=composition,
        factor=factor,
    )
    
    # Convert to pymatgen Structure
    pmg_structure = crystal.to_pymatgen()
    
    # Apply requested lattice parameters if provided
    # This ensures scientific correctness for specific materials (e.g. Si a=5.43)
    if a is not None:
        target_a = a
        target_b = b if b is not None else a
        target_c = c if c is not None else a
        
        # Keep original angles unless specified (TODO: add angle support)
        # For now, we assume the generated angles are correct for the SG
        new_lattice = Lattice.from_parameters(
            a=target_a,
            b=target_b,
            c=target_c,
            alpha=pmg_structure.lattice.alpha,
            beta=pmg_structure.lattice.beta,
            gamma=pmg_structure.lattice.gamma
        )
        
        # Create new structure with target lattice but same fractional coords
        pmg_structure = Structure(
            new_lattice,
            pmg_structure.species,
            pmg_structure.frac_coords
        )
    
    # Get space group info
    sg_info = COMMON_SPACEGROUPS.get(spacegroup, {})
    
    return {
        "success": True,
        "spacegroup_number": spacegroup,
        "spacegroup_symbol": crystal.group.symbol,
        "crystal_system": crystal_system,
        "elements": elements,
        "composition": composition,
        "formula": pmg_structure.formula,
        "n_atoms": len(pmg_structure),
        "lattice_parameters": {
            "a": pmg_structure.lattice.a,
            "b": pmg_structure.lattice.b,
            "c": pmg_structure.lattice.c,
            "alpha": pmg_structure.lattice.alpha,
            "beta": pmg_structure.lattice.beta,
            "gamma": pmg_structure.lattice.gamma,
        },
        "wyckoff_positions": [str(site.wp) for site in crystal.atom_sites] if crystal.atom_sites else [],
        "structure": structure_to_dict(pmg_structure)
    }


def generate_from_spacegroup_symbol(
    symbol: str,
    elements: List[str],
    composition: Optional[List[int]] = None,
    factor: float = 1.0
) -> Dict[str, Any]:
    """
    Generate structure from space group symbol (e.g., "Fm-3m", "P6_3/mmc").
    
    Args:
        symbol: Space group symbol
        elements: List of element symbols
        composition: Number of each element
        factor: Volume scaling factor
    
    Returns:
        Structure in the specified space group
    """
    if not PYXTAL_AVAILABLE:
        return {
            "success": False,
            "error": {"code": "PYXTAL_NOT_INSTALLED", 
                      "message": "PyXtal is required. Install with: pip install pyxtal"}
        }
    
    # Find space group number from symbol
    sg_number = None
    for num, info in COMMON_SPACEGROUPS.items():
        if info["symbol"] == symbol:
            sg_number = num
            break
    
    if sg_number is None:
        is_valid, error, sg_number = validate_space_group(symbol)
        if not is_valid:
            return {
                "success": False,
                "error": {"code": "INVALID_SYMBOL",
                          "message": error}
            }
    
    return generate_from_spacegroup(sg_number, elements, composition, factor=factor)


def list_spacegroups_for_system(crystal_system: str) -> Dict[str, Any]:
    """
    List all space groups for a crystal system.
    
    Args:
        crystal_system: One of triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic
    
    Returns:
        List of space groups in the system
    """
    if crystal_system.lower() not in CRYSTAL_SYSTEMS:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", 
                      "message": f"Unknown crystal system: {crystal_system}",
                      "available": list(CRYSTAL_SYSTEMS.keys())}
        }
    
    info = CRYSTAL_SYSTEMS[crystal_system.lower()]
    min_sg, max_sg = info["min_sg"], info["max_sg"]
    
    spacegroups = []
    for sg in range(min_sg, max_sg + 1):
        sg_info = COMMON_SPACEGROUPS.get(sg, {})
        spacegroups.append({
            "number": sg,
            "symbol": sg_info.get("symbol", ""),
            "examples": sg_info.get("examples", [])
        })
    
    return {
        "success": True,
        "crystal_system": crystal_system,
        "n_spacegroups": max_sg - min_sg + 1,
        "range": [min_sg, max_sg],
        "common_spacegroups": [sg for sg in spacegroups if sg["symbol"]]
    }


def get_spacegroup_info(spacegroup: int) -> Dict[str, Any]:
    """
    Get detailed information about a space group.
    
    Args:
        spacegroup: Space group number (1-230)
    
    Returns:
        Space group information
    """
    if not 1 <= spacegroup <= 230:
        return {
            "success": False,
            "error": {"code": "INVALID_SPACEGROUP", "message": f"Invalid space group: {spacegroup}"}
        }
    
    # Determine crystal system
    crystal_system = None
    for system, info in CRYSTAL_SYSTEMS.items():
        if info["min_sg"] <= spacegroup <= info["max_sg"]:
            crystal_system = system
            break
    
    sg_info = COMMON_SPACEGROUPS.get(spacegroup, {})
    
    # Try to get more info from PyXtal
    wyckoff_positions = []
    if PYXTAL_AVAILABLE:
        group = Group(spacegroup)
        wyckoff_positions = [
            {"letter": wp.letter, "multiplicity": wp.multiplicity}
            for wp in group.Wyckoff_positions
        ]
    
    return {
        "success": True,
        "spacegroup_number": spacegroup,
        "symbol": sg_info.get("symbol", f"SG {spacegroup}"),
        "crystal_system": crystal_system,
        "example_materials": sg_info.get("examples", []),
        "wyckoff_positions": wyckoff_positions[:10],  # First 10
        "pyxtal_available": PYXTAL_AVAILABLE
    }


def generate_all_cubic_prototypes(element: str = "Si") -> Dict[str, Any]:
    """
    Generate all common cubic structure types for an element.
    
    Args:
        element: Element symbol
    
    Returns:
        Dictionary of cubic structures
    """
    if not PYXTAL_AVAILABLE:
        return {"success": False, "error": {"code": "PYXTAL_NOT_INSTALLED", "message": "PyXtal required"}}
    
    prototypes = {
        "FCC": {"sg": 225, "composition": [4]},
        "BCC": {"sg": 229, "composition": [2]},
        "SC": {"sg": 221, "composition": [1]},
        "diamond": {"sg": 227, "composition": [8]},
    }
    
    results = {}
    for name, info in prototypes.items():
        results[name] = generate_from_spacegroup(info["sg"], [element], info["composition"])
    
    return {
        "success": True,
        "element": element,
        "prototypes": results
    }


# Compatibility aliases
def generate_prototype(
    prototype: str,
    elements: Dict[str, str],
    a: float = 5.0,
    c: Optional[float] = None
) -> Dict[str, Any]:
    """
    Generate common prototype structure.
    
    Args:
        prototype: Prototype name (rocksalt, perovskite, zincblende, etc.)
        elements: Mapping of sites to elements
        a: Lattice parameter
        c: c lattice parameter (for non-cubic)
    
    Returns:
        Prototype structure
    """
    prototypes = {
        "rocksalt": {"sg": 225, "elements_order": ["cation", "anion"], "composition": [4, 4]},
        "perovskite": {"sg": 221, "elements_order": ["A", "B", "O"], "composition": [1, 1, 3]},
        "zincblende": {"sg": 216, "elements_order": ["cation", "anion"], "composition": [4, 4]},
        "wurtzite": {"sg": 186, "elements_order": ["cation", "anion"], "composition": [2, 2]},
        "rutile": {"sg": 136, "elements_order": ["M", "O"], "composition": [2, 4]},
        "fluorite": {"sg": 225, "elements_order": ["cation", "anion"], "composition": [4, 8]},
        "spinel": {"sg": 227, "elements_order": ["A", "B", "O"], "composition": [8, 16, 32]},
        "bcc": {"sg": 229, "elements_order": ["element"], "composition": [2]},
        "fcc": {"sg": 225, "elements_order": ["element"], "composition": [4]},
        "hcp": {"sg": 194, "elements_order": ["element"], "composition": [2]},
        "diamond": {"sg": 227, "elements_order": ["element"], "composition": [8]},
    }
    
    if prototype.lower() not in prototypes:
        return {
            "success": False,
            "error": {"code": "INVALID_PROTOTYPE", 
                      "message": f"Unknown prototype: {prototype}",
                      "available": list(prototypes.keys())}
        }
    
    proto_info = prototypes[prototype.lower()]
    
    # Map elements
    element_list = []
    for key in proto_info["elements_order"]:
        if key in elements:
            element_list.append(elements[key])
        elif len(elements) == 1:
            element_list.append(list(elements.values())[0])
        else:
            return {
                "success": False,
                "error": {"code": "MISSING_ELEMENT", 
                          "message": f"Missing element for site '{key}'"}
            }
    
    result = generate_from_spacegroup(
        proto_info["sg"], 
        element_list, 
        proto_info["composition"]
    )
    
    if result["success"]:
        result["prototype"] = prototype
    
    return result
