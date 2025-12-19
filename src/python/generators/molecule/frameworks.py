"""
molecule/frameworks.py - Metal-Organic and Covalent-Organic Frameworks

Comprehensive framework generation per structure_catalogue.md Category 5:
(iii) MOFs with reticular chemistry (MOF-5, ZIF-8, HKUST-1, UiO-66)
(iv) COFs – 2D and 3D with linker library
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure, Lattice


# MOF Database
MOF_DATABASE = {
    # IRMOFs (Isoreticular MOFs)
    "MOF-5": {
        "formula": "Zn4O(BDC)3", "topology": "pcu", "sbu": "Zn4O",
        "linker": "BDC", "a": 25.83, "pore_diameter_A": 11.2,
        "surface_area_m2g": 3800, "void_fraction": 0.79
    },
    "IRMOF-3": {
        "formula": "Zn4O(NH2-BDC)3", "topology": "pcu", "sbu": "Zn4O",
        "linker": "NH2-BDC", "a": 25.76, "pore_diameter_A": 10.5,
        "surface_area_m2g": 2446
    },
    "IRMOF-8": {
        "formula": "Zn4O(NDC)3", "topology": "pcu", "sbu": "Zn4O",
        "linker": "NDC", "a": 30.09, "pore_diameter_A": 18.0
    },
    
    # ZIFs (Zeolitic Imidazolate Frameworks)
    "ZIF-8": {
        "formula": "Zn(mIM)2", "topology": "sod", "sbu": "Zn",
        "linker": "2-methylimidazolate", "a": 16.99, "pore_diameter_A": 11.6,
        "surface_area_m2g": 1810, "aperture_A": 3.4
    },
    "ZIF-67": {
        "formula": "Co(mIM)2", "topology": "sod", "sbu": "Co",
        "linker": "2-methylimidazolate", "a": 16.96, "magnetic": True
    },
    "ZIF-7": {
        "formula": "Zn(bIM)2", "topology": "sod", "sbu": "Zn",
        "linker": "benzimidazolate", "flexible": True
    },
    
    # HKUST-1 (Cu-BTC)
    "HKUST-1": {
        "formula": "Cu3(BTC)2", "topology": "tbo", "sbu": "Cu2_paddle",
        "linker": "BTC", "a": 26.34, "pore_diameter_A": 13.2,
        "surface_area_m2g": 1850, "open_metal_sites": True
    },
    
    # UiO series (Zr-based, high stability)
    "UiO-66": {
        "formula": "Zr6O4(OH)4(BDC)6", "topology": "fcu", "sbu": "Zr6O8",
        "linker": "BDC", "a": 20.75, "pore_diameter_A": 6.0,
        "surface_area_m2g": 1187, "stability": "high"
    },
    "UiO-67": {
        "formula": "Zr6O4(OH)4(BPDC)6", "topology": "fcu", "sbu": "Zr6O8",
        "linker": "BPDC", "a": 26.88, "pore_diameter_A": 11.0,
        "surface_area_m2g": 2064
    },
    "UiO-68": {
        "formula": "Zr6O4(OH)4(TPDC)6", "topology": "fcu", "sbu": "Zr6O8",
        "linker": "TPDC", "a": 32.80, "pore_diameter_A": 16.0
    },
    
    # MIL series (Fe/Al-based)
    "MIL-53": {
        "formula": "Al(OH)(BDC)", "topology": "wine_rack", "sbu": "Al_chain",
        "linker": "BDC", "breathing": True, "flexible": True
    },
    "MIL-101": {
        "formula": "Cr3O(BDC)3", "topology": "mtn", "sbu": "Cr3O",
        "linker": "BDC", "a": 88.87, "pore_diameter_A": 34.0,
        "surface_area_m2g": 5900
    },
    "MIL-125": {
        "formula": "Ti8O8(OH)4(BDC)6", "topology": "bcu", "sbu": "Ti8O8",
        "linker": "BDC", "photocatalytic": True
    },
    
    # Other notable MOFs
    "PCN-222": {
        "formula": "Zr6(OH)8(TCPP)2", "topology": "csq", "sbu": "Zr6O8",
        "linker": "TCPP_porphyrin", "catalytic": True
    },
    "NU-1000": {
        "formula": "Zr6(OH)16(TBAPy)2", "topology": "csq", "sbu": "Zr6O8",
        "linker": "TBAPy_pyrene", "mesoporous": True
    },
}


# COF Database
COF_DATABASE = {
    # 2D COFs (layered)
    "COF-1": {
        "formula": "C3H2BO", "topology": "hcb", "dimension": "2D",
        "linkage": "boroxine", "a": 7.2, "pore_diameter_A": 9.0,
        "surface_area_m2g": 711
    },
    "COF-5": {
        "formula": "C9H4BO2", "topology": "hcb", "dimension": "2D",
        "linkage": "boronate_ester", "a": 16.9, "pore_diameter_A": 27.0,
        "surface_area_m2g": 1590
    },
    "CTF-1": {
        "formula": "C6N3", "topology": "hcb", "dimension": "2D",
        "linkage": "triazine", "a": 12.7, "semiconducting": True
    },
    "TpPa-1": {
        "formula": "C21H15N3O3", "topology": "sql", "dimension": "2D",
        "linkage": "beta_ketoenamine", "stable": True, "crystalline": True
    },
    
    # 3D COFs
    "COF-102": {
        "formula": "C25H24B4O8", "topology": "ctn", "dimension": "3D",
        "linkage": "boronate", "a": 27.2, "surface_area_m2g": 3472
    },
    "COF-103": {
        "formula": "C21H12B4O8Si", "topology": "ctn", "dimension": "3D",
        "linkage": "boronate", "a": 28.1, "surface_area_m2g": 4210
    },
    "COF-300": {
        "formula": "C45H24N6", "topology": "dia", "dimension": "3D",
        "linkage": "imine", "interpenetrated": True
    },
}


# Linker database
LINKER_DATABASE = {
    "BDC": {"name": "1,4-benzenedicarboxylic acid", "length_A": 7.0, "arms": 2},
    "BTC": {"name": "1,3,5-benzenetricarboxylic acid", "arms": 3},
    "NDC": {"name": "2,6-naphthalenedicarboxylic acid", "length_A": 9.5, "arms": 2},
    "BPDC": {"name": "4,4'-biphenyldicarboxylic acid", "length_A": 11.0, "arms": 2},
    "mIM": {"name": "2-methylimidazolate", "angle_deg": 145},
    "TCPP": {"name": "tetrakis(4-carboxyphenyl)porphyrin", "arms": 4, "metallocenter": True},
}


def structure_to_dict(structure: Structure) -> Dict[str, Any]:
    lattice = structure.lattice
    return {
        "lattice": {"a": lattice.a, "b": lattice.b, "c": lattice.c,
                    "matrix": lattice.matrix.tolist()},
        "atoms": [{"element": str(s.specie), "coords": list(s.frac_coords)} for s in structure],
        "metadata": {"formula": structure.formula, "n_atoms": len(structure)}
    }


def generate_mof(
    mof_name: str = "MOF-5",
    supercell: List[int] = [1, 1, 1],
    defective: bool = False,
    defect_fraction: float = 0.1
) -> Dict[str, Any]:
    """
    Generate MOF structure.
    
    Args:
        mof_name: MOF from database
        supercell: Supercell dimensions
        defective: Include missing linker defects
        defect_fraction: Fraction of linkers removed
    
    Returns:
        MOF structure
    """
    if mof_name not in MOF_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_MOF", "message": f"Unknown MOF '{mof_name}'",
                      "available": list(MOF_DATABASE.keys())}
        }
    
    info = MOF_DATABASE[mof_name]
    a = info.get("a", 20.0)
    
    lattice = Lattice.cubic(a)
    
    species = []
    coords = []
    
    # Generate simplified MOF structure based on topology
    if info["topology"] == "pcu":
        # MOF-5 type: Zn4O clusters at corners, linkers on edges
        sbu_positions = [[0, 0, 0], [0.5, 0.5, 0.5]]
        
        for pos in sbu_positions:
            # Zn4O cluster
            species.extend(["Zn", "Zn", "Zn", "Zn", "O"])
            coords.extend([
                [pos[0] + 0.02, pos[1], pos[2]],
                [pos[0] - 0.02, pos[1], pos[2]],
                [pos[0], pos[1] + 0.02, pos[2]],
                [pos[0], pos[1] - 0.02, pos[2]],
                pos
            ])
        
        # Simplified linkers
        for i in range(3):
            species.extend(["C", "C", "O", "O"])
            coords.extend([
                [0.25, 0, 0], [0.35, 0, 0], [0.2, 0, 0], [0.4, 0, 0]
            ])
            
    elif info["topology"] == "sod":
        # ZIF-8 type: sodalite topology
        n_zn = 24
        np.random.seed(42)
        for i in range(n_zn):
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
            r = 0.3
            x = 0.5 + r * np.sin(phi) * np.cos(theta)
            y = 0.5 + r * np.sin(phi) * np.sin(theta)
            z = 0.5 + r * np.cos(phi)
            species.append("Zn")
            coords.append([x % 1, y % 1, z % 1])
            
            # Imidazolate linker
            species.extend(["C", "N", "N"])
            coords.extend([
                [(x + 0.05) % 1, y % 1, z % 1],
                [(x + 0.07) % 1, (y + 0.02) % 1, z % 1],
                [(x + 0.07) % 1, (y - 0.02) % 1, z % 1]
            ])
    else:
        # Generic framework
        species = ["Zn", "O", "C", "C", "O"]
        coords = [[0, 0, 0], [0.1, 0.1, 0.1], [0.25, 0, 0], [0.75, 0, 0], [0.5, 0.5, 0.5]]
    
    structure = Structure(lattice, species, coords)
    
    if supercell != [1, 1, 1]:
        structure.make_supercell(supercell)
    
    return {
        "success": True,
        "mof_name": mof_name,
        "formula": info["formula"],
        "topology": info["topology"],
        "sbu": info["sbu"],
        "linker": info["linker"],
        "pore_diameter_A": info.get("pore_diameter_A", 0),
        "surface_area_m2g": info.get("surface_area_m2g", 0),
        "has_open_metal_sites": info.get("open_metal_sites", False),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }


def generate_cof(
    cof_name: str = "COF-5",
    n_layers: int = 3,
    stacking: str = "AA"
) -> Dict[str, Any]:
    """
    Generate COF structure.
    
    Args:
        cof_name: COF from database
        n_layers: Number of layers (for 2D COFs)
        stacking: 'AA', 'AB', or 'ABC' stacking
    
    Returns:
        COF structure
    """
    if cof_name not in COF_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_COF", "message": f"Unknown COF '{cof_name}'",
                      "available": list(COF_DATABASE.keys())}
        }
    
    info = COF_DATABASE[cof_name]
    a = info.get("a", 15.0)
    c = 3.4 * n_layers  # Interlayer spacing
    
    if info["dimension"] == "2D":
        lattice = Lattice.hexagonal(a, c)
    else:
        lattice = Lattice.cubic(a)
    
    species = []
    coords = []
    
    if info["topology"] == "hcb":
        # Honeycomb network for 2D COFs
        for layer in range(n_layers):
            z = (layer + 0.5) / n_layers
            
            # Offset based on stacking
            if stacking == "AA":
                offset = [0, 0]
            elif stacking == "AB":
                offset = [1/3 * (layer % 2), 1/3 * (layer % 2)]
            else:  # ABC
                offset = [1/3 * (layer % 3), 1/3 * (layer % 3)]
            
            # Hexagonal pore
            for node in range(6):
                angle = node * np.pi / 3 + np.pi / 6
                r = 0.2
                x = 0.5 + r * np.cos(angle) + offset[0]
                y = 0.5 + r * np.sin(angle) + offset[1]
                
                if "B" in info["linkage"]:
                    species.extend(["B", "O", "C"])
                else:
                    species.extend(["C", "N", "C"])
                coords.extend([
                    [x % 1, y % 1, z],
                    [(x + 0.05) % 1, y % 1, z],
                    [(x + 0.05) % 1, (y + 0.03) % 1, z]
                ])
    else:
        # 3D or other topology
        species = ["C", "B", "O", "C", "C"]
        coords = [[0, 0, 0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.3, 0.3, 0.3], [0.4, 0.4, 0.4]]
    
    structure = Structure(lattice, species, coords)
    
    return {
        "success": True,
        "cof_name": cof_name,
        "formula": info["formula"],
        "topology": info["topology"],
        "dimension": info["dimension"],
        "linkage": info["linkage"],
        "n_layers": n_layers if info["dimension"] == "2D" else 1,
        "stacking": stacking if info["dimension"] == "2D" else "N/A",
        "pore_diameter_A": info.get("pore_diameter_A", 0),
        "n_atoms": len(structure),
        "structure": structure_to_dict(structure)
    }
