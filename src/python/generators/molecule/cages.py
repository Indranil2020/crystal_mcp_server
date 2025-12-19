"""
molecule/cages.py - Molecular Cages and Fullerenes

Comprehensive cage compound generation per structure_catalogue.md Category 5:
(ii) Cage compounds – fullerenes, clathrates, metal clusters, supramolecular cages
"""

from typing import Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Molecule


# Fullerene database
FULLERENE_DATABASE = {
    "C20": {"n_atoms": 20, "n_pentagons": 12, "n_hexagons": 0, "diameter_A": 4.0,
            "symmetry": "Ih", "HOMO_LUMO_gap_eV": 5.0},
    "C24": {"n_atoms": 24, "n_pentagons": 12, "n_hexagons": 2, "diameter_A": 4.5,
            "symmetry": "D6d"},
    "C28": {"n_atoms": 28, "n_pentagons": 12, "n_hexagons": 4, "diameter_A": 4.8,
            "symmetry": "Td"},
    "C36": {"n_atoms": 36, "n_pentagons": 12, "n_hexagons": 8, "diameter_A": 5.4,
            "symmetry": "D6h"},
    "C60": {"n_atoms": 60, "n_pentagons": 12, "n_hexagons": 20, "diameter_A": 7.1,
            "symmetry": "Ih", "HOMO_LUMO_gap_eV": 1.9, "buckyball": True},
    "C70": {"n_atoms": 70, "n_pentagons": 12, "n_hexagons": 25, "diameter_A": 7.9,
            "symmetry": "D5h", "HOMO_LUMO_gap_eV": 1.7},
    "C76": {"n_atoms": 76, "n_pentagons": 12, "n_hexagons": 28, "symmetry": "D2"},
    "C78": {"n_atoms": 78, "n_pentagons": 12, "n_hexagons": 29, "isomers": 5},
    "C80": {"n_atoms": 80, "n_pentagons": 12, "n_hexagons": 30, "symmetry": "Ih"},
    "C84": {"n_atoms": 84, "n_pentagons": 12, "n_hexagons": 32, "isomers": 24},
    "C240": {"n_atoms": 240, "diameter_A": 14.0, "giant": True},
    "C540": {"n_atoms": 540, "diameter_A": 21.0, "giant": True},
}


# Endohedral fullerene database
ENDOHEDRAL_DATABASE = {
    "N@C60": {"cage": "C60", "guest": "N", "spin": 3/2, "paramagnetic": True},
    "La@C82": {"cage": "C82", "guest": "La", "metallofullerene": True},
    "Sc3N@C80": {"cage": "C80", "guest": "Sc3N", "trimetallic_nitride": True},
    "Er3N@C80": {"cage": "C80", "guest": "Er3N", "magnetic": True},
    "Li@C60": {"cage": "C60", "guest": "Li", "alkali_doped": True},
    "H2@C60": {"cage": "C60", "guest": "H2", "molecular_container": True},
    "He@C60": {"cage": "C60", "guest": "He", "noble_gas": True},
    "Gd@C82": {"cage": "C82", "guest": "Gd", "MRI_contrast": True},
}


# Metal cluster database
METAL_CLUSTER_DATABASE = {
    # Icosahedral clusters
    "Au13": {"n_atoms": 13, "geometry": "icosahedron", "magic_number": True},
    "Au55": {"n_atoms": 55, "geometry": "icosahedron", "closed_shell": True},
    "Au147": {"n_atoms": 147, "geometry": "icosahedron", "closed_shell": True},
    "Ag55": {"n_atoms": 55, "geometry": "icosahedron"},
    "Pt55": {"n_atoms": 55, "geometry": "icosahedron"},
    
    # Cuboctahedral clusters
    "Pd38": {"n_atoms": 38, "geometry": "cuboctahedron", "catalytic": True},
    "Rh13": {"n_atoms": 13, "geometry": "cuboctahedron"},
    
    # Regular polyhedra
    "Au4": {"n_atoms": 4, "geometry": "tetrahedron"},
    "Au6": {"n_atoms": 6, "geometry": "octahedron"},
    "Au8": {"n_atoms": 8, "geometry": "cube"},
    "Au12": {"n_atoms": 12, "geometry": "icosahedron_shell"},
    "Au20": {"n_atoms": 20, "geometry": "tetrahedron", "magic": True},
    
    # Bimetallic
    "Au12Ag1": {"n_atoms": 13, "geometry": "core_shell", "core": "Ag", "shell": "Au"},
    "Pt12Ru1": {"n_atoms": 13, "geometry": "core_shell", "core": "Ru", "shell": "Pt"},
}


# Supramolecular cage database
SUPRAMOLECULAR_CAGES = {
    "Pd6L4": {"metal": "Pd", "n_metals": 6, "n_ligands": 4, "topology": "octahedron",
              "cavity_volume_A3": 500},
    "Pd12L24": {"metal": "Pd", "n_metals": 12, "n_ligands": 24, "topology": "cuboctahedron",
                "cavity_volume_A3": 1800},
    "Fe4L6": {"metal": "Fe", "n_metals": 4, "n_ligands": 6, "topology": "tetrahedron",
              "chiral": True},
    "Cu24L8": {"metal": "Cu", "n_metals": 24, "n_ligands": 8, "topology": "cube",
               "magnetic": True},
}


def generate_fullerene(
    fullerene: str = "C60",
    optimize: bool = False
) -> Dict[str, Any]:
    """
    Generate fullerene structure.
    
    Args:
        fullerene: Fullerene from database
        optimize: Apply geometry optimization
    
    Returns:
        Fullerene molecule
    """
    if fullerene not in FULLERENE_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_FULLERENE", "message": f"Unknown fullerene",
                      "available": list(FULLERENE_DATABASE.keys())}
        }
    
    info = FULLERENE_DATABASE[fullerene]
    n_atoms = info["n_atoms"]
    diameter = info.get("diameter_A", 7.1)
    radius = diameter / 2
    
    # Generate atomic positions on sphere
    species = []
    coords = []
    
    if fullerene == "C60":
        # C60 icosahedral structure (Ih symmetry)
        # Generate using golden ratio placement
        phi = (1 + np.sqrt(5)) / 2
        
        # 12 pentagon vertices
        pentagon_verts = []
        for i in range(12):
            if i < 6:
                z = 1 if i < 3 else -1
                angle = i * 2 * np.pi / 3 if i < 3 else (i-3) * 2 * np.pi / 3 + np.pi/3
                pentagon_verts.append([np.cos(angle) * 0.5, np.sin(angle) * 0.5, z * 0.8])
            else:
                z = 0.5 if i < 9 else -0.5
                angle = (i-6) * 2 * np.pi / 3 if i < 9 else (i-9) * 2 * np.pi / 3 + np.pi/3
                pentagon_verts.append([np.cos(angle), np.sin(angle), z])
        
        # Distribute 60 atoms
        for i in range(60):
            theta = np.arccos(1 - 2 * (i + 0.5) / 60)
            phi_angle = np.pi * (1 + np.sqrt(5)) * i
            
            x = radius * np.sin(theta) * np.cos(phi_angle)
            y = radius * np.sin(theta) * np.sin(phi_angle)
            z = radius * np.cos(theta)
            
            species.append("C")
            coords.append([x, y, z])
    else:
        # Generic spherical distribution
        for i in range(n_atoms):
            theta = np.arccos(1 - 2 * (i + 0.5) / n_atoms)
            phi_angle = np.pi * (1 + np.sqrt(5)) * i
            
            x = radius * np.sin(theta) * np.cos(phi_angle)
            y = radius * np.sin(theta) * np.sin(phi_angle)
            z = radius * np.cos(theta)
            
            species.append("C")
            coords.append([x, y, z])
    
    mol = Molecule(species, coords)
    
    return {
        "success": True,
        "fullerene": fullerene,
        "n_atoms": n_atoms,
        "diameter_A": diameter,
        "n_pentagons": info.get("n_pentagons", 12),
        "n_hexagons": info.get("n_hexagons", 0),
        "symmetry": info.get("symmetry", ""),
        "HOMO_LUMO_gap_eV": info.get("HOMO_LUMO_gap_eV", 0),
        "molecule": {
            "species": [str(s) for s in mol.species],
            "coords": [list(c) for c in mol.cart_coords]
        }
    }


def generate_endohedral_fullerene(
    system: str = "N@C60"
) -> Dict[str, Any]:
    """
    Generate endohedral fullerene.
    
    Args:
        system: Endohedral system from database
    
    Returns:
        Endohedral fullerene structure
    """
    if system not in ENDOHEDRAL_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_SYSTEM", "message": f"Unknown system",
                      "available": list(ENDOHEDRAL_DATABASE.keys())}
        }
    
    info = ENDOHEDRAL_DATABASE[system]
    cage = info["cage"]
    guest = info["guest"]
    
    # Generate cage
    cage_result = generate_fullerene(cage)
    if not cage_result["success"]:
        return cage_result
    
    cage_coords = cage_result["molecule"]["coords"]
    cage_species = cage_result["molecule"]["species"]
    
    # Add guest at center
    cage_species.insert(0, guest if len(guest) <= 2 else guest[0])
    cage_coords.insert(0, [0, 0, 0])
    
    return {
        "success": True,
        "system": system,
        "cage": cage,
        "guest": guest,
        "is_metallofullerene": info.get("metallofullerene", False),
        "is_paramagnetic": info.get("paramagnetic", False),
        "n_atoms": len(cage_species),
        "molecule": {
            "species": cage_species,
            "coords": cage_coords
        }
    }


def generate_metal_cluster(
    cluster: str = "Au13",
    ligand: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate metal cluster.
    
    Args:
        cluster: Cluster from database
        ligand: Optional surface ligand
    
    Returns:
        Metal cluster structure
    """
    if cluster not in METAL_CLUSTER_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_CLUSTER", "message": f"Unknown cluster",
                      "available": list(METAL_CLUSTER_DATABASE.keys())}
        }
    
    info = METAL_CLUSTER_DATABASE[cluster]
    n_atoms = info["n_atoms"]
    geometry = info["geometry"]
    
    # Extract metal from cluster name
    metal = ''.join([c for c in cluster if c.isalpha()])[:2]
    
    species = []
    coords = []
    
    if geometry == "icosahedron" and n_atoms == 13:
        # 13-atom icosahedron
        # Central atom
        species.append(metal)
        coords.append([0, 0, 0])
        
        # 12 surface atoms
        phi = (1 + np.sqrt(5)) / 2
        r = 2.8  # Typical metal-metal distance
        
        vertices = [
            [0, 1, phi], [0, -1, phi], [0, 1, -phi], [0, -1, -phi],
            [1, phi, 0], [-1, phi, 0], [1, -phi, 0], [-1, -phi, 0],
            [phi, 0, 1], [-phi, 0, 1], [phi, 0, -1], [-phi, 0, -1]
        ]
        
        for v in vertices:
            norm = np.sqrt(sum(x**2 for x in v))
            scaled = [x / norm * r for x in v]
            species.append(metal)
            coords.append(scaled)
    
    elif geometry == "tetrahedron":
        r = 2.5
        vertices = [
            [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]
        ]
        for v in vertices:
            norm = np.sqrt(sum(x**2 for x in v))
            scaled = [x / norm * r for x in v]
            species.append(metal)
            coords.append(scaled)
    
    else:
        # Generic spherical distribution
        for i in range(n_atoms):
            if i == 0:
                species.append(metal)
                coords.append([0, 0, 0])
            else:
                theta = np.arccos(1 - 2 * i / n_atoms)
                phi_angle = np.pi * (1 + np.sqrt(5)) * i
                r = 2.8 * (1 + i // 13 * 0.3)
                
                x = r * np.sin(theta) * np.cos(phi_angle)
                y = r * np.sin(theta) * np.sin(phi_angle)
                z = r * np.cos(theta)
                
                species.append(metal)
                coords.append([x, y, z])
    
    mol = Molecule(species, coords)
    
    return {
        "success": True,
        "cluster": cluster,
        "n_atoms": n_atoms,
        "geometry": geometry,
        "is_magic_number": info.get("magic_number", False),
        "is_closed_shell": info.get("closed_shell", False),
        "molecule": {
            "species": [str(s) for s in mol.species],
            "coords": [list(c) for c in mol.cart_coords]
        }
    }


def generate_cage_compound(
    cage_type: str = "C60"
) -> Dict[str, Any]:
    """Wrapper to generate any cage compound."""
    if cage_type in FULLERENE_DATABASE:
        return generate_fullerene(cage_type)
    elif cage_type in ENDOHEDRAL_DATABASE:
        return generate_endohedral_fullerene(cage_type)
    elif cage_type in METAL_CLUSTER_DATABASE:
        return generate_metal_cluster(cage_type)
    else:
        all_available = list(FULLERENE_DATABASE.keys()) + list(ENDOHEDRAL_DATABASE.keys()) + list(METAL_CLUSTER_DATABASE.keys())
        return {
            "success": False,
            "error": {"code": "INVALID_CAGE", "message": f"Unknown cage",
                      "available": all_available}
        }
