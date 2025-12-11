"""
biomolecular/lipids.py - Lipid Bilayers and Membranes

Generates lipid membrane structures:
- POPC, DPPC, DOPC bilayers
- Cholesterol incorporation
- Membrane proteins
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Lipid database
LIPID_DATABASE = {
    "POPC": {
        "name": "1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine",
        "n_carbons": [16, 18], "area_per_lipid": 68.3,
        "bilayer_thickness": 39.1, "tail_length": 15.0
    },
    "DPPC": {
        "name": "1,2-dipalmitoyl-sn-glycero-3-phosphocholine",
        "n_carbons": [16, 16], "area_per_lipid": 64.0,
        "bilayer_thickness": 38.3, "tail_length": 14.5
    },
    "DOPC": {
        "name": "1,2-dioleoyl-sn-glycero-3-phosphocholine",
        "n_carbons": [18, 18], "area_per_lipid": 72.4,
        "bilayer_thickness": 36.9, "tail_length": 14.0
    },
    "DLPC": {
        "name": "1,2-dilauroyl-sn-glycero-3-phosphocholine",
        "n_carbons": [12, 12], "area_per_lipid": 63.2,
        "bilayer_thickness": 30.8, "tail_length": 10.5
    },
    "DMPC": {
        "name": "1,2-dimyristoyl-sn-glycero-3-phosphocholine",
        "n_carbons": [14, 14], "area_per_lipid": 60.6,
        "bilayer_thickness": 35.3, "tail_length": 12.0
    },
    "POPE": {
        "name": "1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine",
        "n_carbons": [16, 18], "area_per_lipid": 56.0,
        "bilayer_thickness": 40.2, "tail_length": 15.5
    },
    "POPS": {
        "name": "1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-L-serine",
        "n_carbons": [16, 18], "area_per_lipid": 55.1,
        "bilayer_thickness": 39.0, "charge": -1
    },
    "cholesterol": {
        "name": "Cholesterol",
        "area": 38.0, "length": 20.0
    },
    "cardiolipin": {
        "name": "Cardiolipin",
        "n_tails": 4, "area_per_lipid": 130.0, "charge": -2
    },
    "GM1": {
        "name": "Ganglioside GM1",
        "area_per_lipid": 85.0, "has_sugar": True
    },
}


def generate_lipid_bilayer(
    lipid: str = "POPC",
    size: List[int] = [10, 10],
    cholesterol_fraction: float = 0.0,
    water_buffer: float = 20.0,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate lipid bilayer structure.
    
    Args:
        lipid: Lipid type from database
        size: Number of lipids [nx, ny]
        cholesterol_fraction: Fraction of cholesterol (0-0.5)
        water_buffer: Water layer thickness on each side
        seed: Random seed
    
    Returns:
        Lipid bilayer structure
    """
    if lipid not in LIPID_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_LIPID", "message": f"Unknown lipid '{lipid}'",
                      "available": list(LIPID_DATABASE.keys())}
        }
    
    np.random.seed(seed)
    info = LIPID_DATABASE[lipid]
    
    area = info.get("area_per_lipid", 65.0)
    thickness = info.get("bilayer_thickness", 38.0)
    tail_length = info.get("tail_length", 14.0)
    
    nx, ny = size
    spacing = np.sqrt(area)
    
    lx = nx * spacing
    ly = ny * spacing
    lz = thickness + 2 * water_buffer
    
    atoms = []
    n_lipids = 0
    n_cholesterol = 0
    
    for i in range(nx):
        for j in range(ny):
            x = (i + 0.5) * spacing
            y = (j + 0.5) * spacing
            
            for leaflet in ["upper", "lower"]:
                # Determine if this is cholesterol
                is_chol = np.random.random() < cholesterol_fraction
                
                if leaflet == "upper":
                    z_head = lz / 2 + thickness / 2
                    z_tail = lz / 2
                else:
                    z_head = lz / 2 - thickness / 2
                    z_tail = lz / 2
                
                # Add some randomness
                x += np.random.uniform(-1, 1)
                y += np.random.uniform(-1, 1)
                
                if is_chol:
                    # Simplified cholesterol (4 rings + tail)
                    atoms.append({"element": "O", "cartesian": [x, y, z_head], "molecule": "CHOL", "leaflet": leaflet})
                    for k in range(4):
                        z_ring = z_head - k * 3 if leaflet == "upper" else z_head + k * 3
                        atoms.append({"element": "C", "cartesian": [x, y, z_ring], "molecule": "CHOL"})
                    n_cholesterol += 1
                else:
                    # Phospholipid head group
                    atoms.append({"element": "P", "cartesian": [x, y, z_head], "molecule": lipid, "leaflet": leaflet})
                    atoms.append({"element": "N", "cartesian": [x + 1.5, y, z_head + 2], "molecule": lipid})
                    atoms.append({"element": "O", "cartesian": [x - 1, y, z_head], "molecule": lipid})
                    atoms.append({"element": "O", "cartesian": [x + 1, y, z_head], "molecule": lipid})
                    
                    # Glycerol
                    z_glyc = z_head - 3 if leaflet == "upper" else z_head + 3
                    atoms.append({"element": "C", "cartesian": [x, y, z_glyc], "molecule": lipid})
                    atoms.append({"element": "O", "cartesian": [x, y + 1, z_glyc], "molecule": lipid})
                    atoms.append({"element": "O", "cartesian": [x, y - 1, z_glyc], "molecule": lipid})
                    
                    # Two tails (simplified as multiple C atoms)
                    n_carbons = info.get("n_carbons", [16, 16])
                    for tail, n_c in enumerate(n_carbons):
                        y_offset = (tail - 0.5) * 4
                        for c in range(min(n_c // 2, 8)):  # Simplified
                            if leaflet == "upper":
                                z_c = z_glyc - (c + 1) * 1.54
                            else:
                                z_c = z_glyc + (c + 1) * 1.54
                            atoms.append({"element": "C", "cartesian": [x, y + y_offset, z_c], "molecule": lipid})
                    
                    n_lipids += 1
    
    return {
        "success": True,
        "lipid": lipid,
        "description": info["name"],
        "size": size,
        "n_lipids": n_lipids,
        "n_cholesterol": n_cholesterol,
        "cholesterol_fraction": cholesterol_fraction,
        "area_per_lipid_A2": area,
        "bilayer_thickness_A": thickness,
        "box_dimensions_A": [round(lx, 1), round(ly, 1), round(lz, 1)],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_membrane_with_protein(
    lipid: str = "POPC",
    membrane_size: List[int] = [20, 20],
    protein_type: str = "GPCR",
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate membrane with embedded protein.
    
    Args:
        lipid: Lipid type
        membrane_size: Membrane size
        protein_type: Protein type (GPCR, channel, porin)
        seed: Random seed
    
    Returns:
        Membrane-protein structure
    """
    # First generate membrane
    membrane = generate_lipid_bilayer(lipid, membrane_size, seed=seed)
    
    if not membrane["success"]:
        return membrane
    
    atoms = membrane["structure"]["atoms"]
    
    # Protein parameters
    protein_params = {
        "GPCR": {"n_helices": 7, "radius": 15, "height": 40},
        "channel": {"n_helices": 4, "radius": 10, "height": 50},
        "porin": {"n_barrels": 16, "radius": 12, "height": 35},
    }
    
    if protein_type not in protein_params:
        return {"success": False, "error": {"code": "INVALID_PROTEIN", "message": f"Unknown protein"}}
    
    params = protein_params[protein_type]
    lz = membrane["box_dimensions_A"][2]
    
    # Remove lipids in protein region
    protein_center = [membrane["box_dimensions_A"][0] / 2, membrane["box_dimensions_A"][1] / 2]
    protein_radius = params["radius"]
    
    filtered_atoms = []
    for atom in atoms:
        pos = atom["cartesian"]
        dist = np.sqrt((pos[0] - protein_center[0])**2 + (pos[1] - protein_center[1])**2)
        if dist > protein_radius:
            filtered_atoms.append(atom)
    
    # Add protein helices
    n_helices = params.get("n_helices", params.get("n_barrels", 4))
    helix_radius = protein_radius * 0.7
    
    for h in range(n_helices):
        angle = 2 * np.pi * h / n_helices
        center_x = protein_center[0] + helix_radius * np.cos(angle)
        center_y = protein_center[1] + helix_radius * np.sin(angle)
        
        # Helix atoms
        height = params["height"]
        n_residues = int(height / 1.5)
        
        for r in range(n_residues):
            z = lz / 2 - height / 2 + r * 1.5
            helix_angle = r * np.radians(100)  # Alpha helix rotation
            
            x = center_x + 2.3 * np.cos(helix_angle)
            y = center_y + 2.3 * np.sin(helix_angle)
            
            # Backbone
            filtered_atoms.append({"element": "N", "cartesian": [x, y, z], "molecule": "protein", "helix": h})
            filtered_atoms.append({"element": "C", "cartesian": [x + 0.5, y + 0.5, z + 0.5], "molecule": "protein"})
            filtered_atoms.append({"element": "C", "cartesian": [x + 1, y, z + 1], "molecule": "protein"})
            filtered_atoms.append({"element": "O", "cartesian": [x + 2, y, z + 1], "molecule": "protein"})
    
    return {
        "success": True,
        "lipid": lipid,
        "protein_type": protein_type,
        "n_helices": n_helices,
        "protein_radius_A": protein_radius,
        "n_lipids": membrane["n_lipids"],
        "n_atoms": len(filtered_atoms),
        "structure": {"atoms": filtered_atoms}
    }
