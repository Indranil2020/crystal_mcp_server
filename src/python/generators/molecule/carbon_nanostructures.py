"""
molecule/carbon_nanostructures.py - Carbon Nanostructures

Generates advanced carbon nanostructures:
- Carbon nanohorns
- Peapods (CNT filled with fullerenes)
- Nanobuds (CNT with fullerene attachments)
- Schwarzite structures
- Graphene quantum dots
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


# Nanohorn types
NANOHORN_DATABASE = {
    "dahlia": {"cone_angle": 20, "length": 50, "tip_radius": 1.0, "aggregate": "dahlia-like"},
    "bud": {"cone_angle": 25, "length": 40, "tip_radius": 1.5, "aggregate": "bud-like"},
    "seed": {"cone_angle": 30, "length": 30, "tip_radius": 2.0, "aggregate": "seed-like"},
    "pointed": {"cone_angle": 15, "length": 60, "tip_radius": 0.5, "aggregate": "pointed"},
}


# Peapod configurations
PEAPOD_DATABASE = {
    "C60@(10,10)": {"tube": [10, 10], "fullerene": "C60", "filling": 1.0},
    "C70@(11,11)": {"tube": [11, 11], "fullerene": "C70", "filling": 1.0},
    "C60@(9,9)": {"tube": [9, 9], "fullerene": "C60", "filling": 1.0},
    "Gd@C82@(13,13)": {"tube": [13, 13], "fullerene": "Gd@C82", "filling": 0.8},
    "C60@DWCNT": {"tube": [15, 15], "inner": [10, 10], "fullerene": "C60", "filling": 1.0},
}


def generate_nanohorn(
    horn_type: str = "dahlia",
    n_atoms: int = 500,
    closed_tip: bool = True
) -> Dict[str, Any]:
    """
    Generate carbon nanohorn structure.
    
    Args:
        horn_type: Type of nanohorn
        n_atoms: Approximate number of atoms
        closed_tip: Whether tip is closed
    
    Returns:
        Nanohorn structure
    """
    if horn_type not in NANOHORN_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_TYPE", "message": f"Unknown type '{horn_type}'",
                      "available": list(NANOHORN_DATABASE.keys())}
        }
    
    info = NANOHORN_DATABASE[horn_type]
    cone_angle = np.radians(info["cone_angle"])
    length = info["length"]
    tip_radius = info["tip_radius"]
    
    atoms = []
    
    # Generate cone surface with graphene-like arrangement
    c_c = 1.42  # C-C bond length
    
    # Number of rings along length
    n_rings = int(length / (c_c * np.sqrt(3) / 2))
    atoms_per_ring_base = 30  # At base
    
    for i in range(n_rings):
        z = i * c_c * np.sqrt(3) / 2
        # Radius decreases toward tip
        r = (length - z) * np.tan(cone_angle) + tip_radius
        
        if r < tip_radius:
            r = tip_radius
        
        # Number of atoms in this ring (proportional to circumference)
        n_in_ring = max(6, int(2 * np.pi * r / c_c))
        
        for j in range(n_in_ring):
            angle = 2 * np.pi * j / n_in_ring
            # Alternate z slightly for graphene-like structure
            z_offset = (j % 2) * 0.4
            
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            
            atoms.append({"element": "C", "cartesian": [x, y, z + z_offset]})
            
            if len(atoms) >= n_atoms:
                break
        
        if len(atoms) >= n_atoms:
            break
    
    # Close tip with pentagons
    if closed_tip and len(atoms) < n_atoms:
        for i in range(5):
            angle = 2 * np.pi * i / 5
            atoms.append({"element": "C", "cartesian": [
                tip_radius * 0.5 * np.cos(angle),
                tip_radius * 0.5 * np.sin(angle),
                length
            ]})
        atoms.append({"element": "C", "cartesian": [0, 0, length + 0.5]})
    
    return {
        "success": True,
        "horn_type": horn_type,
        "cone_angle_deg": info["cone_angle"],
        "length_angstrom": length,
        "tip_radius_angstrom": tip_radius,
        "aggregate_type": info["aggregate"],
        "closed_tip": closed_tip,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_peapod(
    configuration: str = "C60@(10,10)",
    n_fullerenes: int = 5,
    spacing_factor: float = 1.0
) -> Dict[str, Any]:
    """
    Generate carbon nanotube peapod (fullerenes inside CNT).
    
    Args:
        configuration: Peapod configuration from database
        n_fullerenes: Number of fullerenes
        spacing_factor: Spacing adjustment (1.0 = touching)
    
    Returns:
        Peapod structure
    """
    if configuration not in PEAPOD_DATABASE:
        return {
            "success": False,
            "error": {"code": "INVALID_CONFIG", "message": f"Unknown configuration '{configuration}'",
                      "available": list(PEAPOD_DATABASE.keys())}
        }
    
    info = PEAPOD_DATABASE[configuration]
    n, m = info["tube"]
    fullerene_type = info["fullerene"]
    
    atoms = []
    
    # CNT parameters
    c_c = 1.42
    a = 2.46  # Graphene lattice constant
    
    # CNT radius
    r_cnt = a * np.sqrt(n**2 + m**2 + n*m) / (2 * np.pi)
    
    # Fullerene radius (approximate)
    if "C60" in fullerene_type:
        r_full = 3.55
        n_c_full = 60
    elif "C70" in fullerene_type:
        r_full = 3.8
        n_c_full = 70
    else:
        r_full = 4.0
        n_c_full = 82
    
    # Fullerene spacing
    fullerene_spacing = (2 * r_full + 0.5) * spacing_factor
    
    # CNT length needed
    cnt_length = (n_fullerenes + 1) * fullerene_spacing
    
    # Generate CNT shell
    n_rings = int(cnt_length / (c_c * np.sqrt(3) / 2))
    n_per_ring = int(2 * np.pi * r_cnt / c_c)
    
    for ring in range(n_rings):
        z = ring * c_c * np.sqrt(3) / 2
        for i in range(n_per_ring):
            angle = 2 * np.pi * i / n_per_ring
            x = r_cnt * np.cos(angle)
            y = r_cnt * np.sin(angle)
            atoms.append({"element": "C", "cartesian": [x, y, z]})
    
    # Generate fullerenes inside
    for f_idx in range(n_fullerenes):
        z_center = (f_idx + 0.5) * fullerene_spacing
        
        # Simplified fullerene as spherical shell
        # Icosahedral arrangement of 60 atoms
        golden = (1 + np.sqrt(5)) / 2
        
        # Approximate uniform distribution on sphere
        for i in range(n_c_full):
            phi = np.arccos(1 - 2 * (i + 0.5) / n_c_full)
            theta = np.pi * (1 + np.sqrt(5)) * i
            
            x = r_full * np.sin(phi) * np.cos(theta)
            y = r_full * np.sin(phi) * np.sin(theta)
            z = r_full * np.cos(phi) + z_center
            
            atoms.append({"element": "C", "cartesian": [x, y, z]})
    
    return {
        "success": True,
        "configuration": configuration,
        "tube_chirality": [n, m],
        "tube_radius_angstrom": round(r_cnt, 2),
        "fullerene_type": fullerene_type,
        "n_fullerenes": n_fullerenes,
        "fullerene_spacing_angstrom": round(fullerene_spacing, 2),
        "total_length_angstrom": round(cnt_length, 1),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_nanobud(
    tube_chirality: Tuple[int, int] = (10, 10),
    n_buds: int = 3,
    bud_type: str = "C60",
    attachment_pattern: str = "spiral"
) -> Dict[str, Any]:
    """
    Generate nanobud (CNT with fullerene attachments).
    
    Args:
        tube_chirality: (n, m) indices
        n_buds: Number of attached fullerenes
        bud_type: Fullerene type (C60, C70)
        attachment_pattern: 'spiral', 'line', 'random'
    
    Returns:
        Nanobud structure
    """
    n, m = tube_chirality
    
    atoms = []
    
    # CNT parameters
    c_c = 1.42
    a = 2.46
    r_cnt = a * np.sqrt(n**2 + m**2 + n*m) / (2 * np.pi)
    
    # Fullerene radius
    r_full = 3.55 if bud_type == "C60" else 3.8
    n_c_full = 60 if bud_type == "C60" else 70
    
    cnt_length = 40  # Angstrom
    
    # Generate CNT
    n_rings = int(cnt_length / (c_c * np.sqrt(3) / 2))
    n_per_ring = int(2 * np.pi * r_cnt / c_c)
    
    for ring in range(n_rings):
        z = ring * c_c * np.sqrt(3) / 2
        for i in range(n_per_ring):
            angle = 2 * np.pi * i / n_per_ring
            x = r_cnt * np.cos(angle)
            y = r_cnt * np.sin(angle)
            atoms.append({"element": "C", "cartesian": [x, y, z]})
    
    # Attach fullerenes on outside
    bud_distance = r_cnt + r_full + 1.5  # Gap between surfaces
    
    for b_idx in range(n_buds):
        if attachment_pattern == "spiral":
            z_pos = cnt_length * (b_idx + 0.5) / n_buds
            angle = 2 * np.pi * b_idx / n_buds
        elif attachment_pattern == "line":
            z_pos = cnt_length * (b_idx + 0.5) / n_buds
            angle = 0
        else:
            z_pos = np.random.uniform(5, cnt_length - 5)
            angle = np.random.uniform(0, 2 * np.pi)
        
        # Fullerene center
        cx = bud_distance * np.cos(angle)
        cy = bud_distance * np.sin(angle)
        
        # Generate fullerene
        for i in range(n_c_full):
            phi = np.arccos(1 - 2 * (i + 0.5) / n_c_full)
            theta = np.pi * (1 + np.sqrt(5)) * i
            
            x = r_full * np.sin(phi) * np.cos(theta) + cx
            y = r_full * np.sin(phi) * np.sin(theta) + cy
            z = r_full * np.cos(phi) + z_pos
            
            atoms.append({"element": "C", "cartesian": [x, y, z]})
    
    return {
        "success": True,
        "tube_chirality": list(tube_chirality),
        "tube_radius_angstrom": round(r_cnt, 2),
        "bud_type": bud_type,
        "n_buds": n_buds,
        "attachment_pattern": attachment_pattern,
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_graphene_quantum_dot(
    shape: str = "hexagonal",
    size_nm: float = 2.0,
    edge_type: str = "zigzag",
    edge_passivation: str = "H"
) -> Dict[str, Any]:
    """
    Generate graphene quantum dot.
    
    Args:
        shape: 'hexagonal', 'triangular', 'circular'
        size_nm: Approximate diameter in nm
        edge_type: 'zigzag', 'armchair'
        edge_passivation: Edge atoms (H, F, OH, etc.)
    
    Returns:
        Graphene quantum dot structure
    """
    atoms = []
    
    c_c = 1.42
    size_ang = size_nm * 10  # Convert to Angstrom
    
    # Generate hexagonal lattice
    a1 = np.array([c_c * np.sqrt(3), 0])
    a2 = np.array([c_c * np.sqrt(3) / 2, c_c * 1.5])
    
    # Basis atoms
    basis = [np.array([0, 0]), np.array([c_c * np.sqrt(3) / 2, c_c / 2])]
    
    n_cells = int(size_ang / (c_c * np.sqrt(3))) + 2
    center = size_ang / 2
    
    carbon_positions = []
    
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            for b in basis:
                pos = i * a1 + j * a2 + b
                
                # Apply shape cutoff
                dist = np.linalg.norm(pos)
                
                include = False
                if shape == "hexagonal":
                    # Hexagonal boundary
                    include = dist < size_ang / 2
                elif shape == "triangular":
                    # Triangular boundary
                    include = (pos[0] > 0 and pos[1] > 0 and 
                              pos[0] + pos[1] * np.sqrt(3) < size_ang)
                else:  # circular
                    include = dist < size_ang / 2
                
                if include:
                    carbon_positions.append(pos)
                    atoms.append({"element": "C", "cartesian": [pos[0], pos[1], 0]})
    
    # Find edge atoms and add passivation
    carbon_set = set(tuple(p.round(2)) for p in carbon_positions)
    
    for pos in carbon_positions:
        # Check neighbors
        neighbors = 0
        for dp in [a1, -a1, a2, -a2, a1-a2, a2-a1]:
            neighbor_pos = tuple((pos + dp * 0.5).round(2))
            if neighbor_pos in carbon_set:
                neighbors += 1
        
        # Edge atoms have fewer than 3 neighbors
        if neighbors < 3:
            # Add passivation
            h_distance = 1.09  # C-H bond
            for angle in np.linspace(0, 2*np.pi, 6, endpoint=False):
                h_pos = pos + np.array([h_distance * np.cos(angle), h_distance * np.sin(angle)])
                if tuple(h_pos.round(2)) not in carbon_set:
                    atoms.append({"element": edge_passivation, "cartesian": [h_pos[0], h_pos[1], 0]})
                    break
    
    return {
        "success": True,
        "shape": shape,
        "size_nm": size_nm,
        "edge_type": edge_type,
        "edge_passivation": edge_passivation,
        "n_carbon": sum(1 for a in atoms if a["element"] == "C"),
        "n_edge_atoms": sum(1 for a in atoms if a["element"] != "C"),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }
