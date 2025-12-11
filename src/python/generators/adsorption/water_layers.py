"""
adsorption/water_layers.py - Water and Ice Adlayers

Generates water/ice structures on surfaces:
- Bilayer ice structures
- Half-dissociated water layers
- H-up/H-down orientations
- Proton ordering patterns
"""

from typing import Dict, Any, List, Optional
import numpy as np


# Water/ice structure database
WATER_LAYER_DATABASE = {
    "bilayer-ice": {
        "structure": "Ih",
        "n_molecules_per_cell": 2,
        "height_low": 2.5,
        "height_high": 3.2,
        "description": "Bilayer ice Ih structure"
    },
    "monolayer-water": {
        "structure": "flat",
        "n_molecules_per_cell": 1,
        "height": 2.3,
        "description": "Flat water monolayer"
    },
    "half-dissociated": {
        "structure": "mixed",
        "n_molecules_per_cell": 2,
        "description": "H2O + OH + H mixed layer"
    },
    "ice-Ih": {
        "structure": "hexagonal",
        "n_molecules_per_cell": 4,
        "description": "Hexagonal ice"
    },
    "ice-Ic": {
        "structure": "cubic",
        "n_molecules_per_cell": 8,
        "description": "Cubic ice"
    },
}


def generate_water_bilayer(
    surface: Dict[str, Any],
    orientation: str = "H-down",
    coverage: float = 1.0,
    dissociation_fraction: float = 0.0
) -> Dict[str, Any]:
    """
    Generate water bilayer on surface.
    
    Args:
        surface: Surface structure
        orientation: 'H-down', 'H-up', 'random'
        coverage: Coverage relative to full bilayer
        dissociation_fraction: Fraction of dissociated water (0-1)
    
    Returns:
        Surface with water bilayer
    """
    # Water molecule geometry
    r_OH = 0.96  # O-H bond length
    angle_HOH = np.radians(104.5)
    
    # Get surface info
    surface_atoms = surface.get("structure", {}).get("atoms", [])
    z_coords = [a.get("cartesian", [0,0,0])[2] for a in surface_atoms]
    z_max = max(z_coords) if z_coords else 0
    
    x_coords = [a.get("cartesian", [0,0,0])[0] for a in surface_atoms]
    y_coords = [a.get("cartesian", [0,0,0])[1] for a in surface_atoms]
    
    x_min, x_max = min(x_coords) if x_coords else 0, max(x_coords) if x_coords else 10
    y_min, y_max = min(y_coords) if y_coords else 0, max(y_coords) if y_coords else 10
    
    # Bilayer spacing
    height_low = 2.5
    height_high = 3.2
    lateral_spacing = 4.5  # Approximate ice-like spacing
    
    new_atoms = list(surface_atoms)
    n_water = 0
    n_dissociated = 0
    
    np.random.seed(42)
    
    x = x_min
    while x < x_max:
        y = y_min
        row = 0
        while y < y_max:
            # Lower water molecule
            z_o = z_max + height_low
            
            # Check if this water dissociates
            if np.random.random() < dissociation_fraction:
                # Add OH*
                new_atoms.append({"element": "O", "cartesian": [x, y, z_o], "adsorbate": True, "molecule": "OH"})
                # One H pointing up
                new_atoms.append({"element": "H", "cartesian": [x, y, z_o + r_OH], "adsorbate": True})
                # Dissociated H on surface
                new_atoms.append({"element": "H", "cartesian": [x + 1.5, y, z_max + 1.0], "adsorbate": True, "molecule": "H*"})
                n_dissociated += 1
            else:
                # Add intact water
                new_atoms.append({"element": "O", "cartesian": [x, y, z_o], "adsorbate": True, "molecule": "H2O"})
                
                if orientation == "H-down":
                    h1_z = z_o - r_OH * np.cos(angle_HOH/2)
                    h2_z = z_o - r_OH * np.cos(angle_HOH/2)
                elif orientation == "H-up":
                    h1_z = z_o + r_OH * np.cos(angle_HOH/2)
                    h2_z = z_o + r_OH * np.cos(angle_HOH/2)
                else:
                    h1_z = z_o + r_OH * np.cos(angle_HOH/2) * (1 if np.random.random() > 0.5 else -1)
                    h2_z = h1_z
                
                h_offset = r_OH * np.sin(angle_HOH/2)
                new_atoms.append({"element": "H", "cartesian": [x - h_offset, y, h1_z], "adsorbate": True})
                new_atoms.append({"element": "H", "cartesian": [x + h_offset, y, h2_z], "adsorbate": True})
                n_water += 1
            
            # Upper water molecule (staggered)
            if row % 2 == 0:
                z_o_high = z_max + height_high
                x_high = x + lateral_spacing / 2
                
                if x_high < x_max and np.random.random() < coverage:
                    new_atoms.append({"element": "O", "cartesian": [x_high, y + lateral_spacing/2, z_o_high], "adsorbate": True, "molecule": "H2O"})
                    new_atoms.append({"element": "H", "cartesian": [x_high - 0.5, y + lateral_spacing/2, z_o_high + 0.8], "adsorbate": True})
                    new_atoms.append({"element": "H", "cartesian": [x_high + 0.5, y + lateral_spacing/2, z_o_high + 0.8], "adsorbate": True})
                    n_water += 1
            
            y += lateral_spacing
            row += 1
        x += lateral_spacing
    
    return {
        "success": True,
        "layer_type": "water_bilayer",
        "orientation": orientation,
        "coverage": coverage,
        "n_water_molecules": n_water,
        "n_dissociated": n_dissociated,
        "dissociation_fraction": dissociation_fraction,
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms}
    }


def generate_ice_layer(
    surface: Dict[str, Any],
    ice_type: str = "Ih",
    n_bilayers: int = 1,
    proton_ordering: str = "random"
) -> Dict[str, Any]:
    """
    Generate ice layer on surface.
    
    Args:
        surface: Surface structure
        ice_type: 'Ih' (hexagonal), 'Ic' (cubic)
        n_bilayers: Number of ice bilayers
        proton_ordering: 'random', 'ordered', 'antiferroelectric'
    
    Returns:
        Surface with ice layer
    """
    # Get surface info
    surface_atoms = surface.get("structure", {}).get("atoms", [])
    z_coords = [a.get("cartesian", [0,0,0])[2] for a in surface_atoms]
    z_max = max(z_coords) if z_coords else 0
    
    x_coords = [a.get("cartesian", [0,0,0])[0] for a in surface_atoms]
    y_coords = [a.get("cartesian", [0,0,0])[1] for a in surface_atoms]
    x_span = max(x_coords) - min(x_coords) if x_coords else 10
    y_span = max(y_coords) - min(y_coords) if y_coords else 10
    
    # Ice parameters
    a_ice = 4.5  # Approximate ice lattice constant
    c_ice = 7.3 / 2  # Half of c for bilayer
    r_OH = 0.96
    
    new_atoms = list(surface_atoms)
    n_water = 0
    
    np.random.seed(42)
    
    for bilayer in range(n_bilayers):
        z_base = z_max + 2.5 + bilayer * c_ice
        
        nx = max(1, int(x_span / a_ice))
        ny = max(1, int(y_span / (a_ice * np.sqrt(3)/2)))
        
        for i in range(nx):
            for j in range(ny):
                x = min(x_coords) + i * a_ice + (j % 2) * a_ice / 2 if x_coords else i * a_ice
                y = min(y_coords) + j * a_ice * np.sqrt(3) / 2 if y_coords else j * a_ice * np.sqrt(3) / 2
                
                # Lower molecule
                z_o = z_base
                new_atoms.append({"element": "O", "cartesian": [x, y, z_o], "adsorbate": True})
                
                # H orientations based on proton ordering
                if proton_ordering == "ordered":
                    h_angles = [0, np.pi]
                elif proton_ordering == "antiferroelectric":
                    h_angles = [0, np.pi] if (i + j) % 2 == 0 else [np.pi/2, 3*np.pi/2]
                else:
                    h_angles = np.random.uniform(0, 2*np.pi, 2)
                
                for angle in h_angles:
                    hx = x + r_OH * np.cos(angle) * 0.5
                    hy = y + r_OH * np.sin(angle) * 0.5
                    hz = z_o + r_OH * 0.866
                    new_atoms.append({"element": "H", "cartesian": [hx, hy, hz], "adsorbate": True})
                
                n_water += 1
    
    return {
        "success": True,
        "ice_type": ice_type,
        "n_bilayers": n_bilayers,
        "proton_ordering": proton_ordering,
        "n_water_molecules": n_water,
        "n_atoms": len(new_atoms),
        "structure": {"atoms": new_atoms}
    }
