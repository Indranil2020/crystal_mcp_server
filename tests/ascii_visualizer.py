#!/usr/bin/env python3
"""
ascii_visualizer.py - ASCII Ball-Stick Diagram Generator

Generates ASCII representations of molecular structures from JSON output.
Used for visual validation in end-to-end testing.

Example output:
    H           
    │           
H───C───H   
    │           
    C           
   /|\          
  ···           
"""

import numpy as np
from typing import List, Tuple, Dict, Set

# Covalent radii for bond detection (in Angstroms)
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 0.99, 'Br': 1.14,
    'Fe': 1.25, 'Cu': 1.28, 'Zn': 1.22, 'Na': 1.66, 'K': 2.03,
}

# ASCII atom representations
ATOM_CHARS = {
    'C': '●',   # Carbon - filled circle
    'H': '○',   # Hydrogen - small hollow
    'O': '◎',   # Oxygen - double circle
    'N': '◆',   # Nitrogen - filled diamond
    'S': '◇',   # Sulfur - hollow diamond
    'F': '□',   # Fluorine - hollow square
    'Cl': '■',  # Chlorine - filled square
    'default': '●'
}


def get_bond_order(elem1: str, elem2: str, distance: float) -> int:
    """Determine bond order based on distance."""
    r1 = COVALENT_RADII.get(elem1, 1.0)
    r2 = COVALENT_RADII.get(elem2, 1.0)
    expected = r1 + r2
    
    if distance < expected * 0.8:
        return 3  # Triple bond
    elif distance < expected * 0.95:
        return 2  # Double bond
    elif distance < expected * 1.3:
        return 1  # Single bond
    return 0  # No bond


def detect_bonds(atoms: List[str], coords: List[Tuple[float, float, float]], 
                 tolerance: float = 1.3) -> List[Tuple[int, int, int]]:
    """Detect bonds between atoms based on covalent radii."""
    bonds = []
    n = len(atoms)
    
    for i in range(n):
        for j in range(i + 1, n):
            r1 = COVALENT_RADII.get(atoms[i], 1.0)
            r2 = COVALENT_RADII.get(atoms[j], 1.0)
            max_dist = (r1 + r2) * tolerance
            
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            dz = coords[j][2] - coords[i][2]
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < max_dist:
                order = get_bond_order(atoms[i], atoms[j], dist)
                if order > 0:
                    bonds.append((i, j, order))
    
    return bonds



def project_2d(coords: List[Tuple[float, float, float]], 
               plane: str = 'xy') -> List[Tuple[float, float]]:
    """
    Project 3D coordinates to 2D plane.
    
    Args:
        coords: List of (x,y,z)
        plane: 'xy' (Top) or 'xz' (Side)
        
    Returns:
        List of (u, v) 2D coordinates
    """
    result = []
    for x, y, z in coords:
        if plane == 'xz':
            # Side view: X represents lateral, Z represents vertical
            result.append((x, z))
        else:
            # Top view: X, Y
            result.append((x, y))
    return result


def draw_line(canvas: List[List[str]], x1: int, y1: int, x2: int, y2: int, order: int = 1):
    """Draw a bond line on the canvas using Bresenham's algorithm."""
    if not canvas: return
    height = len(canvas)
    width = len(canvas[0])
    
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x2 > x1 else -1
    sy = 1 if y2 > y1 else -1
    
    char = '─'
    if dx > dy:
        # Horizontal-ish bond
        char = '═' if order == 2 else '≡' if order == 3 else '─'
    else:
        # Vertical-ish bond
        char = '║' if order == 2 else '┃' if order == 3 else '│'
    
    # Diagonal check
    if dx > 0 and dy > 0:
        # Use simple diagonal chars
        if (x2 - x1) * (y2 - y1) > 0:
            char = '╲'
        else:
            char = '╱'
    
    err = dx - dy
    x, y = x1, y1
    
    while True:
        if 0 <= x < width and 0 <= y < height:
            if canvas[y][x] == ' ':
                canvas[y][x] = char
        
        if x == x2 and y == y2:
            break
            
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x += sx
        if e2 < dx:
            err += dx
            y += sy



def render_ascii_structure(atoms: List[str], 
                           coords: List[Tuple[float, float, float]],
                           width: int = 40,
                           height: int = 20,
                           view_plane: str = 'xy',
                           title: str = None) -> str:
    """
    Render atoms and bonds as ASCII diagram for a specific plane.
    """
    if not atoms or not coords:
        return "  (No structure data)"
    
    # Project
    coords_2d = project_2d(coords, plane=view_plane)
    
    xs = [c[0] for c in coords_2d]
    ys = [c[1] for c in coords_2d]
    
    if not xs: return "  (Empty)"
    
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    
    # Scale to canvas
    # Aspect ratio: Characters are ~2x higher than wide.
    # To preserve squareness: 1 unit x = 1 char, 1 unit y = 0.5 char (visually)
    # Scale factors: sx = W / range_x, sy = H / range_y * 2.0
    
    range_x = max(max_x - min_x, 1e-3)
    range_y = max(max_y - min_y, 1e-3)
    
    pad = 2
    avail_w = max(width - 2*pad, 1)
    avail_h = max(height - 2*pad, 1)
    
    scale_x = avail_w / range_x
    scale_y = avail_h / range_y * 1.8 # Correction for char aspect
    
    scale = min(scale_x, scale_y)
    
    # Canvas
    canvas = [[' ' for _ in range(width)] for _ in range(height)]
    
    def to_pixel(u, v):
        # Center the structure
        cx = int((u - min_x) * scale + (avail_w - range_x * scale)/2) + pad
        # Flip Y (screen coords go down, physical y/z goes up)
        cy = int((max_y - v) * scale / 1.8 + (avail_h - range_y * scale / 1.8)/2) + pad
        return cx, cy

    # Bonds
    bonds = detect_bonds(atoms, coords)
    for i, j, order in bonds:
        x1, y1 = to_pixel(coords_2d[i][0], coords_2d[i][1])
        x2, y2 = to_pixel(coords_2d[j][0], coords_2d[j][1])
        draw_line(canvas, x1, y1, x2, y2, order)
        
    # Atoms
    # Sort by depth (perpendicular axis)
    # If xy view, depth is z. If xz view, depth is y (going into page).
    depths = []
    for k, (_, _, z_orig) in enumerate(coords):
        if view_plane == 'xz':
            depth = coords[k][1] # Y is depth
        else:
            depth = z_orig       # Z is depth
        depths.append((depth, k))
        
    depths.sort() # Back to front
    
    for _, idx in depths:
        elem = atoms[idx]
        u, v = coords_2d[idx]
        px, py = to_pixel(u, v)
        
        if 0 <= px < width and 0 <= py < height:
            char = ATOM_CHARS.get(elem, ATOM_CHARS['default'])
            canvas[py][px] = char
            # Label
            if px + 1 < width:
                canvas[py][px] = elem[0]
                if len(elem) > 1 and px + 2 < width:
                    canvas[py][px + 1] = elem[1].lower()

    # Border
    color_code = "\033[1m" # Bold
    reset = "\033[0m"
    
    # Just standard ASCII border
    canvas[0] = ['-'] * width
    canvas[-1] = ['-'] * width
    for r in range(height):
        canvas[r][0] = '|'
        canvas[r][-1] = '|'
    canvas[0][0] = '+'
    canvas[0][-1] = '+'
    canvas[-1][0] = '+'
    canvas[-1][-1] = '+'
    
    if title:
        t_start = max(1, (width - len(title)) // 2)
        for i, char in enumerate(title):
            if t_start + i < width - 1:
                canvas[0][t_start + i] = char

    return '\n'.join(''.join(row) for row in canvas)


def render_structure_from_json(structure_json: dict, max_atoms: int = 50) -> str:
    """
    Render Dual-View ASCII diagram (Top + Side).
    """
    atoms_data = structure_json.get('atoms', structure_json.get('sites', []))
    if not atoms_data: return "  (No atom data)\n"
    
    atoms = []
    coords = []
    for atom in atoms_data[:max_atoms]:
        elem = atom.get('element', 'X')
        cart_raw = atom.get('cartesian', atom.get('coords', [0, 0, 0]))
        if isinstance(cart_raw, dict):
            cart = [cart_raw.get('x',0), cart_raw.get('y',0), cart_raw.get('z',0)]
        else:
            cart = list(cart_raw)
        atoms.append(elem)
        coords.append(tuple(cart))

    top_view = render_ascii_structure(atoms, coords, width=38, height=18, view_plane='xy', title="TOP VIEW (XY)")
    side_view = render_ascii_structure(atoms, coords, width=38, height=18, view_plane='xz', title="SIDE VIEW (Stacking)")
    
    # Combine side-by-side
    top_lines = top_view.split('\n')
    side_lines = side_view.split('\n')
    
    combined = []
    for t, s in zip(top_lines, side_lines):
        combined.append(f"{t}  {s}")
        
    return '\n'.join(combined)


def format_test_result(nl_input: str, atoms: List[str], coords: List[Tuple], 
                       formula: str = None, stacking: str = None) -> str:
    """Format test result."""
    
    # Generate Dual View
    top_view = render_ascii_structure(atoms, coords, width=44, height=20, view_plane='xy', title="TOP VIEW (XY)")
    side_view = render_ascii_structure(atoms, coords, width=44, height=20, view_plane='xz', title="SIDE VIEW (Stack/Height)")
    
    top_lines = top_view.split('\n')
    side_lines = side_view.split('\n')
    
    output = []
    output.append("=" * 92)
    output.append(f" NL INPUT: {nl_input}")
    if formula:
        output.append(f" RESULT: Formula: {formula} | Atoms: {len(atoms)} | Stacking: {stacking or 'N/A'}")
    else:
        output.append(f" RESULT: Atoms: {len(atoms)}")
    output.append("=" * 92)
    
    for t, s in zip_longest(top_lines, side_lines, fillvalue=" "*44):
        output.append(f" {t}   {s}")
        
    output.append("=" * 92 + "\n")
    return '\n'.join(output)

from itertools import zip_longest

if __name__ == "__main__":
    # Test benzene stacking
    # Create two benzene rings stacked in Z
    atoms = []
    coords = []
    # Ring 1 at z=0
    for i in range(6):
        angle = i * 60 * 3.14159 / 180
        atoms.append('C')
        coords.append((1.4*np.cos(angle), 1.4*np.sin(angle), 0.0))
        atoms.append('H')
        coords.append((2.4*np.cos(angle), 2.4*np.sin(angle), 0.0))
        
    # Ring 2 at z=3.5
    for i in range(6):
        angle = i * 60 * 3.14159 / 180
        atoms.append('C')
        coords.append((1.4*np.cos(angle), 1.4*np.sin(angle), 3.5))
        atoms.append('H')
        coords.append((2.4*np.cos(angle), 2.4*np.sin(angle), 3.5))
        
    print(format_test_result("Test Stacking", atoms, coords, "stack example", "pi_pi"))

