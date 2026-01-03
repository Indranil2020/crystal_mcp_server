#!/usr/bin/env python3
"""
Visualization Tool
Generates interactive HTML (3Dmol.js) and static PNG (ASE) visualizations.
"""

import sys
import json
import numpy as np
import io
from typing import Dict, Any, Optional
from ase import Atoms
from ase.io import write
from ase.build import make_supercell
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

def dict_to_atoms(structure_dict: Dict[str, Any]) -> Optional[Atoms]:
    """Convert structure dictionary to ASE Atoms object."""
    if not structure_dict or "lattice" not in structure_dict or "atoms" not in structure_dict:
        return None
    
    lattice = structure_dict["lattice"]
    if "matrix" not in lattice:
        return None
    
    cell = np.array(lattice["matrix"])
    atoms_data = structure_dict["atoms"]
    
    symbols = []
    positions = []
    
    for atom in atoms_data:
        if "element" not in atom or "cartesian" not in atom:
            return None
        symbols.append(atom["element"])
        positions.append(atom["cartesian"])
    
    if not symbols:
        return None
        
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

def generate_html(atoms: Atoms, output_file: str, style: str = "ball-stick"):
    """Generate interactive HTML using 3Dmol.js."""
    # Convert to XYZ format string for embedding
    xyz_io = io.StringIO()
    write(xyz_io, atoms, format='xyz')
    xyz_data = xyz_io.getvalue().replace('`', '\\`') # Escape backticks if any
    
    # 3Dmol.js template using vanilla JavaScript for compatibility
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Crystal Structure Visualization</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; overflow: hidden; }}
        #container {{ width: 100vw; height: 100vh; position: relative; }}
    </style>
</head>
<body>
    <div id="container"></div>
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            let element = document.getElementById('container');
            let config = {{ backgroundColor: "white" }};
            let viewer = $3Dmol.createViewer(element, config);

            let data = `{xyz_data}`;
            viewer.addModel(data, "xyz");

            viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.25}}}});
            viewer.addUnitCell(viewer.getModel());

            viewer.zoomTo();
            viewer.render();
        }});
    </script>
</body>
</html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    return {"format": "html", "path": output_file}

def generate_image(atoms: Atoms, output_file: str, rotation: str = '10x,10y,10z'):
    """Generate static PNG using ASE/Matplotlib."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Clean up axes
    ax.axis('off')
    
    # Plot
    plot_atoms(atoms, ax, radii=0.3, rotation=rotation)
    
    # Save
    fig.savefig(output_file, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    return {"format": "png", "path": output_file}

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "error": "Usage: python visualization.py <input.json>"}))
        sys.exit(1)
        
    input_file = sys.argv[1]
    
    with open(input_file, 'r') as f:
        params = json.load(f)
        
    structure = params.get("structure")
    output_file = params.get("output_file", "visualization.html")
    fmt = params.get("format", "html").lower()
    
    atoms = dict_to_atoms(structure)
    if atoms is None:
        print(json.dumps({"success": False, "error": "Invalid structure data"}))
        sys.exit(1)
        
    if fmt == "html":
        result = generate_html(atoms, output_file)
    elif fmt == "png":
        result = generate_image(atoms, output_file)
    else:
        print(json.dumps({"success": False, "error": f"Unknown format: {fmt}"}))
        sys.exit(1)
        
    print(json.dumps({"success": True, "result": result}))

if __name__ == "__main__":
    main()
