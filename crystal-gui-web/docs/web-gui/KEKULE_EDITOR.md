# Kekule 2D Molecular Editor

The web GUI includes a professional-grade 2D molecular editor powered by [Kekule.js](http://partridgejiang.github.io/Kekule.js/).

## Features

### Drawing Tools

| Tool | Description |
|------|-------------|
| **Manipulate** | Select and move atoms/bonds |
| **Erase** | Delete selected objects |
| **Bond** | Draw single, double, triple bonds |
| **Atom & Formula** | Input atom symbols or formulas |
| **Ring** | Insert ring structures (3-6 membered, benzene) |
| **Charge** | Add +/- charges to atoms |
| **Glyph** | Reaction arrows and symbols |
| **Text & Image** | Add annotations |

### Templates

Pre-built molecular templates available from the dropdown:

- **Benzene** (C₆H₆) - Aromatic ring
- **Cyclohexane** (C₆H₁₂) - Saturated 6-membered ring
- **Pyridine** (C₅H₅N) - Heteroaromatic with nitrogen
- **Naphthalene** (C₁₀H₈) - Fused aromatic rings
- **Phenol** (C₆H₅OH) - Benzene with hydroxyl group
- **Aniline** (C₆H₅NH₂) - Benzene with amino group

### File Operations

| Button | Action |
|--------|--------|
| **Clear** | Reset the editor canvas |
| **Export** | Download current molecule as `.mol` file |
| **Import** | Load a `.mol` or `.sdf` file |
| **→ 3D** | Push molecule to the 3D Mol* viewer |

## Usage Guide

### Drawing a Molecule

1. Select the **Bond** tool from the vertical toolbar.
2. Click and drag on the canvas to draw bonds.
3. Atoms default to Carbon; use **Atom & Formula** to change elements.
4. Use the **Ring** tool to quickly insert ring structures.

### Using Templates

1. Select a molecule from the **Templates** dropdown.
2. The molecule will appear in the 2D editor.
3. Edit as needed, then click **→ 3D** to visualize.

### Exporting Molecules

1. Draw or load your molecule.
2. Click **Export** to download as MDL MOL V2000 format.
3. The file can be opened in other chemistry software.

### Importing Molecules

1. Click **Import** and select a `.mol` or `.sdf` file.
2. The molecule will load into the 2D editor.
3. Make edits if needed, then push to 3D.

### 3D Visualization

1. Draw or load your molecule in the 2D editor.
2. Click **→ 3D** to transfer to the Mol* viewer.
3. The 3D viewer will display:
   - Ball-and-stick representation
   - Aromatic bonds as dashed lines
   - Correct 3D geometry

## Technical Notes

- Molecules are exported in **MDL MOL V2000** format for maximum compatibility.
- The 2D→3D pipeline preserves:
  - Bond orders (single, double, triple, aromatic)
  - Atom types and charges
  - Stereochemistry (when present)
- Mol* automatically interprets the MOL data and generates 3D coordinates.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Editor not loading | Check browser console for Kekule.js errors |
| Template not appearing | Ensure dropdown value changed; try clicking again |
| 3D push fails | Check console for "[KekuleEditor]" messages |
| Bonds not visible in 3D | Verify MOL export includes "V2000" in console |

## Related Documentation

- [Web GUI Architecture](./ARCHITECTURE.md)
- [Components Overview](./COMPONENTS.md)
- [State Management](./STATE.md)
