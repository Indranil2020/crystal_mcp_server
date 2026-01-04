# Crystal MCP Server - Tool Registry

This document provides a complete inventory of all tools exposed by the Crystal MCP Server.

## Overview

The MCP server exposes **26 tools** across 5 categories:
- Generation (10 tools)
- Transformation (7 tools)
- Analysis (3 tools)
- Export (2 tools)
- Optimization (4 tools)

All tools are defined in `src/types/tools.ts` and their schemas are automatically converted to JSON Schema format for MCP protocol compliance.

---

## Generation Tools

### 1. generate_crystal
Generate a crystal structure with specified composition and space group using PyXtal.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| composition | array | Yes | Chemical composition as array of element symbols |
| space_group | number/string | Yes | Space group number (1-230) or Hermann-Mauguin symbol |
| num_atoms | number | No | Total number of atoms in the unit cell |
| lattice_params | object | No | Lattice parameters (a, b, c, alpha, beta, gamma) |
| volume_factor | number | No | Relative volume factor (default: 1.0) |
| min_distance | object | No | Minimum distances between element pairs |
| wyckoff_positions | array | No | Specify exact Wyckoff positions |
| seed | number | No | Random seed for reproducibility |
| dimensionality | number | No | 0=cluster, 1=rod, 2=slab, 3=bulk |

### 2. generate_prototype
Generate common prototype structures (rocksalt, perovskite, zincblende, etc.).

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| prototype | string | Yes | Prototype name: rocksalt, zincblende, wurtzite, fluorite, perovskite, spinel, heusler, rutile, diamond |
| elements | object | Yes | Element mapping, e.g., {"A": "Na", "B": "Cl"} |
| lattice_constant | number | No | Lattice constant in Angstroms |
| c_over_a | number | No | c/a ratio for non-cubic systems |

### 3. generate_nanostructure
Generate nanostructures (nanotubes, graphene, ribbons, fullerenes).

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| type | string | Yes | Type: nanotube, graphene, nanoribbon, fullerene, mos2, nanowire |
| params | object | Yes | Type-specific parameters |

**Type-specific params:**
- **nanotube**: n, m, length, element
- **graphene**: size [x, y, z], vacuum
- **nanoribbon**: width, length, edge_type (zigzag/armchair)
- **fullerene**: n_atoms (60, 70, 80)

### 4. generate_2d_material
Generate 2D materials (hBN, MoS2, WS2, phosphorene, silicene, MXene).

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| material | string | Yes | Material type: hBN, MoS2, WS2, phosphorene, silicene, MXene |
| size | array | No | Supercell size [a, b, c] |
| vacuum | number | No | Vacuum thickness in Angstroms |
| extra_params | object | No | Material-specific parameters |

### 5. generate_molecular_crystal
Generate a molecular crystal structure using PyXtal.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| molecule | string | Yes | Molecule identifier (SMILES, name, or formula) |
| space_group | number/string | Yes | Space group |
| num_molecules | number | No | Number of molecules in unit cell |

### 6. build_molecule
Build individual molecules from SMILES or common names.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| identifier | string | Yes | SMILES string or molecule name |
| optimize | boolean | No | Whether to optimize geometry |

### 7. generate_twisted_bilayer
Generate twisted bilayer/multilayer structures.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| material | string | Yes | Base material: graphene, MoS2, hBN |
| twist_angle | number | Yes | Twist angle in degrees |
| layers | number | No | Number of layers (default: 2) |

### 8. generate_high_entropy_alloy
Generate high-entropy alloy structures with 4+ elements.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| elements | array | Yes | Array of 4+ element symbols |
| lattice_type | string | Yes | fcc, bcc, or hcp |
| lattice_constant | number | No | Lattice constant in Angstroms |
| supercell | array | No | Supercell dimensions |

### 9. generate_mof
Generate metal-organic framework structures.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| mof_type | string | Yes | MOF type: MOF-5, HKUST-1, UiO-66, ZIF-8 |
| metal | string | No | Metal center element |
| supercell | array | No | Supercell dimensions |

### 10. generate_cage
Generate cage structures (fullerenes, clathrates).

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| cage_type | string | Yes | Type: fullerene, clathrate |
| n_atoms | number | No | Number of atoms (for fullerene: 60, 70, 80) |
| guest | string | No | Guest atom/molecule for clathrate |

---

## Transformation Tools

### 1. generate_slab
Generate a surface slab from a bulk structure.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure (or prototype definition) |
| miller_indices | array | Yes | Miller indices [h, k, l] |
| thickness | number | Yes | Number of atomic layers |
| vacuum | number | No | Vacuum thickness in Angstroms |
| symmetric | boolean | No | Whether to create symmetric slab |
| center_slab | boolean | No | Center slab in cell |

### 2. make_supercell
Create a supercell from an existing structure.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| scaling_matrix | array | Yes | 3x3 matrix or [a, b, c] scaling |

### 3. create_defect
Create point defects in a crystal structure.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| defect_type | string | Yes | vacancy, substitution, or interstitial |
| site_index | number | No | Atom index for defect |
| element | string | No | Element for substitution/interstitial |

### 4. create_alloy
Create a substitutional alloy with random mixing.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| substitutions | object | Yes | Element substitution map |
| concentration | number | No | Substitution concentration |

### 5. create_heterostructure
Create a vertical heterostructure by stacking layers.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| layer1 | object | Yes | First layer structure |
| layer2 | object | Yes | Second layer structure |
| distance | number | No | Interlayer distance |
| vacuum | number | No | Vacuum thickness |

### 6. add_adsorbate
Add an adsorbate molecule to a surface.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Surface structure |
| adsorbate | string | Yes | Adsorbate molecule (SMILES or name) |
| site | string | No | Adsorption site: ontop, bridge, hollow |
| height | number | No | Adsorbate height above surface |

### 7. apply_strain
Apply strain tensor to a structure.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| strain | array | Yes | 3x3 strain tensor or 6 Voigt components |

---

## Analysis Tools

### 1. analyze_symmetry
Analyze and detect symmetry properties.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| symprec | number | No | Symmetry precision (default: 0.01) |

**Returns:** Space group, point group, Wyckoff positions, symmetry operations

### 2. validate_structure
Validate a crystal structure for physical/chemical correctness.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| check_distances | boolean | No | Check interatomic distances |
| check_coordination | boolean | No | Check coordination numbers |

### 3. analyze_symmetry_relations
Explore group-subgroup relations.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| space_group | number | Yes | Space group number |
| relation_type | string | No | subgroups or supergroups |

---

## Export Tools

### 1. export_structure
Export crystal structure to multiple file formats.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| formats | array | Yes | Output formats: cif, poscar, xyz, json |
| output_directory | string | No | Directory to save files |

### 2. generate_visualization
Generate interactive HTML or static PNG visualization.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| format | string | No | Output format: html, png |
| style | string | No | Visualization style: ball_stick, spacefill |

---

## Optimization Tools

### 1. optimize_structure_mlff
Optimize crystal structure using machine learning force fields.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| model | string | No | MLFF model: chgnet, m3gnet, mace |
| fmax | number | No | Force convergence criterion |
| steps | number | No | Maximum optimization steps |

### 2. calculate_energy_mlff
Calculate energy, forces, and stress using MLFF.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| structure | object | Yes | Input structure |
| model | string | No | MLFF model: chgnet, m3gnet, mace |

### 3. ground_state_search
Search for ground state across multiple space groups.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| composition | array | Yes | Chemical composition |
| space_groups | array | No | Space groups to search |
| model | string | No | MLFF model for energy ranking |

### 4. generate_space_group_scan
Generate structures across multiple space groups.

**Parameters:**
| Name | Type | Required | Description |
|------|------|----------|-------------|
| composition | array | Yes | Chemical composition |
| space_groups | array | Yes | List of space group numbers |
| num_atoms | number | No | Target number of atoms |

---

## How Tools Are Exposed

### MCP Protocol

When a client connects and calls `tools/list`, the server returns:

```json
{
  "tools": [
    {
      "name": "generate_prototype",
      "description": "Generate common prototype structures...",
      "inputSchema": {
        "type": "object",
        "properties": {
          "prototype": {
            "type": "string",
            "description": "Prototype name..."
          },
          ...
        },
        "required": ["prototype", "elements"]
      }
    },
    ...
  ]
}
```

### Schema Conversion

Tool schemas are defined using Zod in TypeScript and automatically converted to JSON Schema using `zod-to-json-schema`. This ensures:

1. Type safety in the server code
2. Standard JSON Schema for MCP protocol
3. Automatic validation of inputs

---

## Version

- Document Version: 1.0
- Last Updated: 2026-01-04
- Total Tools: 26
