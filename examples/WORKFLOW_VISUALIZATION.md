# Crystal MCP Server - Complete Workflow Visualization

## Project Overview

The **Crystal MCP Server** is a production-ready Model Context Protocol (MCP) server that enables AI assistants (like Claude) to generate, analyze, and export crystal structures for computational materials science.

## Architecture Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           USER INTERACTION                                   │
│                                                                              │
│  Natural Language Query: "Generate a silicon crystal with diamond structure"│
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         LLM (Claude) + MCP                                  │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  Intent Recognition:                                                 │   │
│  │  • Tool: generate_crystal                                            │   │
│  │  • Composition: ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"]    │   │
│  │  • Space Group: 227 (Fd-3m, diamond cubic)                          │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      MCP SERVER (TypeScript/Node.js)                        │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  • stdio transport (stdin/stdout)                                    │   │
│  │  • Zod input validation                                              │   │
│  │  • Tool routing: 27+ MCP tools                                       │   │
│  │  • JSON serialization to temp file                                   │   │
│  │  • Python subprocess invocation                                      │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     PYTHON BACKEND (Scientific Core)                        │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  PyXtal: Crystal structure generation with all 230 space groups     │   │
│  │  Pymatgen: Materials analysis and structure manipulation            │   │
│  │  Spglib: Symmetry detection and analysis                            │   │
│  │  ASE: Atomic Simulation Environment for I/O and MLFF                │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                              │
│  Processing Steps:                                                           │
│  1. Validate composition against space group Wyckoff positions              │
│  2. Initialize PyXtal with species and numIons                              │
│  3. Generate structure with from_random()                                   │
│  4. Validate distances, symmetry, stoichiometry                             │
│  5. Extract lattice, atoms, and metadata                                    │
│  6. Return structured JSON response                                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          OUTPUT FORMATS                                      │
│                                                                              │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐  ┌────────────────────────┐ │
│  │    CIF     │  │   POSCAR   │  │    XYZ     │  │     JSON Complete      │ │
│  │            │  │   (VASP)   │  │            │  │                        │ │
│  │ • Lattice  │  │ • Lattice  │  │ • Atoms    │  │ • Full structure       │ │
│  │ • Atoms    │  │ • Species  │  │ • Coords   │  │ • Symmetry analysis    │ │
│  │ • Symmetry │  │ • Coords   │  │            │  │ • Validation           │ │
│  └────────────┘  └────────────┘  └────────────┘  │ • MCP request          │ │
│                                                   └────────────────────────┘ │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                    HTML VISUALIZATION (3Dmol.js)                       │ │
│  │  • Interactive 3D viewer                                               │ │
│  │  • Unit cell display                                                   │ │
│  │  • Multiple render styles (ball-stick, space-filling, stick)          │ │
│  │  • Crystal property sidebar                                            │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Available MCP Tools (27 Total)

### Generation Tools
| Tool | Description |
|------|-------------|
| `generate_crystal` | Generate crystal with composition + space group |
| `comprehensive_generate` | Unified entry for 50+ operations across 18 categories |
| `generate_molecular_crystal` | Molecular crystal structures |
| `generate_nanostructure` | Nanotubes, graphene, ribbons, fullerenes |
| `generate_prototype` | Common prototypes (rocksalt, perovskite, etc.) |
| `build_molecule` | Isolated molecular structures |
| `generate_space_group_scan` | Scan multiple space groups |

### Transformation Tools
| Tool | Description |
|------|-------------|
| `make_supercell` | Create supercells with scaling matrix |
| `generate_slab` | Surface slabs with Miller indices |
| `create_defect` | Point defects (vacancy, substitution) |
| `create_alloy` | Substitutional alloys |
| `create_heterostructure` | Stack layers vertically |
| `add_adsorbate` | Adsorbate molecules on surfaces |
| `apply_strain` | Strain tensors |
| `generate_twisted_bilayer` | Twisted bilayers (graphene, MoS2) |

### Analysis & Optimization
| Tool | Description |
|------|-------------|
| `analyze_symmetry` | Space group, point group detection |
| `validate_structure` | Comprehensive validation |
| `optimize_structure_mlff` | MLFF optimization (CHGNet, M3GNet, MACE) |
| `calculate_energy_mlff` | Energy, forces, stress calculation |
| `ground_state_search` | Search across space groups |

### Export & Visualization
| Tool | Description |
|------|-------------|
| `export_structure` | CIF, POSCAR, XYZ, JSON |
| `generate_visualization` | Interactive HTML or static PNG |

## Real Demo Results (8 Structures Generated)

| Structure | Formula | Space Group | Crystal System | Volume (Å³) | Density (g/cm³) |
|-----------|---------|-------------|----------------|-------------|-----------------|
| Silicon Diamond | Si8 | Fd-3m (#227) | cubic | 185.46 | 2.012 |
| NaCl Rocksalt | Na4Cl4 | Fm-3m (#225) | cubic | 146.41 | 2.651 |
| BaTiO3 Perovskite | Ba1Ti1O3 | Pm-3m (#221) | cubic | 93.76 | 4.130 |
| Graphite | C4 | P63/mmc (#194) | hexagonal | 63.13 | 1.264 |
| hBN | B2N2 | P63/mmc (#194) | hexagonal | 56.07 | 1.470 |
| Iron BCC | Fe2 | Im-3m (#229) | cubic | 46.21 | 4.014 |
| ZnO Wurtzite | Zn2O2 | P63mc (#186) | hexagonal | 22.09 | 12.237 |
| TiO2 Rutile | Ti2O4 | P42/mnm (#136) | tetragonal | 54.09 | 4.904 |

## Output File Structure

```
examples/real_outputs/
├── silicon_diamond.cif          # Crystallographic Information File
├── silicon_diamond.vasp         # VASP POSCAR format
├── silicon_diamond.xyz          # XYZ coordinates
├── silicon_diamond_complete.json # Full JSON with all metadata
├── silicon_diamond_viz.html     # Interactive 3D visualization
├── nacl_rocksalt.cif
├── nacl_rocksalt.vasp
├── nacl_rocksalt.xyz
├── nacl_rocksalt_complete.json
├── nacl_rocksalt_viz.html
├── batio3_perovskite.cif
├── ... (5 more structures × 5 files each)
└── tio2_rutile_viz.html
```

## Usage with DFT Codes

| DFT Code | File Format | Notes |
|----------|-------------|-------|
| VASP | `*.vasp` | Direct POSCAR input |
| Quantum ESPRESSO | `*.cif` | Import via `cif2cell` |
| LAMMPS | `*.xyz` | Convert via `atomsk` |
| VESTA | `*.cif` | Visualization |
| Materials Studio | `*.cif` | Full import |
| Avogadro | `*.xyz` | Visualization |
| OVITO | `*.vasp`, `*.xyz` | Visualization & analysis |

## Project Ready for Publication

The Crystal MCP Server is production-ready with:

- **27+ MCP tools** covering the full crystal structure workflow
- **All 230 space groups** supported via PyXtal
- **228+ generator operations** across 18 material categories
- **Multi-format export** (CIF, POSCAR, XYZ, JSON, HTML)
- **Interactive visualization** with 3Dmol.js
- **Comprehensive validation** including distance checks, symmetry verification
- **MLFF optimization** support (CHGNet, M3GNet, MACE)
- **Defensive programming** without try/except for control flow
- **Full type safety** via TypeScript + Zod validation

To use with Claude Code or other MCP-compatible LLMs, simply configure the server in your MCP settings.
