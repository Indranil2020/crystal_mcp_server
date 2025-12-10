# Crystal Structure Generation MCP Server - Ultrafine Implementation Plan

## Executive Summary

A comprehensive MCP server for generating accurate crystal structures for DFT, quantum chemistry, and condensed matter physics calculations using PyXtal as the core engine, with support for all 230 space groups, complex structure transformations (supercells, slabs, surfaces), symmetry-constrained optimization workflows, and integration with MLFF (Machine Learning Force Fields) for ground state energy calculations.

---

## Phase 1: Deep Research and Requirements Analysis (Week 1-2)

### 1.1 Core Research Areas

#### A. PyXtal Capabilities Deep Dive
- **Study PyXtal API documentation thoroughly**
  - Crystal generation from scratch (all 230 space groups)
  - Wyckoff position handling (general and special positions)
  - Molecular crystal generation capabilities
  - Sub-periodic systems (0D clusters, 1D rods, 2D layers, 3D bulk)
  - Structure manipulation functions
  - Symmetry analysis tools
  - Export formats (CIF, POSCAR, XYZ, etc.)

- **Understand PyXtal limitations**
  - Maximum system sizes
  - Computational constraints
  - Edge cases in space group generation
  - Molecule placement constraints

#### B. Space Group Theory & Crystallography
- **International Tables for Crystallography**
  - All 230 space group specifications
  - Wyckoff positions for each space group
  - Symmetry operations and generators
  - Hermann-Mauguin notation
  - Schoenflies notation
  - Crystal systems and lattice parameters

- **Wyckoff Position Theory**
  - Site symmetry
  - Multiplicity
  - Special vs general positions
  - Compatible occupancies for different elements

#### C. Structure Generation Constraints
- **Lattice parameter constraints**
  - Crystal system-specific constraints (cubic, tetragonal, etc.)
  - Reasonable ranges for different materials
  - Volume per atom considerations
  - Cell optimization strategies

- **Atomic distance constraints**
  - Minimum interatomic distances (element-dependent)
  - Coordination number constraints
  - Bond length reasonableness checks

#### D. Modern Crystal Generation Methods
- **Study DiffCSP++ paper (arXiv:2402.03992)**
  - Space group constrained diffusion models
  - Lattice matrix invariant logarithmic space
  - Wyckoff position constraints in generative models
  - How to integrate ML-based generation with PyXtal

- **Related Methods**
  - CDVAE (Crystal Diffusion VAE)
  - Flow-based models (CrystalFlow)
  - Autoregressive LLM approaches
  - Symmetry-preserving generation

### 1.2 MCP Architecture Research

#### A. Load MCP Documentation
```
- https://modelcontextprotocol.io/sitemap.xml
- TypeScript SDK README from GitHub
- MCP Best Practices guide
```

#### B. Tool Design Philosophy
- **Comprehensive API coverage** vs **High-level workflows**
  - Basic tools: Generate single structure
  - Workflow tools: Generate full space group scan
  - Analysis tools: Symmetry detection, validation
  - Transformation tools: Supercell, slab, surface generation
  - Export tools: Multiple format support

#### C. Transport Selection
- **Recommendation: Streamable HTTP**
  - Stateless design (important for batch jobs)
  - Easier scaling for compute-intensive structure generation
  - Better for remote computational resources
- **Alternative: stdio for local use**

### 1.3 Integration Requirements

#### A. File Format Support
- **Input formats:**
  - CIF (Crystallographic Information File)
  - POSCAR/CONTCAR (VASP)
  - XYZ (simple Cartesian coordinates)
  - PDB (protein data bank, for molecules)
  - JSON (custom internal format)

- **Output formats:**
  - CIF (standard crystallography)
  - POSCAR/CONTCAR (VASP DFT)
  - XYZ (visualization)
  - extXYZ (extended with lattice)
  - ASE atoms format
  - JSON (with full metadata)
  - PDB (for molecular crystals)

#### B. DFT Software Integration Paths
- **VASP:** POSCAR/POTCAR generation
- **Quantum ESPRESSO:** PWscf input format
- **CP2K:** Coordinate and cell input
- **CASTEP:** Cell and param files
- **LAMMPS:** Data file format (for MLFF)

#### C. MLFF Framework Compatibility
- **CHGNet, M3GNet, MACE, SevenNet, etc.**
  - ASE calculator interface
  - PyTorch Geometric compatibility
  - Structure relaxation workflows
  - Energy/force calculation integration

---

## Phase 2: Core Tool Architecture Design (Week 2-3)

### 2.1 Tool Categories & Hierarchy

#### Category 1: Basic Structure Generation Tools

**Tool: `generate_crystal`**
```typescript
Input Schema:
{
  composition: string[],          // e.g., ["Si", "Si"] or ["Na", "Cl"]
  space_group: number | string,   // 1-230 or Hermann-Mauguin symbol
  num_atoms?: number,             // Optional constraint on total atoms
  lattice_params?: {
    a?: number, b?: number, c?: number,
    alpha?: number, beta?: number, gamma?: number
  },
  volume_factor?: number,         // Relative volume (default: 1.0)
  min_distance?: Record<string, number>, // Min distances between elements
  wyckoff_positions?: Array<{
    element: string,
    wyckoff: string,              // e.g., "4a", "8c"
    coords?: number[]             // Optional fractional coordinates
  }>,
  seed?: number,                  // Random seed for reproducibility
  max_attempts?: number           // Max attempts to generate valid structure
}

Output:
{
  structure: {
    lattice: {
      a: number, b: number, c: number,
      alpha: number, beta: number, gamma: number,
      matrix: number[][]          // 3x3 lattice vectors
    },
    atoms: Array<{
      element: string,
      coords: number[],           // Fractional coordinates
      cartesian: number[],        // Cartesian coordinates
      wyckoff: string,
      multiplicity: number,
      site_symmetry: string
    }>,
    space_group: {
      number: number,
      symbol: string,
      hall_symbol: string,
      point_group: string,
      crystal_system: string
    },
    metadata: {
      formula: string,
      natoms: number,
      volume: number,
      density: number,
      packing_fraction: number
    }
  },
  files: {
    cif: string,
    poscar: string,
    xyz: string
  },
  validation: {
    valid: boolean,
    issues: string[],
    warnings: string[]
  }
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true (with same seed)
- openWorldHint: true
```

**Tool: `generate_space_group_scan`**
```typescript
Input Schema:
{
  composition: string[],
  space_groups?: number[],        // Specific groups, or all 230 if null
  space_group_range?: [number, number], // Range like [1, 230]
  crystal_systems?: string[],     // Filter by system: cubic, hexagonal, etc.
  num_atoms?: number,
  volume_factor?: number,
  parallel?: boolean,             // Enable parallel generation
  output_directory?: string,
  naming_scheme?: string          // Template for filenames
}

Output:
{
  generated_structures: Array<{
    space_group: number,
    success: boolean,
    structure?: {...},            // Same as generate_crystal output
    file_path?: string,
    error_message?: string,
    generation_time: number
  }>,
  summary: {
    total_attempted: number,
    successful: number,
    failed: number,
    total_time: number,
    statistics: {
      by_crystal_system: Record<string, number>,
      avg_atoms: number,
      avg_volume: number
    }
  }
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: true
```

#### Category 2: Structure Transformation Tools

**Tool: `make_supercell`**
```typescript
Input Schema:
{
  structure: string | object,     // CIF file, JSON structure, or structure object
  scaling_matrix: number[][] | number[], // 3x3 matrix or [nx, ny, nz]
  wrap_atoms?: boolean,           // Wrap atoms back into cell
  preserve_symmetry?: boolean,    // Try to preserve space group if possible
  min_distance_check?: boolean    // Check for overlapping atoms
}

Output:
{
  supercell: {...},               // Structure in same format
  transformation_matrix: number[][],
  original_space_group: number,
  new_space_group?: number,       // May be P1 if symmetry lost
  volume_multiplier: number,
  files: {...}
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `generate_slab`**
```typescript
Input Schema:
{
  structure: string | object,
  miller_indices: number[],       // [h, k, l]
  thickness: number,              // Number of layers
  vacuum: number,                 // Vacuum thickness in Angstroms
  center_slab?: boolean,
  fix_bottom_layers?: number,     // For DFT relaxation
  fix_atoms?: number[],           // Specific atom indices to fix
  orthogonalize?: boolean,        // Create orthogonal cell
  min_slab_size?: number,         // Minimum slab size in Angstroms
  symmetric?: boolean             // Make symmetric slab
}

Output:
{
  slab: {...},
  surface_area: number,
  slab_thickness: number,
  vacuum_thickness: number,
  termination: string,
  miller_indices: number[],
  fixed_atoms: number[],
  metadata: {
    original_space_group: number,
    slab_space_group: number,     // Usually P1
    n_original_atoms: number,
    n_slab_atoms: number
  },
  files: {...}
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `generate_surface_cell`**
```typescript
Input Schema:
{
  structure: string | object,
  miller_indices: number[],
  layers: number,                 // Number of layers to include
  surface_reconstruction?: string, // Type of reconstruction if any
  adsorbate?: {
    molecule: string,             // SMILES or molecule name
    sites: string[],              // Adsorption sites: "top", "bridge", "hollow"
    coverage?: number             // ML coverage (0-1)
  }
}

Output: {...}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `create_heterostructure`**
```typescript
Input Schema:
{
  substrate: string | object,
  overlayer: string | object,
  interface_distance: number,     // Distance between layers in Angstroms
  strain_tolerance?: number,      // Maximum acceptable strain
  rotation_angles?: number[],     // Try different rotations
  max_area_mismatch?: number      // Maximum area difference percentage
}

Output: {...}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

#### Category 3: Symmetry Analysis Tools

**Tool: `analyze_symmetry`**
```typescript
Input Schema:
{
  structure: string | object,
  symprec?: number,               // Symmetry precision (default: 1e-3)
  angle_tolerance?: number,       // Angle tolerance in degrees
  detect_primitive?: boolean,     // Find primitive cell
  standardize?: boolean           // Standardize to conventional cell
}

Output:
{
  space_group: {
    number: number,
    symbol: string,
    hall_symbol: string,
    point_group: string,
    crystal_system: string,
    international_number: number
  },
  symmetry_operations: {
    rotations: number[][][],
    translations: number[][],
    n_operations: number
  },
  wyckoff_positions: Array<{
    element: string,
    wyckoff: string,
    multiplicity: number,
    site_symmetry: string,
    coords: number[]
  }>,
  primitive_cell?: {...},
  conventional_cell?: {...},
  is_standardized: boolean
}

Annotations:
- readOnlyHint: true
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `validate_structure`**
```typescript
Input Schema:
{
  structure: string | object,
  checks: string[],               // ["distances", "symmetry", "stoichiometry", ...]
  min_distance?: number,          // Minimum allowed distance
  max_coordination?: number,
  expected_space_group?: number
}

Output:
{
  valid: boolean,
  errors: Array<{
    type: string,
    severity: "error" | "warning",
    message: string,
    details: any
  }>,
  metrics: {
    min_distance: number,
    avg_distance: number,
    density: number,
    packing_fraction: number
  },
  suggestions: string[]
}

Annotations:
- readOnlyHint: true
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

#### Category 4: Molecular Crystal Tools

**Tool: `generate_molecular_crystal`**
```typescript
Input Schema:
{
  molecules: Array<string>,       // SMILES strings or molecule names
  space_group: number | string,
  num_molecules: number[],        // Number of each molecule type
  volume_factor?: number,
  allow_special_positions?: boolean,
  orientational_constraints?: Array<{
    molecule_index: number,
    constraints: string           // e.g., "C2v symmetry"
  }>
}

Output: {...}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: true
- openWorldHint: true
```

#### Category 5: Workflow & Batch Tools

**Tool: `optimize_structure_mlff`**
```typescript
Input Schema:
{
  structure: string | object,
  mlff_model: string,             // "chgnet", "m3gnet", "mace", etc.
  optimizer: string,              // "BFGS", "FIRE", "LBFGS"
  fmax: number,                   // Force convergence criterion
  steps: number,                  // Maximum optimization steps
  constrain_symmetry?: boolean,   // Preserve space group during optimization
  fix_lattice?: boolean,          // Fix lattice, only relax positions
  fix_volume?: boolean,           // Fix volume, allow shape changes
  pressure?: number,              // External pressure in GPa
  trajectory_file?: string        // Save optimization trajectory
}

Output:
{
  optimized_structure: {...},
  initial_energy: number,
  final_energy: number,
  energy_change: number,
  max_force_initial: number,
  max_force_final: number,
  n_steps: number,
  converged: boolean,
  preserved_symmetry: boolean,
  final_space_group?: number,
  trajectory?: Array<{
    step: number,
    energy: number,
    max_force: number,
    structure: {...}
  }>,
  timing: {
    setup_time: number,
    optimization_time: number,
    total_time: number
  }
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: false
- openWorldHint: false
```

**Tool: `calculate_energy_mlff`**
```typescript
Input Schema:
{
  structure: string | object,
  mlff_model: string,
  calculate_forces?: boolean,
  calculate_stress?: boolean
}

Output:
{
  energy: number,
  energy_per_atom: number,
  forces?: number[][],
  stress?: number[][],
  metadata: {
    model: string,
    model_version?: string,
    calculation_time: number
  }
}

Annotations:
- readOnlyHint: true
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `ground_state_search`**
```typescript
Input Schema:
{
  composition: string[],
  space_groups?: number[],        // Search specific groups or all 230
  num_structures_per_group: number, // Multiple attempts per space group
  mlff_model: string,
  optimization_settings: {
    optimizer: string,
    fmax: number,
    steps: number,
    constrain_symmetry: boolean
  },
  parallel: boolean,
  n_workers?: number,
  output_directory: string,
  save_trajectories?: boolean
}

Output:
{
  ground_state: {
    structure: {...},
    energy: number,
    space_group: number,
    formation_energy?: number
  },
  all_results: Array<{
    space_group: number,
    attempt: number,
    initial_structure: {...},
    optimized_structure: {...},
    initial_energy: number,
    final_energy: number,
    converged: boolean,
    preserved_symmetry: boolean
  }>,
  energy_ranking: Array<{
    rank: number,
    space_group: number,
    energy: number,
    energy_above_ground_state: number
  }>,
  statistics: {
    total_structures: number,
    converged_structures: number,
    unique_space_groups: number,
    energy_range: [number, number],
    most_stable_crystal_system: string
  },
  convex_hull?: {...}              // If doing alloy or compound
}

Annotations:
- readOnlyHint: false
- destructiveHint: false
- idempotentHint: false
- openWorldHint: false
```

#### Category 6: Advanced Structure Manipulation

**Tool: `apply_strain`**
```typescript
Input Schema:
{
  structure: string | object,
  strain_tensor: number[][],      // 3x3 strain tensor
  strain_type?: string,           // "biaxial", "uniaxial", "hydrostatic", "shear"
  strain_value?: number           // If using predefined strain type
}

Output: {...}
```

**Tool: `create_defect`**
```typescript
Input Schema:
{
  structure: string | object,
  defect_type: string,            // "vacancy", "substitution", "interstitial"
  defect_species?: string,
  defect_site?: number | number[], // Atom index or fractional coords
  concentration?: number,         // For multiple defects
  charge_state?: number
}

Output: {...}
```

**Tool: `create_alloy`**
```typescript
Input Schema:
{
  base_structure: string | object,
  substitutions: Record<string, {
    element: string,
    concentration: number
  }>,
  sqs?: boolean,                  // Use special quasi-random structure
  ordering_type?: string          // "random", "clustered", "ordered"
}

Output: {...}
```

#### Category 7: Export & Visualization Tools

**Tool: `export_structure`**
```typescript
Input Schema:
{
  structure: string | object,
  formats: string[],              // ["cif", "poscar", "xyz", "json", ...]
  include_metadata?: boolean,
  dft_software?: string,          // "vasp", "qe", "cp2k", etc.
  output_directory?: string
}

Output:
{
  files: Record<string, string>,  // format -> file_content mapping
  file_paths: Record<string, string>
}

Annotations:
- readOnlyHint: true
- destructiveHint: false
- idempotentHint: true
- openWorldHint: false
```

**Tool: `generate_visualization`**
```typescript
Input Schema:
{
  structure: string | object,
  style: string,                  // "ball-stick", "polyhedral", "wireframe"
  view_direction?: number[],      // Miller indices
  show_unit_cell?: boolean,
  show_axes?: boolean,
  atom_colors?: Record<string, string>,
  output_format?: string          // "png", "svg", "html"
}

Output:
{
  visualization: string,          // Base64 encoded or HTML
  file_path?: string
}
```

### 2.2 Error Handling Strategy

#### Error Categories
1. **Input Validation Errors**
   - Invalid space group number
   - Incompatible composition with space group
   - Invalid lattice parameters
   - Unrealistic constraints

2. **Generation Failures**
   - Cannot satisfy distance constraints
   - Cannot place atoms in Wyckoff positions
   - Maximum attempts exceeded
   - Numerical instabilities

3. **File I/O Errors**
   - Cannot read input file
   - Cannot write output file
   - Invalid file format
   - Corrupted structure data

4. **Computational Errors**
   - MLFF model not available
   - Optimization convergence failure
   - Out of memory errors
   - Timeout errors

#### Error Response Format
```typescript
{
  error: {
    code: string,                 // e.g., "INVALID_SPACE_GROUP"
    message: string,              // User-friendly message
    details: any,                 // Technical details
    suggestions: string[],        // How to fix
    recoverable: boolean
  }
}
```

#### Actionable Error Messages Examples
- ❌ "Generation failed"
- ✅ "Structure generation failed after 100 attempts. The minimum distance constraint of 1.5 Å is too restrictive for space group 225 with 8 Si atoms. Try: (1) Increase volume_factor to 1.2-1.5, (2) Reduce min_distance to 1.2 Å, or (3) Use fewer atoms (4-6)."

---

## Phase 3: Implementation Details (Week 3-6)

### 3.1 Technology Stack

#### Primary Language: TypeScript
**Rationale:**
- Best MCP SDK support
- Strong typing for crystallographic data
- Good package ecosystem
- Easy integration with scientific libraries via child processes

#### Core Dependencies
```json
{
  "@modelcontextprotocol/sdk": "latest",
  "zod": "^3.x",                  // Schema validation
  "axios": "^1.x",                // HTTP client if needed
  "uuid": "^9.x"                  // Generate unique IDs
}
```

#### Python Backend Integration
```json
{
  // Use child_process to call Python scripts
  "dependencies": {
    "python-shell": "^5.x"        // Bridge to Python
  }
}
```

#### Python Requirements (backend)
```python
# Core crystallography
pyxtal>=1.0.0
pymatgen>=2024.1.1
spglib>=2.0.0
ase>=3.22.0

# MLFF models
chgnet>=0.3.0
matgl>=1.0.0                      # M3GNet
mace-torch>=0.3.0

# Structure manipulation
numpy>=1.24.0
scipy>=1.10.0
pandas>=2.0.0

# File I/O
openbabel-wheel>=3.1.0           # For molecule handling
```

### 3.2 Project Structure

```
crystal-mcp-server/
├── src/
│   ├── index.ts                 # MCP server entry point
│   ├── server.ts                # Server initialization
│   ├── types/
│   │   ├── crystal.ts           # Crystal structure types
│   │   ├── tools.ts             # Tool input/output schemas
│   │   └── errors.ts            # Error types
│   ├── tools/
│   │   ├── generation/
│   │   │   ├── generate-crystal.ts
│   │   │   ├── space-group-scan.ts
│   │   │   └── molecular-crystal.ts
│   │   ├── transformation/
│   │   │   ├── supercell.ts
│   │   │   ├── slab.ts
│   │   │   └── surface.ts
│   │   ├── analysis/
│   │   │   ├── symmetry.ts
│   │   │   └── validation.ts
│   │   ├── optimization/
│   │   │   ├── mlff-optimize.ts
│   │   │   └── ground-state-search.ts
│   │   └── export/
│   │       ├── export-structure.ts
│   │       └── visualization.ts
│   ├── python/
│   │   ├── crystal_generator.py  # PyXtal wrapper
│   │   ├── symmetry_analyzer.py
│   │   ├── mlff_calculator.py    # MLFF integration
│   │   ├── structure_tools.py
│   │   └── validators.py
│   ├── utils/
│   │   ├── file-io.ts
│   │   ├── python-bridge.ts     # Python subprocess management
│   │   ├── validation.ts
│   │   └── formatting.ts
│   └── config/
│       ├── space-groups.json    # Space group metadata
│       ├── elements.json        # Element properties
│       └── mlff-models.json     # Available MLFF models
├── tests/
│   ├── unit/
│   ├── integration/
│   └── fixtures/
├── docs/
│   ├── API.md
│   ├── examples/
│   └── troubleshooting.md
├── scripts/
│   ├── install-deps.sh
│   └── test-python-env.sh
├── package.json
├── tsconfig.json
├── pyproject.toml              # Python dependencies
└── README.md
```

### 3.3 Core Python Backend Implementation

#### crystal_generator.py
```python
"""
PyXtal wrapper for structure generation with enhanced error handling
"""
from pyxtal import pyxtal
from pyxtal.symmetry import Group
import numpy as np
from typing import Dict, List, Optional, Tuple
import json

class CrystalGenerator:
    """Wrapper for PyXtal with JSON I/O for MCP integration"""
    
    def __init__(self):
        self.crystal = pyxtal()
    
    def generate_crystal(
        self,
        composition: List[str],
        space_group: int,
        num_atoms: Optional[int] = None,
        lattice_params: Optional[Dict] = None,
        volume_factor: float = 1.0,
        min_distance: Optional[Dict[str, float]] = None,
        wyckoff_positions: Optional[List[Dict]] = None,
        seed: Optional[int] = None,
        max_attempts: int = 100
    ) -> Dict:
        """Generate crystal structure with PyXtal"""
        
        if seed is not None:
            np.random.seed(seed)
        
        # Validate space group
        if not (1 <= space_group <= 230):
            raise ValueError(f"Invalid space group: {space_group}")
        
        # Validate composition compatibility
        group = Group(space_group, dim=3)
        if not self._check_composition_compatibility(
            composition, num_atoms, group
        ):
            raise ValueError(
                f"Composition {composition} incompatible with "
                f"space group {space_group}"
            )
        
        # Generate structure
        success = False
        attempt = 0
        
        while attempt < max_attempts and not success:
            try:
                self.crystal.from_random(
                    dim=3,
                    group=space_group,
                    species=composition,
                    numIons=num_atoms if num_atoms else None,
                    factor=volume_factor,
                    lattice=lattice_params,
                    sites=wyckoff_positions
                )
                success = True
            except Exception as e:
                attempt += 1
                if attempt >= max_attempts:
                    raise RuntimeError(
                        f"Failed to generate structure after {max_attempts} "
                        f"attempts. Last error: {str(e)}"
                    )
        
        # Extract structure data
        structure_data = self._extract_structure_data()
        
        # Validate structure
        validation = self._validate_structure(min_distance)
        
        return {
            "structure": structure_data,
            "validation": validation,
            "metadata": {
                "attempts": attempt + 1,
                "seed": seed
            }
        }
    
    def _extract_structure_data(self) -> Dict:
        """Extract structure data in standard format"""
        lattice = self.crystal.lattice
        
        structure_data = {
            "lattice": {
                "a": lattice.a,
                "b": lattice.b,
                "c": lattice.c,
                "alpha": lattice.alpha,
                "beta": lattice.beta,
                "gamma": lattice.gamma,
                "matrix": lattice.matrix.tolist(),
                "volume": lattice.volume
            },
            "atoms": [],
            "space_group": {
                "number": self.crystal.group.number,
                "symbol": self.crystal.group.symbol,
                "hall_symbol": self.crystal.group.hall_symbol,
                "point_group": self.crystal.group.point_group,
                "crystal_system": self.crystal.group.crystal_system
            },
            "metadata": {
                "formula": self.crystal.formula,
                "natoms": len(self.crystal.atom_sites),
                "volume": lattice.volume,
                "density": self._calculate_density()
            }
        }
        
        # Extract atomic positions
        for site in self.crystal.atom_sites:
            for atom in site.coords:
                structure_data["atoms"].append({
                    "element": site.specie,
                    "coords": atom.tolist(),  # Fractional
                    "cartesian": lattice.get_cartesian_coords(atom).tolist(),
                    "wyckoff": site.wp.letter,
                    "multiplicity": site.wp.multiplicity,
                    "site_symmetry": site.wp.site_symm
                })
        
        return structure_data
    
    def _validate_structure(
        self, 
        min_distance: Optional[Dict[str, float]] = None
    ) -> Dict:
        """Validate generated structure"""
        issues = []
        warnings = []
        
        # Check minimum distances
        if min_distance:
            distances = self._calculate_min_distances()
            for pair, dist in distances.items():
                if pair in min_distance and dist < min_distance[pair]:
                    issues.append(
                        f"Distance between {pair} is {dist:.3f} Å, "
                        f"less than minimum {min_distance[pair]:.3f} Å"
                    )
        
        # Check lattice reasonableness
        lattice_warnings = self._check_lattice_parameters()
        warnings.extend(lattice_warnings)
        
        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "warnings": warnings
        }
    
    def _calculate_min_distances(self) -> Dict[str, float]:
        """Calculate minimum interatomic distances"""
        from pymatgen.core import Structure
        
        pmg_structure = self.crystal.to_pymatgen()
        distances = {}
        
        for i in range(len(pmg_structure)):
            for j in range(i+1, len(pmg_structure)):
                site_i = pmg_structure[i]
                site_j = pmg_structure[j]
                dist = site_i.distance(site_j)
                pair = tuple(sorted([site_i.specie.symbol, site_j.specie.symbol]))
                pair_key = f"{pair[0]}-{pair[1]}"
                
                if pair_key not in distances or dist < distances[pair_key]:
                    distances[pair_key] = dist
        
        return distances
    
    def _calculate_density(self) -> float:
        """Calculate crystal density in g/cm³"""
        from pymatgen.core import Element
        
        mass = 0
        for site in self.crystal.atom_sites:
            element = Element(site.specie)
            mass += element.atomic_mass * len(site.coords)
        
        # Convert to g/cm³
        volume_cm3 = self.crystal.lattice.volume * 1e-24  # Å³ to cm³
        avogadro = 6.022e23
        density = (mass / avogadro) / volume_cm3
        
        return density
    
    def _check_composition_compatibility(
        self,
        composition: List[str],
        num_atoms: Optional[int],
        group: Group
    ) -> bool:
        """Check if composition is compatible with space group"""
        # Get available Wyckoff positions
        wyckoffs = group.Wyckoff_positions
        
        # Check if we can distribute atoms across Wyckoff positions
        if num_atoms:
            multiplicities = [wp.multiplicity for wp in wyckoffs]
            # Try to find valid combination
            # This is a simplified check
            return any(
                num_atoms % mult == 0 
                for mult in multiplicities
            )
        
        return True  # If no constraint, assume compatible
    
    def _check_lattice_parameters(self) -> List[str]:
        """Check if lattice parameters are reasonable"""
        warnings = []
        lattice = self.crystal.lattice
        
        # Check for very small or very large parameters
        params = [lattice.a, lattice.b, lattice.c]
        if any(p < 2.0 for p in params):
            warnings.append(
                "Some lattice parameters are very small (< 2 Å)"
            )
        if any(p > 50.0 for p in params):
            warnings.append(
                "Some lattice parameters are very large (> 50 Å)"
            )
        
        # Check for very unusual angles
        angles = [lattice.alpha, lattice.beta, lattice.gamma]
        if any(a < 30 or a > 150 for a in angles):
            warnings.append(
                "Some lattice angles are unusual (< 30° or > 150°)"
            )
        
        return warnings
    
    def export_to_format(self, format_type: str) -> str:
        """Export structure to various formats"""
        if format_type == "cif":
            return self.crystal.to_file()
        elif format_type == "poscar":
            return self.crystal.to_file(fmt="poscar")
        elif format_type == "xyz":
            pmg_structure = self.crystal.to_pymatgen()
            from pymatgen.io.xyz import XYZ
            xyz = XYZ(pmg_structure)
            return str(xyz)
        elif format_type == "json":
            return json.dumps(self._extract_structure_data(), indent=2)
        else:
            raise ValueError(f"Unsupported format: {format_type}")


def main():
    """CLI interface for testing"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python crystal_generator.py <input_json>")
        sys.exit(1)
    
    # Read input JSON
    with open(sys.argv[1], 'r') as f:
        params = json.load(f)
    
    # Generate crystal
    generator = CrystalGenerator()
    result = generator.generate_crystal(**params)
    
    # Output result
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
```

#### mlff_calculator.py
```python
"""
MLFF calculator wrapper for structure optimization and energy calculations
"""
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
from ase.constraints import FixAtoms, UnitCellFilter
import numpy as np
from typing import Dict, List, Optional, Tuple
import json
import time

class MLFFCalculator:
    """Wrapper for MLFF calculators with symmetry constraints"""
    
    AVAILABLE_MODELS = {
        "chgnet": "CHGNet",
        "m3gnet": "M3GNet",
        "mace": "MACE"
    }
    
    def __init__(self, model_name: str):
        self.model_name = model_name
        self.calculator = self._load_calculator(model_name)
    
    def _load_calculator(self, model_name: str):
        """Load the requested MLFF calculator"""
        if model_name == "chgnet":
            from chgnet.model import CHGNet
            model = CHGNet.load()
            return model.get_ase_calculator()
        
        elif model_name == "m3gnet":
            import matgl
            from matgl.ext.ase import PESCalculator
            pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            return PESCalculator(pot)
        
        elif model_name == "mace":
            from mace.calculators import MACECalculator
            return MACECalculator(
                model_path="path/to/mace/model",
                device="cpu"
            )
        
        else:
            raise ValueError(
                f"Unknown model: {model_name}. "
                f"Available: {list(self.AVAILABLE_MODELS.keys())}"
            )
    
    def optimize_structure(
        self,
        structure_dict: Dict,
        optimizer: str = "BFGS",
        fmax: float = 0.01,
        steps: int = 500,
        constrain_symmetry: bool = False,
        fix_lattice: bool = False,
        fix_atoms: Optional[List[int]] = None,
        trajectory_file: Optional[str] = None
    ) -> Dict:
        """Optimize structure using MLFF"""
        
        # Convert to ASE Atoms
        atoms = self._dict_to_atoms(structure_dict)
        
        # Attach calculator
        atoms.calc = self.calculator
        
        # Apply constraints
        if fix_atoms:
            c = FixAtoms(indices=fix_atoms)
            atoms.set_constraint(c)
        
        # Track initial state
        initial_energy = atoms.get_potential_energy()
        initial_forces = atoms.get_forces()
        max_force_initial = np.max(np.linalg.norm(initial_forces, axis=1))
        
        # Set up optimizer
        if constrain_symmetry:
            # Use symmetry-preserving optimization
            # This is a simplified version - full implementation needs
            # spglib integration to preserve space group
            raise NotImplementedError(
                "Symmetry-constrained optimization not yet implemented"
            )
        
        if fix_lattice:
            opt_atoms = atoms
        else:
            # Optimize both cell and positions
            opt_atoms = UnitCellFilter(atoms)
        
        # Choose optimizer
        if optimizer == "BFGS":
            opt = BFGS(opt_atoms, trajectory=trajectory_file)
        elif optimizer == "FIRE":
            opt = FIRE(opt_atoms, trajectory=trajectory_file)
        elif optimizer == "LBFGS":
            opt = LBFGS(opt_atoms, trajectory=trajectory_file)
        else:
            raise ValueError(f"Unknown optimizer: {optimizer}")
        
        # Run optimization
        start_time = time.time()
        opt.run(fmax=fmax, steps=steps)
        optimization_time = time.time() - start_time
        
        # Get final state
        final_energy = atoms.get_potential_energy()
        final_forces = atoms.get_forces()
        max_force_final = np.max(np.linalg.norm(final_forces, axis=1))
        converged = max_force_final < fmax
        
        # Convert back to dict
        optimized_structure = self._atoms_to_dict(atoms)
        
        # Check if symmetry was preserved (if requested)
        preserved_symmetry = False
        final_space_group = None
        if constrain_symmetry:
            from symmetry_analyzer import SymmetryAnalyzer
            analyzer = SymmetryAnalyzer()
            final_sym = analyzer.analyze(optimized_structure)
            original_spg = structure_dict["space_group"]["number"]
            final_space_group = final_sym["space_group"]["number"]
            preserved_symmetry = (original_spg == final_space_group)
        
        return {
            "optimized_structure": optimized_structure,
            "initial_energy": float(initial_energy),
            "final_energy": float(final_energy),
            "energy_change": float(final_energy - initial_energy),
            "max_force_initial": float(max_force_initial),
            "max_force_final": float(max_force_final),
            "n_steps": opt.get_number_of_steps(),
            "converged": converged,
            "preserved_symmetry": preserved_symmetry,
            "final_space_group": final_space_group,
            "timing": {
                "optimization_time": optimization_time
            }
        }
    
    def calculate_energy(
        self,
        structure_dict: Dict,
        calculate_forces: bool = True,
        calculate_stress: bool = True
    ) -> Dict:
        """Calculate energy (and optionally forces/stress)"""
        
        atoms = self._dict_to_atoms(structure_dict)
        atoms.calc = self.calculator
        
        start_time = time.time()
        energy = atoms.get_potential_energy()
        calculation_time = time.time() - start_time
        
        result = {
            "energy": float(energy),
            "energy_per_atom": float(energy / len(atoms)),
            "metadata": {
                "model": self.model_name,
                "calculation_time": calculation_time
            }
        }
        
        if calculate_forces:
            forces = atoms.get_forces()
            result["forces"] = forces.tolist()
        
        if calculate_stress:
            stress = atoms.get_stress()
            result["stress"] = stress.tolist()
        
        return result
    
    def _dict_to_atoms(self, structure_dict: Dict) -> Atoms:
        """Convert structure dictionary to ASE Atoms"""
        lattice = structure_dict["lattice"]
        cell = lattice["matrix"]
        
        positions = []
        symbols = []
        
        for atom in structure_dict["atoms"]:
            symbols.append(atom["element"])
            # Use Cartesian coordinates
            positions.append(atom["cartesian"])
        
        atoms = Atoms(
            symbols=symbols,
            positions=positions,
            cell=cell,
            pbc=True
        )
        
        return atoms
    
    def _atoms_to_dict(self, atoms: Atoms) -> Dict:
        """Convert ASE Atoms to structure dictionary"""
        cell = atoms.get_cell()
        positions_cartesian = atoms.get_positions()
        positions_fractional = atoms.get_scaled_positions()
        
        # Get lattice parameters
        a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
        
        structure = {
            "lattice": {
                "a": float(a),
                "b": float(b),
                "c": float(c),
                "alpha": float(alpha),
                "beta": float(beta),
                "gamma": float(gamma),
                "matrix": cell.tolist(),
                "volume": float(atoms.get_volume())
            },
            "atoms": [],
            "metadata": {
                "natoms": len(atoms)
            }
        }
        
        for i, (symbol, frac_pos, cart_pos) in enumerate(
            zip(atoms.get_chemical_symbols(), 
                positions_fractional, 
                positions_cartesian)
        ):
            structure["atoms"].append({
                "element": symbol,
                "coords": frac_pos.tolist(),
                "cartesian": cart_pos.tolist()
            })
        
        return structure
```

### 3.4 TypeScript Tool Implementation Example

#### generate-crystal.ts
```typescript
import { z } from "zod";
import { PythonShell } from "python-shell";
import { Tool } from "@modelcontextprotocol/sdk/types.js";
import * as path from "path";
import * as fs from "fs/promises";

// Input schema
const GenerateCrystalSchema = z.object({
  composition: z.array(z.string())
    .describe("Chemical composition as array of element symbols, e.g., ['Si', 'O', 'O']"),
  space_group: z.union([z.number().int().min(1).max(230), z.string()])
    .describe("Space group number (1-230) or Hermann-Mauguin symbol"),
  num_atoms: z.number().int().positive().optional()
    .describe("Total number of atoms in the unit cell"),
  lattice_params: z.object({
    a: z.number().positive().optional(),
    b: z.number().positive().optional(),
    c: z.number().positive().optional(),
    alpha: z.number().min(0).max(180).optional(),
    beta: z.number().min(0).max(180).optional(),
    gamma: z.number().min(0).max(180).optional()
  }).optional()
    .describe("Optional lattice parameters in Angstroms and degrees"),
  volume_factor: z.number().positive().default(1.0)
    .describe("Relative volume factor (1.0 = standard)"),
  min_distance: z.record(z.string(), z.number()).optional()
    .describe("Minimum distances between element pairs, e.g., {'Si-Si': 2.0}"),
  wyckoff_positions: z.array(z.object({
    element: z.string(),
    wyckoff: z.string(),
    coords: z.array(z.number()).optional()
  })).optional()
    .describe("Specify exact Wyckoff positions for atoms"),
  seed: z.number().int().optional()
    .describe("Random seed for reproducibility"),
  max_attempts: z.number().int().positive().default(100)
    .describe("Maximum attempts to generate valid structure")
});

type GenerateCrystalInput = z.infer<typeof GenerateCrystalSchema>;

export async function generateCrystal(
  input: GenerateCrystalInput
): Promise<any> {
  // Validate input
  const validated = GenerateCrystalSchema.parse(input);
  
  // Prepare Python script call
  const pythonDir = path.join(__dirname, "..", "python");
  const scriptPath = path.join(pythonDir, "crystal_generator.py");
  
  // Write input to temp file
  const tempInputPath = path.join("/tmp", `crystal_input_${Date.now()}.json`);
  await fs.writeFile(tempInputPath, JSON.stringify(validated, null, 2));
  
  try {
    // Call Python script
    const results = await PythonShell.run(
      scriptPath,
      {
        args: [tempInputPath],
        pythonPath: "python3",  // or path to venv python
        mode: "text"
      }
    );
    
    // Parse output
    const output = JSON.parse(results.join("\n"));
    
    // Clean up temp file
    await fs.unlink(tempInputPath).catch(() => {});
    
    // Generate file exports
    const fileExports = await generateFileExports(output.structure);
    output.files = fileExports;
    
    return {
      content: [
        {
          type: "text",
          text: formatStructureOutput(output)
        }
      ]
    };
    
  } catch (error) {
    // Clean up on error
    await fs.unlink(tempInputPath).catch(() => {});
    
    throw new Error(
      `Crystal generation failed: ${error.message}\n` +
      `Suggestions: ${getSuggestions(error, validated)}`
    );
  }
}

function formatStructureOutput(result: any): string {
  const struct = result.structure;
  const spg = struct.space_group;
  const lattice = struct.lattice;
  
  let output = `## Generated Crystal Structure\n\n`;
  output += `**Formula:** ${struct.metadata.formula}\n`;
  output += `**Space Group:** ${spg.number} (${spg.symbol})\n`;
  output += `**Crystal System:** ${spg.crystal_system}\n`;
  output += `**Lattice Parameters:**\n`;
  output += `  a = ${lattice.a.toFixed(4)} Å\n`;
  output += `  b = ${lattice.b.toFixed(4)} Å\n`;
  output += `  c = ${lattice.c.toFixed(4)} Å\n`;
  output += `  α = ${lattice.alpha.toFixed(2)}°\n`;
  output += `  β = ${lattice.beta.toFixed(2)}°\n`;
  output += `  γ = ${lattice.gamma.toFixed(2)}°\n`;
  output += `**Volume:** ${lattice.volume.toFixed(3)} Ų\n`;
  output += `**Density:** ${struct.metadata.density.toFixed(3)} g/cm³\n`;
  output += `**Number of Atoms:** ${struct.metadata.natoms}\n\n`;
  
  if (result.validation && !result.validation.valid) {
    output += `### ⚠️ Validation Issues\n`;
    for (const issue of result.validation.issues) {
      output += `- ${issue}\n`;
    }
    output += `\n`;
  }
  
  if (result.validation && result.validation.warnings.length > 0) {
    output += `### ⚠️ Warnings\n`;
    for (const warning of result.validation.warnings) {
      output += `- ${warning}\n`;
    }
    output += `\n`;
  }
  
  output += `### Atomic Positions\n\n`;
  output += `| Element | Wyckoff | Fractional Coords | Cartesian Coords (Å) |\n`;
  output += `|---------|---------|-------------------|----------------------|\n`;
  
  for (const atom of struct.atoms) {
    const frac = atom.coords.map((x: number) => x.toFixed(4)).join(", ");
    const cart = atom.cartesian.map((x: number) => x.toFixed(3)).join(", ");
    output += `| ${atom.element} | ${atom.wyckoff} | (${frac}) | (${cart}) |\n`;
  }
  
  return output;
}

function getSuggestions(error: any, input: GenerateCrystalInput): string {
  let suggestions: string[] = [];
  
  if (error.message.includes("max_attempts")) {
    suggestions.push(
      "1. Increase volume_factor (try 1.2-1.5)"
    );
    suggestions.push(
      "2. Relax minimum distance constraints"
    );
    suggestions.push(
      "3. Try a different space group with similar symmetry"
    );
    suggestions.push(
      "4. Reduce the number of atoms if specified"
    );
  }
  
  if (error.message.includes("incompatible")) {
    suggestions.push(
      "1. Check if the composition can fit in this space group"
    );
    suggestions.push(
      "2. Try adjusting the number of atoms to match Wyckoff multiplicities"
    );
    suggestions.push(
      "3. Consult the International Tables for compatible compositions"
    );
  }
  
  return suggestions.join("\n");
}

async function generateFileExports(structure: any): Promise<any> {
  const files: any = {};
  
  // Generate CIF
  files.cif = await exportToCIF(structure);
  
  // Generate POSCAR
  files.poscar = await exportToPOSCAR(structure);
  
  // Generate XYZ
  files.xyz = await exportToXYZ(structure);
  
  // Generate JSON
  files.json = JSON.stringify(structure, null, 2);
  
  return files;
}

// Export functions implementation...
```

### 3.5 Server Configuration

#### server.ts
```typescript
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import {
  CallToolRequestSchema,
  ListToolsRequestSchema
} from "@modelcontextprotocol/sdk/types.js";

// Import all tools
import { generateCrystal } from "./tools/generation/generate-crystal.js";
import { spaceGroupScan } from "./tools/generation/space-group-scan.js";
// ... import other tools

const server = new Server(
  {
    name: "crystal-structure-generator",
    version: "1.0.0"
  },
  {
    capabilities: {
      tools: {}
    }
  }
);

// Register all tools
server.setRequestHandler(ListToolsRequestSchema, async () => {
  return {
    tools: [
      {
        name: "generate_crystal",
        description: "Generate a crystal structure with specified composition and space group",
        inputSchema: {
          // Full Zod schema converted to JSON Schema
          type: "object",
          properties: {
            composition: {
              type: "array",
              items: { type: "string" },
              description: "Chemical composition as array of element symbols"
            },
            space_group: {
              oneOf: [
                { type: "number", minimum: 1, maximum: 230 },
                { type: "string" }
              ],
              description: "Space group number (1-230) or Hermann-Mauguin symbol"
            },
            // ... rest of schema
          },
          required: ["composition", "space_group"]
        },
        annotations: {
          readOnlyHint: false,
          destructiveHint: false,
          idempotentHint: true,
          openWorldHint: true
        }
      },
      // ... other tools
    ]
  };
});

server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;
  
  switch (name) {
    case "generate_crystal":
      return await generateCrystal(args);
    
    case "generate_space_group_scan":
      return await spaceGroupScan(args);
    
    // ... other tools
    
    default:
      throw new Error(`Unknown tool: ${name}`);
  }
});

// Start server
async function main() {
  const transport = new StdioServerTransport();
  await server.connect(transport);
  console.error("Crystal Structure Generator MCP Server running on stdio");
}

main().catch((error) => {
  console.error("Fatal error:", error);
  process.exit(1);
});
```

---

## Phase 4: Testing & Quality Assurance (Week 6-7)

### 4.1 Unit Tests

#### Test Coverage Requirements
- **Tool Input Validation:** 100%
- **Python Bridge:** 100%
- **Structure Generation:** >90%
- **File I/O:** >95%
- **Error Handling:** 100%

#### Key Test Cases

```typescript
// tests/unit/generate-crystal.test.ts
describe("generate_crystal tool", () => {
  test("generates cubic Si structure", async () => {
    const result = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,  // Fd-3m (diamond)
      seed: 42
    });
    
    expect(result.structure.space_group.number).toBe(227);
    expect(result.structure.metadata.formula).toBe("Si2");
    expect(result.validation.valid).toBe(true);
  });
  
  test("handles invalid space group gracefully", async () => {
    await expect(
      generateCrystal({
        composition: ["Si"],
        space_group: 999
      })
    ).rejects.toThrow("Invalid space group");
  });
  
  test("respects minimum distance constraints", async () => {
    const result = await generateCrystal({
      composition: ["Na", "Cl"],
      space_group: 225,  // Fm-3m (rock salt)
      min_distance: { "Na-Cl": 2.5 },
      seed: 42
    });
    
    // Check all Na-Cl distances >= 2.5 Å
    // Implementation...
  });
});
```

### 4.2 Integration Tests

```python
# tests/integration/test_workflows.py
def test_ground_state_search_si():
    """Test full ground state search workflow for Si"""
    params = {
        "composition": ["Si", "Si"],
        "space_groups": [227, 141, 194],  # Diamond, I41/amd, P63/mmc
        "num_structures_per_group": 5,
        "mlff_model": "chgnet",
        "optimization_settings": {
            "optimizer": "BFGS",
            "fmax": 0.01,
            "steps": 500,
            "constrain_symmetry": True
        },
        "parallel": False
    }
    
    result = call_mcp_tool("ground_state_search", params)
    
    assert result["ground_state"]["space_group"] == 227  # Diamond is ground state
    assert result["ground_state"]["energy"] < -5.0  # eV/atom
    assert result["statistics"]["converged_structures"] > 0
```

### 4.3 Performance Benchmarks

Target performance metrics:
- **Basic structure generation:** <1s per structure
- **Space group scan (230 groups):** <10 minutes (parallel)
- **MLFF optimization:** <30s per structure
- **Ground state search (230 groups, 5 attempts each):** <4 hours

---

## Phase 5: Documentation & Deployment (Week 7-8)

### 5.1 Documentation Structure

#### README.md
```markdown
# Crystal Structure Generator MCP Server

Generate accurate crystal structures for DFT and quantum chemistry calculations.

## Features
- All 230 space groups
- MLFF-based optimization
- Supercell, slab, and surface generation
- Ground state searches
- Multiple export formats

## Installation
[Installation instructions]

## Quick Start
[Examples]

## Available Tools
[Tool reference]
```

#### API.md
Full API documentation with:
- All tool schemas
- Input/output examples
- Error codes and solutions
- Best practices

#### Examples Directory
```
examples/
├── basic/
│   ├── generate_simple_crystal.json
│   ├── space_group_scan.json
│   └── supercell.json
├── advanced/
│   ├── ground_state_search.json
│   ├── slab_generation.json
│   └── heterostructure.json
└── workflows/
    ├── dft_preparation.json
    └── mlff_training_data.json
```

### 5.2 Deployment

#### Package Structure
```json
{
  "name": "@your-org/crystal-mcp-server",
  "version": "1.0.0",
  "main": "dist/index.js",
  "bin": {
    "crystal-mcp": "dist/index.js"
  },
  "scripts": {
    "build": "tsc",
    "test": "jest",
    "lint": "eslint src/**/*.ts",
    "install-python": "pip install -r requirements.txt"
  }
}
```

#### Docker Support
```dockerfile
FROM node:20-slim

RUN apt-get update && apt-get install -y python3 python3-pip

WORKDIR /app
COPY package*.json ./
RUN npm install

COPY requirements.txt ./
RUN pip3 install -r requirements.txt

COPY . .
RUN npm run build

CMD ["node", "dist/index.js"]
```

---

## Phase 6: Advanced Features (Week 8-10)

### 6.1 Symmetry-Constrained Optimization

Implement true symmetry-preserving optimization:
- Use spglib to detect symmetry operations
- Apply constraints to maintain space group during optimization
- Project forces onto symmetry-allowed directions
- Constrain lattice parameter changes to preserve crystal system

### 6.2 Machine Learning Integration

#### Structure Fingerprints
- Implement SOAP descriptors
- Crystal graph neural network features
- Use for structure similarity detection

#### Generative Models
- Integrate DiffCSP++ for ML-based structure generation
- VAE-based crystal generation
- Compare ML-generated vs PyXtal structures

### 6.3 High-Throughput Workflows

- Database integration (MongoDB/PostgreSQL)
- Job queue management (Celery/Bull)
- Parallel execution with proper resource management
- Checkpoint/restart capabilities

### 6.4 Visualization Tools

- 3D structure viewers (integrate with NGLViewer, 3Dmol.js)
- Crystal symmetry visualization
- Wyckoff position highlights
- Interactive space group exploration

---

## Appendix A: Critical Implementation Notes

### A.1 Space Group Constraints

**Lattice Parameter Constraints by Crystal System:**

1. **Cubic:** a = b = c, α = β = γ = 90°
2. **Tetragonal:** a = b ≠ c, α = β = γ = 90°
3. **Orthorhombic:** a ≠ b ≠ c, α = β = γ = 90°
4. **Hexagonal:** a = b ≠ c, α = β = 90°, γ = 120°
5. **Trigonal:** a = b = c, α = β = γ ≠ 90° (rhombohedral) OR a = b ≠ c, α = β = 90°, γ = 120° (hexagonal axes)
6. **Monoclinic:** a ≠ b ≠ c, α = γ = 90° ≠ β
7. **Triclinic:** a ≠ b ≠ c, α ≠ β ≠ γ

### A.2 Wyckoff Position Multiplicity Table

[Full table of multiplicities for all 230 space groups - reference International Tables]

### A.3 Common Structure Types

Database of common crystal structures to validate against:
- Rock salt (NaCl): Fm-3m (225)
- Diamond (C, Si): Fd-3m (227)
- Wurtzite (ZnS): P63mc (186)
- Perovskite (CaTiO3): Pm-3m (221)
- Spinel (MgAl2O4): Fd-3m (227)
- Rutile (TiO2): P42/mnm (136)

### A.4 MLFF Model Comparison

| Model | Accuracy | Speed | Training Data | Best For |
|-------|----------|-------|---------------|----------|
| CHGNet | High | Fast | Materials Project | General materials |
| M3GNet | High | Medium | Materials Project | General materials |
| MACE | Very High | Slow | Custom | High accuracy needed |
| SevenNet | High | Fast | Multiple sources | Large systems |

### A.5 Recommended Element Radii (Å)

For minimum distance calculations:
```python
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Ti": 1.60,
    # ... complete table
}

def get_min_distance(elem1, elem2, scale=0.8):
    """Reasonable minimum distance = scale * sum of covalent radii"""
    return scale * (COVALENT_RADII[elem1] + COVALENT_RADII[elem2])
```

---

## Appendix B: Error Codes Reference

```typescript
enum CrystalErrorCode {
  // Input validation
  INVALID_SPACE_GROUP = "E1001",
  INVALID_COMPOSITION = "E1002",
  INVALID_LATTICE_PARAMS = "E1003",
  
  // Generation failures
  GENERATION_TIMEOUT = "E2001",
  MAX_ATTEMPTS_EXCEEDED = "E2002",
  INCOMPATIBLE_COMPOSITION = "E2003",
  INVALID_WYCKOFF_ASSIGNMENT = "E2004",
  
  // Structure validation
  ATOMS_TOO_CLOSE = "E3001",
  INVALID_SYMMETRY = "E3002",
  UNREALISTIC_DENSITY = "E3003",
  
  // MLFF errors
  MODEL_NOT_AVAILABLE = "E4001",
  OPTIMIZATION_FAILED = "E4002",
  ENERGY_CALCULATION_FAILED = "E4003",
  
  // File I/O
  FILE_READ_ERROR = "E5001",
  FILE_WRITE_ERROR = "E5002",
  INVALID_FORMAT = "E5003"
}
```

---

## Appendix C: Testing Checklist

### Before Release
- [ ] All unit tests passing
- [ ] Integration tests for main workflows
- [ ] Performance benchmarks meet targets
- [ ] Documentation complete and accurate
- [ ] Example files work correctly
- [ ] Error messages are helpful
- [ ] Python environment setup automated
- [ ] Docker image builds successfully
- [ ] Security audit passed
- [ ] License compliance verified

### Space Group Coverage Testing
- [ ] Test cubic space groups (195-230)
- [ ] Test hexagonal space groups (168-194)
- [ ] Test trigonal space groups (143-167)
- [ ] Test tetragonal space groups (75-142)
- [ ] Test orthorhombic space groups (16-74)
- [ ] Test monoclinic space groups (3-15)
- [ ] Test triclinic space groups (1-2)

---

## Summary & Next Steps

This comprehensive plan provides:

1. **Complete tool architecture** covering all requirements
2. **Detailed implementation guidance** for TypeScript/Python integration
3. **MLFF workflow integration** for ground state searches
4. **Robust error handling** with actionable suggestions
5. **Production-ready structure** with testing and deployment

**Recommended Development Timeline:**
- **Weeks 1-2:** Research & Planning
- **Weeks 3-6:** Core Implementation
- **Weeks 6-7:** Testing
- **Weeks 7-8:** Documentation
- **Weeks 8-10:** Advanced Features

**Critical Success Factors:**
1. PyXtal mastery for all 230 space groups
2. Robust symmetry preservation during optimization
3. Efficient MLFF integration
4. Clear error messages and user guidance
5. Comprehensive testing across space groups

**Ready to proceed?** Start with Phase 1: Load MCP documentation, study PyXtal examples, and test space group generation for common materials (Si, NaCl, TiO2) to validate the approach.
