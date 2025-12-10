# Crystal Structure MCP Server - System Architecture

## High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        MCP Client (Claude)                       │
│                  (Requesting crystal structures)                 │
└────────────────────────────┬────────────────────────────────────┘
                             │ MCP Protocol
                             │ (stdio / HTTP)
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                   Crystal MCP Server (TypeScript)                │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │                    Tool Registry                            │ │
│  │  • generate_crystal        • make_supercell                │ │
│  │  • space_group_scan        • generate_slab                 │ │
│  │  • analyze_symmetry        • optimize_structure_mlff       │ │
│  │  • validate_structure      • ground_state_search          │ │
│  └────────────────────────────────────────────────────────────┘ │
│                             │                                    │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │                  Python Bridge Layer                        │ │
│  │              (python-shell / child_process)                 │ │
│  │  • JSON serialization     • Error handling                 │ │
│  │  • Process management     • Result parsing                 │ │
│  └────────────────────────────────────────────────────────────┘ │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│               Python Backend (Scientific Computing)              │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │                   Crystal Generator                         │ │
│  │  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐  │ │
│  │  │  PyXtal  │  │ Pymatgen │  │  Spglib  │  │   ASE    │  │ │
│  │  │          │  │          │  │          │  │          │  │ │
│  │  │ Generate │  │ Convert  │  │ Symmetry │  │Structure │  │ │
│  │  │ Crystals │  │ Formats  │  │ Analysis │  │  I/O     │  │ │
│  │  └──────────┘  └──────────┘  └──────────┘  └──────────┘  │ │
│  └────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │              MLFF Optimization Engine                       │ │
│  │  ┌──────────┐  ┌──────────┐  ┌──────────┐                 │ │
│  │  │  CHGNet  │  │ M3GNet   │  │   MACE   │                 │ │
│  │  │          │  │          │  │          │                 │ │
│  │  │ Energy & │  │ Energy & │  │ Energy & │                 │ │
│  │  │  Forces  │  │  Forces  │  │  Forces  │                 │ │
│  │  └──────────┘  └──────────┘  └──────────┘                 │ │
│  │                                                              │ │
│  │  ┌────────────────────────────────────────────────┐        │ │
│  │  │          ASE Optimizers                        │        │ │
│  │  │  BFGS | FIRE | LBFGS | GPMin                  │        │ │
│  │  └────────────────────────────────────────────────┘        │ │
│  └────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │                Structure Transformations                    │ │
│  │  • Supercell builder     • Slab generator                  │ │
│  │  • Defect creator        • Strain application              │ │
│  │  • Alloy generator       • Surface builder                 │ │
│  └────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                        Output Formats                            │
│  CIF │ POSCAR │ XYZ │ JSON │ PDB │ PWscf │ CASTEP │ LAMMPS      │
└─────────────────────────────────────────────────────────────────┘
```

## Data Flow: Generate Crystal Structure

```
1. Client Request
   └──> generate_crystal({composition: ["Si", "Si"], space_group: 227})
         │
         ▼
2. TypeScript Tool (Validation)
   └──> Validate input schema with Zod
         │ Valid? ✓
         ▼
3. Python Bridge
   └──> Serialize to JSON
   └──> Call crystal_generator.py
         │
         ▼
4. PyXtal Generation
   └──> Initialize pyxtal()
   └──> from_random(dim=3, group=227, species=["Si", "Si"])
   └──> Check constraints (distances, volume, etc.)
   └──> Retry if needed (up to max_attempts)
         │
         ▼
5. Structure Extraction
   └──> Extract lattice parameters
   └──> Extract atomic positions (fractional + Cartesian)
   └──> Extract Wyckoff positions
   └──> Extract space group info
         │
         ▼
6. Validation
   └──> Check minimum distances
   └──> Verify lattice parameters reasonable
   └──> Check symmetry consistency
         │
         ▼
7. Format Conversion
   └──> Generate CIF
   └──> Generate POSCAR
   └──> Generate XYZ
   └──> Create JSON representation
         │
         ▼
8. Return to Client
   └──> Structure data + files + validation results
```

## Data Flow: Ground State Search

```
1. Client Request
   └──> ground_state_search({
          composition: ["Si", "Si"],
          space_groups: [227, 141, 194],
          mlff_model: "chgnet"
        })
         │
         ▼
2. Parallel Generation
   ┌────────────┬────────────┬────────────┐
   │ SG 227     │ SG 141     │ SG 194     │
   │ (Diamond)  │ (I41/amd)  │ (P63/mmc)  │
   └─────┬──────┴─────┬──────┴─────┬──────┘
         │            │            │
         ▼            ▼            ▼
3. Generate Multiple Attempts per Space Group
   ├──> Attempt 1 (different random seed)
   ├──> Attempt 2
   ├──> Attempt 3
   └──> Attempt n
         │
         ▼
4. MLFF Optimization (per structure)
   └──> Load CHGNet calculator
   └──> Attach to ASE Atoms
   └──> Run BFGS optimization (fmax=0.01)
   └──> Calculate final energy
         │
         ▼
5. Energy Ranking
   └──> Sort all structures by energy
   └──> Identify ground state (lowest energy)
   └──> Calculate energy above ground state for each
         │
         ▼
6. Statistical Analysis
   └──> Group by crystal system
   └──> Calculate average energies
   └──> Identify metastable structures
         │
         ▼
7. Return Results
   └──> Ground state structure
   └──> Complete energy ranking
   └──> Statistics and analysis
```

## Component Interaction Matrix

```
┌─────────────────┬──────┬──────┬──────┬──────┬──────┬──────┐
│  Component      │PyXtal│Pymat.│Spglib│ ASE  │CHGNet│Export│
├─────────────────┼──────┼──────┼──────┼──────┼──────┼──────┤
│generate_crystal │  ●   │  ○   │  ○   │  ○   │      │  ●   │
│analyze_symmetry │      │  ○   │  ●   │      │      │      │
│validate_struct  │      │  ●   │  ●   │      │      │      │
│make_supercell   │      │  ●   │  ○   │  ●   │      │  ●   │
│generate_slab    │      │  ●   │  ○   │  ●   │      │  ●   │
│optimize_mlff    │      │  ○   │  ○   │  ●   │  ●   │  ●   │
│ground_state_sea │  ●   │  ○   │  ○   │  ●   │  ●   │  ●   │
└─────────────────┴──────┴──────┴──────┴──────┴──────┴──────┘

Legend: ● = Primary dependency, ○ = Secondary dependency
```

## Module Dependencies

```
TypeScript Layer
├── @modelcontextprotocol/sdk ─┐
├── zod                        │─> Core MCP functionality
├── python-shell               │
└── file I/O utilities          ─┘

Python Backend
├── Core Crystallography
│   ├── pyxtal ───────────────┐
│   ├── pymatgen              │─> Structure generation & analysis
│   └── spglib                ─┘
│
├── Machine Learning
│   ├── chgnet ───────────────┐
│   ├── matgl (M3GNet)        │─> MLFF calculations
│   └── mace-torch             ─┘
│
├── Structure Manipulation
│   ├── ase ──────────────────┐
│   ├── numpy                 │─> Optimization & transforms
│   └── scipy                 ─┘
│
└── Utilities
    ├── json
    ├── pathlib
    └── typing
```

## Tool Categorization & Dependencies

```
┌──────────────────────────────────────────────────────┐
│              Generation Tools                         │
│  ┌────────────────┐  ┌────────────────┐             │
│  │generate_crystal│  │space_group_scan│             │
│  │   [PyXtal]     │  │   [PyXtal]     │             │
│  └────────────────┘  └────────────────┘             │
│  ┌────────────────┐                                  │
│  │molecular_crystal│                                  │
│  │[PyXtal+OpenBabel]                                 │
│  └────────────────┘                                  │
└──────────────────────────────────────────────────────┘
                    │
                    ▼
┌──────────────────────────────────────────────────────┐
│           Transformation Tools                        │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐    │
│  │make_super- │  │generate_   │  │create_     │    │
│  │cell        │  │slab        │  │defect      │    │
│  │[Pymatgen]  │  │[Pymatgen]  │  │[Pymatgen]  │    │
│  └────────────┘  └────────────┘  └────────────┘    │
└──────────────────────────────────────────────────────┘
                    │
                    ▼
┌──────────────────────────────────────────────────────┐
│             Analysis Tools                            │
│  ┌────────────┐  ┌────────────┐                     │
│  │analyze_    │  │validate_   │                     │
│  │symmetry    │  │structure   │                     │
│  │[Spglib]    │  │[Spglib+PM] │                     │
│  └────────────┘  └────────────┘                     │
└──────────────────────────────────────────────────────┘
                    │
                    ▼
┌──────────────────────────────────────────────────────┐
│          Optimization Tools                           │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐    │
│  │optimize_   │  │calculate_  │  │ground_state│    │
│  │mlff        │  │energy_mlff │  │_search     │    │
│  │[ASE+MLFF]  │  │[ASE+MLFF]  │  │[All above] │    │
│  └────────────┘  └────────────┘  └────────────┘    │
└──────────────────────────────────────────────────────┘
                    │
                    ▼
┌──────────────────────────────────────────────────────┐
│             Export Tools                              │
│  ┌────────────┐  ┌────────────┐                     │
│  │export_     │  │generate_   │                     │
│  │structure   │  │visualization│                     │
│  │[Multiple]  │  │[3Dmol.js]  │                     │
│  └────────────┘  └────────────┘                     │
└──────────────────────────────────────────────────────┘
```

## File System Layout

```
crystal-mcp-server/
│
├── src/                          # TypeScript source
│   ├── index.ts                  # Entry point
│   ├── server.ts                 # MCP server setup
│   │
│   ├── types/                    # Type definitions
│   │   ├── crystal.ts            # Structure types
│   │   ├── tools.ts              # Tool schemas
│   │   └── errors.ts             # Error types
│   │
│   ├── tools/                    # Tool implementations
│   │   ├── generation/
│   │   │   ├── generate-crystal.ts
│   │   │   ├── space-group-scan.ts
│   │   │   └── molecular-crystal.ts
│   │   │
│   │   ├── transformation/
│   │   │   ├── supercell.ts
│   │   │   ├── slab.ts
│   │   │   └── defect.ts
│   │   │
│   │   ├── analysis/
│   │   │   ├── symmetry.ts
│   │   │   └── validation.ts
│   │   │
│   │   ├── optimization/
│   │   │   ├── mlff-optimize.ts
│   │   │   └── ground-state-search.ts
│   │   │
│   │   └── export/
│   │       ├── export-structure.ts
│   │       └── visualization.ts
│   │
│   ├── python/                   # Python backend
│   │   ├── crystal_generator.py  # PyXtal wrapper
│   │   ├── symmetry_analyzer.py  # Spglib wrapper
│   │   ├── mlff_calculator.py    # MLFF integration
│   │   ├── structure_tools.py    # Transformations
│   │   └── validators.py         # Validation logic
│   │
│   ├── utils/                    # Utilities
│   │   ├── python-bridge.ts      # Python subprocess
│   │   ├── file-io.ts            # File operations
│   │   ├── validation.ts         # Input validation
│   │   └── formatting.ts         # Output formatting
│   │
│   └── config/                   # Configuration
│       ├── space-groups.json     # Space group metadata
│       ├── elements.json         # Element properties
│       └── mlff-models.json      # MLFF model configs
│
├── tests/                        # Test suites
│   ├── unit/                     # Unit tests
│   ├── integration/              # Integration tests
│   └── fixtures/                 # Test data
│
├── dist/                         # Compiled JS (generated)
├── node_modules/                 # Dependencies (generated)
│
├── docs/                         # Documentation
│   ├── API.md                    # API reference
│   ├── examples/                 # Usage examples
│   └── troubleshooting.md        # Common issues
│
├── scripts/                      # Utility scripts
│   ├── install-deps.sh           # Setup script
│   └── test-python-env.sh        # Verify Python env
│
├── package.json                  # Node.js dependencies
├── tsconfig.json                 # TypeScript config
├── requirements.txt              # Python dependencies
├── pyproject.toml               # Python project config
├── .gitignore
└── README.md
```

## Communication Flow

```
MCP Client (Claude)
        │
        │ 1. List Tools Request
        │────────────────────────────>
        │                              MCP Server (TS)
        │ 2. Available Tools          │
        │<────────────────────────────│
        │                              │
        │ 3. Call Tool                 │
        │    (generate_crystal)        │
        │────────────────────────────> │
        │                              │ 4. Validate Input
        │                              │    (Zod schema)
        │                              │
        │                              │ 5. Prepare Python Call
        │                              │    (Serialize JSON)
        │                              │
        │                              ▼
        │                         Python Backend
        │                              │
        │                              │ 6. Generate Crystal
        │                              │    (PyXtal)
        │                              │
        │                              │ 7. Validate Result
        │                              │    (Distances, etc.)
        │                              │
        │                              │ 8. Format Output
        │                              │    (CIF, POSCAR, etc.)
        │                              │
        │                              │ 9. Return JSON
        │                              ▼
        │                         MCP Server (TS)
        │                              │
        │                              │ 10. Parse Result
        │                              │     (Error handling)
        │                              │
        │ 11. Tool Response            │
        │     (Structure + Files)      │
        │<────────────────────────────│
        │
        ▼
   Client receives structure
```

## Error Handling Flow

```
                  ┌──────────────┐
                  │  User Input  │
                  └──────┬───────┘
                         │
                         ▼
            ┌────────────────────────┐
            │  Input Validation (TS) │
            └────┬──────────────┬────┘
                 │              │
           Valid │              │ Invalid
                 │              └─────────> Return error immediately
                 ▼                          (E1001-E1003)
        ┌─────────────────┐
        │ Python Call     │
        └────┬────────┬───┘
             │        │
       OK    │        │ Process Error
             │        └────────────────────> Log & format error
             ▼                               (E5001-E5003)
    ┌──────────────────┐
    │ PyXtal Generation│
    └────┬────────┬────┘
         │        │
   OK    │        │ Generation Failed
         │        └────────────────────────> Retry with adjusted params
         │                                   or return E2001-E2004
         ▼
┌──────────────────┐
│  Validation      │
└────┬────────┬────┘
     │        │
 OK  │        │ Validation Failed
     │        └────────────────────────────> Return warnings
     │                                       (E3001-E3003)
     ▼
┌──────────────────┐
│ Format & Return  │
└──────────────────┘
```

## Performance Considerations

```
Operation                      Time Complexity    Space Complexity
───────────────────────────────────────────────────────────────────
generate_crystal               O(n²)              O(n)
  - Wyckoff assignment         O(n)               O(1)
  - Distance checks            O(n²)              O(n²)
  - Symmetry operations        O(n × s)           O(s)
    where s = # symmetry ops

space_group_scan (serial)      O(230 × n²)        O(230 × n)
space_group_scan (parallel)    O((230/p) × n²)    O(230 × n)
  where p = # parallel workers

make_supercell                 O(n × m)           O(n × m)
  where m = multiplication

generate_slab                  O(n × l)           O(n × l)
  where l = # layers

optimize_mlff                  O(steps × n³)      O(n)
  - Energy evaluation          O(n³)              O(n)
  - Force calculation          O(n³)              O(n)

ground_state_search            O(sg × att × steps × n³)
  where:
    sg = # space groups
    att = attempts per group
    steps = optimization steps
    n = atoms per structure

Memory Requirements:
─────────────────────
Single structure:     ~10 KB
MLFF model:          ~100 MB - 1 GB
Optimization cache:   ~100 MB per structure
Parallel workers:     ~(model_size + cache) × n_workers
```

## Scalability Strategy

```
                    Small Scale             Medium Scale           Large Scale
                    (Single user)           (Team)                 (HPC Cluster)
─────────────────────────────────────────────────────────────────────────────────
Architecture        Single process          Multi-process          Distributed

Concurrency         Sequential              Thread pool            Job queue
                                           (4-8 workers)          (100+ workers)

MLFF Backend        CPU                     CPU + GPU              GPU cluster

Structure DB        File system             SQLite                 PostgreSQL

Caching             None                    In-memory              Redis

Rate Limiting       None                    Per-user               Per-organization

Monitoring          Logs                    Metrics                Full observability
─────────────────────────────────────────────────────────────────────────────────
```

This architecture supports scaling from personal research to high-throughput materials discovery campaigns.
