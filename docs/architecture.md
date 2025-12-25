# Crystal Structure Generator MCP Server - Architecture

## System Overview

The Crystal Structure Generator is an MCP (Model Context Protocol) server that provides comprehensive crystal structure generation, manipulation, and analysis capabilities. It bridges TypeScript/JavaScript front-end with Python crystallography libraries.

```mermaid
graph TB
    subgraph "Client Layer"
        LLM[LLM Client<br/>Claude/GPT]
        USER[User Interface]
    end

    subgraph "MCP Protocol Layer"
        MCP[MCP Server<br/>index.ts]
        TOOLS[Tool Registry<br/>Zod Schemas]
    end

    subgraph "Python Bridge"
        DISPATCHER[Python Dispatcher<br/>crystal_generator.py]
        ADAPTER[Schema Adapters]
    end

    subgraph "Generator Modules"
        BULK[Bulk Generators<br/>Space Groups]
        2D[2D Materials<br/>Graphene, TMDs]
        SURFACE[Surface Tools<br/>Slabs, Adsorbates]
        DEFECT[Defect Generation<br/>Vacancies, Dopants]
        WORKFLOW[Workflow Editing<br/>Iterative Ops]
    end

    subgraph "Core Libraries"
        PYXTAL[PyXtal<br/>Random Structures]
        PYMATGEN[PyMatGen<br/>Structure Manipulation]
        ASE[ASE<br/>Molecular Structures]
        MLFF[MLFF Models<br/>CHGNet, M3GNet]
    end

    LLM -->|JSON-RPC 2.0| MCP
    USER -->|JSON-RPC 2.0| MCP
    MCP --> TOOLS
    TOOLS -->|Spawn Python| DISPATCHER
    DISPATCHER --> ADAPTER
    ADAPTER --> BULK
    ADAPTER --> 2D
    ADAPTER --> SURFACE
    ADAPTER --> DEFECT
    ADAPTER --> WORKFLOW

    BULK --> PYXTAL
    BULK --> PYMATGEN
    2D --> PYMATGEN
    SURFACE --> PYMATGEN
    DEFECT --> PYMATGEN
    WORKFLOW --> PYMATGEN
    WORKFLOW --> ASE
    BULK --> MLFF

    style MCP fill:#e1f5ff
    style DISPATCHER fill:#fff4e1
    style PYXTAL fill:#e8f5e9
    style PYMATGEN fill:#e8f5e9
    style MLFF fill:#fce4ec
```

## System Architecture Layers

### 1. Client Layer
- **LLM Clients**: Claude, GPT, or other AI assistants
- **Direct API Clients**: Custom applications using MCP protocol
- **Communication**: JSON-RPC 2.0 over stdio/HTTP

### 2. MCP Protocol Layer (TypeScript)
**Location**: `src/index.ts`, `src/types/`

Responsibilities:
- Protocol compliance (JSON-RPC 2.0, MCP specification)
- Tool discovery and schema validation
- Request routing
- Error handling and formatting

```
┌─────────────────────────────────────────┐
│      MCP Server Entry Point             │
│          (index.ts)                     │
├─────────────────────────────────────────┤
│ • Initialize MCP server                 │
│ • Register tool handlers                │
│ • Validate requests (Zod)               │
│ • Format responses                      │
└──────────────┬──────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────┐
│      Tool Schema Definitions            │
│         (src/types/tools.ts)            │
├─────────────────────────────────────────┤
│ • GenerateCrystalSchema                 │
│ • ComprehensiveGenerateSchema           │
│ • MakeSupercellSchema                   │
│ • GenerateSlabSchema                    │
│ • ExportStructureSchema                 │
│ • ... (50+ tool schemas)                │
└─────────────────────────────────────────┘
```

### 3. Python Bridge Layer
**Location**: `src/python/crystal_generator.py`

```mermaid
sequenceDiagram
    participant MCP as MCP Server (TS)
    participant Spawn as Python Spawner
    participant CG as crystal_generator.py
    participant Gen as Generator Module

    MCP->>Spawn: spawn(crystal_generator.py, args)
    Spawn->>CG: sys.argv (JSON params)
    CG->>CG: Parse & validate params
    CG->>Gen: Call appropriate generator
    Gen->>Gen: Generate/transform structure
    Gen-->>CG: Return result dict
    CG->>CG: JSON.stringify(result)
    CG-->>Spawn: stdout: JSON result
    Spawn-->>MCP: Return JSON
    MCP-->>MCP: Format MCP response
```

**Responsibilities**:
- Parse command-line arguments
- Route to appropriate generator module
- Centralize error handling
- Standardize response format

### 4. Generator Modules (Python)

```
src/python/generators/
├── bulk/
│   ├── spacegroups.py          # PyXtal space group generation
│   ├── prototypes.py           # Common structure prototypes
│   └── clathrates.py           # Clathrate frameworks
│
├── two_d/
│   ├── base.py                 # Generic 2D material generator
│   ├── graphene.py             # Graphene, nanotubes, ribbons
│   ├── transition_metal.py     # TMDs (MoS2, WS2, etc.)
│   ├── janus.py                # Janus 2D structures
│   └── magnetic_2d.py          # Magnetic 2D materials
│
├── surface/
│   └── (uses structure_tools.py)
│
├── defect/
│   └── point_defects.py        # Vacancies, interstitials, substitutions
│
├── workflow/
│   ├── iterative_editing.py   # StructureEditor class
│   ├── parametric_scans.py    # Z-scans, parameter sweeps
│   └── schema_adapter.py      # Workflow ↔ Core schema conversion
│
├── quality_control/
│   └── symmetry.py            # Symmetry analysis & validation
│
└── output_formats/
    └── converters.py          # VASP, QE, XYZ, CIF exporters
```

## Module Interaction Diagram

```mermaid
graph LR
    subgraph Input
        A[User Request<br/>Space Group + Composition]
    end

    subgraph Type Safety
        B[Zod Schema<br/>Validation]
    end

    subgraph Routing
        C[crystal_generator.py<br/>Dispatcher]
    end

    subgraph Generation
        D[spacegroups.py<br/>PyXtal Integration]
        E[Structure Dict<br/>Generation]
    end

    subgraph Validation
        F[Symmetry Analyzer<br/>spglib Check]
        G[Min Distance<br/>Validator]
    end

    subgraph Output
        H[Export Formats<br/>CIF/POSCAR/XYZ]
        I[JSON Response<br/>MCP Format]
    end

    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    F --> G
    G -->|Valid| H
    G -->|Invalid| C
    H --> I

    style D fill:#fff4e1
    style F fill:#e1f5ff
    style I fill:#e8f5e9
```

## Critical Data Structures

### Core Structure Schema

```python
{
    "atoms": [
        {
            "element": "Si",           # Element symbol
            "coords": [0.0, 0.0, 0.0], # Fractional coordinates
            "cartesian": [x, y, z]     # Cartesian (optional)
        }
    ],
    "lattice": {
        "matrix": [[ax,ay,az], [bx,by,bz], [cx,cy,cz]],  # 3x3 matrix
        "a": 5.43, "b": 5.43, "c": 5.43,                 # Lengths (Å)
        "alpha": 90.0, "beta": 90.0, "gamma": 90.0,       # Angles (°)
        "volume": 160.1                                   # Volume (ų)
    },
    "space_group": {
        "number": 227,        # International Tables number (1-230)
        "symbol": "Fd-3m",    # Hermann-Mauguin symbol
        "hall_symbol": "F 4d 2 3 -1d",
        "point_group": "m-3m",
        "crystal_system": "cubic"
    },
    "metadata": {
        "formula": "Si2",
        "natoms": 2,
        "source": "pyxtal",
        "composition_warnings": []  # Wyckoff adjustments
    }
}
```

### Workflow Editing Schema

```python
{
    "species": ["Si", "Si"],                    # Element list
    "positions": [[0,0,0], [2.715,2.715,2.715]], # Cartesian coords
    "cell": [[5.43,0,0], [0,5.43,0], [0,0,5.43]], # Cell matrix
    "n_atoms": 2
}
```

**Schema Adapter**: `generators/workflow/schema_adapter.py`
- `workflow_to_core()`: Converts Cartesian → Fractional, adds lattice params
- `core_to_workflow()`: Converts Fractional → Cartesian, extracts cell matrix

## Key Design Decisions

### 1. TypeScript + Python Hybrid Architecture

**Rationale**:
- TypeScript: Protocol layer benefits from strong typing (Zod schemas)
- Python: Rich ecosystem for crystallography (PyXtal, pymatgen, ASE)
- Bridge: Clean separation of concerns, type-safe boundaries

### 2. Schema Validation at Multiple Layers

```
Layer 1: Zod (TypeScript)  →  Type safety, basic validation
Layer 2: Python Dispatcher  →  Parameter ranges, element validity
Layer 3: Generator Modules  →  Scientific constraints (min_distance, Wyckoff)
Layer 4: Post-Generation    →  Symmetry verification, structure validation
```

### 3. Fractional vs. Cartesian Coordinates

**Core System**: Fractional coordinates
- Reason: Crystallographic convention, symmetry operations are simpler
- Conversion: Always done at boundaries (import/export, workflow editing)

**Workflow Editing**: Cartesian coordinates
- Reason: Intuitive for direct manipulation (translation, rotation)
- Trade-off: Schema adapter adds complexity but enables better UX

### 4. Error Handling Philosophy

**No try/except in generators** (user preference):
```python
# Instead of:
try:
    crystal.from_random(...)
except Exception as e:
    return error_dict

# Use:
if not crystal.valid:
    return {"success": False, "error": {...}}
```

**Explicit validation checks**:
- Pre-flight parameter validation
- Post-generation structure verification
- Detailed error codes and messages

### 5. Composition Transparency

Always report:
- `requested_composition`: What user asked for
- `actual_composition`: What was generated
- `composition_warnings`: Wyckoff multiplicity adjustments

Example:
```python
{
    "requested_composition": [1, 2],  # User wanted Si1O2
    "actual_composition": [2, 4],     # Generated Si2O4
    "composition_adjusted": true,
    "composition_warnings": ["Si: requested 1, adjusted to 2 (Wyckoff 8a requires multiples of 2)"]
}
```

## Performance Characteristics

### Structure Generation Times (Typical)

```
Bulk (simple cubic, 10 atoms):        0.1-0.5s
Bulk (complex, low symmetry, 50 atoms): 2-10s
2D Material (graphene sheet):         0.2-1s
Slab from bulk (3 layers):            0.5-2s
Supercell (2×2×2):                    0.1-0.5s
MLFF optimization (50 atoms):         5-30s
```

### Bottlenecks

1. **PyXtal Random Generation**: Space groups with low Wyckoff multiplicities
   - Mitigation: Retry loop with crystal.valid check

2. **MLFF Optimization**: Energy/force calculations
   - Mitigation: Optional, user-controlled

3. **Symmetry Analysis (spglib)**: Large structures (>1000 atoms)
   - Mitigation: Adjustable symprec tolerance

## Extension Points

### Adding New Generators

1. **Create generator module** in appropriate category
2. **Register in `__init__.py`** for module exports
3. **Add Zod schema** in `src/types/tools.ts`
4. **Add dispatcher route** in `crystal_generator.py`
5. **Add tests** in `tests/testsuit/`

### Adding New Export Formats

1. **Implement converter** in `generators/output_formats/converters.py`
2. **Add format to schema** in `ExportStructureSchema`
3. **Update dispatcher** in `crystal_generator.py` export handler
4. **Test roundtrip**: Structure → Export → Import → Validate

### Adding MLFF Models

1. **Add import** in `crystal_generator.py` (conditional)
2. **Implement calculator** wrapper if needed
3. **Add to `mlff_model` enum** in schema
4. **Test on representative structures**

## Security Considerations

### Input Validation

- **Element symbols**: Validated against periodic table (Z=1-118)
- **Space groups**: Range check (1-230)
- **Numeric parameters**: Range checks, NaN detection
- **File paths**: Not exposed (all operations in-memory)

### Resource Limits

- **Structure size**: Implicit limit via generation timeouts
- **MLFF steps**: Capped at reasonable values (default: 500)
- **Supercell size**: No hard limit, but warns for large expansions

### Sandboxing

- Python subprocess isolation
- No filesystem access in generators
- Temporary files cleaned up automatically

## Dependencies and Versions

```
Critical Dependencies (Crystallography):
├── pyxtal >= 1.0.0         # Structure generation
├── pymatgen >= 2024.1.1    # Structure manipulation
├── spglib >= 2.5.0         # Symmetry analysis
└── ase >= 3.22.0           # Molecular structures

MLFF Models (Optional):
├── chgnet >= 0.3.0         # M3GNet-based
├── matgl >= 1.0.0          # Graph neural network
└── mace-torch >= 0.3.0     # MACE models

Numerical/Scientific:
├── numpy >= 1.20.0
└── scipy >= 1.10.0         # For rotations (workflow editing)

File I/O:
└── openbabel-wheel >= 3.1.0  # Molecule parsing
```

## Testing Strategy

```
tests/
├── testsuit/
│   ├── test_mcp_comprehensive.py    # 56 MCP protocol tests
│   ├── test_mcp_e2e.py             # 18 end-to-end tests
│   ├── test_scientific_accuracy.py  # Scientific validation
│   └── test_all_operations.py       # Coverage for all generators
│
├── test_schema_adapter.py           # Workflow schema conversion
└── test_adsorbate_placement.py      # Surface normal placement

Coverage Target: >80% for generator modules
```

### Test Categories

1. **Protocol Compliance**: MCP JSON-RPC 2.0 spec
2. **Scientific Correctness**: Lattice preservation, composition accuracy
3. **Error Handling**: Invalid inputs, edge cases
4. **Performance**: Generation times, large structures
5. **Integration**: Multi-step workflows

## Monitoring and Observability

### Error Codes

```python
Error Code Categories:
├── INVALID_INPUT          # Bad parameters
├── GENERATION_FAILED      # PyXtal/structure creation failed
├── MIN_DISTANCE_VIOLATION # Atoms too close
├── COMPOSITION_MISMATCH   # Stoichiometry impossible
└── EXPORT_ERROR          # Format conversion failed
```

### Logging Strategy

- **TypeScript Layer**: Tool calls, validation failures
- **Python Layer**: Generation attempts, warnings, errors
- **Format**: Structured JSON for parsing

## Future Architecture Improvements

1. **Caching Layer**:
   - Cache generated structures by (space_group, composition, seed)
   - TTL-based eviction
   - Estimated 50-80% hit rate for common materials

2. **Async Python Bridge**:
   - Use async subprocess spawning
   - Handle multiple concurrent requests
   - Requires protocol-level changes

3. **Database Backend**:
   - Store generated structures
   - Enable search by composition, space group, properties
   - Materials Project / OQMD integration

4. **GPU Acceleration**:
   - MLFF calculations on GPU
   - Requires CUDA-enabled MLFF models

5. **Distributed Generation**:
   - Split space group scan across workers
   - Redis/RabbitMQ for job queue
   - Aggregate results

---

**Document Version**: 2.0
**Last Updated**: 2025-12-25
**Maintainer**: Crystal Structure Generator Team

---

## COMPREHENSIVE ARCHITECTURAL DIAGRAMS

This section provides detailed visual representations of all system components, their interactions, and data flows for scientific accuracy and developer clarity.

---

## Complete Component Dependency Graph

This diagram shows all 20 generator categories and their dependencies on core libraries:

```mermaid
graph TB
    subgraph "MCP Layer"
        MCP[MCP Server<br/>index.ts]
    end

    subgraph "Router Layer"
        CR[Comprehensive Router<br/>comprehensive_structures.py]
        CG[Crystal Generator<br/>crystal_generator.py]
    end

    subgraph "Generator Categories (20)"
        BULK[Bulk<br/>spacegroups, prototypes,<br/>perovskites, clathrates]
        TWO_D[2D Materials<br/>graphene, TMDs, xenes,<br/>Janus, ribbons]
        SURFACE[Surface<br/>slabs, adsorbates,<br/>interfaces, reconstructions]
        DEFECT[Defect<br/>point, extended,<br/>grain boundaries]
        TWIST[Twist<br/>bilayers, moiré,<br/>stacking, intercalation]
        MOLECULE[Molecule<br/>crystals, frameworks,<br/>MOFs, COFs]
        NANO[Nanotube<br/>CNTs, BNNTs,<br/>nanowires]
        ELEC[Electronic<br/>semiconductors,<br/>topological, Weyl]
        MAG[Magnetic<br/>ferro, antiferro,<br/>spin spirals]
        BATTERY[Battery<br/>cathodes, anodes,<br/>solid electrolytes]
        CATALYST[Catalyst<br/>metal oxides,<br/>single atoms, clusters]
        THERMO[Thermoelectric<br/>Half-Heusler, Zintl,<br/>skutterudites]
        PHOTONIC[Photonic<br/>crystals, metamaterials,<br/>inverse opals]
        QUANTUM[Quantum<br/>dots, wells,<br/>superlattices]
        HIGH_P[High Pressure<br/>diamond anvil,<br/>phase transitions]
        EXTERNAL[External Fields<br/>electric, magnetic,<br/>strain responses]
        QC[Quality Control<br/>symmetry, validation,<br/>optimization]
        OUTPUT[Output Formats<br/>VASP, QE, LAMMPS,<br/>CIF, XYZ]
        WORKFLOW[Workflow<br/>iterative editing,<br/>parametric scans]
        TOOLS[Structure Tools<br/>utilities, transformations,<br/>analysis]
    end

    subgraph "Core Libraries"
        PYXTAL[PyXtal<br/>Random crystal<br/>generation]
        PYMATGEN[PyMatGen<br/>Structure manipulation<br/>and analysis]
        ASE_LIB[ASE<br/>Molecular structures<br/>and dynamics]
        SPGLIB[Spglib<br/>Symmetry analysis<br/>and detection]
        MLFF[MLFF Models<br/>CHGNet, M3GNet,<br/>MACE]
    end

    %% MCP to Router
    MCP --> CR
    MCP --> CG

    %% Router to Generators
    CR --> BULK
    CR --> TWO_D
    CR --> SURFACE
    CR --> DEFECT
    CR --> TWIST
    CR --> MOLECULE
    CR --> NANO
    CR --> ELEC
    CR --> MAG
    CR --> BATTERY
    CR --> CATALYST
    CR --> THERMO
    CR --> PHOTONIC
    CR --> QUANTUM
    CR --> HIGH_P
    CR --> EXTERNAL
    CR --> QC
    CR --> OUTPUT
    CR --> WORKFLOW

    %% All use tools
    BULK --> TOOLS
    TWO_D --> TOOLS
    SURFACE --> TOOLS
    DEFECT --> TOOLS
    TWIST --> TOOLS

    %% Generators to Core Libraries
    BULK --> PYXTAL
    BULK --> PYMATGEN
    BULK --> SPGLIB
    
    TWO_D --> PYMATGEN
    TWO_D --> ASE_LIB
    
    SURFACE --> PYMATGEN
    SURFACE --> TOOLS
    
    DEFECT --> PYMATGEN
    
    TWIST --> PYMATGEN
    TWIST --> ASE_LIB
    
    MOLECULE --> PYMATGEN
    MOLECULE --> ASE_LIB
    
    NANO --> ASE_LIB
    
    ELEC --> PYMATGEN
    
    MAG --> PYMATGEN
    
    BATTERY --> PYMATGEN
    
    CATALYST --> PYMATGEN
    
    THERMO --> PYMATGEN
    
    PHOTONIC --> PYMATGEN
    
    QUANTUM --> PYMATGEN
    
    HIGH_P --> PYMATGEN
    
    EXTERNAL --> PYMATGEN
    
    QC --> PYMATGEN
    QC --> SPGLIB
    QC --> MLFF
    
    OUTPUT --> PYMATGEN
    
    WORKFLOW --> PYMATGEN
    WORKFLOW --> ASE_LIB

    style MCP fill:#e1f5ff
    style CR fill:#fff4e1
    style PYXTAL fill:#e8f5e9
    style PYMATGEN fill:#e8f5e9
    style MLFF fill:#fce4ec
    style TOOLS fill:#f3e5f5
```

---

## Detailed Request Lifecycle with Timing

This sequence diagram shows a complete request from Claude Desktop through all layers with typical timing information:

```mermaid
sequenceDiagram
    autonumber
    actor User
    participant Claude as Claude Desktop
    participant MCP as MCP Server (TS)<br/>index.ts
    participant Zod as Zod Validator
    participant Spawn as Python Spawner
    participant CG as crystal_generator.py
    participant Router as comprehensive_structures.py
    participant Gen as Generator Module
    participant PyXtal as PyXtal Library
    participant PM as PyMatGen Library
    
    Note over User,PM: Typical Timeline: 100-500ms total
    
    User->>Claude: "Generate Si diamond structure"
    Note over Claude: ~10ms: Parse intent
    
    Claude->>MCP: tools/call: comprehensive_generate
    Note over Claude,MCP: {operation: "generate_from_spacegroup",<br/>spacegroup: 227, elements: ["Si"], composition: [8]}
    
    MCP->>Zod: Validate request schema
    Note over Zod: ~1ms: Type checking
    Zod-->>MCP: ✓ Valid
    
    MCP->>Spawn: spawn("python3", ["crystal_generator.py", "--json", "..."])
    Note over Spawn: ~20ms: Process startup
    
    Spawn->>CG: Execute with sys.argv
    
    CG->>CG: Parse command-line JSON
    Note over CG: ~2ms: JSON parsing
    
    CG->>Router: Route to operation handler
    Note over Router: ~1ms: Function lookup
    
    Router->>Gen: generate_from_spacegroup(227, ["Si"], [8])
    
    Gen->>Gen: Validate inputs
    Note over Gen: - Check space group range (1-230)<br/>- Validate elements<br/>- Check composition ratios
    
    Gen->>Gen: Adjust composition for Wyckoff
    Note over Gen: Si: 2 atoms → 8 atoms (4a multiplicity)
    
    Gen->>PyXtal: crystal.from_random(...)
    Note over PyXtal: ~50-200ms: Random structure<br/>generation with retry logic
    
    PyXtal-->>Gen: Structure object
    
    Gen->>PM: Symmetry verification
    Note over PM: ~10ms: spglib analysis
    PM-->>Gen: Space group 227 confirmed
    
    Gen->>Gen: Build result dictionary
    Note over Gen: - Convert to standard schema<br/>- Add metadata<br/>- Include warnings
    
    Gen-->>Router: Return success dict
    Router-->>CG: Return result
    
    CG->>CG: JSON.stringify(result)
    Note over CG: ~5ms: Serialization
    
    CG-->>Spawn: stdout: JSON string
    Spawn-->>MCP: Capture stdout
    
    MCP->>MCP: Format MCP response
    Note over MCP: Wrap in MCP protocol:<br/>{content: [{type: "text", text: "..."}]}
    
    MCP-->>Claude: Return tool result
    Note over Claude: ~10ms: Process response
    
    Claude-->>User: "Generated Si8 diamond structure (Fd-3m)"
    
    Note over User,PM: Total Time: ~100-250ms typical<br/>Complex structures: up to 2-3 seconds
```

---

## Schema Transformation Flow

This diagram illustrates how structures are converted between different schemas as they flow through the system:

```mermaid
graph TB
    subgraph "Input Schemas"
        USER_JSON["User JSON Input<br/>{spacegroup: 227,<br/>composition: ['Si', 8]}"]
        MCP_REQUEST["MCP Request<br/>{method: 'tools/call',<br/>params: {...}}"]
    end

    subgraph "Validation Layer"
        ZOD_SCHEMA["Zod Schema<br/>Type: number, string[], etc.<br/>Constraints: min, max, enum"]
    end

    subgraph "Internal Processing"
        PYTHON_ARGS["Python Arguments<br/>Parsed from sys.argv<br/>JSON deserialized"]
        
        PYMATGEN_INPUT["PyMatGen Input<br/>Structure object<br/>Fractional coords"]
    end

    subgraph "Generation"
        PYXTAL_STRUCT["PyXtal Structure<br/>Internal representation<br/>Symmetry operations"]
    end

    subgraph "Standard Schema"
        CORE_SCHEMA["Core Structure Schema<br/>{atoms: [{element, coords}],<br/>lattice: {matrix, a, b, c, ...},<br/>metadata: {...}}"]
    end

    subgraph "Workflow Schema (Optional)"
        WORKFLOW_SCHEMA["Workflow Schema<br/>{species: [],<br/>positions: [] (Cartesian),<br/>cell: 3x3 matrix}"]
    end

    subgraph "Output Formats"
        CIF["CIF File<br/>Crystallographic<br/>Information File"]
        POSCAR["VASP POSCAR<br/>Direct or Cartesian"]
        XYZ["XYZ Format<br/>Element + Cartesian"]
        QE["Quantum Espresso<br/>pw.x input"]
    end

    subgraph "Response"
        MCP_RESPONSE["MCP Response<br/>{content: [{type: 'text',<br/>text: JSON.stringify(...)}]}"]
    end

    USER_JSON --> MCP_REQUEST
    MCP_REQUEST --> ZOD_SCHEMA
    ZOD_SCHEMA --> PYTHON_ARGS
    PYTHON_ARGS --> PYMATGEN_INPUT
    PYMATGEN_INPUT --> PYXTAL_STRUCT
    PYXTAL_STRUCT --> CORE_SCHEMA
    CORE_SCHEMA --> WORKFLOW_SCHEMA
    CORE_SCHEMA --> CIF
    CORE_SCHEMA --> POSCAR
    CORE_SCHEMA --> XYZ
    CORE_SCHEMA --> QE
    CIF --> MCP_RESPONSE
    POSCAR --> MCP_RESPONSE
    XYZ --> MCP_RESPONSE
    QE --> MCP_RESPONSE
    CORE_SCHEMA --> MCP_RESPONSE

    style USER_JSON fill:#e3f2fd
    style CORE_SCHEMA fill:#c8e6c9
    style WORKFLOW_SCHEMA fill:#fff9c4
    style MCP_RESPONSE fill:#f3e5f5
```

---

## Complete Generator Category Map (All 228 Operations)

This hierarchical diagram shows all 20 categories and representative operations:

```mermaid
graph TB
    ROOT[Crystal MCP Server<br/>228 Operations]

    ROOT --> CAT1[Bulk Structures<br/>35 operations]
    ROOT --> CAT2[2D Materials<br/>25 operations]
    ROOT --> CAT3[Surface<br/>18 operations]
    ROOT --> CAT4[Defect<br/>15 operations]
    ROOT --> CAT5[Twist<br/>12 operations]
    ROOT --> CAT6[Molecule<br/>10 operations]
    ROOT --> CAT7[Nanotube<br/>8 operations]
    ROOT --> CAT8[Electronic<br/>15 operations]
    ROOT --> CAT9[Magnetic<br/>12 operations]
    ROOT --> CAT10[Battery<br/>10 operations]
    ROOT --> CAT11[Catalyst<br/>10 operations]
    ROOT --> CAT12[Thermoelectric<br/>8 operations]
    ROOT --> CAT13[Photonic<br/>8 operations]
    ROOT --> CAT14[Quantum<br/>10 operations]
    ROOT --> CAT15[High Pressure<br/>8 operations]
    ROOT --> CAT16[External Fields<br/>8 operations]
    ROOT --> CAT17[Quality Control<br/>12 operations]
    ROOT --> CAT18[Output Formats<br/>15 operations]
    ROOT --> CAT19[Workflow<br/>6 operations]
    ROOT --> CAT20[Adsorption<br/>3 operations]

    CAT1 --> BULK_OPS["• generate_from_spacegroup<br/>• generate_from_prototype<br/>• generate_perovskite<br/>• generate_clathrate<br/>• generate_zeolite<br/>• make_supercell<br/>• apply_strain<br/>• generate_polytype<br/>...27 more"]

    CAT2 --> TWO_D_OPS["• generate_2d_material<br/>• generate_graphene<br/>• generate_tmd<br/>• generate_xene<br/>• generate_janus<br/>• generate_nanoribbon<br/>• generate_2d_magnetic<br/>...18 more"]

    CAT3 --> SURFACE_OPS["• generate_slab<br/>• add_adsorbate<br/>• create_heterostructure<br/>• generate_interface<br/>• generate_stepped_surface<br/>• generate_reconstruction<br/>...12 more"]

    CAT4 --> DEFECT_OPS["• generate_vacancy<br/>• generate_interstitial<br/>• generate_substitution<br/>• create_defect<br/>• create_alloy<br/>• generate_grain_boundary<br/>• generate_dislocation<br/>...8 more"]

    CAT5 --> TWIST_OPS["• generate_twisted_bilayer<br/>• generate_moire_superlattice<br/>• generate_stacking_variant<br/>• generate_intercalated<br/>• generate_ferroelectric_twist<br/>...7 more"]

    CAT6 --> MOLECULE_OPS["• build_molecule<br/>• generate_molecular_crystal<br/>• generate_mof<br/>• generate_cof<br/>• generate_molecular_cluster<br/>...5 more"]

    CAT7 --> NANO_OPS["• generate_nanotube<br/>• generate_nanowire<br/>• generate_nanoribbon<br/>• generate_fullerene<br/>...4 more"]

    CAT8 --> ELEC_OPS["• generate_semiconductor<br/>• generate_topological<br/>• generate_weyl_semimetal<br/>• generate_quantum_spin_hall<br/>...11 more"]

    CAT17 --> QC_OPS["• analyze_symmetry<br/>• validate_structure<br/>• optimize_structure_mlff<br/>• calculate_energy_mlff<br/>• ground_state_search<br/>...7 more"]

    CAT18 --> OUTPUT_OPS["• export_vasp<br/>• export_qe<br/>• export_lammps<br/>• export_cif<br/>• export_xyz<br/>• generate_visualization<br/>...9 more"]

    style ROOT fill:#e1f5ff
    style CAT1 fill:#c8e6c9
    style CAT2 fill:#c8e6c9
    style CAT17 fill:#fff9c4
    style CAT18 fill:#f3e5f5
```

---

## Error Handling and Propagation Flow

This diagram shows how errors are detected, propagated, and formatted through all layers:

```mermaid
graph TB
    subgraph "Error Sources"
        USER_ERR[User Input Error<br/>Invalid parameters]
        VALIDATION_ERR[Validation Error<br/>Schema mismatch]
        GENERATION_ERR[Generation Error<br/>PyXtal failure]
        SCIENTIFIC_ERR[Scientific Error<br/>Min distance violation]
        SYSTEM_ERR[System Error<br/>Library exception]
    end

    subgraph "Detection Layer"
        ZOD_DETECT[Zod Validator<br/>Type/constraint errors]
        PARAM_DETECT[Parameter Validator<br/>Range/value checks]
        GENERATOR_DETECT[Generator Validation<br/>Scientific constraints]
        EXCEPTION_DETECT[Exception Handler<br/>Catch-all]
    end

    subgraph "Error Formatting"
        ERROR_DICT["Error Dictionary<br/>{success: false,<br/>error: {code, message, details}}"]
    end

    subgraph "Response Path"
        JSON_RESPONSE[JSON Response<br/>Serialized error]
        MCP_RESPONSE[MCP Response<br/>Wrapped in content block]
        USER_MESSAGE[User-Friendly Message<br/>LLM interprets error]
    end

    USER_ERR --> ZOD_DETECT
    VALIDATION_ERR --> ZOD_DETECT
    VALIDATION_ERR --> PARAM_DETECT
    
    GENERATION_ERR --> GENERATOR_DETECT
    SCIENTIFIC_ERR --> GENERATOR_DETECT
    SYSTEM_ERR --> EXCEPTION_DETECT

    ZOD_DETECT --> ERROR_DICT
    PARAM_DETECT --> ERROR_DICT
    GENERATOR_DETECT --> ERROR_DICT
    EXCEPTION_DETECT --> ERROR_DICT

    ERROR_DICT --> JSON_RESPONSE
    JSON_RESPONSE --> MCP_RESPONSE
    MCP_RESPONSE --> USER_MESSAGE

    style USER_ERR fill:#ffebee
    style VALIDATION_ERR fill:#ffebee
    style GENERATION_ERR fill:#ffebee
    style ERROR_DICT fill:#fff3e0
    style USER_MESSAGE fill:#e8f5e9
```

### Error Code Reference

| Error Code | Layer | Description | User Action |
|------------|-------|-------------|-------------|
| `INVALID_PARAMETER` | Validation | Parameter out of range or wrong type | Check API documentation for valid values |
| `INVALID_SPACE_GROUP` | Validation | Space group not in range 1-230 | Use valid space group number |
| `INVALID_COMPOSITION` | Validation | Element symbols invalid or composition impossible | Check element symbols and stoichiometry |
| `GENERATION_FAILED` | Generation | PyXtal could not generate structure | Try different parameters or seed |
| `MIN_DISTANCE_VIOLATION` | Scientific | Atoms too close (< min_distance) | Increase min_distance or adjust composition |
| `WYCKOFF_MISMATCH` | Scientific | Composition incompatible with Wyckoff positions | Adjust atom counts or use different space group |
| `SYMMETRY_VERIFICATION_FAILED` | Validation | Generated structure doesn't match requested space group | Report as bug - should not occur |
| `EXPORT_ERROR` | Export | Structure cannot be converted to target format | Check structure validity first |
| `SCHEMA_CONVERSION_ERROR` | Workflow | Cannot convert between schemas | Check that structure has required fields |

---

## Memory and Performance Characteristics

### Typical Resource Usage

```mermaid
graph LR
    subgraph "Operation Complexity"
        SIMPLE[Simple Generation<br/>Single unit cell<br/>< 50 atoms]
        MEDIUM[Medium Generation<br/>Supercell 2x2x2<br/>50-500 atoms]
        COMPLEX[Complex Generation<br/>Large supercell<br/>500-5000 atoms]
        MLFF_OP[MLFF Optimization<br/>Energy minimization<br/>Any size]
    end

    subgraph "Resource Requirements"
        R_SIMPLE[Memory: ~50 MB<br/>Time: 50-200 ms<br/>CPU: 1 core]
        R_MEDIUM[Memory: ~100 MB<br/>Time: 200-500 ms<br/>CPU: 1 core]
        R_COMPLEX[Memory: ~500 MB<br/>Time: 1-5 seconds<br/>CPU: 1 core]
        R_MLFF[Memory: ~1-2 GB<br/>Time: 5-30 seconds<br/>CPU/GPU: intensive]
    end

    SIMPLE --> R_SIMPLE
    MEDIUM --> R_MEDIUM
    COMPLEX --> R_COMPLEX
    MLFF_OP --> R_MLFF

    style SIMPLE fill:#c8e6c9
    style MEDIUM fill:#fff9c4
    style COMPLEX fill:#ffccbc
    style MLFF_OP fill:#f8bbd0
```

### Performance Optimization Strategies

1. **Caching (Future)**
   - Cache generated structures by hash(space_group, composition, seed)
   - Estimated hit rate: 50-80% for common materials
   - TTL: 1 hour for development, 24 hours for production

2. **Parallel Generation**
   - Space group scans can be parallelized (each space group independent)
   - Python multiprocessing for CPU-bound tasks
   - Requires protocol-level changes (async support)

3. **GPU Acceleration**
   - MLFF models support GPU (CUDA)
   - 5-10x speedup for optimization
   - Requires: `chgnet[cuda]`, CUDA toolkit

4. **Memory Management**
   - Clear PyXtal/PyMatGen objects after generation
   - Limit supercell size (max 10,000 atoms recommended)
   - Stream large structure lists instead of holding in memory

---

## Security Considerations

### Input Validation
```mermaid
graph TB
    INPUT[User Input] --> L1[Layer 1: Zod Type Check]
    L1 --> L2[Layer 2: Range Validation]
    L2 --> L3[Layer 3: Scientific Constraints]
    L3 --> L4[Layer 4: Resource Limits]
    L4 --> SAFE[Safe Execution]

    L1 -->|Type Error| REJECT1[Reject: Invalid Type]
    L2 -->|Out of Range| REJECT2[Reject: Invalid Value]
    L3 -->|Unphysical| REJECT3[Reject: Scientific Error]
    L4 -->|Too Large| REJECT4[Reject: Resource Limit]

    style SAFE fill:#c8e6c9
    style REJECT1 fill:#ffebee
    style REJECT2 fill:#ffebee
    style REJECT3 fill:#ffebee
    style REJECT4 fill:#ffebee
```

### Resource Limits
- **Max atoms per structure**: 10,000 (configurable)
- **Max supercell expansion**: 10x10x10 = 1000x volume
- **Max generation attempts**: 100 (PyXtal retries)
- **Python subprocess timeout**: 60 seconds (configurable)

### Sandboxing
- Python subprocess runs in isolated environment
- No file system access except `/tmp` for IPC
- No network access required
- Environment variables sanitized

---

**Document Version**: 3.0 (Enhanced with Comprehensive Diagrams)  
**Last Updated**: 2025-12-25  
**Maintainer**: Crystal Structure Generator Team

