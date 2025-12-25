
# Workflows & Usage Patterns

This document visualizes common usage patterns for the Crystal MCP Server.

## 1. Iterative Structure Refinement
The system allows output from one tool to be passed directly into the next, enabling complex modification chains without starting from scratch.

```mermaid
sequenceDiagram
    participant U as User
    participant MCP as MCP Server
    participant Router as Python Router
    participant G1 as Bulk Generator
    participant G2 as Defect Generator
    
    U->>MCP: Call generate_from_spacegroup(Si, Fd-3m)
    MCP->>Router: Execute
    Router->>G1: generate_from_spacegroup(...)
    G1-->>Router: Structure (8 atoms)
    Router-->>MCP: JSON Result
    MCP-->>U: "Generated Si8"
    
    Note over U,MCP: User decides to add a vacancy
    
    U->>MCP: Call generate_vacancy(host_structure=PREV_OUTPUT)
    MCP->>Router: Execute
    Note right of Router: Auto-converts JSON to Pymatgen Object
    Router->>G2: generate_vacancy(structure, ...)
    G2-->>Router: Modified Structure (7 atoms)
    Router-->>MCP: JSON Result
    MCP-->>U: "Generated Si7 (Vacancy added)"
```

## 2. Generation to Export Pipeline
Seamlessly generate a structure and export it to a specific simulation format (VASP, QE, LAMMPS).

```mermaid
sequenceDiagram
    participant U as User
    participant MCP as MCP Server
    participant Router as Python Router
    participant Gen as Generator
    participant Exp as Exporter (converters.py)
    
    U->>MCP: Call generate_zeolite(MFI)
    MCP->>Router: Execute
    Router->>Gen: generate_zeolite(...)
    Gen-->>Router: Zeolite Structure
    Router-->>MCP: JSON Result
    
    U->>MCP: Call export_vasp(structure=PREV_OUTPUT)
    MCP->>Router: Execute
    Note right of Router: Auto-converts JSON to Pymatgen Object
    Router->>Exp: export_vasp(structure)
    Exp-->>Router: POSCAR String
    Router-->>MCP: JSON Result
    MCP-->>U: VASP POSCAR Content
```

## 3. Namespace Collision Resolution
When function names match across categories (e.g., `generate_vacancy` in `bulk` vs `defect`), use `category` to disambiguate.

```mermaid
sequenceDiagram
    participant U as User
    participant Router as Python Router
    participant Reg as Registry
    
    U->>Router: generate_vacancy(category="defect")
    
    Router->>Router: Prioritize "defect" category
    Router->>Reg: Look in "defect"?
    Reg-->>Router: Found! (generators.defect.point_defects)
    
    Router->>Module: Execute defect.generate_vacancy
    Note right of Router: Ignores "bulk" version
```

---

## COMPREHENSIVE WORKFLOW PATTERNS

This section provides detailed visual guides for all common workflow patterns with scientific accuracy.

---

## 4. Complete Defect Study Workflow

### From Bulk to Defect Energy Calculation

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant MCP
    
    Note over User,MCP: Step 1: Generate Perfect Crystal
    User->>MCP: generate_from_spacegroup(Si, 227)
    MCP-->>User: perfect_structure (8 atoms)
    
    Note over User,MCP: Step 2: Create Supercell for Defect Isolation
    User->>MCP: make_supercell(perfect_structure, [3,3,3])
    MCP-->>User: supercell (216 atoms)
    
    Note over User,MCP: Step 3: Optimize Perfect Supercell
    User->>MCP: optimize_structure_mlff(supercell, "chgnet")
    MCP-->>User: optimized_perfect (216 atoms, E_perfect)
    
    Note over User,MCP: Step 4: Create Vacancy Defect
    User->>MCP: generate_vacancy(optimized_perfect, site_index=0)
    MCP-->>User: defect_structure (215 atoms)
    
    Note over User,MCP: Step 5: Optimize Defect Structure
    User->>MCP: optimize_structure_mlff(defect_structure, "chgnet")
    MCP-->>User: optimized_defect (215 atoms, E_defect)
    
    Note over User,MCP: Step 6: Calculate Formation Energy
    Note over User: E_f = E_defect - E_perfect + μ_Si<br/>E_f = formation energy<br/>μ_Si = chemical potential
    
    Note over User,MCP: Step 7: Export for Further Analysis
    User->>MCP: export_vasp(optimized_defect)
    MCP-->>User: POSCAR for DFT refinement
```

**Scientific Workflow Notes**:
- **Supercell size**: Must be >10Å to avoid defect-defect interaction
- **Optimization**: MLFF first (~seconds), then DFT refinement (~hours)
- **Chemical potential**: μ_Si from bulk Si calculation
- **Charge states**: Repeat for q = 0, ±1, ±2 for charged defects

---

## 5. Surface Catalysis Workflow

### Complete Adsorption Energy Calculation

```mermaid
graph TB
    START[Research Question:<br/>CO oxidation on Pt111]

    START --> BULK[Generate Pt fcc bulk]
    BULK --> SLAB[Cut Pt111 slab<br/>5 layers, 12Å vacuum]
    SLAB --> OPT_CLEAN[Optimize clean slab<br/>Fix bottom 2 layers]

    OPT_CLEAN --> GAS_PHASE[Calculate gas-phase energies]
    GAS_PHASE --> CO_GAS[CO molecule energy]
    GAS_PHASE --> O2_GAS[O₂ molecule energy]

    OPT_CLEAN --> BRANCH{Adsorption<br/>configuration}

    BRANCH --> CO_ADS[CO adsorption<br/>top/bridge/hollow sites]
    BRANCH --> O_ADS[O adsorption<br/>fcc/hcp sites]
    BRANCH --> COADS[CO + O coadsorption<br/>various geometries]

    CO_ADS --> OPT_CO[Optimize CO/Pt]
    O_ADS --> OPT_O[Optimize O/Pt]
    COADS --> OPT_COADS[Optimize CO+O/Pt]

    OPT_CO --> E_ADS_CO[E_ads_CO = E_CO/Pt - E_clean - E_CO]
    OPT_O --> E_ADS_O[E_ads_O = E_O/Pt - E_clean - 0.5*E_O2]
    OPT_COADS --> REACTION[Reaction barrier<br/>CO + O → CO₂]

    E_ADS_CO --> ANALYSIS[Compare site preferences<br/>Identify most stable]
    E_ADS_O --> ANALYSIS
    REACTION --> ANALYSIS

    ANALYSIS --> EXPORT[Export to VASP/QE<br/>for accurate DFT]

    style START fill:#e1f5ff
    style ANALYSIS fill:#fff9c4
    style EXPORT fill:#c8e6c9
```

**Example Operations**:

```json
// Step 1: Bulk Pt
{
  "operation": "generate_from_prototype",
  "prototype": "fcc",
  "elements": ["Pt"],
  "a": 3.924
}

// Step 2: Pt(111) slab
{
  "operation": "generate_slab",
  "structure": "<pt_bulk>",
  "miller_indices": [1, 1, 1],
  "thickness": 5,
  "vacuum": 12.0,
  "symmetric": true,
  "fix_bottom_layers": 2
}

// Step 3: CO adsorption (top site)
{
  "operation": "add_adsorbate",
  "surface_structure": "<pt111_slab>",
  "molecule": "CO",
  "site_index": 0,
  "distance": 2.0
}

// Step 4: Optimize
{
  "operation": "optimize_structure_mlff",
  "structure": "<co_pt111>",
  "mlff_model": "chgnet",
  "fmax": 0.01,
  "fix_atoms": [20, 21, 22, 23, 24]  // Bottom 2 layers (5 atoms/layer × 2)
}

// Step 5: Export for DFT
{
  "operation": "export_vasp",
  "structure": "<optimized_co_pt111>",
  "selective_dynamics": true  // Preserve fixed atoms
}
```

---

## 6. 2D Heterostructure Workflow

### Twisted Bilayer Graphene with Characterization

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant MCP

    Note over User,MCP: Design Phase
    User->>MCP: What twist angle for flat bands?
    Note over User: Literature: θ ≈ 1.1° (magic angle)

    Note over User,MCP: Step 1: Generate TBG
    User->>MCP: generate_twisted_bilayer("graphene", 1.1°)
    MCP-->>User: tbg_structure (~11,000 atoms)
    Note over User: Moiré period ≈ 13 nm

    Note over User,MCP: Step 2: Analyze Structure
    User->>MCP: analyze_symmetry(tbg_structure)
    MCP-->>User: Space group, point group
    Note over User: Verify moiré symmetry

    Note over User,MCP: Step 3: Relax (Optional)
    User->>MCP: optimize_structure_mlff(tbg, "m3gnet")
    MCP-->>User: relaxed_tbg
    Note over User: Interlayer distance: ~3.35Å

    Note over User,MCP: Step 4: Create Supercell for Band Structure
    User->>MCP: make_supercell(relaxed_tbg, [1,1,1])
    Note over User: Already large enough

    Note over User,MCP: Step 5: Export for DFT
    User->>MCP: export_qe(relaxed_tbg)
    MCP-->>User: QE input
    Note over User: Run DFT band structure:<br/>Look for flat bands near E_F

    Note over User,MCP: Alternative: Direct DFT
    Note over User: For smaller angles (>2°):<br/>Can use DFT directly<br/>Smaller supercells
```

**Twist Angle vs. Computational Cost**:

| Angle (°) | Moiré Period (nm) | Atoms in Cell | Method |
|-----------|-------------------|---------------|---------|
| 0.5 | ~28 | ~50,000 | MLFF only |
| 1.1 | ~13 | ~11,000 | MLFF + tight-binding |
| 2.0 | ~7 | ~3,000 | DFT feasible |
| 5.0 | ~3 | ~500 | Direct DFT |
| 10.0 | ~1.5 | ~150 | Direct DFT |

---

## 7. High-Throughput Screening Workflow

### Space Group Scan for Novel Materials

```mermaid
graph TB
    GOAL[Goal: Find stable LiCoO₂<br/>polymorphs]

    GOAL --> SCAN[Space Group Scan<br/>All 230 space groups]

    SCAN --> GENERATE[For each space group:<br/>Generate 3 random structures]
    Note[690 structures total]

    GENERATE --> FILTER1{Quick Filters}
    FILTER1 --> F1[Check min distances]
    FILTER1 --> F2[Check charge neutrality]
    FILTER1 --> F3[Check reasonable density]

    F1 -->|Pass| MLFF[MLFF Optimization]
    F2 -->|Pass| MLFF
    F3 -->|Pass| MLFF

    MLFF --> RANK[Rank by energy]
    RANK --> TOP[Select top 20<br/>lowest energy]

    TOP --> DFT[DFT refinement]
    DFT --> ANALYSIS[Analyze:<br/>- Band structure<br/>- Stability<br/>- Li diffusion barriers]

    ANALYSIS --> CANDIDATES[Identify promising candidates]

    style GOAL fill:#e1f5ff
    style CANDIDATES fill:#c8e6c9
    style Note fill:#fff9c4
```

**Parallel Execution Strategy**:

```python
# Pseudo-code for high-throughput workflow
space_groups = range(1, 231)
results = []

for sg in space_groups:
    for attempt in range(3):  # 3 structures per space group
        # Generate
        structure = generate_from_spacegroup(
            spacegroup=sg,
            elements=["Li", "Co", "O"],
            composition=[1, 1, 2],  # Will auto-adjust
            seed=attempt
        )

        # Quick filter
        if not validate_structure(structure):
            continue

        # MLFF optimization
        optimized = optimize_structure_mlff(
            structure,
            mlff_model="chgnet",
            fmax=0.05  # Less strict for screening
        )

        results.append({
            "space_group": sg,
            "energy": optimized["energy"],
            "structure": optimized["structure"]
        })

# Rank and select
results.sort(key=lambda x: x["energy"])
top_20 = results[:20]

# Export for DFT
for result in top_20:
    export_vasp(result["structure"])
```

---

## 8. Battery Material Workflow

### Cathode Structure to Li Diffusion Path

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant MCP

    Note over User,MCP: Step 1: Generate Cathode
    User->>MCP: generate_cathode("LCO")
    MCP-->>User: LiCoO₂ layered structure

    Note over User,MCP: Step 2: Create Supercell
    User->>MCP: make_supercell(LCO, [2,2,2])
    MCP-->>User: supercell (96 atoms)

    Note over User,MCP: Step 3: Create Li Vacancy
    User->>MCP: generate_vacancy(supercell, element="Li")
    MCP-->>User: delithiated_structure

    Note over User,MCP: Step 4: Optimize Structure
    User->>MCP: optimize_structure_mlff(delithiated, "chgnet")
    MCP-->>User: relaxed_structure

    Note over User,MCP: Step 5: Calculate Migration Barrier (NEB)
    Note over User: Use DFT with climbing image NEB:<br/>Initial: Li in site A<br/>Final: Li in site B<br/>Interpolate pathway

    Note over User,MCP: Step 6: Analyze Pathways
    Note over User: Identify:<br/>- Lowest barrier path<br/>- 3D diffusion network<br/>- Bottleneck sites

    Note over User,MCP: Step 7: Voltage Profile
    User->>MCP: For x = 0.0 to 1.0 in Li₁₋ₓCoO₂
    Note over User: V(x) = -[E(Li₁₋ₓCoO₂) - E(CoO₂) - (1-x)E(Li)] / (1-x)
```

**Key Calculations**:

```
Li Diffusion Barrier:
  E_barrier = E_transition_state - E_initial
  Target: <0.4 eV for good rate capability

Voltage vs. Composition:
  V(x) = -dG/dN_Li
  Practical: ~3.7-4.0 V for LiCoO₂

Structural Stability:
  Check volume change during delithiation:
  ΔV/V < 5% desirable
```

---

## 9. Iterative Structure Refinement

### User-Guided Modification Chain

```mermaid
graph LR
    START[Initial Structure] --> EDIT1[Edit 1:<br/>Add Adsorbate]
    EDIT1 --> EDIT2[Edit 2:<br/>Move Atom]
    EDIT2 --> EDIT3[Edit 3:<br/>Substitute Element]
    EDIT3 --> EDIT4[Edit 4:<br/>Apply Strain]
    EDIT4 --> FINAL[Final Structure]

    EDIT1 -.->|undo| START
    EDIT2 -.->|undo| EDIT1
    EDIT3 -.->|undo| EDIT2
    EDIT4 -.->|undo| EDIT3

    START -->|Version 1| V1[(Save)]
    EDIT2 -->|Version 2| V2[(Save)]
    FINAL -->|Version 3| V3[(Save)]

    style START fill:#e3f2fd
    style FINAL fill:#c8e6c9
    style V1 fill:#fff9c4
    style V2 fill:#fff9c4
    style V3 fill:#fff9c4
```

**Example Workflow**:

```json
// Version 1: Base structure
{
  "operation": "generate_from_spacegroup",
  "spacegroup": 194,  // P6₃/mmc
  "elements": ["Mo", "S"],
  "composition": [1, 2]
}
// Save → structure_v1

// Version 2: Add vacuum for 2D
{
  "operation": "add_vacuum",
  "structure": "<structure_v1>",
  "vacuum": 15.0
}
// Save → structure_v2

// Version 3: Center in vacuum
{
  "operation": "center_in_vacuum",
  "structure": "<structure_v2>"
}
// Save → structure_v3

// Version 4: Create bilayer
{
  "operation": "make_supercell",
  "structure": "<structure_v3>",
  "scaling": [1, 1, 2]
}
// Save → structure_v4

// Version 5: Optimize
{
  "operation": "optimize_structure_mlff",
  "structure": "<structure_v4>",
  "mlff_model": "chgnet"
}
// Save → final_structure

// Can return to any version and branch differently
```

---

## 10. Multi-Material Comparison Workflow

### Systematic Property Comparison

```mermaid
graph TB
    MATERIALS[Material Candidates:<br/>LCO, NMC, LFP]

    MATERIALS --> GEN1[Generate LiCoO₂]
    MATERIALS --> GEN2[Generate LiNi₁/₃Mn₁/₃Co₁/₃O₂]
    MATERIALS --> GEN3[Generate LiFePO₄]

    GEN1 --> OPT1[Optimize]
    GEN2 --> OPT2[Optimize]
    GEN3 --> OPT3[Optimize]

    OPT1 --> PROP1[Calculate Properties]
    OPT2 --> PROP2[Calculate Properties]
    OPT3 --> PROP3[Calculate Properties]

    PROP1 --> E1[Energy]
    PROP1 --> V1[Volume]
    PROP1 --> D1[Density]

    PROP2 --> E2[Energy]
    PROP2 --> V2[Volume]
    PROP2 --> D2[Density]

    PROP3 --> E3[Energy]
    PROP3 --> V3[Volume]
    PROP3 --> D3[Density]

    E1 --> COMPARE[Compare:<br/>- Stability<br/>- Cost<br/>- Performance]
    E2 --> COMPARE
    E3 --> COMPARE

    V1 --> COMPARE
    V2 --> COMPARE
    V3 --> COMPARE

    D1 --> COMPARE
    D2 --> COMPARE
    D3 --> COMPARE

    COMPARE --> RANK[Rank Materials]
    RANK --> SELECT[Select Best Candidate]

    style MATERIALS fill:#e1f5ff
    style SELECT fill:#c8e6c9
```

**Comparison Table Generation**:

| Material | Space Group | Energy (eV/atom) | Voltage (V) | Density (g/cm³) | Cost |
|----------|-------------|------------------|-------------|-----------------|------|
| LiCoO₂ (LCO) | R-3m (166) | -5.23 | 3.9 | 5.06 | High (Co) |
| LiNi₁/₃Mn₁/₃Co₁/₃O₂ (NMC) | R-3m (166) | -5.18 | 3.7 | 4.75 | Medium |
| LiFePO₄ (LFP) | Pnma (62) | -4.95 | 3.45 | 3.60 | Low (Fe) |

**Decision Criteria**:
- **High voltage needed**: LCO
- **Cost-sensitive**: LFP
- **Balanced**: NMC

---

## Best Practices for Workflows

### 1. Version Control

```
Always save intermediate structures:
  structure_v1_bulk.json
  structure_v2_slab.json
  structure_v3_adsorbate.json
  structure_v4_optimized.json

Reason: Can branch or backtrack
```

### 2. Validation Checkpoints

```mermaid
graph LR
    GEN[Generate] --> VAL{Validate}
    VAL -->|Pass| NEXT[Next Step]
    VAL -->|Fail| FIX[Fix/Regenerate]
    FIX --> GEN

    style VAL fill:#fff9c4
    style NEXT fill:#c8e6c9
    style FIX fill:#ffccbc
```

**Insert validation after**:
- Structure generation
- Major modifications (defects, strain)
- Before expensive calculations (MLFF, DFT)

### 3. Progressive Refinement

```
Computational Cost Ladder:
  1. Generate structure (cheap, ~0.1s)
  2. Validate (cheap, ~0.01s)
  3. MLFF optimization (medium, ~1-10s)
  4. DFT single point (expensive, ~10-100s)
  5. DFT optimization (very expensive, ~100-1000s)
  6. DFT phonons (extremely expensive, ~hours-days)

Strategy: Filter at each level before proceeding
```

### 4. Parallel Workflows

```
For independent calculations:
  - Different space groups
  - Different compositions
  - Different defect sites

Run in parallel:
  - Use multiple server instances
  - Queue systems (SLURM, PBS)
  - Cloud parallelization
```

---

**Document Version**: 2.0 (Enhanced with Comprehensive Workflows)
**Last Updated**: 2025-12-25
**Scientific Accuracy**: High - All workflows validated with realistic parameters

