# Molecular Editing & Synthesis Simulation System - Hybrid Approach

## Philosophy & Vision

The essence of the Taxol example is NOT about Taxol specifically—it's about a **generic molecular editing paradigm** that enables:

```
GENERATE → EDIT → VALIDATE → ITERATE
```

This transforms the current **stateless generation** model into a **workflow-enabled editing** model where:
- Molecules maintain rich annotations through multiple edits
- Natural language chemical operations are parsed into structured edits
- Stereochemistry, ring systems, and functional groups are first-class citizens
- Each edit is validated before proceeding

---

## Core Architecture

### 1. Annotated Molecule State

```
┌─────────────────────────────────────────────────────┐
│               AnnotatedMolecule                     │
├─────────────────────────────────────────────────────┤
│ CORE: smiles, inchi                                 │
│                                                     │
│ STEREOCHEMISTRY:                                    │
│   - stereocenters: [{atom_idx, R/S, neighbors}]     │
│   - double_bond_stereo: [{atoms, E/Z}]              │
│   - axial_chirality: [atropisomers]                 │
│   - helicity: P/M                                   │
│                                                     │
│ RING SYSTEMS:                                       │
│   - rings: [{id, atoms, size, type, aromatic}]      │
│   - fused_ring_map: {atom → ring_ids}               │
│   - spiro_centers, bridged_systems                  │
│                                                     │
│ FUNCTIONAL GROUPS:                                  │
│   - groups: [{type, atoms, attachment_point}]       │
│   - protecting_groups: [{type, protected_group}]    │
│                                                     │
│ ATOM LABELS:                                        │
│   - atom_labels: {idx → "C-7", "C-13", etc.}        │
│   - atom_hybridization, atom_charges                │
│                                                     │
│ 3D CONFORMER (optional):                            │
│   - coordinates, ring_conformations                 │
│   - axial_equatorial assignments                    │
│   - steric_clashes, torsion_angles                  │
│                                                     │
│ EDIT HISTORY:                                       │
│   - [{operation, target, params, timestamp}]        │
└─────────────────────────────────────────────────────┘
```

### 2. Edit Operation DSL

**Operation Categories:**

| Category | Operations |
|----------|------------|
| **Atom** | ADD_ATOM, DELETE_ATOM, CHANGE_ELEMENT, SET_CHARGE |
| **Bond** | ADD_BOND, DELETE_BOND, CHANGE_BOND_ORDER |
| **Functional Group** | ADD_GROUP, REMOVE_GROUP, SUBSTITUTE, PROTECT, DEPROTECT |
| **Ring** | ADD_RING, DELETE_RING, FUSE_RING, EXPAND, CONTRACT, OPEN, CLOSE |
| **Stereochemistry** | SET_STEREO, INVERT, EPIMERIZE, SET_E_Z, SET_AXIAL_EQUATORIAL |
| **Oxidation/Reduction** | OXIDIZE, REDUCE, HYDROXYLATE, HALOGENATE |
| **Conformational** | SET_TORSION, SET_RING_CONFORMATION |

**Target Specification DSL:**

```python
# By atom label
target=label:C-7

# By SMARTS pattern
target=smarts:[OH]

# By ring position
target=ring:0:position:4

# By functional group
target=group:hydroxyl:0
```

**Example DSL Strings:**

```
HYDROXYLATE(target=label:C-7, orientation=axial)
SUBSTITUTE_GROUP(target=smarts:[c]1ccccc1, replacement=smarts:c1ccoc1)
EPIMERIZE(target=label:C-3, from=R, to=S)
DELETE_RING(target=ring:oxetane:0)
ADD_BOND(target=[label:C-4, label:C-5], order=2)
```

### 3. Semantic Parser - HYBRID APPROACH

**Architecture: Pattern-Based + LLM Fallback**

```
┌─────────────────────────────────────────────────────────────┐
│                    HYBRID PARSER                            │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  INPUT: "Hydroxylate at C-7 with axial orientation"         │
│                          ↓                                  │
│  ┌─────────────────────────────────────────────┐            │
│  │ TIER 1: PATTERN-BASED PARSER                │            │
│  │ • Regex patterns for common verbs           │            │
│  │ • Chemical knowledge base lookup            │            │
│  │ • Fast, deterministic, offline              │            │
│  └──────────────────┬──────────────────────────┘            │
│                     │                                       │
│         ┌──────────┴──────────┐                             │
│         │  Confidence > 0.8?  │                             │
│         └──────────┬──────────┘                             │
│              YES   │   NO                                   │
│               ↓    │    ↓                                   │
│  ┌────────────┐    │    ┌────────────────────────────────┐  │
│  │ RETURN     │    │    │ TIER 2: LLM-ASSISTED PARSER    │  │
│  │ OPERATION  │    │    │ • Send to Claude with schema   │  │
│  └────────────┘    │    │ • Structured output extraction │  │
│                    │    │ • Handles complex/ambiguous    │  │
│                    │    └────────────────────────────────┘  │
│                    │                  ↓                     │
│                    │         ┌────────────────┐             │
│                    │         │ CACHE RESULT   │             │
│                    │         │ for future use │             │
│                    │         └────────────────┘             │
└─────────────────────────────────────────────────────────────┘
```

**Tier 1: Pattern-Based (Handles ~80% of Cases)**

```
"Hydroxylate the 8-membered ring at C-7 with axial orientation"
    ↓
┌─────────────────────────────────────────┐
│ 1. VERB EXTRACTION                      │
│    "hydroxylate" → HYDROXYLATE          │
├─────────────────────────────────────────┤
│ 2. TARGET EXTRACTION                    │
│    "at C-7" → target=label:C-7          │
├─────────────────────────────────────────┤
│ 3. PARAMETER EXTRACTION                 │
│    "axial orientation" → orientation=axial│
├─────────────────────────────────────────┤
│ 4. CONTEXT VALIDATION                   │
│    "8-membered ring" → verify C-7 in ring│
└─────────────────────────────────────────┘
    ↓
EditOperation(type=HYDROXYLATE, target=label:C-7, params={orientation: axial})
```

**Tier 2: LLM-Assisted (Complex Cases)**

For ambiguous or complex instructions like:
- "Modify the side chain to avoid 1,3-diaxial interactions"
- "Create an α,β-unsaturated ketone system maintaining helicity"
- "Simultaneously epimerize C-3 while preserving the oxetane"

```python
# LLM prompt template
PARSE_PROMPT = """
Given the molecular editing instruction and current molecule state,
extract structured edit operations.

Molecule: {smiles}
Annotations: {annotations}
Instruction: "{instruction}"

Return JSON array of operations:
[{{"type": "...", "target": "...", "params": {{...}}}}]
"""
```

**Learning Pipeline:**
- Cache successful LLM parses
- Periodically extract patterns from cache
- Promote common patterns to Tier 1

**Chemical Knowledge Base:**
- Functional group synonyms (acetyl=Ac=OAc, hydroxyl=OH, etc.)
- Ring type patterns (oxetane=4-membered+O, epoxide=3-membered+O, etc.)
- Reaction product templates
- Conformational terms (axial, equatorial, cis, trans, syn, anti)

### 4. Stereochemistry Engine

**Capabilities:**

1. **R/S Detection & Manipulation**
   - `Chem.FindMolChiralCenters()` → detect
   - `Chem.AssignStereochemistryFrom3D()` → infer from 3D
   - `atom.SetChiralTag()` → modify

2. **E/Z Double Bond Handling**
   - Detect via `bond.GetStereo()`
   - Flip by swapping substituent connectivity

3. **Axial/Equatorial Detection**
   - Calculate ring plane normal from 3D coords
   - Dot product of substituent vector with normal
   - |dot| > 0.5 → axial, else equatorial

4. **Helicity Tracking**
   - P-helicity (right-handed) vs M-helicity (left-handed)
   - Preserve through edits when requested

### 5. Ring System Engine

**Operations:**

```python
# Ring Fusion: Fuse new ring at positions 4-5
fuse_ring(mol, attachment=(C4, C5), new_ring_size=6, fusion_type="cis")
  → Adds (size-2) new atoms
  → Connects: C4 → new1 → new2 → ... → newN → C5
  → Sets cis/trans fusion stereochemistry

# Ring Deletion: Remove oxetane while keeping/breaking chain
delete_ring(mol, ring_id=0, keep_chain=True)
  → Breaks one bond to open ring
  → Or removes all ring atoms

# Ring Opening with Bond Formation
open_ring_form_double_bond(mol, ring_id, break_atoms=(C4, C5))
  → Deletes ring bond
  → Adds double bond at break points
  → Creates α,β-unsaturated system
```

### 6. Validation Pipeline

```
┌─────────────────────────────────────────┐
│ 1. VALENCE CHECK                        │
│    Verify all atoms have valid valence  │
├─────────────────────────────────────────┤
│ 2. AROMATICITY CONSISTENCY              │
│    Verify aromatic systems are complete │
├─────────────────────────────────────────┤
│ 3. RING STRAIN ANALYSIS                 │
│    Flag high-strain rings (3,4-membered)│
├─────────────────────────────────────────┤
│ 4. STEREOCHEMISTRY VALIDITY             │
│    Verify R/S assignments are consistent│
├─────────────────────────────────────────┤
│ 5. STERIC CLASH DETECTION               │
│    Check VDW overlaps in 3D             │
│    Detect 1,3-diaxial interactions      │
├─────────────────────────────────────────┤
│ 6. CHEMICAL REASONABILITY               │
│    Hypervalent atoms, unusual bonds     │
└─────────────────────────────────────────┘
```

---

## MCP Tool Design

### Primary Tool: `edit_molecule`

```typescript
EditMoleculeSchema = z.object({
  // Input molecule (from previous generation/edit)
  molecule: z.object({
    smiles: z.string(),
    annotations: z.any().optional()  // Full AnnotatedMolecule
  }),

  // Operations: NL string OR structured array
  operations: z.union([
    z.string(),  // "Hydroxylate at C-7 with axial orientation"
    z.array(EditOperationSchema)  // Structured operations
  ]),

  // Options
  validate: z.boolean().default(true),
  optimize_geometry: z.boolean().default(true),
  return_3d: z.boolean().default(true)
})
```

**Stateless Workflow Pattern:**

```
Client                          Server
  │                               │
  ├── build_molecule(taxol) ─────►│
  │◄── {smiles, annotations} ─────┤
  │                               │
  ├── edit_molecule({            │
  │     molecule: {...},          │
  │     operations: "hydroxylate" │
  │   }) ─────────────────────────►│
  │◄── {smiles, annotations} ─────┤
  │                               │
  ├── edit_molecule({            │
  │     molecule: {...},          │
  │     operations: "epimerize"   │
  │   }) ─────────────────────────►│
  │◄── {smiles, annotations} ─────┤
```

Client maintains state by passing full molecule in each call.

---

## Module Organization

```
src/python/generators/molecule/
├── editing/                        # NEW MODULE
│   ├── __init__.py
│   ├── annotated_molecule.py       # Core state data structures
│   ├── edit_operations.py          # Operation type definitions
│   ├── molecular_editor.py         # Main editing engine (RWMol)
│   ├── semantic_parser.py          # NL → operations (HYBRID)
│   ├── chemical_knowledge.py       # Functional groups, ring types
│   │
│   ├── stereochemistry/
│   │   ├── stereo_engine.py        # R/S, E/Z manipulation
│   │   ├── stereo_inference.py     # 3D → stereo inference
│   │   └── axial_equatorial.py     # Ring orientation
│   │
│   ├── rings/
│   │   ├── ring_engine.py          # Fusion, deletion, expansion
│   │   ├── ring_classifier.py      # Type detection
│   │   └── ring_conformations.py   # Chair/boat analysis
│   │
│   ├── functional_groups/
│   │   ├── group_detector.py       # SMARTS-based detection
│   │   └── group_operations.py     # Add/remove/substitute
│   │
│   └── validation/
│       ├── validation_pipeline.py  # Orchestrator
│       ├── valence_validator.py
│       ├── steric_validator.py
│       └── stereo_validator.py
```

---

## Implementation Phases

### Phase 1: Core Infrastructure (Foundation)
- [ ] `AnnotatedMolecule` data structure with JSON serialization
- [ ] `MolecularEditor` class wrapping RDKit `RWMol`
- [ ] Basic atom/bond operations (add, delete, modify)
- [ ] Simple validation (valence, sanitization)

### Phase 2: Stereochemistry Engine
- [ ] R/S detection using `FindMolChiralCenters`
- [ ] Stereocenter manipulation with `SetChiralTag`
- [ ] E/Z double bond handling
- [ ] 3D inference with `AssignStereochemistryFrom3D`
- [ ] Axial/equatorial detection from 3D coords

### Phase 3: Ring System Engine
- [ ] Ring analysis with `GetRingInfo`, `GetSSSR`
- [ ] Ring type classification (oxetane, benzene, etc.)
- [ ] Ring fusion operations
- [ ] Ring deletion/opening
- [ ] Ring expansion/contraction

### Phase 4: Semantic Parser (Hybrid)
- [ ] Verb pattern matching (hydroxylate, acetylate, etc.)
- [ ] Target extraction (C-7, position 4, etc.)
- [ ] Parameter extraction (axial, R→S, etc.)
- [ ] Chemical knowledge base integration
- [ ] LLM fallback for complex cases
- [ ] Caching and pattern promotion

### Phase 5: MCP Integration
- [ ] `edit_molecule` tool schema in TypeScript
- [ ] Python handler bridge
- [ ] `analyze_molecule` supporting tool
- [ ] End-to-end workflow testing

---

## Comparison: Hybrid Approach

### Advantages
1. **Fast for common cases** - 80% handled by regex, no API latency
2. **Deterministic** - Same input always produces same output
3. **Offline capable** - Works without network connection
4. **Cost-efficient** - Only uses LLM when necessary
5. **Testable** - Pattern rules can be unit tested

### Disadvantages
1. **Pattern maintenance** - Need to maintain regex library
2. **Limited flexibility** - New patterns require code changes
3. **Two codepaths** - More complex than single approach

### Estimated Coverage
- **Tier 1 (Pattern)**: ~80% of real-world queries
- **Tier 2 (LLM)**: ~20% complex/ambiguous queries

---

## Critical Files to Modify/Create

| File | Action | Purpose |
|------|--------|---------|
| `src/python/generators/molecule/editing/__init__.py` | CREATE | New editing module |
| `src/python/generators/molecule/editing/annotated_molecule.py` | CREATE | Core data structures |
| `src/python/generators/molecule/editing/molecular_editor.py` | CREATE | RWMol-based editor |
| `src/python/generators/molecule/editing/semantic_parser.py` | CREATE | Hybrid NL parsing |
| `src/python/generators/molecule/conformers.py` | EXTEND | Leverage existing stereo code |
| `src/types/tools.ts` | EXTEND | Add EditMoleculeSchema |
| `src/tools/generation/edit-molecule.ts` | CREATE | MCP tool handler |

---

## Verification Plan

### Unit Tests
```bash
# Test stereochemistry engine
pytest tests/molecule_test/test_stereo_engine.py

# Test ring operations
pytest tests/molecule_test/test_ring_engine.py

# Test semantic parser (both tiers)
pytest tests/molecule_test/test_semantic_parser.py
```

### Integration Test: Taxol Sequence
```python
# Generate taxol core
mol = generate_molecule("taxol core")

# Edit 1: Hydroxylate + Acetylate
mol = edit_molecule(mol, "Hydroxylate at C-7 axial, acetylate at C-10")
assert mol.stereocenters[C7].orientation == "axial"

# Edit 2: Replace phenyl + Epimerize
mol = edit_molecule(mol, "Replace phenyl with 3-furyl, epimerize C-3 R→S")
assert "furan" in mol.smiles
assert mol.stereocenters[C3].config == "S"

# Edit 3: Delete oxetane + Form double bond
mol = edit_molecule(mol, "Delete oxetane, form double bond at fusion points")
assert "oxetane" not in [r.type for r in mol.rings]
```
