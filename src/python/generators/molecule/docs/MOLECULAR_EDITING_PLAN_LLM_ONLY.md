# Molecular Editing & Synthesis Simulation System - LLM-Only Approach

## Philosophy & Vision

Same core paradigm as the hybrid approach:

```
GENERATE → EDIT → VALIDATE → ITERATE
```

**Key Difference:** ALL natural language parsing is handled by Claude, eliminating the pattern-based parser entirely. This creates a simpler, more flexible system at the cost of API dependency.

---

## Core Architecture

### 1. Annotated Molecule State

*(Identical to Hybrid approach)*

```
┌─────────────────────────────────────────────────────┐
│               AnnotatedMolecule                     │
├─────────────────────────────────────────────────────┤
│ CORE: smiles, inchi                                 │
│ STEREOCHEMISTRY: stereocenters, double_bond_stereo  │
│ RING SYSTEMS: rings, fused_ring_map                 │
│ FUNCTIONAL GROUPS: groups, protecting_groups        │
│ ATOM LABELS: atom_labels, atom_hybridization        │
│ 3D CONFORMER: coordinates, ring_conformations       │
│ EDIT HISTORY: [{operation, target, params, ts}]     │
└─────────────────────────────────────────────────────┘
```

### 2. Edit Operation DSL

*(Identical to Hybrid approach - same structured output format)*

| Category | Operations |
|----------|------------|
| **Atom** | ADD_ATOM, DELETE_ATOM, CHANGE_ELEMENT, SET_CHARGE |
| **Bond** | ADD_BOND, DELETE_BOND, CHANGE_BOND_ORDER |
| **Functional Group** | ADD_GROUP, REMOVE_GROUP, SUBSTITUTE, PROTECT, DEPROTECT |
| **Ring** | ADD_RING, DELETE_RING, FUSE_RING, EXPAND, CONTRACT, OPEN, CLOSE |
| **Stereochemistry** | SET_STEREO, INVERT, EPIMERIZE, SET_E_Z, SET_AXIAL_EQUATORIAL |
| **Oxidation/Reduction** | OXIDIZE, REDUCE, HYDROXYLATE, HALOGENATE |
| **Conformational** | SET_TORSION, SET_RING_CONFORMATION |

### 3. Semantic Parser - LLM-ONLY APPROACH

**Architecture: Direct Claude API Integration**

```
┌─────────────────────────────────────────────────────────────┐
│                    LLM-ONLY PARSER                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  INPUT: "Hydroxylate at C-7 with axial orientation"         │
│                          ↓                                  │
│  ┌─────────────────────────────────────────────┐            │
│  │ MOLECULE CONTEXT BUILDER                    │            │
│  │ • Serialize AnnotatedMolecule to JSON       │            │
│  │ • Include current stereocenters             │            │
│  │ • Include ring systems and atom labels      │            │
│  └──────────────────┬──────────────────────────┘            │
│                     ↓                                       │
│  ┌─────────────────────────────────────────────┐            │
│  │ CLAUDE API CALL                             │            │
│  │ • Tool use with structured output schema    │            │
│  │ • Returns EditOperation[] directly          │            │
│  │ • Handles ALL complexity levels             │            │
│  └──────────────────┬──────────────────────────┘            │
│                     ↓                                       │
│  ┌─────────────────────────────────────────────┐            │
│  │ RESPONSE VALIDATION                         │            │
│  │ • Verify operations reference valid atoms   │            │
│  │ • Check target resolution                   │            │
│  │ • Validate parameter types                  │            │
│  └─────────────────────────────────────────────┘            │
└─────────────────────────────────────────────────────────────┘
```

**Claude Tool Definition:**

```python
PARSE_TOOL = {
    "name": "extract_molecular_edits",
    "description": "Extract structured molecular editing operations from natural language",
    "input_schema": {
        "type": "object",
        "properties": {
            "operations": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "type": {
                            "type": "string",
                            "enum": ["HYDROXYLATE", "ACETYLATE", "EPIMERIZE",
                                     "DELETE_RING", "ADD_BOND", "SUBSTITUTE", ...]
                        },
                        "target": {
                            "type": "string",
                            "description": "Target specification (label:C-7, smarts:[OH], ring:oxetane:0)"
                        },
                        "params": {
                            "type": "object",
                            "description": "Operation-specific parameters"
                        }
                    },
                    "required": ["type", "target"]
                }
            },
            "reasoning": {
                "type": "string",
                "description": "Explanation of interpretation decisions"
            }
        },
        "required": ["operations"]
    }
}
```

**System Prompt for Parsing:**

```python
SYSTEM_PROMPT = """
You are a molecular editing interpreter. Given a molecule's current state and a
natural language instruction, extract structured editing operations.

MOLECULE STATE:
{molecule_json}

OPERATION TYPES:
- HYDROXYLATE: Add -OH group (params: orientation=axial|equatorial)
- ACETYLATE: Add acetyl group (params: orientation)
- EPIMERIZE: Invert stereocenter (params: from=R|S, to=R|S)
- DELETE_RING: Remove a ring system (params: keep_chain=true|false)
- ADD_BOND: Add bond between atoms (params: order=1|2|3)
- SUBSTITUTE: Replace one group with another (params: replacement)
... [full operation reference]

TARGET SPECIFICATION:
- label:X → Atom with label X (e.g., label:C-7)
- smarts:PATTERN → First match of SMARTS pattern
- ring:TYPE:N → Nth ring of TYPE (e.g., ring:oxetane:0)
- group:TYPE:N → Nth functional group of TYPE

Rules:
1. Resolve ambiguous references using molecule context
2. If multiple interpretations exist, choose the most chemically reasonable
3. Include reasoning for complex decisions
4. Reject instructions that would create invalid chemistry
"""
```

**API Call Flow:**

```python
async def parse_edit_instruction(
    molecule: AnnotatedMolecule,
    instruction: str
) -> List[EditOperation]:
    """Parse NL instruction using Claude API."""

    # Build context
    context = {
        "smiles": molecule.smiles,
        "stereocenters": molecule.stereocenters,
        "rings": molecule.rings,
        "atom_labels": molecule.atom_labels,
        "functional_groups": molecule.functional_groups
    }

    # Call Claude with tool use
    response = await anthropic.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=1024,
        system=SYSTEM_PROMPT.format(molecule_json=json.dumps(context)),
        messages=[{
            "role": "user",
            "content": instruction
        }],
        tools=[PARSE_TOOL],
        tool_choice={"type": "tool", "name": "extract_molecular_edits"}
    )

    # Extract structured output
    tool_use = next(
        block for block in response.content
        if block.type == "tool_use"
    )

    # Validate and return
    operations = validate_operations(tool_use.input["operations"], molecule)
    return operations
```

### 4. Caching Layer (Performance Optimization)

To mitigate latency and API costs, add a semantic cache:

```
┌─────────────────────────────────────────────────────────────┐
│                    SEMANTIC CACHE                           │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  CACHE KEY:                                                 │
│    hash(instruction_normalized + molecule_fingerprint)      │
│                                                             │
│  NORMALIZATION:                                             │
│    "Hydroxylate at C-7 axial" → "hydroxylate c-7 axial"     │
│                                                             │
│  FINGERPRINT:                                               │
│    Only include structurally-relevant fields                │
│    (stereocenters, rings, labels - NOT coordinates)         │
│                                                             │
│  EXPIRATION:                                                │
│    TTL: 24 hours (instructions don't change meaning)        │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### 5. Batch Processing Mode

For multi-step workflows, batch instructions in a single API call:

```python
async def parse_batch_instructions(
    molecule: AnnotatedMolecule,
    instructions: List[str]
) -> List[List[EditOperation]]:
    """Parse multiple instructions in single API call."""

    combined_prompt = """
    Parse each of these instructions as separate edit operations.
    Return an array of operation arrays, one per instruction.

    Instructions:
    1. {instructions[0]}
    2. {instructions[1]}
    ...
    """

    # Single API call for all instructions
    response = await anthropic.messages.create(...)

    return response.tool_use.input["instruction_operations"]
```

### 6. Engines (Same as Hybrid)

**Stereochemistry Engine**, **Ring System Engine**, **Validation Pipeline** - all identical to the hybrid approach. These are pure RDKit operations that don't involve NL parsing.

---

## Module Organization (Simplified)

```
src/python/generators/molecule/
├── editing/                        # NEW MODULE
│   ├── __init__.py
│   ├── annotated_molecule.py       # Core state data structures
│   ├── edit_operations.py          # Operation type definitions
│   ├── molecular_editor.py         # Main editing engine (RWMol)
│   ├── llm_parser.py               # Claude API integration (SIMPLIFIED)
│   ├── semantic_cache.py           # Caching layer
│   │
│   ├── stereochemistry/            # Same as hybrid
│   ├── rings/                      # Same as hybrid
│   ├── functional_groups/          # Same as hybrid
│   └── validation/                 # Same as hybrid
```

**Key Difference:** `semantic_parser.py` (complex hybrid) replaced with `llm_parser.py` (simple API wrapper)

---

## Implementation Phases

### Phase 1: Core Infrastructure
- [ ] `AnnotatedMolecule` data structure
- [ ] `MolecularEditor` class with RWMol
- [ ] Basic atom/bond operations
- [ ] Simple validation

### Phase 2: LLM Parser
- [ ] Claude API integration
- [ ] Tool schema definition
- [ ] System prompt engineering
- [ ] Response validation
- [ ] Error handling and retries

### Phase 3: Caching & Optimization
- [ ] Semantic cache implementation
- [ ] Instruction normalization
- [ ] Molecule fingerprinting
- [ ] Batch processing mode

### Phase 4: Chemistry Engines
- [ ] Stereochemistry engine (RDKit)
- [ ] Ring system engine (RWMol)
- [ ] Validation pipeline

### Phase 5: MCP Integration
- [ ] `edit_molecule` tool schema
- [ ] Python handler bridge
- [ ] End-to-end testing

---

## Comparison: LLM-Only Approach

### Advantages

| Advantage | Description |
|-----------|-------------|
| **Simpler codebase** | No regex patterns, no chemical knowledge base maintenance |
| **Zero-shot flexibility** | Handles novel instructions without code changes |
| **Better disambiguation** | LLM understands context and chemistry semantics |
| **Self-documenting** | Reasoning field explains interpretation decisions |
| **Future-proof** | Model improvements automatically enhance parsing |
| **Fewer edge cases** | No pattern-matching bugs or corner cases |

### Disadvantages

| Disadvantage | Description | Mitigation |
|--------------|-------------|------------|
| **API dependency** | Requires network connection | Caching layer |
| **Latency** | ~500ms per parse call | Batch processing |
| **Cost** | API usage costs | Caching, use Haiku for simple cases |
| **Non-deterministic** | Same input may produce different output | Structured output schema, temperature=0 |
| **Privacy** | Molecule data sent to API | Local deployment option |

### Cost Analysis

```
Assumptions:
- 100 edits/day per user
- Average instruction: 50 tokens
- Average molecule context: 200 tokens
- Average response: 100 tokens
- Using Claude Haiku (cheapest)

Cost calculation:
- Input: (50 + 200) × 100 = 25,000 tokens/day
- Output: 100 × 100 = 10,000 tokens/day
- Haiku pricing: $0.25/M input, $1.25/M output
- Daily cost: $0.006 + $0.013 = $0.02/user/day
- Monthly cost: ~$0.60/user/month

With 80% cache hit rate:
- Actual API calls: 20 edits/day
- Monthly cost: ~$0.12/user/month
```

### Latency Analysis

```
LLM parse call: ~500ms (median)
- 300ms: Network + API overhead
- 200ms: Model inference

With caching (80% hit rate):
- Cache hit: ~5ms (local lookup)
- Weighted average: 0.8×5 + 0.2×500 = 104ms

With batch processing (5 instructions):
- Single API call: 600ms
- Per instruction: 120ms
```

---

## Critical Design Decisions

### 1. Model Selection Strategy

```python
def select_model(instruction: str, molecule: AnnotatedMolecule) -> str:
    """Choose model based on complexity."""

    # Simple cases → Haiku (fast, cheap)
    simple_patterns = ["hydroxylate", "acetylate", "methylate"]
    if any(p in instruction.lower() for p in simple_patterns):
        if len(molecule.rings) <= 3:
            return "claude-3-haiku-20240307"

    # Complex cases → Sonnet (better reasoning)
    return "claude-sonnet-4-20250514"
```

### 2. Structured Output Enforcement

```python
# Force tool use - no free-form text
tool_choice={"type": "tool", "name": "extract_molecular_edits"}
```

### 3. Temperature Control

```python
# Deterministic output for reproducibility
temperature=0
```

### 4. Retry Strategy

```python
async def parse_with_retry(instruction: str, max_retries: int = 3):
    for attempt in range(max_retries):
        try:
            result = await parse_edit_instruction(molecule, instruction)
            if validate_operations(result):
                return result
        except (APIError, ValidationError) as e:
            if attempt == max_retries - 1:
                raise
            await asyncio.sleep(2 ** attempt)  # Exponential backoff
```

---

## MCP Tool Design

*(Same as hybrid approach)*

```typescript
EditMoleculeSchema = z.object({
  molecule: z.object({
    smiles: z.string(),
    annotations: z.any().optional()
  }),
  operations: z.union([
    z.string(),  // NL always goes to LLM
    z.array(EditOperationSchema)  // Structured skips LLM
  ]),
  validate: z.boolean().default(true),
  optimize_geometry: z.boolean().default(true),
  return_3d: z.boolean().default(true)
})
```

**Important:** When operations is already a structured array, skip the LLM entirely and execute directly.

---

## Verification Plan

### Unit Tests
```bash
# Test LLM parser with mocked API
pytest tests/molecule_test/test_llm_parser.py

# Test caching layer
pytest tests/molecule_test/test_semantic_cache.py

# Test full editing pipeline
pytest tests/molecule_test/test_molecular_editor.py
```

### Integration Test: Taxol Sequence
```python
# Same test as hybrid - system behavior should be identical
mol = generate_molecule("taxol core")
mol = edit_molecule(mol, "Hydroxylate at C-7 axial, acetylate at C-10")
mol = edit_molecule(mol, "Replace phenyl with 3-furyl, epimerize C-3 R→S")
mol = edit_molecule(mol, "Delete oxetane, form double bond at fusion points")
```

---

## When to Choose LLM-Only

**Choose LLM-Only when:**
- Development speed is priority
- User instructions are highly variable
- Chemistry domain is broad (not specialized)
- API costs are acceptable
- Network connectivity is reliable

**Choose Hybrid when:**
- Offline operation required
- Strict latency requirements (<100ms)
- Very high volume (millions of requests/day)
- Need complete determinism
- Cost is critical concern

---

## Summary Comparison

| Aspect | Hybrid | LLM-Only |
|--------|--------|----------|
| **Codebase size** | ~2000 LOC parser | ~200 LOC parser |
| **Maintenance** | Pattern updates needed | Prompt refinement only |
| **Latency (cold)** | ~10ms | ~500ms |
| **Latency (cached)** | N/A | ~5ms |
| **Flexibility** | Constrained by patterns | Unlimited |
| **Determinism** | 100% | 95%+ with temp=0 |
| **Offline** | Yes | No |
| **Cost/month** | $0 | ~$0.12-0.60/user |
