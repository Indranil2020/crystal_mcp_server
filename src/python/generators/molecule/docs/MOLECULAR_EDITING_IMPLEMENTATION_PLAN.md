# Molecular Editing System - Implementation Plan

## Executive Summary

Add a new `edit_molecule` MCP tool that enables iterative molecular modification workflows. This integrates seamlessly with existing `build_molecule` and `build_molecular_cluster` tools following established architecture patterns.

**Approach:** LLM-Only semantic parsing (simpler, more flexible, matches existing MCP architecture)

---

## Part 1: Current Architecture Context

### Existing Molecular Tools
| Tool | Entry Point (Python) | Handler (TypeScript) | Timeout |
|------|---------------------|---------------------|---------|
| `build_molecule` | `molecule_generator.py` | `build-molecule.ts` | 30s |
| `build_molecular_cluster` | `molecular_cluster_generator.py` | `build-molecular-cluster.ts` | 600s |
| `suggest_molecules` | `molecule_suggester.py` | `suggest-molecules.ts` | 15s |

### Data Flow Pattern (MUST FOLLOW)
```
MCP Request
    ↓
server.ts (CallToolRequestSchema switch)
    ↓
src/tools/generation/<tool>.ts (handler)
    ↓
python-bridge.ts::executePythonWithJSON()
    ↓
Write JSON → /tmp/crystal_mcp_*.json
    ↓
Spawn: python <script>.py /tmp/...json
    ↓
Python reads argv[0], processes, prints JSON to stdout
    ↓
Parse stdout JSON, validate {success: bool, ...}
    ↓
Format response with <json-data> wrapper
    ↓
MCP Response
```

### Key Pattern: Stateless Tool Calls
Each MCP tool call is independent. State (molecule data) must be passed in the request and returned in the response. The client maintains workflow state.

---

## Part 2: New Tool Design

### Tool: `edit_molecule`

**Purpose:** Modify an existing molecule using natural language or structured operations.

**Schema (TypeScript - src/types/tools.ts):**
```typescript
export const EditMoleculeSchema = z.object({
  // Input molecule (from previous build_molecule or edit_molecule call)
  molecule: z.object({
    smiles: z.string().describe("SMILES string of the molecule"),
    annotations: z.any().optional().describe("Rich annotations from previous call")
  }).describe("The molecule to edit, obtained from build_molecule or previous edit"),

  // Operations: Natural language OR structured array
  operations: z.union([
    z.string().describe("Natural language edit instruction, e.g., 'Hydroxylate at C-7 with axial orientation'"),
    z.array(z.object({
      type: z.string().describe("Operation type: HYDROXYLATE, ACETYLATE, EPIMERIZE, DELETE_RING, etc."),
      target: z.string().describe("Target specification: label:C-7, smarts:[OH], ring:oxetane:0"),
      params: z.record(z.any()).optional().describe("Operation-specific parameters")
    })).describe("Structured edit operations array")
  ]).describe("Edit instruction(s) - either natural language or structured operations"),

  // Options
  validate: z.boolean().default(true).optional()
    .describe("Run validation pipeline after edit"),
  optimize_geometry: z.boolean().default(true).optional()
    .describe("Re-optimize 3D geometry after edit"),
  return_3d: z.boolean().default(true).optional()
    .describe("Include 3D coordinates in response")
});
```

**Return Format (matches existing build_molecule pattern):**
```python
{
  "success": True,
  "structure": {
    "lattice": {...},
    "atoms": [...],
    "metadata": {
      "formula": "C9H8O4",
      "natoms": 21,
      "smiles": "...",
      "source": "edit"
    }
  },
  "annotations": {
    "stereocenters": [...],
    "rings": [...],
    "functional_groups": [...],
    "atom_labels": {...},
    "edit_history": [...]
  },
  "edit_result": {
    "operations_applied": [...],
    "validation": {
      "valid": True,
      "warnings": [],
      "errors": []
    }
  }
}
```

---

## Part 3: Implementation Steps

### Step 1: Python Core Module
**Location:** `src/python/generators/molecule/editing/`

```
editing/
├── __init__.py
├── annotated_molecule.py    # AnnotatedMolecule dataclass
├── molecular_editor.py      # RWMol-based editing engine
├── edit_operations.py       # Operation types & execution
├── semantic_parser.py       # NL → operations (LLM-based)
└── validation.py            # Post-edit validation
```

**annotated_molecule.py:**
```python
@dataclass
class AnnotatedMolecule:
    smiles: str
    mol: Chem.RWMol  # RDKit editable molecule

    # Rich annotations
    stereocenters: List[StereocenterInfo]
    rings: List[RingInfo]
    functional_groups: List[FunctionalGroupInfo]
    atom_labels: Dict[int, str]

    # Edit tracking
    edit_history: List[EditOperation]

    def to_dict(self) -> Dict:
        """Serialize for MCP response"""

    @classmethod
    def from_dict(cls, data: Dict) -> 'AnnotatedMolecule':
        """Deserialize from MCP request"""

    @classmethod
    def from_smiles(cls, smiles: str) -> 'AnnotatedMolecule':
        """Create from SMILES with auto-annotation"""
```

**molecular_editor.py:**
```python
class MolecularEditor:
    def __init__(self, molecule: AnnotatedMolecule):
        self.molecule = molecule
        self.rwmol = Chem.RWMol(molecule.mol)

    def apply_operation(self, op: EditOperation) -> EditResult:
        """Execute single edit operation"""

    def apply_operations(self, ops: List[EditOperation]) -> List[EditResult]:
        """Execute multiple operations in sequence"""

    def finalize(self, optimize: bool = True) -> AnnotatedMolecule:
        """Sanitize, validate, optionally optimize, return edited molecule"""
```

**semantic_parser.py (LLM-based):**
```python
def parse_nl_instruction(
    instruction: str,
    molecule: AnnotatedMolecule
) -> List[EditOperation]:
    """
    Parse natural language to structured operations.

    Uses the calling LLM (Claude) via the MCP response mechanism.
    Since this is an MCP tool, we can return a structured prompt
    asking for clarification if needed.

    For v1: Simple pattern matching for common operations.
    For v2: Full LLM integration via external API.
    """
```

### Step 2: Python Entry Point Script
**Location:** `src/python/molecule_editor.py`

```python
#!/usr/bin/env python3
"""
MCP entry point for edit_molecule tool.
Reads JSON from temp file (argv[0]), executes edit, prints JSON to stdout.
"""
import sys
import json
from generators.molecule.editing import (
    AnnotatedMolecule,
    MolecularEditor,
    parse_nl_instruction
)

def main():
    # Read input from temp file
    with open(sys.argv[1], 'r') as f:
        input_data = json.load(f)

    try:
        # Parse input molecule
        mol_data = input_data['molecule']
        if 'annotations' in mol_data:
            molecule = AnnotatedMolecule.from_dict(mol_data)
        else:
            molecule = AnnotatedMolecule.from_smiles(mol_data['smiles'])

        # Parse operations
        operations = input_data['operations']
        if isinstance(operations, str):
            # Natural language → structured operations
            ops = parse_nl_instruction(operations, molecule)
        else:
            ops = [EditOperation.from_dict(op) for op in operations]

        # Apply edits
        editor = MolecularEditor(molecule)
        results = editor.apply_operations(ops)

        # Finalize
        edited = editor.finalize(
            optimize=input_data.get('optimize_geometry', True)
        )

        # Build response
        response = build_response(edited, results, input_data)
        print(json.dumps(response))

    except Exception as e:
        print(json.dumps({
            "success": False,
            "error": {
                "code": "EDIT_FAILED",
                "message": str(e)
            }
        }))
        sys.exit(1)

if __name__ == "__main__":
    main()
```

### Step 3: TypeScript Handler
**Location:** `src/tools/generation/edit-molecule.ts`

```typescript
import { z } from "zod";
import { EditMoleculeSchema } from "../../types/tools.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { createSuccess, createFailure, Result } from "../../utils/result.js";

export async function editMolecule(
  input: unknown
): Promise<Result<any>> {
  // 1. Parse input with Zod schema
  const parsed = EditMoleculeSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure({
      code: "INVALID_INPUT",
      message: "Invalid edit_molecule input",
      details: parsed.error.flatten()
    });
  }

  // 2. Call Python with JSON bridge (60s timeout for complex edits)
  const result = await executePythonWithJSON<typeof parsed.data, any>(
    "molecule_editor.py",
    parsed.data,
    { timeout: 60000 }
  );

  // 3. Error handling
  if (!result.success) {
    return createFailure(result.error);
  }

  return createSuccess(result.data);
}

export async function handleEditMolecule(args: unknown): Promise<any> {
  const result = await editMolecule(args);

  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `Error: ${result.error.message}`
      }],
      isError: true
    };
  }

  const data = result.data;
  const metadata = data.structure?.metadata || {};

  return {
    content: [
      {
        type: "text",
        text: `### Molecule Edited: ${metadata.formula || 'Unknown'}\n` +
              `- Atoms: ${metadata.natoms || 'N/A'}\n` +
              `- Operations applied: ${data.edit_result?.operations_applied?.length || 0}\n` +
              `- Validation: ${data.edit_result?.validation?.valid ? '✓ Valid' : '✗ Issues found'}`
      },
      {
        type: "text",
        text: `<json-data>\n${JSON.stringify(data, null, 2)}\n</json-data>`
      }
    ]
  };
}
```

### Step 4: Register Tool in Server
**Location:** `src/server.ts` (add to switch statement)

```typescript
case "edit_molecule":
  return await handleEditMolecule(args);
```

**Location:** `src/types/tools.ts` (add to TOOL_DEFINITIONS)

```typescript
{
  name: "edit_molecule",
  description: "Modify an existing molecule using natural language or structured operations. " +
    "Supports: hydroxylation, acetylation, epimerization, ring operations, functional group changes. " +
    "Input a molecule from build_molecule and apply edits iteratively.",
  inputSchema: EditMoleculeSchema,
  annotations: {
    readOnlyHint: false,
    destructiveHint: false,
    idempotentHint: true,
    openWorldHint: false
  }
}
```

---

## Part 4: Supported Operations (Phase 1)

| Operation | Target Types | Parameters | Example |
|-----------|-------------|------------|---------|
| `HYDROXYLATE` | label, smarts, atom_idx | orientation | "Hydroxylate at C-7 axial" |
| `ACETYLATE` | label, smarts, atom_idx | orientation | "Acetylate the primary hydroxyl" |
| `METHYLATE` | label, smarts, atom_idx | - | "Methylate the nitrogen" |
| `HALOGENATE` | label, smarts, atom_idx | element (F,Cl,Br,I) | "Brominate the benzene ring" |
| `EPIMERIZE` | label, atom_idx | from_config, to_config | "Epimerize C-3 from R to S" |
| `DELETE_ATOM` | label, smarts, atom_idx | - | "Remove the methyl group" |
| `ADD_BOND` | atom_pair | order (1,2,3) | "Form double bond between C4 and C5" |
| `DELETE_BOND` | atom_pair | - | "Break the ring at C4-C5" |
| `SUBSTITUTE` | smarts | replacement_smarts | "Replace phenyl with furyl" |

---

## Part 5: Files to Create/Modify

### NEW FILES
| File | Purpose |
|------|---------|
| `src/python/molecule_editor.py` | MCP entry point |
| `src/python/generators/molecule/editing/__init__.py` | Module init |
| `src/python/generators/molecule/editing/annotated_molecule.py` | Core data structure |
| `src/python/generators/molecule/editing/molecular_editor.py` | Edit engine |
| `src/python/generators/molecule/editing/edit_operations.py` | Operation definitions |
| `src/python/generators/molecule/editing/semantic_parser.py` | NL parsing |
| `src/python/generators/molecule/editing/validation.py` | Post-edit validation |
| `src/tools/generation/edit-molecule.ts` | TypeScript handler |
| `tests/molecule_test/test_molecular_editor.py` | Unit tests |

### MODIFY FILES
| File | Change |
|------|--------|
| `src/types/tools.ts` | Add EditMoleculeSchema + TOOL_DEFINITIONS entry |
| `src/server.ts` | Add case to switch statement |

---

## Part 6: Verification Plan

### Unit Tests (Python)
```bash
# Test annotated molecule creation and serialization
python -m pytest tests/molecule_test/test_annotated_molecule.py -v

# Test edit operations
python -m pytest tests/molecule_test/test_edit_operations.py -v

# Test molecular editor
python -m pytest tests/molecule_test/test_molecular_editor.py -v
```

### Integration Test (MCP)
```bash
# Build TypeScript
npm run build

# Test via direct MCP call
echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"edit_molecule","arguments":{"molecule":{"smiles":"c1ccccc1"},"operations":"hydroxylate the benzene ring"}}}' | node dist/index.js
```

### End-to-End Workflow Test
```python
# 1. Generate base molecule
result1 = call_mcp("build_molecule", {"name": "aspirin"})
molecule = result1["structure"]

# 2. Edit: Add hydroxyl group
result2 = call_mcp("edit_molecule", {
    "molecule": {"smiles": molecule["metadata"]["smiles"]},
    "operations": "Add hydroxyl group to the benzene ring"
})

# 3. Verify edit succeeded
assert result2["success"] == True
assert "O" in result2["structure"]["metadata"]["formula"]
```

---

## Part 7: Implementation Order

1. **Core Data Structures** (Day 1)
   - `annotated_molecule.py` - AnnotatedMolecule class
   - `edit_operations.py` - EditOperation enum and classes

2. **Molecular Editor** (Day 2)
   - `molecular_editor.py` - RWMol-based editing
   - Basic operations: HYDROXYLATE, METHYLATE, DELETE_ATOM

3. **MCP Integration** (Day 3)
   - `molecule_editor.py` - Python entry point
   - `edit-molecule.ts` - TypeScript handler
   - `tools.ts` and `server.ts` modifications

4. **Semantic Parser** (Day 4)
   - `semantic_parser.py` - Pattern-based NL parsing
   - Support common chemistry verbs

5. **Validation & Testing** (Day 5)
   - `validation.py` - Post-edit validation
   - Unit tests and integration tests

---

## Part 8: Non-Breaking Integration

### Guarantees
1. **No changes to existing tools** - `build_molecule`, `build_molecular_cluster`, `suggest_molecules` unchanged
2. **Additive only** - New files, new switch case, new schema entry
3. **Same response format** - Follows `{success, structure, ...}` pattern
4. **Same Python bridge** - Uses `executePythonWithJSON()` like other tools
5. **Backward compatible annotations** - If no annotations provided, auto-generate from SMILES

### Risk Mitigation
- All new code in isolated `editing/` module
- New entry point script (not modifying existing generators)
- TypeScript handler follows exact pattern of `build-molecule.ts`
- Comprehensive tests before merge

---

## Summary

This plan adds molecular editing capabilities by:
1. Creating a new `editing/` Python module with RDKit-based editing
2. Adding `edit_molecule` MCP tool following existing patterns
3. Supporting both natural language and structured operations
4. Preserving rich annotations through edit workflows
5. NOT modifying any existing functionality
