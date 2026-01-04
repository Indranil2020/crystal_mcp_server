# Crystal MCP Server - Architecture Documentation

## Overview

This document describes the architecture of the Crystal MCP Server system, which enables natural language control of crystal structure generation through the Model Context Protocol (MCP).

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              USER INPUT                                      │
│                     "Create a Ge (100) slab with 4 layers"                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           LLM CLIENT (GUI)                                   │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 1. CONNECT TO MCP SERVER                                            │    │
│  │    - Establish stdio connection to MCP server process               │    │
│  │    - Send: initialize request                                       │    │
│  │    - Receive: server capabilities                                   │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 2. FETCH TOOL REGISTRY                                              │    │
│  │    - Send: tools/list request                                       │    │
│  │    - Receive: Array of tool definitions                             │    │
│  │      - name: string                                                 │    │
│  │      - description: string                                          │    │
│  │      - inputSchema: JSON Schema object                              │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 3. BUILD LLM PROMPT                                                 │    │
│  │    - System message with:                                           │    │
│  │      - Role definition                                              │    │
│  │      - Tool descriptions (from MCP)                                 │    │
│  │      - Input schemas (from MCP)                                     │    │
│  │      - Output format instructions                                   │    │
│  │    - User message (natural language request)                        │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 4. SEND TO LLM (Ollama)                                             │    │
│  │    - POST to Ollama API                                             │    │
│  │    - Model processes prompt                                         │    │
│  │    - Returns structured JSON response                               │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           LLM OUTPUT                                         │
│                                                                              │
│  {                                                                           │
│    "tool": "generate_slab",                                                  │
│    "params": {                                                               │
│      "structure": {"prototype": "diamond", "elements": {"A": "Ge"}, ...},   │
│      "miller_indices": [1, 0, 0],                                           │
│      "thickness": 4,                                                         │
│      "vacuum": 15.0                                                          │
│    }                                                                         │
│  }                                                                           │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           LLM CLIENT (GUI)                                   │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 5. PARSE LLM OUTPUT                                                 │    │
│  │    - Extract tool name                                              │    │
│  │    - Extract parameters                                             │    │
│  │    - Validate tool exists in registry                               │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 6. CALL MCP TOOL                                                    │    │
│  │    - Send: tools/call request                                       │    │
│  │      - name: "generate_slab"                                        │    │
│  │      - arguments: {structure: ..., miller_indices: ...}             │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           MCP SERVER (Node.js)                               │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 7. ROUTE TO TOOL HANDLER                                            │    │
│  │    - Match tool name to handler                                     │    │
│  │    - Validate input against schema                                  │    │
│  │    - Prepare parameters for Python backend                          │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 8. EXECUTE PYTHON BACKEND                                           │    │
│  │    - Spawn Python process                                           │    │
│  │    - Pass JSON input via stdin/file                                 │    │
│  │    - Receive JSON output                                            │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
│                                     │                                        │
│                                     ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 9. FORMAT RESPONSE                                                  │    │
│  │    - Human-readable text (Markdown)                                 │    │
│  │    - Machine-readable JSON in <json-data> block                     │    │
│  │    - Return via MCP protocol                                        │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           LLM CLIENT (GUI)                                   │
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │ 10. PROCESS RESULT                                                  │    │
│  │     - Display text output in chat                                   │    │
│  │     - Extract structure from <json-data> block                      │    │
│  │     - Render structure in 3D viewer                                 │    │
│  └─────────────────────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           USER SEES                                          │
│                                                                              │
│  - Chat response with structure details                                      │
│  - 3D rendered crystal structure                                             │
│  - Option to iterate or save                                                 │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Component Details

### 1. MCP Server (TypeScript/Node.js)

**Location:** `/home/niel/git/crystal-mcp-server/src/`

**Responsibilities:**
- Expose tools via MCP protocol
- Handle `tools/list` requests (return tool registry)
- Handle `tools/call` requests (execute tools)
- Validate input schemas
- Format output responses

**Key Files:**
- `server.ts` - Main MCP server entry point
- `types/tools.ts` - Tool definitions and schemas
- `tools/**/*.ts` - Individual tool handlers

### 2. Python Backend

**Location:** `/home/niel/git/crystal-mcp-server/src/python/`

**Responsibilities:**
- Execute actual crystal structure computations
- Use pymatgen, ASE, and other scientific libraries
- Return structured JSON output

**Key Files:**
- `advanced_structures.py` - Prototype generation
- `structure_tools.py` - Slab generation, transformations
- `nanostructure_generator.py` - Nanostructures
- `generators/**/*.py` - Specialized generators

### 3. LLM Client (Rust GUI)

**Location:** `/home/niel/git/crystal-mcp-server/crystal-gui/`

**Responsibilities:**
- Connect to MCP server
- Fetch tool registry dynamically
- Build LLM prompts from tool schemas
- Parse LLM output
- Call MCP tools
- Render results

**Key Files:**
- `src/app.rs` - Main application logic
- `src/mcp_client.rs` - MCP protocol client
- `src/llm_client.rs` - Ollama LLM client
- `src/crystal_viewer.rs` - 3D structure renderer

---

## Tool Registry

The tool registry is the **single source of truth** for all available tools. It is defined in the MCP server and exposed via the `tools/list` method.

### Registry Location

**Primary:** `src/types/tools.ts` - Contains `TOOL_DEFINITIONS` array

### Registry Structure

Each tool in the registry has:

```typescript
{
  name: string,           // Unique tool identifier
  description: string,    // Human-readable description
  inputSchema: {          // JSON Schema for input validation
    type: "object",
    properties: {...},
    required: [...]
  }
}
```

### How the Registry is Used

1. **MCP Server** reads `TOOL_DEFINITIONS` at startup
2. **Client** calls `tools/list` to get the full registry
3. **LLM prompt** is constructed from registry data
4. **Tool calls** are validated against registry schemas

---

## Data Flow

### Input Flow
```
Natural Language → LLM → JSON Tool Call → MCP Server → Python → Computation
```

### Output Flow
```
Computation → Python JSON → MCP Response → GUI Parser → 3D Renderer
```

---

## Key Design Principles

1. **Single Source of Truth**
   - Tool definitions exist ONLY in MCP server
   - GUI fetches definitions dynamically
   - No hardcoded tool knowledge in client

2. **Separation of Concerns**
   - LLM: Natural language understanding, planning
   - MCP Server: Tool routing, schema validation
   - Python: Deterministic computation

3. **Protocol Compliance**
   - Standard MCP protocol
   - Any MCP-compatible client can use this server
   - Not tied to specific LLM or GUI

4. **Structured Output**
   - All tools return consistent JSON format
   - `<json-data>` blocks for machine parsing
   - Markdown text for human reading

---

## File Structure

```
crystal-mcp-server/
├── src/
│   ├── server.ts              # MCP server entry
│   ├── index.ts               # Main export
│   ├── types/
│   │   ├── tools.ts           # TOOL REGISTRY (schemas)
│   │   ├── crystal.ts         # Crystal structure types
│   │   └── errors.ts          # Error types
│   ├── tools/
│   │   ├── generation/        # Structure generation tools
│   │   ├── transformation/    # Structure modification tools
│   │   ├── analysis/          # Analysis tools
│   │   ├── export/            # Export tools
│   │   └── optimization/      # MLFF optimization tools
│   ├── utils/
│   │   └── formatting.ts      # Output formatting
│   └── python/
│       ├── advanced_structures.py
│       ├── structure_tools.py
│       ├── nanostructure_generator.py
│       └── generators/        # Specialized generators
│
├── crystal-gui/
│   └── src/
│       ├── app.rs             # Main GUI application
│       ├── mcp_client.rs      # MCP protocol client
│       ├── llm_client.rs      # Ollama LLM client
│       └── crystal_viewer.rs  # 3D renderer
│
└── docs/
    └── MCP_ARCHITECTURE.md    # This document
```

---

## Tool Categories

### Generation Tools
- `generate_crystal` - Arbitrary space group structures
- `generate_prototype` - Common prototype structures
- `generate_nanostructure` - Nanotubes, graphene, etc.
- `generate_2d_material` - MXene, hBN, MoS2
- `generate_molecular_crystal` - Molecular crystals
- `build_molecule` - Individual molecules

### Transformation Tools
- `generate_slab` - Surface slabs
- `make_supercell` - Supercell expansion
- `create_defect` - Point defects
- `create_alloy` - Alloy structures
- `add_adsorbate` - Surface adsorbates
- `apply_strain` - Lattice strain
- `create_heterostructure` - Layered heterostructures

### Analysis Tools
- `analyze_symmetry` - Space group analysis
- `analyze_symmetry_relations` - Subgroup relations
- `validate_structure` - Structure validation

### Export Tools
- `export_structure` - CIF, VASP, XYZ formats
- `generate_visualization` - 3D visualizations

### Optimization Tools
- `calculate_energy_mlff` - MLFF energy calculation
- `mlff_optimize` - Structure optimization
- `ground_state_search` - Ground state finding

---

## Version

- Document Version: 1.0
- Last Updated: 2026-01-04
- MCP Protocol Version: 2024-11-05
