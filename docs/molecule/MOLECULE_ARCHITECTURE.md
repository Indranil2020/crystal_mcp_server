# Molecule Tools Architecture

> **Quick Reference for Debugging** – All molecule-related files and data flows.

---

## 🗂️ File Map

```
crystal-mcp-server/
│
├── src/
│   ├── server.ts                          ← MCP Server (tool routing)
│   │
│   ├── tools/generation/
│   │   └── build-molecule.ts              ← TypeScript tool handler
│   │
│   ├── types/
│   │   └── tools.ts                       ← Schemas (BuildMoleculeSchema, BuildMolecularClusterSchema)
│   │
│   ├── utils/
│   │   └── python-bridge.ts               ← Python subprocess executor
│   │
│   ├── python/
│   │   ├── molecule_generator.py          ← Python entry point (build_molecule)
│   │   ├── molecular_cluster_generator.py ← Python entry point (build_molecular_cluster)
│   │   │
│   │   └── generators/molecule/           ← Core generation logic
│   │       ├── universal_molecule.py      🔑 Main resolver (130M+ molecules)
│   │       ├── molecular_cluster.py       🔑 Cluster arrangements
│   │       ├── molecule_database.py       ← SQLite database manager
│   │       ├── small_molecules.py         ← Built-in molecule coords
│   │       ├── biomolecules.py            ← Amino acids, DNA bases
│   │       ├── conformers.py              ← Conformer generation
│   │       ├── cages.py                   ← Cage molecules
│   │       ├── carbon_nanostructures.py   ← C60, nanotubes
│   │       ├── frameworks.py              ← MOFs, COFs
│   │       ├── organometallics.py         ← Metal complexes
│   │       └── porphyrins.py              ← Porphyrin rings
│   │
│   └── data/molecule/
│       └── molecules.db                   ← SQLite cache
│
├── crystal-gui/src/                       ← Rust GUI
│   ├── app.rs                             ← Main app (chat, tool calls)
│   ├── llm_client.rs                      ← Ollama integration
│   ├── mcp_client.rs                      ← MCP communication
│   └── crystal_viewer.rs                  ← 3D visualization
│
└── tests/
    └── end_to_end_test.py                 ← E2E testing
```

---

## 🔄 Data Flow: Single Molecule (`build_molecule`)

```mermaid
sequenceDiagram
    participant GUI as Crystal GUI (Rust)
    participant LLM as Ollama LLM
    participant MCP as MCP Server (TS)
    participant PY as Python Backend
    participant DB as molecule_database
    participant EXT as External APIs

    GUI->>LLM: "Generate benzene"
    LLM-->>GUI: {tool: "build_molecule", params: {name: "benzene"}}
    
    GUI->>MCP: call_tool("build_molecule", {name: "benzene"})
    MCP->>MCP: BuildMoleculeSchema.parse()
    MCP->>PY: executePythonWithJSON("molecule_generator.py")
    
    PY->>PY: generate_molecule("benzene")
    
    alt Found in local DB
        PY->>DB: lookup_molecule("benzene")
        DB-->>PY: {smiles: "c1ccccc1", ...}
    else Not in local DB
        PY->>EXT: fetch_from_pubchem("benzene")
        EXT-->>PY: {cid: 241, smiles: "c1ccccc1"}
        PY->>DB: cache molecule
    end
    
    PY->>PY: smiles_to_3d_structure() [RDKit]
    PY->>PY: MMFF94/UFF optimization
    PY->>PY: molecule_to_structure_dict()
    
    PY-->>MCP: {success: true, structure: {...}}
    MCP-->>GUI: JSON response with <json-data>
    
    GUI->>GUI: Parse structure
    GUI->>GUI: CrystalViewer.render()
```

---

## 🔄 Data Flow: Molecular Cluster (`build_molecular_cluster`)

```mermaid
sequenceDiagram
    participant GUI as Crystal GUI
    participant MCP as MCP Server
    participant PY as molecular_cluster.py
    participant UNI as universal_molecule.py

    GUI->>MCP: build_molecular_cluster({molecules: [{identifier: "benzene", count: 2}], stacking: "pi_pi_parallel"})
    
    MCP->>PY: molecular_cluster_generator.py
    
    loop For each molecule
        PY->>UNI: generate_molecule_universal("benzene")
        UNI-->>PY: {atoms: [...], coords: [...]}
    end
    
    PY->>PY: classify_molecule() → is_aromatic?
    PY->>PY: auto_select_stacking() or use specified
    
    alt Stacking Type
        PY->>PY: arrange_stacked() [π-π]
    else
        PY->>PY: arrange_linear()
    else
        PY->>PY: arrange_circular()
    else
        PY->>PY: arrange_t_shaped()
    else
        PY->>PY: arrange_custom()
    end
    
    PY->>PY: combine_molecules()
    PY->>PY: add_vacuum_box()
    
    PY-->>MCP: {success: true, structure: {...}}
    MCP-->>GUI: Response
```

---

## 📊 Molecule Resolution Priority

```mermaid
flowchart TB
    subgraph INPUT["📥 Input"]
        ID[/"name: benzene"/]
    end
    
    subgraph DETECT["🔍 detect_input_type()"]
        D1{Is SMILES?}
        D2{Is InChI?}
        D3{Is CID?}
        D4{Is IUPAC?}
        D5[Assume Name]
    end
    
    subgraph RESOLVE["🎯 Resolution Priority"]
        P1["1️⃣ Local Database<br/>(small_molecules.py)"]
        P2["2️⃣ SQLite Cache<br/>(molecules.db)"]
        P3["3️⃣ Aliases<br/>(MOLECULE_ALIASES)"]
        P4["4️⃣ RDKit SMILES→3D"]
        P5["5️⃣ PubChem API<br/>(130M+ molecules)"]
        P6["6️⃣ OPSIN<br/>(IUPAC names)"]
    end
    
    subgraph OUTPUT["📤 Output"]
        OUT["{atoms, coords, formula}"]
    end
    
    ID --> D1
    D1 -->|Yes| P4
    D1 -->|No| D2
    D2 -->|Yes| P4
    D2 -->|No| D3
    D3 -->|Yes| P5
    D3 -->|No| D4
    D4 -->|Yes| P6
    D4 -->|No| D5
    D5 --> P1
    
    P1 -->|Found| OUT
    P1 -->|Not Found| P2
    P2 -->|Found| OUT
    P2 -->|Not Found| P3
    P3 -->|Found| P4
    P3 -->|Not Found| P5
    P4 --> OUT
    P5 -->|Found| OUT
    P5 -->|Not Found| P6
    P6 --> OUT
```

---

## 🏗️ Cluster Stacking Types

```mermaid
graph LR
    subgraph STACKING["Stacking Arrangements"]
        direction TB
        A["π-π Parallel<br/>3.4 Å"]
        B["π-π Antiparallel<br/>180° rotated"]
        C["π-π Offset<br/>slip-stacked"]
        D["T-Shaped<br/>edge-to-face"]
        E["Herringbone<br/>60° tilt"]
        F["H-Bonded<br/>2.8 Å"]
        G["Linear<br/>along axis"]
        H["Circular<br/>ring"]
        I["Swastika<br/>4-fold cross"]
    end
    
    AUTO["auto_select_stacking()"] --> |Aromatic| A
    AUTO --> |H-bond capable| F
    AUTO --> |Mixed| D
    AUTO --> |Default| VDW["Van der Waals<br/>3.5 Å"]
```

---

## 🔧 Key Functions Reference

| File | Function | Purpose |
|------|----------|---------|
| `molecule_generator.py` | `generate_molecule()` | Entry point for single molecules |
| `universal_molecule.py` | `generate_molecule_universal()` | Main resolver (all sources) |
| `universal_molecule.py` | `detect_input_type()` | Auto-detect SMILES/name/CID |
| `universal_molecule.py` | `smiles_to_3d_structure()` | RDKit 3D generation |
| `universal_molecule.py` | `fetch_from_pubchem()` | External API lookup |
| `molecular_cluster.py` | `generate_molecular_cluster()` | Cluster assembly |
| `molecular_cluster.py` | `arrange_stacked()` | π-stacking arrangement |
| `molecular_cluster.py` | `rotate_molecule()` | 3D rotation utility |
| `molecule_database.py` | `lookup_molecule()` | SQLite cache lookup |

---

## 🐛 Debug Checkpoints

```
┌─────────────────────────────────────────────────────────────────┐
│  1. GUI → LLM                                                   │
│     Check: app.rs:send_chat_message() line 133                  │
│     Debug: [DEBUG] LLM REQUEST printed to stderr                │
├─────────────────────────────────────────────────────────────────┤
│  2. LLM → Tool Parse                                            │
│     Check: app.rs:parse_tool_call() line 266                    │
│     Issue: LLM returns text instead of tool_calls array         │
├─────────────────────────────────────────────────────────────────┤
│  3. Tool → MCP                                                  │
│     Check: app.rs:call_tool() line 322                          │
│     Check: mcp_client.rs:call_tool()                            │
├─────────────────────────────────────────────────────────────────┤
│  4. MCP → Python                                                │
│     Check: server.ts case "build_molecule" line 157             │
│     Check: python-bridge.ts:executePythonWithJSON()             │
├─────────────────────────────────────────────────────────────────┤
│  5. Python Resolution                                           │
│     Check: molecule_generator.py:generate_molecule()            │
│     Check: universal_molecule.py:generate_molecule_universal()  │
│     Debug: Add logging.info() statements                        │
├─────────────────────────────────────────────────────────────────┤
│  6. Response → GUI                                              │
│     Check: <json-data> tag in response                          │
│     Check: crystal_viewer.rs structure parsing                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 📁 Quick Links

| What | Where |
|------|-------|
| MCP Tool Schema | [tools.ts#L626](file:///home/niel/git/crystal-mcp-server/src/types/tools.ts#L626) |
| TS Handler | [build-molecule.ts](file:///home/niel/git/crystal-mcp-server/src/tools/generation/build-molecule.ts) |
| Python Entry | [molecule_generator.py](file:///home/niel/git/crystal-mcp-server/src/python/molecule_generator.py) |
| Universal Resolver | [universal_molecule.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/universal_molecule.py) |
| Cluster Generator | [molecular_cluster.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/molecular_cluster.py) |
| GUI App | [app.rs](file:///home/niel/git/crystal-mcp-server/crystal-gui/src/app.rs) |
| LLM Client | [llm_client.rs](file:///home/niel/git/crystal-mcp-server/crystal-gui/src/llm_client.rs) |
