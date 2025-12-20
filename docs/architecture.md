
# System Architecture

## Overview
The Crystal MCP Server is a hybrid TypeScript/Python system designed to provide comprehensive crystal structure generation capabilities to LLMs via the Model Context Protocol (MCP).

```mermaid
graph TD
    User[LLM / Client] -->|JSON-RPC| TS[TypeScript MCP Server]
    TS -->|Bridge Protocol| Py[Python Implementation]
    
    subgraph "TypeScript Layer"
        SRV[server.ts]
        HND[comprehensive-generator.ts]
        BRG[python-bridge.ts]
        
        SRV --> HND
        HND -->|Exec| BRG
    end
    
    subgraph "Python Layer"
        RTR[comprehensive_structures.py]
        REG[generators/__init__.py]
        MODS[Generator Modules]
        
        BRG -->|Stdin/File| RTR
        RTR -->|Lookup| REG
        RTR -->|Dynamic Import| MODS
        
        MODS -->|pymatgen| CORE[Pymatgen Core]
        MODS -->|pyxtal| XTAL[PyXtal]
    end
```

## Unified Router Logic
The heart of the system is the **Unified Router** (`comprehensive_structures.py`), which enables a single MCP tool to access 228+ operations.

```mermaid
flowchart TD
    A[Input JSON] --> B{Has Structure?}
    B -- Yes --> C[Auto-Detect & Convert to Pymatgen Object]
    B -- No --> D[Check 'category' Arg]
    
    C --> D
    D -- Provided --> E[Prioritize Category Search]
    D -- Empty --> F[Standard Search]
    
    E --> G{Find Operation}
    F --> G
    
    G -- Found --> H[Dynamic Import Module]
    G -- Not Found --> I[Return Error + Hints]
    
    H --> J[Execute Function]
    J --> K[Return Result JSON]
```

## Directory Structure
The codebase is organized to separate concerns while maintaining a flat, discoverable registry for generators.

```mermaid
graph LR
    Root --> src
    src --> types[TypeScript Types]
    src --> tools[Tool Handlers]
    src --> python[Python Backend]
    
    python --> generators
    python --> scripts
    
    generators --> bulk
    generators --> surface
    generators --> defect
    generators --> twist
    generators --> output_formats
    generators --> ...[15+ More]
    
    scripts --> update_registry[Registry Scanner]
    scripts --> generate_docs[Doc Generator]
```
