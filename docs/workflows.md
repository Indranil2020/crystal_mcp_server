
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
