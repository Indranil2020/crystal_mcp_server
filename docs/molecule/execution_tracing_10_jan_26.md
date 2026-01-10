# Molecular Execution Traces & Dependency Analysis
**Date:** 10 Jan 2026
**Status:** Comprehensive Path Analysis (Visual + Forensic Text)

## 1. System Routing Logic (The "Brain")
The `arrangement_adapter.py` acts as the Facade, implementing a bifurcated routing strategy that preserves legacy stability while enabling new capabilities.

```mermaid
flowchart TD
    Start([User Request]) --> Adapter[arrangement_adapter.py]
    Adapter --> Decision{Constraints / NL / Formulas?}
    
    Decision -- No --> LegacyRoute
    Decision -- Yes --> NewRoute
    
    subgraph LegacyRoute [Legacy Path]
        direction TB
        L_Call[_legacy_generate_cluster]
        L_File[molecular_cluster.py]
        L_Call --> L_File
    end
    
    subgraph NewRoute [New Advanced Path]
        direction TB
        N_Call[_arrange_with_new_engine]
        N_File[molecular_arrangement.py]
        N_Solver[ConstraintSolver]
        N_Call --> N_File
        N_File --> N_Solver
    end
    
    style LegacyRoute fill:#f9f,stroke:#333
    style NewRoute fill:#bbf,stroke:#333
```

**Key Logic**:
- If the request contains **Natural Language**, **Constraints**, or **Formulas**, `_needs_new_engine()` returns `True`, activating the `NewRoute`.
- If the request defines simple parameters (e.g. `stacking="custom"`) without advanced features, it defaults to the `LegacyRoute` (using `molecular_cluster.py`).

---

## 2. Trace: Polycatenane (Level 1)
**Status:** ❌ FAILED (Legacy Fallback + Name Error)

```mermaid
sequenceDiagram
    autonumber
    participant GUI as Frontend GUI
    participant Adapter as Adapter
    participant Legacy as molecular_cluster.py
    participant UnivMol as universal_molecule.py

    Note over GUI: Prompt: "3 Interlocked Macrocycles..."
    GUI->>Adapter: generate_cluster(stacking="custom")
    Note right of GUI: Missing 'constraints' param
    
    Adapter->>Adapter: _needs_new_engine() -> False
    Adapter->>Legacy: _legacy_generate_cluster()
    
    Legacy->>UnivMol: generate("1,4-phenylene-2,3-dicarboximide")
    activate UnivMol
    UnivMol--xLegacy: ERROR: Name Resolution Failed
    deactivate UnivMol
    
    Legacy-->>Adapter: Return Failure
    Adapter-->>GUI: Error Message
```

### Forensic Analysis
1.  **Frontend Failure**: The LLM parsed the prompt but failed to translate the complex topological requirement ("interlocked") into mathematical `constraints`. It passed `stacking="custom"` without any coordinate formulas.
2.  **Routing**: Because no advanced parameters were detected, `arrangement_adapter` fell back to the **Legacy Route**.
3.  **Backend Failure**: The legacy engine attempted to generate the base molecule. The name `"1,4-phenylene-2,3-dicarboximide"` was not recognized by `universal_molecule` (which relies on internal catalogs or basic RDKit generic names).
4.  **Dependencies Involved**:
    - `molecular_cluster.py` (Active)
    - `universal_molecule.py` (Active)

---

## 3. Trace: HOF Diamondoid (Level 2)
**Status:** ❌ FAILED (Legacy Fallback + Data Missing)

```mermaid
sequenceDiagram
    autonumber
    participant GUI as Frontend GUI
    participant Adapter as Adapter
    participant Legacy as molecular_cluster.py

    Note over GUI: Prompt: "HOF with Guanidinium & Sulfonate..."
    
    Note right of GUI: PARSING ERROR: Dropped 'Sulfonate'
    GUI->>Adapter: generate_cluster(molecules=[Guanidinium], stacking="auto")
    
    Adapter->>Adapter: _needs_new_engine() -> False
    Adapter->>Legacy: _legacy_generate_cluster()
    
    Legacy->>Legacy: Auto-detect Pattern -> "linear"
    Legacy->>Legacy: Place 18 Guanidiniums
    
    Legacy-->>Adapter: Return Structure (Partial)
    Adapter-->>GUI: Render (Blob at Origin)
```

### Forensic Analysis
1.  **Frontend Failure**: The LLM completely **dropped** the second component ("18 sulfonate anions") from the JSON payload effectively silencing half the prompt.
2.  **Routing**: The prompt was passed with `stacking="auto"` and no constraints (the "diamondoid network" requirement was lost). This triggered the **Legacy Route**.
3.  **Backend Behavior**: The legacy engine received a single list of atoms. It auto-detected a "linear" default pattern (or clump) and placed them.
4.  **Dependencies Involved**:
    - `molecular_cluster.py` (Active)

---

## 4. Trace: Lipid Raft (Level 5)
**Status:** ❌ TIMEOUT (New Engine Path)

```mermaid
sequenceDiagram
    autonumber
    participant GUI as Frontend GUI
    participant Adapter as Adapter
    participant Engine as molecular_arrangement.py
    participant Solver as ConstraintSolver

    Note over GUI: Prompt: "Lipid Raft... Bilayer..."
    GUI->>Adapter: generate_cluster(natural_language="bilayer...")
    
    Note right of Adapter: Log: use_new_engine: True
    Adapter->>Engine: arrange_molecules()
    
    rect rgb(240, 255, 240)
    Engine->>Engine: Resolve POPC (Success)
    end
    rect rgb(255, 220, 220)
    Engine->>Engine: Resolve Sphingomyelin (FAIL)
    end
    
    Engine->>Solver: solve(70 molecules)
    activate Solver
    Note over Solver: Optimizing ~7000 Atoms (O(N^2))
    
    Solver--xAdapter: TIMEOUT (>60s)
    deactivate Solver
    
    Adapter-->>GUI: Error: Execution Timed Out
```

### Forensic Analysis
1.  **Frontend Success**: The LLM successfully passed the raw prompt in the `natural_language` field.
2.  **Routing**: `arrangement_adapter` detected the `natural_language` parameter and correctly routed to the **New Advanced Route**.
3.  **Backend Performance Failure**: The system attempted to optimize the arrangement of 70 large biomolecules (POPC is ~130 atoms).
    - The `ConstraintSolver` performs pairwise clash detection.
    - O(N_atoms²) for ~7000 atoms is computationally expensive in Python.
    - The operation exceeded the 60-second HTTP timeout.
4.  **Dependencies Involved**:
    - `molecular_arrangement.py` (Active)
    - `ConstraintSolver` Class

---

## 5. Dependency Map
Files currently participating in the Unified System:

| Component | Status | File Path |
|-----------|--------|-----------|
| **Facade** | Active | `src/python/generators/molecule/arrangement_adapter.py` |
| **New Engine** | Active | `src/python/generators/molecule/molecular_arrangement.py` |
| **Legacy Engine** | **Active (Fallback)** | `src/python/generators/molecule/molecular_cluster.py` |
| **Molecule Source** | Active | `src/python/generators/molecule/universal_molecule.py` |
