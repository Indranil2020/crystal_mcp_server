# Module Dependencies

> Visual graph of Python imports and dependencies for molecule generation.

---

## 🔗 Import Dependency Graph

```mermaid
graph TB
    subgraph TS["TypeScript Layer"]
        SERVER["server.ts"]
        HANDLER["build-molecule.ts"]
        BRIDGE["python-bridge.ts"]
        SCHEMA["tools.ts"]
    end
    
    subgraph PY_ENTRY["Python Entry Points"]
        MOL_GEN["molecule_generator.py"]
        CLUSTER_GEN["molecular_cluster_generator.py"]
    end
    
    subgraph PY_CORE["Core Generators"]
        UNIVERSAL["universal_molecule.py"]
        CLUSTER["molecular_cluster.py"]
        DB["molecule_database.py"]
        SMALL["small_molecules.py"]
    end
    
    subgraph PY_SPECIAL["Specialized Generators"]
        BIO["biomolecules.py"]
        CONF["conformers.py"]
        CAGE["cages.py"]
        CARBON["carbon_nanostructures.py"]
        FRAME["frameworks.py"]
        ORGANO["organometallics.py"]
        PORPH["porphyrins.py"]
    end
    
    subgraph EXTERNAL["External Libraries"]
        RDKIT["RDKit<br/>(3D generation)"]
        PYMATGEN["Pymatgen<br/>(Molecule class)"]
        NUMPY["NumPy<br/>(coordinates)"]
        SCIPY["SciPy<br/>(rotations)"]
        REQUESTS["Requests<br/>(PubChem API)"]
        ASE["ASE<br/>(g2 database)"]
    end
    
    SERVER --> HANDLER
    HANDLER --> BRIDGE
    HANDLER --> SCHEMA
    BRIDGE --> MOL_GEN
    BRIDGE --> CLUSTER_GEN
    
    MOL_GEN --> UNIVERSAL
    MOL_GEN --> ASE
    
    CLUSTER_GEN --> CLUSTER
    
    CLUSTER --> UNIVERSAL
    CLUSTER --> SCIPY
    CLUSTER --> NUMPY
    CLUSTER --> PYMATGEN
    
    UNIVERSAL --> SMALL
    UNIVERSAL --> DB
    UNIVERSAL --> RDKIT
    UNIVERSAL --> REQUESTS
    UNIVERSAL --> PYMATGEN
    
    DB --> SMALL
    
    UNIVERSAL -.-> BIO
    UNIVERSAL -.-> CONF
    UNIVERSAL -.-> CAGE
    UNIVERSAL -.-> CARBON
    UNIVERSAL -.-> FRAME
    UNIVERSAL -.-> ORGANO
    UNIVERSAL -.-> PORPH
    
    style UNIVERSAL fill:#f9f,stroke:#333,stroke-width:3px
    style CLUSTER fill:#f9f,stroke:#333,stroke-width:3px
    style MOL_GEN fill:#bbf,stroke:#333,stroke-width:2px
```

---

## 📦 External Dependencies

| Library | Usage | Required? |
|---------|-------|-----------|
| `rdkit` | SMILES parsing, 3D conformer generation, force field optimization | ✅ Core |
| `pymatgen` | Molecule/Structure objects, coordinate handling | ✅ Core |
| `numpy` | Coordinate arrays, transformations | ✅ Core |
| `scipy` | 3D rotations (Rotation.from_euler) | ⚠️ Optional (fallback exists) |
| `requests` | PubChem API calls | ⚠️ Optional (offline mode) |
| `ase` | G2 molecule database (fallback) | ⚠️ Optional |

---

## 🔄 Runtime Check Pattern

All modules use `importlib.util.find_spec()` for safe imports:

```python
# Pattern used in all molecule generators
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None
SCIPY_AVAILABLE = importlib.util.find_spec("scipy") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem
```

---

## 📊 Call Stack Visualization

```mermaid
graph LR
    subgraph STACK["Call Stack: build_molecule"]
        A["handleBuildMolecule()"]
        B["buildMolecule()"]
        C["executePythonWithJSON()"]
        D["spawn python"]
        E["main()"]
        F["generate_molecule()"]
        G["generate_molecule_universal()"]
        H["smiles_to_3d_structure()"]
        I["molecule_to_structure_dict()"]
    end
    
    A --> B --> C --> D --> E --> F --> G --> H --> I
    
    style A fill:#4af,stroke:#333
    style E fill:#fa4,stroke:#333
    style G fill:#f4a,stroke:#333
```

```mermaid
graph LR
    subgraph STACK2["Call Stack: build_molecular_cluster"]
        A2["handleBuildMolecularCluster()"]
        B2["executePythonWithJSON()"]
        C2["generate_molecular_cluster()"]
        D2["generate_molecule_universal()"]
        E2["classify_molecule()"]
        F2["arrange_stacked()"]
        G2["combine_molecules()"]
        H2["add_vacuum_box()"]
    end
    
    A2 --> B2 --> C2 --> D2
    C2 --> E2 --> F2 --> G2 --> H2
    
    style C2 fill:#f4a,stroke:#333
    style D2 fill:#fa4,stroke:#333
```
