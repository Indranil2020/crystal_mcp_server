# Crystal MCP Server - Backend API Documentation

## Overview
The Crystal MCP Server provides a suite of tools for molecular and crystal generation, analysis, and optimization. It exposes these tools via the Model Context Protocol (MCP) to the AI agent.

## Core Services

### 1. Molecule Generator (`src/python/molecule_generator.py`)
Primary service for generating molecular structures.

**Key Functions:**
- `generate_molecule(name: str, use_cache: bool = True) -> str`
  - Generates a molecule by name.
  - **Strategy:**
    1.  Check **LRU Cache** (256 entries).
    2.  Check internal `MOLECULE_DATABASE` (common molecules).
    3.  Try `rdkit` generation (if valid SMILES/name).
    4.  Try `pubchempy` lookup (requires internet).
    5.  Try `ase.build.molecule` fallback.
    6.  Try `opsin` (for IUPAC names).
  - **Returns:** JSON string (structure data).

- `get_cache_stats() -> dict`
  - Returns debug stats: `{ "hits": int, "misses": int, "size": int }`.

- `clear_cache() -> None`
  - Clears the LRU cache.

### 2. Dependency Management (`src/python/dependency_check.py`)
Handles optional dependencies to ensure graceful degradation.

**Key Functions:**
- `get_capabilities() -> dict`
  - Returns a dictionary of available libraries (`has_rdkit`, `has_ase`, `has_pymatgen`, etc.).
- `check_dependencies() -> None`
  - Prints a human-readable status report to stdout.

### 3. Crystal Generator (`src/python/crystal_generator.py`)
Service for generating crystal structures.

**Key Functions:**
- `generate_crystal(compound: str, lattice: str, ...) -> str`
  - Generates crystal structures based on parameters.
  - Supports `mp_api` (Materials Project), `ase.spacegroup`, and `pyxtal`.

### 4. Geometry Optimization (`src/python/optimization.py`)
Performs partial optimization of molecular geometries.

**Key Functions:**
- `optimize_geometry_rdkit(mol_str: str) -> str`
  - Uses RDKit's MMFF or UFF forcefields to minimize energy and correct bond lengths/angles.
  - Critical for "2D to 3D" workflow from Kekule editor.

## LLM Tools (MCP)
These tools are exposed to the AI agent:

| Tool Name | Description | Parameters |
|-----------|-------------|------------|
| `generate_molecule` | Create a molecule | `name` (str) |
| `generate_crystal_structure` | Create a crystal | `element` (str), `lattice` (str) |
| `optimize_structure` | Relax geometry | `structure_id` (str) |
| `get_capabilities` | Check installed libs | None |

## Architecture
- **Polyglot:** Node.js (MCP Server) <-> Python (Scientific Stack).
- **Communication:** Stdio pipes using JSON-RPC.
- **Data Format:** Structures are passed as JSON objects containing atomic numbers (`atoms`), coordinates (`coords`), and optional unit cell info (`cell`).
