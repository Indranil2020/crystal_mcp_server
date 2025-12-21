
# Server Capabilities & FAQ

## 1. Supported Operations
The server supports **228 operations** organized into 20 categories.
The primary entry point is the tool `comprehensive_generate`.

### Key Categories
- **Bulk**: `generate_from_spacegroup`, `make_supercell`, `apply_strain`
- **Defects**: `generate_vacancy`, `generate_interstitial`, `create_defect`
- **Surfaces**: `generate_slab`, `add_adsorbate`
- **2D / Low Dim**: `generate_2d_material`, `generate_nanostructure`, `generate_twisted_bilayer`
- **Quantum/Electronic**: `generate_quantum_dot`, `generate_semiconductor`
- **Optimization**: `optimize_structure_mlff`, `calculate_energy_mlff` (if ML libraries installed)

See [API Reference](./api_reference.md) for the complete list.

## 2. Resources and Prompts
**Does it expose resources or prompts?**
- **Resources**: No. The server focuses on *generative tools*. It typically does not expose static resources via `resources/list`.
- **Prompts**: No. It does not provide pre-defined prompts via `prompts/list`.

## 3. Domain & Purpose
**Domain**: Materials Science, Crystallography, and Computational Chemistry.
**Purpose**:
- To enable LLMs to generate valid, physical crystal structures represented in standard formats (pymatgen JSON, CIF, POSCAR).
- To bridge the gap between text (LLM) and physics engines (pymatgen, PyXtal).

## 4. File I/O
**How does it handle files?**
- **Input**: Accepts JSON payloads directly via MCP arguments.
- **Output**: Returns generated content as Strings (Text).
- **Export**: The `export_*` operations (e.g., `export_vasp`) return the *file content* as a string. The server does **not** write files to the user's disk arbitrarily. It is the client's responsibility to save the string if desired.
- **Internal**: Uses `/tmp` for secure IPC between TypeScript and Python layers, which is cleaned up automatically.

## 5. State Management
**Does it maintain state?**
- **Stateless**. Each tool call is independent.
- However, the **Design Pattern** supports "Iterative Refinement" by allowing you to pass the Output JSON of one call as the Input argument of the next. Use the client/Model context to hold the state.

## 6. Expected Load
- **Architecture**: Single-tenant process (Stdio).
- **Concurrency**: The Python bridge processes requests sequentially per server instance.
- **Performance**: Most operations take < 0.1s. Complex generations (e.g., large supercells or MLFF optimization) may take seconds.
