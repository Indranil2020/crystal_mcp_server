# Crystal Structure Generator MCP Server

A comprehensive Model Context Protocol (MCP) server for generating accurate crystal structures for DFT, quantum chemistry, and condensed matter physics calculations.

## Features

- ✨ **Complete Coverage**: All 230 space groups supported
- 🔬 **PyXtal Integration**: High-quality structure generation
- ⚡ **MLFF Optimization**: CHGNet, M3GNet, and MACE support
- 🎯 **Ground State Search**: Find lowest energy structures
- 📐 **Structure Transformations**: Supercells, slabs, surfaces
- 🔷 **Symmetry Analysis**: Spglib-powered symmetry detection
- 📁 **Multiple Formats**: CIF, POSCAR, XYZ, JSON export
- ✅ **Comprehensive Validation**: Distance checks, lattice validation

## Installation

### Prerequisites

- Node.js 18+ and npm
- Python 3.8+ and pip
- Git

### Quick Install

```bash
# Clone repository
git clone <repository-url>
cd crystal-mcp-server

# Install Node.js dependencies
npm install

# Install Python dependencies
pip install -r requirements.txt

# Build TypeScript
npm run build

# Run tests
npm test
```

### Detailed Installation

#### 1. Install Node.js Dependencies

```bash
npm install
```

This installs:
- @modelcontextprotocol/sdk
- zod (schema validation)
- python-shell (Python bridge)

#### 2. Install Python Dependencies

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install packages
pip install -r requirements.txt
```

Required Python packages:
- `pyxtal>=1.0.0` - Crystal structure generation
- `pymatgen>=2024.1.1` - Materials analysis
- `spglib>=2.5.0` - Symmetry operations
- `ase>=3.22.0` - Atomic simulation environment
- `chgnet>=0.3.0` - CHGNet MLFF model
- `matgl>=1.0.0` - M3GNet MLFF model
- `mace-torch>=0.3.0` - MACE MLFF model

#### 3. Build and Test

```bash
# Compile TypeScript
npm run build

# Run tests
npm test

# Lint code
npm run lint
```

## Usage

### Starting the Server

```bash
# Start with npm
npm start

# Or run directly
node dist/index.js
```

The server runs on stdio transport for MCP protocol communication.

### Available Tools

#### 1. generate_crystal

Generate a crystal structure with specified composition and space group.

```typescript
{
  composition: ["Si", "Si"],
  space_group: 227,  // Diamond structure
  seed: 42
}
```

#### 2. generate_space_group_scan

Scan multiple space groups to find stable structures.

```typescript
{
  composition: ["Na", "Cl"],
  space_groups: [225, 216, 221],
  parallel: true
}
```

#### 3. make_supercell

Create supercell from existing structure.

```typescript
{
  structure_dict: {...},  // Structure object from generate_crystal
  scaling_matrix: [2, 2, 2]  // 2x2x2 supercell, or [[2,0,0],[0,2,0],[0,0,2]]
}
```

#### 4. generate_slab

Generate surface slab for DFT calculations.

```typescript
{
  structure_dict: {...},  // Structure object from generate_crystal
  miller_indices: [1, 0, 0],
  thickness: 5,  // Number of atomic layers
  vacuum: 15.0   // Vacuum padding in Angstroms
}
```

#### 5. analyze_symmetry

Detect and analyze crystal symmetry.

```typescript
{
  structure_dict: {...},  // Structure object
  symprec: 0.001,         // Symmetry precision tolerance
  detect_primitive: true  // Find primitive cell
}
```

#### 6. validate_structure

Validate crystal structure quality.

```typescript
{
  structure_dict: {...},  // Structure object
  checks: ["distances", "lattice", "stoichiometry"],  // Available checks
  min_distance: 1.5  // Minimum interatomic distance in Angstroms
}
```

#### 7. optimize_structure_mlff

Optimize structure using machine learning force fields (optional - requires MLFF packages).

```typescript
{
  structure_dict: {...},  // Structure object
  mlff_model: "chgnet",   // Options: "chgnet", "m3gnet", "mace"
  optimizer: "BFGS",      // Options: "BFGS", "FIRE", "LBFGS"
  fmax: 0.01,             // Force convergence criterion (eV/A)
  steps: 500              // Maximum optimization steps
}
```

#### 8. ground_state_search

Find ground state across space groups.

```typescript
{
  composition: ["Si", "Si"],
  space_groups: [227, 141, 194],
  mlff_model: "chgnet",
  num_structures_per_group: 5
}
```

#### 9. export_structure

Export structure to multiple formats.

```typescript
{
  structure_dict: {...},  // Structure object
  formats: ["cif", "poscar", "xyz", "json"],  // Output formats
  output_directory: "./output"  // Optional output path
}
```

## Examples

### Example 1: Generate Diamond Si

```typescript
const result = await mcpClient.callTool("generate_crystal", {
  composition: ["Si", "Si"],
  space_group: 227,
  seed: 42
});
```

### Example 2: Create Surface Slab

```typescript
// First generate bulk structure
const bulk = await mcpClient.callTool("generate_crystal", {
  composition: ["Ti", "O", "O"],
  space_group: 136  // Rutile TiO2
});

// Then create (110) surface
const slab = await mcpClient.callTool("generate_slab", {
  structure_dict: bulk.structure,
  miller_indices: [1, 1, 0],
  thickness: 5,
  vacuum: 15.0,
  fix_bottom_layers: 2
});
```

### Example 3: Ground State Search

```typescript
const search = await mcpClient.callTool("ground_state_search", {
  composition: ["Si", "Si"],
  space_groups: [227, 141, 194],  // Diamond, I41/amd, P63/mmc
  num_structures_per_group: 10,
  mlff_model: "chgnet",
  optimization_settings: {
    optimizer: "BFGS",
    fmax: 0.01,
    steps: 500
  }
});
```

## Architecture

```
crystal-mcp-server/
├── src/
│   ├── types/          # TypeScript type definitions
│   ├── utils/          # Utility functions
│   ├── python/         # Python backend scripts
│   ├── tools/          # MCP tool implementations
│   │   ├── generation/
│   │   ├── transformation/
│   │   ├── analysis/
│   │   ├── optimization/
│   │   └── export/
│   ├── server.ts       # MCP server setup
│   └── index.ts        # Entry point
├── tests/              # Test suites
├── docs/               # Documentation
└── package.json
```

## Code Quality Standards

This project follows strict defensive programming practices:

- ✅ **No try/catch blocks** - Use Result<T> pattern
- ✅ **100% type hints** - All functions fully typed
- ✅ **Defensive validation** - Check inputs before processing
- ✅ **Comprehensive tests** - Unit and integration tests
- ✅ **Clear error messages** - Actionable suggestions included
- ✅ **PEP 8 compliant** - Python code follows standards
- ✅ **ESLint compliant** - TypeScript code follows standards

## Performance

### Benchmarks

| Operation | Time | Notes |
|-----------|------|-------|
| Generate crystal | <1s | Single structure |
| Space group scan (230) | ~10 min | Parallel mode |
| MLFF optimization | ~30s | Per structure |
| Ground state search | ~4 hours | 230 groups, 5 attempts each |

### Memory Requirements

- Minimal: 8 GB RAM, 4 cores
- Recommended: 16 GB RAM, 8 cores
- Large-scale: 32+ GB RAM, GPU support

## Troubleshooting

### Python Import Errors

```bash
# Verify Python packages
python -c "import pyxtal; print(pyxtal.__version__)"
python -c "import chgnet; print('CHGNet OK')"
python -c "import spglib; print('Spglib OK')"
```

### MLFF Model Loading Fails

```bash
# Download CHGNet model
python -c "from chgnet.model import CHGNet; CHGNet.load()"

# Check GPU availability (optional)
python -c "import torch; print(torch.cuda.is_available())"
```

### Permission Errors

```bash
# Make Python scripts executable
chmod +x src/python/*.py
```

## Contributing

Contributions welcome! Please ensure:

1. All tests pass: `npm test`
2. Code is linted: `npm run lint`
3. Follow defensive programming guidelines
4. Add tests for new features
5. Update documentation

## License

MIT License - see LICENSE file for details

## Citation

If you use this server in your research, please cite:

```bibtex
@software{crystal_mcp_server,
  title = {Crystal Structure Generator MCP Server},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourname/crystal-mcp-server}
}
```

## Acknowledgments

Built with:
- [PyXtal](https://github.com/qzhu2017/PyXtal) - Crystal structure generation
- [Pymatgen](https://pymatgen.org/) - Materials analysis
- [Spglib](https://spglib.github.io/spglib/) - Symmetry operations
- [CHGNet](https://github.com/CederGroupHub/chgnet) - MLFF model
- [MCP SDK](https://github.com/modelcontextprotocol/sdk) - Protocol implementation

## Support

- 📧 Email: support@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/yourname/crystal-mcp-server/issues)
- 📖 Docs: [Full Documentation](./docs/)

---

**Status:** Production Ready ✅

**Version:** 1.0.0

**Last Updated:** December 2024
