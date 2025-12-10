# Setting Up Crystal MCP Server with Claude Code

## Overview

This guide shows you how to configure the Crystal Structure Generator MCP Server to work with Claude Code, allowing you to generate crystal structures using natural language.

## Prerequisites

✅ Server is built and tested (npm run build completed)
✅ Python environment is set up with all dependencies
✅ Server starts successfully with `npm start`

## Configuration Steps

### 1. Locate Your Claude Code Settings

The Claude Code settings file is located at:
- **Linux/macOS**: `~/.claude/settings.json`
- **Windows**: `%USERPROFILE%\.claude\settings.json`

### 2. Add MCP Server Configuration

Add the following to your `settings.json` under the `mcpServers` section:

```json
{
  "mcpServers": {
    "crystal-structure-generator": {
      "command": "node",
      "args": [
        "/home/niel/git/crystal-mcp-server/dist/index.js"
      ],
      "env": {
        "PYTHON_PATH": "python3"
      }
    }
  }
}
```

**Important**: Update the path `/home/niel/git/crystal-mcp-server/dist/index.js` to match your actual installation directory.

### 3. Restart Claude Code

After saving the settings, restart Claude Code for the changes to take effect.

### 4. Verify Connection

Once Claude Code restarts, the Crystal Structure Generator MCP Server will be available. You can verify it's working by asking Claude:

```
"List the available MCP tools"
```

You should see 10 tools from the crystal-structure-generator server.

## Available Tools

Once configured, you can use these tools through natural language:

### 1. **generate_crystal**
Generate crystal structures with specific space groups.

**Example prompts:**
- "Generate a diamond silicon structure (space group 227)"
- "Create a rock salt NaCl structure"
- "Generate a rutile TiO2 crystal"

### 2. **generate_space_group_scan**
Scan multiple space groups to find stable structures.

**Example prompts:**
- "Scan space groups 225-230 for a NaCl composition"
- "Find all possible structures for silicon in cubic space groups"

### 3. **make_supercell**
Create supercells from existing structures.

**Example prompts:**
- "Create a 2x2x2 supercell of the silicon structure"
- "Make a 3x3x1 supercell for surface calculations"

### 4. **generate_slab**
Generate surface slabs for DFT calculations.

**Example prompts:**
- "Create a (100) surface slab of silicon with 5 layers and 15Å vacuum"
- "Generate a (110) TiO2 surface with 4 layers"

### 5. **analyze_symmetry**
Analyze crystal symmetry.

**Example prompts:**
- "Analyze the symmetry of this structure"
- "Detect the primitive cell of this crystal"

### 6. **validate_structure**
Validate crystal structure quality.

**Example prompts:**
- "Validate this structure for interatomic distances"
- "Check if this crystal structure is physically reasonable"

### 7. **optimize_structure_mlff**
Optimize structures using machine learning force fields.

**Example prompts:**
- "Optimize this structure using CHGNet"
- "Relax the atomic positions with MACE force field"

### 8. **ground_state_search**
Find ground state structures across multiple space groups.

**Example prompts:**
- "Search for the ground state structure of silicon"
- "Find the lowest energy configuration for this composition"

### 9. **export_structure**
Export structures to various formats.

**Example prompts:**
- "Export this structure to CIF and POSCAR formats"
- "Save the structure as XYZ and JSON files"

## Usage Examples

Once configured, you can use natural language like:

```
"Generate a diamond silicon structure and then create a (111) surface slab
with 5 layers and 15 Angstrom vacuum. Optimize it using CHGNet and export
to POSCAR format."
```

Claude will automatically:
1. Use `generate_crystal` to create the Si diamond structure
2. Use `generate_slab` to create the (111) surface
3. Use `optimize_structure_mlff` with CHGNet
4. Use `export_structure` to save as POSCAR

## Troubleshooting

### Server Not Appearing

1. Check that the path in settings.json is correct
2. Verify the server builds successfully: `npm run build`
3. Test the server manually: `npm start`
4. Check Claude Code logs for errors

### Python Errors

1. Ensure Python 3.8+ is installed
2. Verify all dependencies are installed: `pip install -r requirements.txt`
3. Check that python3 is in your PATH

### Permission Errors

1. Ensure the dist/index.js file is readable
2. Check that Python scripts in src/python/ are accessible

## Advanced Configuration

### Using a Virtual Environment

If you use a Python virtual environment:

```json
{
  "mcpServers": {
    "crystal-structure-generator": {
      "command": "node",
      "args": [
        "/home/niel/git/crystal-mcp-server/dist/index.js"
      ],
      "env": {
        "PYTHON_PATH": "/path/to/venv/bin/python"
      }
    }
  }
}
```

### Custom Python Path

To use a specific Python installation:

```json
{
  "mcpServers": {
    "crystal-structure-generator": {
      "command": "node",
      "args": [
        "/home/niel/git/crystal-mcp-server/dist/index.js"
      ],
      "env": {
        "PYTHON_PATH": "/usr/local/bin/python3.11"
      }
    }
  }
}
```

## Next Steps

After configuration, try these workflows:

1. **Materials Discovery**: Generate structures across space groups and find ground states
2. **Surface Science**: Create and optimize surface slabs for catalysis studies
3. **DFT Preparation**: Generate, validate, and export structures for quantum calculations
4. **Structure Analysis**: Analyze symmetry and validate existing structures

## Support

If you encounter issues:

1. Check the server logs: The server prints diagnostic information to stderr
2. Verify Python backend: Run `python3 tests/test_python_backend.py`
3. Test individual tools: See examples in README.md
4. Check TypeScript build: Run `npm run build` and fix any errors

## Features Available Through Claude

With this MCP server, you can ask Claude to:

- Generate any crystal structure from the 230 space groups
- Create supercells and surface slabs
- Optimize structures with ML force fields (CHGNet, M3GNet, MACE)
- Search for ground state configurations
- Validate structures for DFT calculations
- Export to CIF, POSCAR, XYZ, and other formats
- Analyze crystal symmetry
- Create complex workflows combining multiple operations

All through natural language! No need to remember exact tool names or parameter formats.
