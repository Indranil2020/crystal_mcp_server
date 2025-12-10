# Crystal MCP Server - Quick Installation & Usage Guide

## 🚀 Quick Start (5 Minutes)

### Step 1: Prerequisites Check

```bash
# Check Node.js (need 18+)
node --version

# Check Python (need 3.8+)
python3 --version

# Check pip
pip3 --version
```

### Step 2: Extract and Setup

```bash
# Extract the project
cd /path/to/extracted/crystal-mcp-server

# Install Node.js dependencies
npm install

# Install Python dependencies
pip3 install -r requirements.txt
```

### Step 3: Build

```bash
# Compile TypeScript to JavaScript
npm run build
```

### Step 4: Test

```bash
# Test Python backend
python3 tests/test_python_backend.py
```

### Step 5: Run

```bash
# Start the MCP server
npm start
```

---

## 📦 Python Dependencies Installation

### Core Requirements (Essential)

```bash
pip3 install pyxtal>=1.0.0 pymatgen>=2024.1.1 spglib>=2.5.0 ase>=3.22.0 numpy>=1.24.0 scipy>=1.10.0
```

### MLFF Models (Highly Recommended)

```bash
# CHGNet (fast and accurate)
pip3 install chgnet>=0.3.0

# M3GNet (general purpose)
pip3 install matgl>=1.0.0

# MACE (highest accuracy)
pip3 install mace-torch>=0.3.0
```

### Using Virtual Environment (Recommended)

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate  # Linux/Mac
# OR
venv\Scripts\activate  # Windows

# Install everything
pip3 install -r requirements.txt
```

---

## 🧪 Testing the Installation

### Test 1: Python Backend

```bash
python3 tests/test_python_backend.py
```

**Expected Output:**
```
============================================================
CRYSTAL MCP SERVER - PYTHON BACKEND TEST SUITE
============================================================

=== Testing Crystal Generator ===
✓ Valid space group 227
✓ Rejects invalid space group 999
✓ Generated structure with 8 atoms
  Formula: Si8
  Volume: 160.18 ų
  Density: 2.33 g/cm³

✓✓✓ Crystal Generator: ALL TESTS PASSED ✓✓✓

... (more tests) ...

TEST SUMMARY
============================================================
crystal_generator        : ✓ PASSED
symmetry_analyzer        : ✓ PASSED
mlff_calculator          : ✓ PASSED
validators               : ✓ PASSED
structure_tools          : ✓ PASSED
============================================================

🎉 ALL TESTS PASSED! 🎉
```

### Test 2: TypeScript Build

```bash
npm run build
```

**Expected Output:**
```
> crystal-mcp-server@1.0.0 build
> tsc

# No errors = success!
```

---

## 💡 Usage Examples

### Example 1: Generate Diamond Silicon

**Input:**
```json
{
  "composition": ["Si", "Si"],
  "space_group": 227,
  "num_atoms": 8,
  "seed": 42
}
```

**Command (if testing Python directly):**
```bash
echo '{"composition": ["Si", "Si"], "space_group": 227, "seed": 42}' > /tmp/input.json
python3 src/python/crystal_generator.py /tmp/input.json
```

**Output:**
```json
{
  "success": true,
  "structure": {
    "lattice": {
      "a": 5.43,
      "b": 5.43,
      "c": 5.43,
      "alpha": 90.0,
      "beta": 90.0,
      "gamma": 90.0,
      "volume": 160.18
    },
    "space_group": {
      "number": 227,
      "symbol": "Fd-3m",
      "crystal_system": "cubic"
    },
    "metadata": {
      "formula": "Si8",
      "natoms": 8,
      "density": 2.33
    },
    "atoms": [...]
  }
}
```

### Example 2: Analyze Symmetry

**Input:**
```json
{
  "structure_dict": {...},
  "symprec": 0.001,
  "detect_primitive": true
}
```

**Command:**
```bash
python3 src/python/symmetry_analyzer.py input.json
```

### Example 3: Optimize with MLFF

**Input:**
```json
{
  "structure_dict": {...},
  "mlff_model": "chgnet",
  "optimizer": "BFGS",
  "fmax": 0.01,
  "steps": 500
}
```

**Command:**
```bash
python3 src/python/mlff_calculator.py input.json
```

---

## 🔧 Troubleshooting

### Issue 1: Python Packages Not Found

**Error:**
```
ModuleNotFoundError: No module named 'pyxtal'
```

**Solution:**
```bash
# Make sure you're in the correct environment
source venv/bin/activate  # if using venv

# Reinstall packages
pip3 install -r requirements.txt
```

### Issue 2: CHGNet Model Not Loading

**Error:**
```
Model load failed: CHGNet model not downloaded
```

**Solution:**
```bash
# Download CHGNet model (first time only)
python3 -c "from chgnet.model import CHGNet; CHGNet.load()"
```

This will download ~100MB model files.

### Issue 3: TypeScript Compilation Errors

**Error:**
```
Cannot find module '@modelcontextprotocol/sdk'
```

**Solution:**
```bash
# Reinstall Node.js dependencies
rm -rf node_modules package-lock.json
npm install
```

### Issue 4: Structure Generation Fails

**Error:**
```
Generation failed after 100 attempts
```

**Solutions:**
```json
// 1. Increase volume_factor
{"volume_factor": 1.5}

// 2. Relax min_distance constraints
{"min_distance": {"Si-Si": 2.0}}

// 3. Increase max_attempts
{"max_attempts": 500}

// 4. Try different composition
{"num_atoms": 4}  // Fewer atoms
```

### Issue 5: MLFF Optimization Too Slow

**Solutions:**
```json
// 1. Use faster optimizer
{"optimizer": "FIRE"}

// 2. Reduce steps
{"steps": 200}

// 3. Relax convergence
{"fmax": 0.05}

// 4. Use faster model
{"mlff_model": "chgnet"}  // fastest
```

---

## 📊 Common Use Cases

### Use Case 1: Screen Materials Across Space Groups

**Goal:** Generate structures in multiple space groups to find stable phases

**Method:**
```json
{
  "tool": "generate_space_group_scan",
  "input": {
    "composition": ["Ti", "O", "O"],
    "space_groups": [136, 141, 225],
    "num_structures_per_group": 3
  }
}
```

### Use Case 2: Find Ground State Structure

**Goal:** Determine lowest energy configuration

**Method:**
```json
{
  "tool": "ground_state_search",
  "input": {
    "composition": ["Si", "Si"],
    "space_groups": [227, 141, 194, 12, 166],
    "num_structures_per_group": 5,
    "mlff_model": "chgnet",
    "optimization_settings": {
      "optimizer": "BFGS",
      "fmax": 0.01,
      "steps": 500
    }
  }
}
```

**Result:** Finds diamond structure (227) as ground state

### Use Case 3: Create Surface for Catalysis Study

**Goal:** Generate (111) surface slab

**Method:**
```json
{
  "tool": "generate_slab",
  "input": {
    "structure": {...},
    "miller_indices": [1, 1, 1],
    "thickness": 5,
    "vacuum": 15.0,
    "symmetric": true,
    "fix_bottom_layers": 2
  }
}
```

### Use Case 4: Create Supercell for DFT

**Goal:** Make 2×2×2 supercell for phonon calculations

**Method:**
```json
{
  "tool": "make_supercell",
  "input": {
    "structure": {...},
    "scaling_matrix": [2, 2, 2]
  }
}
```

---

## 🎯 Performance Optimization Tips

### Tip 1: Parallel Processing

For space group scans:
```json
{
  "parallel": true,
  "n_workers": 8  // Use 8 CPU cores
}
```

### Tip 2: Model Selection

| Task | Best Model | Reason |
|------|------------|--------|
| Quick screening | CHGNet | Fastest |
| Production calcs | M3GNet | Balanced |
| High accuracy | MACE | Most accurate |

### Tip 3: Optimizer Selection

| System Size | Best Optimizer | Reason |
|-------------|----------------|--------|
| Small (<50 atoms) | BFGS | Fast convergence |
| Medium (50-200) | LBFGS | Memory efficient |
| Large (>200) | FIRE | Handles rough landscapes |

### Tip 4: Memory Management

For large workflows:
```bash
# Increase Node.js memory
export NODE_OPTIONS="--max-old-space-size=8192"
npm start
```

---

## 📝 Configuration Files

### Space Groups Config

Location: `src/config/space-groups.json`

Contains:
- Crystal system definitions
- Lattice constraints
- Common structure types
- Validation rules

### Elements Config

Location: `src/config/elements.json`

Contains:
- Covalent radii (118 elements)
- Atomic masses
- Typical coordination numbers

### MLFF Models Config

Location: `src/config/mlff-models.json`

Contains:
- Available models
- Optimizer settings
- Convergence criteria
- Recommendations

---

## 🐛 Common Error Codes

| Code | Meaning | Solution |
|------|---------|----------|
| E1001 | Invalid space group | Use 1-230 |
| E1002 | Invalid composition | Check element symbols |
| E2002 | Max attempts exceeded | Increase volume_factor |
| E3001 | Atoms too close | Relax min_distance |
| E4001 | Model not available | Install MLFF package |
| E6001 | Python execution failed | Check Python path |

See full list in `src/types/errors.ts`

---

## 📚 Additional Resources

### Documentation
- `README.md` - Complete user guide
- `docs/API.md` - API reference
- `docs/examples/` - Usage examples

### Code
- `src/types/` - Type definitions
- `src/python/` - Python backend
- `src/tools/` - MCP tool implementations

### Configuration
- `package.json` - Node.js dependencies
- `requirements.txt` - Python dependencies
- `tsconfig.json` - TypeScript settings

---

## ✅ Installation Checklist

- [ ] Node.js 18+ installed
- [ ] Python 3.8+ installed
- [ ] pip installed and updated
- [ ] `npm install` completed
- [ ] `pip install -r requirements.txt` completed
- [ ] `npm run build` successful
- [ ] Python tests passing
- [ ] TypeScript compiles without errors
- [ ] Server starts successfully

---

## 🎉 Success Indicators

**Installation is successful when:**

1. ✅ `python3 tests/test_python_backend.py` shows all tests passing
2. ✅ `npm run build` completes without errors
3. ✅ `npm start` launches server
4. ✅ Can generate diamond Si structure
5. ✅ Can detect symmetry correctly

**System is ready for production use!**

---

## 🆘 Getting Help

If you encounter issues:

1. Check this guide's troubleshooting section
2. Review error messages (they include suggestions)
3. Verify all dependencies are installed
4. Check Python and Node.js versions
5. Try the minimal examples first

---

## 📦 Project Archive Location

**Archive:** `/mnt/user-data/outputs/crystal-mcp-server.tar.gz`

**Extract:**
```bash
tar -xzf crystal-mcp-server.tar.gz
cd crystal-mcp-server
```

---

*Quick Guide - Crystal MCP Server v1.0.0*  
*Last Updated: December 9, 2025*
