# Troubleshooting Guide

**Comprehensive Debugging Reference with Decision Trees and Visual Flowcharts**

This guide provides systematic troubleshooting procedures for all common issues encountered when using the Crystal MCP Server.

---

## Table of Contents

1. [Quick Diagnostic Flowchart](#quick-diagnostic-flowchart)
2. [Installation Issues](#installation-issues)
3. [Connection Problems](#connection-problems)
4. [Generation Failures](#generation-failures)
5. [Scientific Accuracy Issues](#scientific-accuracy-issues)
6. [Performance Problems](#performance-problems)
7. [Export/Format Errors](#exportformat-errors)
8. [Error Code Reference](#error-code-reference)

---

## Quick Diagnostic Flowchart

Start here to identify your issue category:

```mermaid
graph TB
    START[Problem Occurred] --> Q1{Can you connect<br/>to the server?}

    Q1 -->|No| CONN[Connection Issues<br/>Section 3]
    Q1 -->|Yes| Q2{Did structure<br/>generation fail?}

    Q2 -->|Yes| Q3{What error<br/>code?}
    Q2 -->|No| Q4{Is output<br/>incorrect?}

    Q3 --> INVALID[INVALID_*<br/>Section 4A]
    Q3 --> GENERATION[GENERATION_FAILED<br/>Section 4B]
    Q3 --> SCIENTIFIC[MIN_DISTANCE_*<br/>WYCKOFF_*<br/>Section 5]

    Q4 -->|Yes| SCI[Scientific Issues<br/>Section 5]
    Q4 -->|No| Q5{Export/format<br/>problem?}

    Q5 -->|Yes| EXPORT[Export Errors<br/>Section 7]
    Q5 -->|No| PERF[Performance<br/>Section 6]

    style START fill:#e1f5ff
    style CONN fill:#ffebee
    style GENERATION fill:#ffebee
    style SCI fill:#fff9c4
    style EXPORT fill:#c8e6c9
```

---

## Installation Issues

### Decision Tree: Installation Debugging

```mermaid
graph TB
    INSTALL[Installation Problem] --> CHECK1{Node.js<br/>installed?}

    CHECK1 -->|No| INSTALL_NODE[Install Node.js ≥18<br/>nodejs.org]
    CHECK1 -->|Yes| CHECK2{Python ≥3.9<br/>installed?}

    CHECK2 -->|No| INSTALL_PY[Install Python 3.9+<br/>python.org]
    CHECK2 -->|Yes| CHECK3{npm install<br/>succeeded?}

    CHECK3 -->|No| NPM_ERROR[Check npm error:<br/>- Network issues?<br/>- Permissions?<br/>- Disk space?]
    CHECK3 -->|Yes| CHECK4{pip install<br/>requirements.txt<br/>succeeded?}

    CHECK4 -->|No| PIP_ERROR[Check pip error]
    CHECK4 -->|Yes| CHECK5{npm run build<br/>succeeded?}

    CHECK5 -->|No| BUILD_ERROR[TypeScript compilation error]
    CHECK5 -->|Yes| SUCCESS[Installation Complete]

    PIP_ERROR --> PIP1{Missing system<br/>libraries?}
    PIP1 -->|Yes| SYSTEM[Install dev packages:<br/>python3-dev<br/>build-essential]
    PIP1 -->|No| PIP2{Version conflict?}
    PIP2 -->|Yes| VENV[Use virtual environment:<br/>python3 -m venv venv]

    style INSTALL fill:#e1f5ff
    style SUCCESS fill:#c8e6c9
    style NPM_ERROR fill:#ffebee
    style PIP_ERROR fill:#ffebee
```

### Common Installation Errors

#### Error: "gyp ERR! find Python"
**Symptom**: npm install fails with Python-related error

**Solution**:
```bash
# macOS
brew install python3

# Ubuntu/Debian
sudo apt install python3-dev python3-pip

# Set Python path
npm config set python /usr/bin/python3
```

#### Error: "Cannot find module 'pymatgen'"
**Symptom**: Python subprocess fails

**Solution**:
```bash
# Ensure correct Python environment
which python3

# Install in correct environment
pip3 install -r requirements.txt

# If using venv:
source venv/bin/activate
pip install -r requirements.txt
```

#### Error: "Permission denied" during pip install
**Solution**:
```bash
# Option 1: User install
pip3 install --user -r requirements.txt

# Option 2: Virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## Connection Problems

### MCP Connection Diagnostic

```mermaid
graph TB
    CONN_START[Connection Failed] --> CHECK_SERVER{Server<br/>process running?}

    CHECK_SERVER -->|No| START_SERVER[Start server:<br/>node dist/index.js]
    CHECK_SERVER -->|Yes| CHECK_STDIO{Correct<br/>stdio transport?}

    CHECK_STDIO -->|No| FIX_CONFIG[Fix Claude config:<br/>Use 'stdio' transport<br/>not 'http']
    CHECK_STDIO -->|Yes| CHECK_PATH{Correct path<br/>in config?}

    CHECK_PATH -->|No| UPDATE_PATH[Update config:<br/>Use absolute path to<br/>dist/index.js]
    CHECK_PATH -->|Yes| CHECK_PERMS{File<br/>permissions OK?}

    CHECK_PERMS -->|No| FIX_PERMS[chmod +x dist/index.js]
    CHECK_PERMS -->|Yes| CHECK_LOGS{Check<br/>server logs}

    CHECK_LOGS --> VIEW_LOGS[View stderr output:<br/>Look for Python errors]

    style CONN_START fill:#e1f5ff
    style VIEW_LOGS fill:#fff9c4
```

### Claude Desktop Connection Issues

#### Symptom: Hammer icon not showing
**Diagnostic Steps**:

1. **Check Config File Location**:
   ```bash
   # macOS
   ~/Library/Application Support/Claude/claude_desktop_config.json

   # Windows
   %APPDATA%\Claude\claude_desktop_config.json

   # Linux
   ~/.config/Claude/claude_desktop_config.json
   ```

2. **Validate JSON**:
   ```json
   {
     "mcpServers": {
       "crystal-gen": {
         "command": "node",
         "args": ["/absolute/path/to/crystal-mcp-server/dist/index.js"]
       }
     }
   }
   ```

3. **Check Logs**:
   ```bash
   # macOS
   tail -f ~/Library/Logs/Claude/mcp*.log

   # Look for:
   # - "Server started" (success)
   # - Python import errors
   # - Path not found errors
   ```

#### Symptom: Server starts then crashes
**Common Causes**:

```mermaid
graph LR
    CRASH[Server Crashes] --> C1[Python Not Found]
    CRASH --> C2[Missing Dependencies]
    CRASH --> C3[Permission Issues]

    C1 --> S1[Install Python 3.9+<br/>Add to PATH]
    C2 --> S2[pip install -r requirements.txt]
    C3 --> S3[chmod +x dist/index.js<br/>Check Python permissions]

    style CRASH fill:#ffebee
    style S1 fill:#c8e6c9
    style S2 fill:#c8e6c9
    style S3 fill:#c8e6c9
```

---

## Generation Failures

### Failure Diagnosis Decision Tree

```mermaid
graph TB
    FAIL[Generation Failed] --> ERROR_CODE{Error Code?}

    ERROR_CODE -->|INVALID_PARAMETER| PARAM[Check Parameters<br/>Section 4A]
    ERROR_CODE -->|GENERATION_FAILED| PYXTAL[PyXtal Issues<br/>Section 4B]
    ERROR_CODE -->|MIN_DISTANCE_VIOLATION| DISTANCE[Distance Issues<br/>Section 5A]
    ERROR_CODE -->|WYCKOFF_MISMATCH| WYCKOFF[Wyckoff Issues<br/>Section 5B]
    ERROR_CODE -->|Other| GENERIC[Generic Debug<br/>Section 4C]

    PARAM --> P1{Space group<br/>1-230?}
    PARAM --> P2{Valid<br/>elements?}
    PARAM --> P3{Positive<br/>parameters?}

    P1 -->|No| FIX_SG[Use valid space group number]
    P2 -->|No| FIX_ELEM[Check element symbols:<br/>Use standard notation<br/>e.g., Si not si]
    P3 -->|No| FIX_PARAM[Use positive values:<br/>lattice_constant > 0<br/>num_atoms > 0]

    style FAIL fill:#ffebee
    style FIX_SG fill:#c8e6c9
    style FIX_ELEM fill:#c8e6c9
```

### 4A. Invalid Parameter Errors

#### Error: "Space group not in range 1-230"
**Cause**: Invalid space group number

**Solution**:
```json
{
  "operation": "generate_from_spacegroup",
  "spacegroup": 227,  // ✓ Valid (1-230)
  // NOT: "spacegroup": 0 or 300 ❌
  "elements": ["Si"],
  "composition": [8]
}
```

**Reference**: See [Crystallography Guide](./crystallography_guide.md) for space group ranges

#### Error: "Invalid element symbol"
**Common Mistakes**:
```json
// ❌ Wrong
"elements": ["si", "Silicon", "SI"]

// ✓ Correct
"elements": ["Si"]

// Element symbol rules:
// - First letter uppercase
// - Second letter lowercase (if present)
// - Use standard IUPAC symbols
```

### 4B. PyXtal Generation Failures

#### Symptom: "Generation failed after 100 attempts"
**Cause**: Impossible or very difficult structure to generate

**Debugging Steps**:

```mermaid
graph TB
    PYXTAL_FAIL[PyXtal Generation Failed] --> CHECK1{Composition<br/>compatible with<br/>space group?}

    CHECK1 -->|No| WYCKOFF_CHECK[Check Wyckoff multiplicities:<br/>Composition must match<br/>available positions]
    CHECK1 -->|Yes| CHECK2{Min distance<br/>too restrictive?}

    CHECK2 -->|Yes| RELAX_DIST[Reduce min_distance:<br/>Try default first,<br/>then adjust]
    CHECK2 -->|No| CHECK3{Volume<br/>too small?}

    CHECK3 -->|Yes| INCREASE_VOL[Increase volume_factor:<br/>Try 1.2-1.5]
    CHECK3 -->|No| CHECK4{Rare/difficult<br/>space group?}

    CHECK4 -->|Yes| ALTERNATIVE[Try alternative:<br/>- Different space group<br/>- Prototype structure<br/>- Increase max_attempts]

    style PYXTAL_FAIL fill:#ffebee
    style ALTERNATIVE fill:#fff9c4
```

**Example Solutions**:

```json
// Problem: Generation fails for Fe8O12 in space group 227
{
  "operation": "generate_from_spacegroup",
  "spacegroup": 227,
  "elements": ["Fe", "O"],
  "composition": [8, 12],  // ❌ May not match Wyckoff positions
  "max_attempts": 200  // Increase attempts
}

// Solution 1: Adjust composition
{
  "spacegroup": 227,
  "elements": ["Fe", "O"],
  "composition": [8, 16],  // ✓ Matches Wyckoff (8a + 16c)
}

// Solution 2: Use prototype
{
  "operation": "generate_from_prototype",
  "prototype": "spinel",  // FeO₄ spinel structure
  "elements": ["Fe", "O"],
  "a": 8.4
}
```

---

## Scientific Accuracy Issues

### 5A. Minimum Distance Violations

**Symptom**: "MIN_DISTANCE_VIOLATION: Atoms too close"

**Diagnostic Flowchart**:

```mermaid
graph TB
    MIN_DIST[Min Distance Violation] --> Q1{User specified<br/>min_distance?}

    Q1 -->|Yes| Q2{Value too<br/>restrictive?}
    Q1 -->|No| Q3{Default too<br/>large for this<br/>system?}

    Q2 -->|Yes| REDUCE[Reduce min_distance:<br/>Check covalent radii:<br/>Si-Si: 2.3Å typical<br/>C-C: 1.5Å typical]
    Q2 -->|No| VOLUME[Increase volume_factor:<br/>Give atoms more space]

    Q3 -->|Yes| CUSTOM[Set custom min_distance:<br/>Based on target structure]
    Q3 -->|No| BUG[Possible bug:<br/>Report with structure]

    style MIN_DIST fill:#ffebee
    style REDUCE fill:#c8e6c9
    style CUSTOM fill:#c8e6c9
```

**Typical Bond Lengths** (for reference):

| Bond | Typical Length (Å) | Min Safe Distance |
|------|-------------------|-------------------|
| C-C | 1.54 | 1.3 |
| Si-Si | 2.35 | 2.0 |
| O-O | 1.48 | 1.2 |
| Metal-O | 1.8-2.2 | 1.5 |
| Au-Au | 2.88 | 2.5 |

**Solution Example**:
```json
{
  "operation": "generate_from_spacegroup",
  "spacegroup": 227,
  "elements": ["Si"],
  "composition": [8],
  "min_distance": {
    "Si-Si": 2.0  // ✓ Reasonable for diamond structure
    // NOT: "Si-Si": 3.0 ❌ Too large, generation will fail
  }
}
```

### 5B. Wyckoff Position Mismatches

**Understanding the Issue**:

```
Space Group 225 (Fm-3m) Wyckoff Positions:
  4a: (0,0,0) + fcc translations → 4 atoms
  4b: (1/2,1/2,1/2) + fcc translations → 4 atoms
  8c: (1/4,1/4,1/4) + fcc translations → 8 atoms
  ...

If you request 6 atoms total:
  ❌ Cannot fit into available multiplicities (4, 8, 12, ...)
  ✓ System adjusts to 8 atoms (nearest valid)
  ⚠️ Warning generated about composition change
```

**Solution Strategies**:

1. **Check Allowed Multiplicities**:
   ```python
   # Use list_category to see requirements
   {
     "operation": "list_category",
     "category": "bulk"
   }
   # Response includes Wyckoff information
   ```

2. **Let System Adjust** (recommended):
   ```json
   {
     "composition": [7]  // System adjusts to [8]
     // Check warnings in response
   }
   ```

3. **Specify Exact Wyckoff Positions**:
   ```json
   {
     "wyckoff_positions": [
       {"element": "Na", "wyckoff": "4a"},
       {"element": "Cl", "wyckoff": "4b"}
     ]
   }
   ```

---

## Performance Problems

### Slow Generation Diagnostic

```mermaid
graph TB
    SLOW[Generation Too Slow] --> SIZE{Structure<br/>size?}

    SIZE -->|Small < 50 atoms| PYXTAL_SLOW[PyXtal retry loop:<br/>Check max_attempts<br/>Difficult structure]
    SIZE -->|Medium 50-500| NORMAL[Normal for<br/>supercells]
    SIZE -->|Large > 500| Q_LARGE{Using MLFF<br/>optimization?}

    Q_LARGE -->|Yes| MLFF_SLOW[MLFF Slow:<br/>Section 6A]
    Q_LARGE -->|No| SUPERCELL_SLOW[Large supercell:<br/>Expected<br/>Consider smaller]

    PYXTAL_SLOW --> SOLUTION1[Solutions:<br/>- Reduce max_attempts<br/>- Use prototype instead<br/>- Simplify constraints]

    MLFF_SLOW --> SOLUTION2[Solutions:<br/>- Reduce fmax<br/>- Reduce max steps<br/>- Use GPU if available<br/>- Try faster model]

    style SLOW fill:#fff9c4
    style SOLUTION1 fill:#c8e6c9
    style SOLUTION2 fill:#c8e6c9
```

### 6A. MLFF Optimization Performance

**Typical Times**:
```
Structure Size → CPU Time (CHGNet)
  10 atoms     → 1-2 seconds
  50 atoms     → 5-10 seconds
  100 atoms    → 10-30 seconds
  500 atoms    → 1-5 minutes
  1000 atoms   → 5-15 minutes
```

**Optimization**:

1. **Use GPU** (if available):
   ```bash
   # Check CUDA available
   python3 -c "import torch; print(torch.cuda.is_available())"

   # Install GPU version
   pip install chgnet[cuda]
   ```

2. **Adjust Convergence**:
   ```json
   {
     "operation": "optimize_structure_mlff",
     "fmax": 0.05,  // Less strict (faster)
     // Instead of: "fmax": 0.001 (very strict, slower)
     "steps": 200   // Limit iterations
   }
   ```

3. **Use Faster Model**:
   ```json
   {
     "mlff_model": "m3gnet"  // Generally faster than CHGNet
     // Trade-off: Slightly less accurate
   }
   ```

---

## Export/Format Errors

### Export Failure Decision Tree

```mermaid
graph TB
    EXPORT_FAIL[Export Failed] --> FORMAT{Which<br/>format?}

    FORMAT -->|CIF| CIF_CHECK{Valid<br/>structure?}
    FORMAT -->|POSCAR| POSCAR_CHECK{Element<br/>order correct?}
    FORMAT -->|QE| QE_CHECK{Pseudopotentials<br/>specified?}
    FORMAT -->|LAMMPS| LAMMPS_CHECK{Atom types<br/>defined?}

    CIF_CHECK -->|No| FIX_STRUCT[Validate structure first:<br/>use validate_structure]
    POSCAR_CHECK -->|No| FIX_ELEM[Elements must be<br/>alphabetically sorted<br/>or by occurrence]
    QE_CHECK -->|No| FIX_PSEUDO[Specify pseudopotential<br/>files in export params]
    LAMMPS_CHECK -->|No| FIX_TYPES[Define atom types<br/>and masses]

    style EXPORT_FAIL fill:#ffebee
    style FIX_STRUCT fill:#c8e6c9
```

### Common Export Issues

#### Issue: POSCAR element order wrong
**Symptom**: VASP reads wrong elements

**Solution**:
```python
# POSCAR requires consistent element ordering
# VASP reads elements in order listed

Correct POSCAR:
Si O        ← Element line
4  8        ← Count line (4 Si, then 8 O)
Direct
<Si coords>  ← First 4 positions are Si
<O coords>   ← Next 8 positions are O
```

#### Issue: CIF export fails with "Invalid symmetry"
**Cause**: Structure doesn't match reported space group

**Solution**:
```json
// Step 1: Verify symmetry
{
  "operation": "analyze_symmetry",
  "structure": "<structure>"
}

// Step 2: If mismatch, either:
// A) Use detected space group
// B) Relax structure to match intended symmetry
{
  "operation": "optimize_structure_mlff",
  "structure": "<structure>",
  "constrain_symmetry": true,  // Preserve symmetry during relax
  "mlff_model": "chgnet"
}
```

---

## Error Code Reference

### Complete Error Code Catalog

| Error Code | Layer | Cause | Solution |
|------------|-------|-------|----------|
| `INVALID_PARAMETER` | Validation | Parameter type/range wrong | Check parameter documentation |
| `INVALID_SPACE_GROUP` | Validation | Space group not 1-230 | Use valid space group number |
| `INVALID_COMPOSITION` | Validation | Elements or counts invalid | Check element symbols, use positive integers |
| `INVALID_INPUT` | Validation | Generic input error | Review all parameters |
| `GENERATION_FAILED` | PyXtal | Could not generate structure | Adjust composition, increase attempts, or use prototype |
| `MIN_DISTANCE_VIOLATION` | Scientific | Atoms too close | Increase volume_factor or reduce min_distance |
| `WYCKOFF_MISMATCH` | Scientific | Composition incompatible | Adjust composition or let system auto-adjust |
| `COMPOSITION_MISMATCH` | Scientific | Stoichiometry impossible | Use valid composition ratios |
| `SYMMETRY_VERIFICATION_FAILED` | Validation | Generated structure doesn't match space group | Should not occur - report as bug |
| `EXPORT_ERROR` | Export | Cannot convert to format | Validate structure first |
| `SCHEMA_CONVERSION_ERROR` | Workflow | Schema incompatible | Check structure has required fields |
| `MISSING_DEPENDENCY` | System | Required library not installed | Install missing Python package |
| `TIMEOUT` | System | Operation took too long | Reduce structure size or complexity |

---

## Debug Mode and Logging

### Enable Debug Output

**For Python Errors**:
```python
# In your generator function, add:
import sys
print(f"DEBUG: parameter_value = {value}", file=sys.stderr)

# stderr output will appear in MCP logs
```

**For MCP Protocol Issues**:
```bash
# Set environment variable
export MCP_DEBUG=1

# Run server
node dist/index.js

# Or in Claude config:
{
  "mcpServers": {
    "crystal-gen": {
      "command": "node",
      "args": ["dist/index.js"],
      "env": {
        "MCP_DEBUG": "1"
      }
    }
  }
}
```

### View Logs

**Claude Desktop Logs**:
```bash
# macOS
tail -f ~/Library/Logs/Claude/mcp-server-crystal-gen.log

# Look for:
# - Server startup messages
# - Python subprocess errors
# - Tool call requests/responses
```

**Server Direct Logs** (when running manually):
```bash
node dist/index.js 2> error.log

# error.log will contain:
# - Python traceback
# - Import errors
# - Runtime exceptions
```

---

## Getting Help

### Before Reporting a Bug

**Checklist**:
- [ ] Reviewed this troubleshooting guide
- [ ] Checked error code in reference table
- [ ] Verified installation (npm install, pip install)
- [ ] Tested with simple example (e.g., Si diamond)
- [ ] Collected error messages and logs

**Information to Include**:
1. **Operation**: Which operation failed
2. **Parameters**: Full JSON parameters used
3. **Error Message**: Complete error text
4. **Error Code**: Error code from response
5. **Environment**:
   - OS (macOS, Linux, Windows)
   - Node.js version (`node --version`)
   - Python version (`python3 --version`)
   - Package versions (`pip list | grep -E "pymatgen|pyxtal"`)
6. **Logs**: Relevant log excerpts

**Report Issue**:
- GitHub Issues: [repository]/issues
- Include minimal reproducible example
- Attach logs if applicable

---

## Quick Fixes Cheat Sheet

```bash
# Server won't start
npm run build
node dist/index.js  # Test directly

# Python import errors
pip3 install -r requirements.txt
which python3  # Check Python path

# Permission errors
chmod +x dist/index.js
pip3 install --user -r requirements.txt

# MCP connection issues
# Check Claude config path is absolute
# Restart Claude Desktop

# Generation failures
# Try with prototype structure first
# Increase max_attempts
# Check space group/composition compatibility

# Slow performance
# Reduce structure size
# Lower fmax for MLFF
# Use GPU if available

# Export errors
# Validate structure first
# Check element ordering
# Verify format requirements
```

---

**Document Version**: 1.0
**Last Updated**: 2025-12-25
**Coverage**: All common issues with visual debugging flowcharts
**Maintainer**: Crystal Structure Generator Team
