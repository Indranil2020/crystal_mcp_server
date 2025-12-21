# Comprehensive Testing Guide for Crystal MCP Server

## Overview

This guide provides a complete testing strategy for your Crystal Structure MCP server with 228 operations. The testing suite ensures full functionality, scientific accuracy, and production-readiness before publishing.

## Test Suite Structure

```
tests/
├── test_mcp_comprehensive.py      # Main comprehensive test suite
├── test_operation_matrix.py       # 228 operation coverage tests
├── test_mcp_e2e.py               # Your original E2E tests
├── test_scientific_accuracy.py   # Scientific validation
└── test_performance_benchmarks.py # Performance testing
```

---

## Quick Start

### Prerequisites

```bash
# 1. Build the server
npm run build

# 2. Install Python testing dependencies
pip install pytest pytest-cov pytest-timeout pytest-xdist

# 3. Verify server is built
ls -l dist/index.js
```

### Run All Tests

```bash
# Full test suite (comprehensive + operation matrix)
pytest test_mcp_comprehensive.py test_operation_matrix.py -v

# With coverage report
pytest test_mcp_comprehensive.py --cov --cov-report=html

# Parallel execution (faster)
pytest test_mcp_comprehensive.py -n 4
```

---

## Test Categories

### 1. Protocol Compliance Tests

**Purpose**: Verify JSON-RPC 2.0 and MCP protocol compliance

**Run**: 
```bash
pytest test_mcp_comprehensive.py::TestProtocolCompliance -v
```

**What it tests**:
- JSON-RPC 2.0 version handling
- MCP protocol version negotiation
- Server initialization and capabilities
- Error handling for malformed requests
- Unknown method handling

**Critical for**: Publishing, as MCP specification compliance is mandatory

---

### 2. Tool Discovery Tests

**Purpose**: Validate tool schema and discoverability

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestToolDiscovery -v
```

**What it tests**:
- `tools/list` returns proper structure
- All tools have required fields (name, description, inputSchema)
- Schemas are valid JSON Schema
- Tool descriptions are informative
- Tool count verification

**Critical for**: LLM discoverability and usability

---

### 3. Core Functionality Tests

**Purpose**: Test actual crystal generation and manipulation

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestCrystalGeneration -v
pytest test_mcp_comprehensive.py::TestTransformations -v
pytest test_mcp_comprehensive.py::TestDefects -v
```

**What it tests**:
- Basic crystal structure generation
- All 230 space groups
- Supercell creation
- Strain application
- Defect generation (vacancies, interstitials)
- Surface slab creation

**Critical for**: Core server functionality

---

### 4. Space Group Coverage Tests

**Purpose**: Verify all 230 crystallographic space groups work

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestSpaceGroups -v
pytest test_operation_matrix.py::TestBulkStructures -v
```

**What it tests**:
- Representative samples from all 7 crystal systems
- Edge cases (space group 1, 230)
- Invalid space groups (0, 999)
- Lattice parameter requirements per system

**Expected runtime**: ~5-15 minutes

---

### 5. Operation Matrix Tests (228 Operations)

**Purpose**: Full coverage of all server capabilities

**Run**:
```bash
pytest test_operation_matrix.py -v
```

**Categories tested**:

#### Bulk Structures
```bash
pytest test_operation_matrix.py::TestBulkStructures -v
```
- All 230 space groups
- Various crystal systems

#### Defects
```bash
pytest test_operation_matrix.py::TestDefectGeneration -v
```
- Point defects (vacancy, interstitial, substitution)
- Complex defects (Frenkel pairs, Schottky defects)

#### Surfaces
```bash
pytest test_operation_matrix.py::TestSurfaceGeneration -v
```
- Miller indices: (100), (110), (111), (210), (211)
- Adsorbate placement
- Surface reconstruction

#### 2D Materials
```bash
pytest test_operation_matrix.py::Test2DMaterials -v
```
- Graphene, hBN, MoS₂, WS₂
- Nanoribbons
- Grain boundaries

#### Twistronics
```bash
pytest test_operation_matrix.py::TestTwistronics -v
```
- Twisted bilayers
- Moiré superlattices
- Magic angle structures

---

### 6. Scientific Correctness Tests

**Purpose**: Validate physical and chemical accuracy

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestScientificCorrectness -v
```

**What it tests**:
- Lattice parameters preserved
- Composition stoichiometry correct
- Coordination numbers accurate
- Symmetry operations valid
- Physical constraints respected

**Critical for**: Scientific credibility and publication

---

### 7. Error Handling Tests

**Purpose**: Ensure graceful failure and helpful error messages

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestErrorHandling -v
```

**What it tests**:
- Missing required parameters
- Wrong parameter types
- Invalid values (negative lattice, invalid elements)
- Empty/malformed structures
- Edge cases

**Critical for**: User experience and robustness

---

### 8. Export Format Tests

**Purpose**: Verify all output formats work correctly

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestExportFormats -v
pytest test_operation_matrix.py::TestAllExportFormats -v
```

**What it tests**:
- VASP POSCAR format
- CIF format
- XYZ format
- Other formats (PDB, LAMMPS, Quantum ESPRESSO, etc.)

---

### 9. Performance Tests

**Purpose**: Ensure operations complete in reasonable time

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestPerformance -v
```

**Benchmarks**:
- Simple generation: < 3 seconds
- Large supercell (3x3x3): < 10 seconds
- Rapid sequential requests: 20 requests without failure
- Space group scan: < 15 seconds per group

---

### 10. Integration Workflow Tests

**Purpose**: Test multi-step workflows

**Run**:
```bash
pytest test_mcp_comprehensive.py::TestIntegrationWorkflows -v
```

**Workflows tested**:
- Generate → Analyze → Export
- Generate → Transform → Validate
- Generate → Defect → Analyze
- Complete research workflows

---

## Test Execution Strategies

### Development (Fast Feedback)

```bash
# Run only protocol and basic functionality
pytest test_mcp_comprehensive.py::TestProtocolCompliance \
       test_mcp_comprehensive.py::TestCrystalGeneration -v

# Stop on first failure
pytest test_mcp_comprehensive.py -x -v
```

### Pre-Commit (Medium Coverage)

```bash
# Representative sample of tests
pytest test_mcp_comprehensive.py -v --tb=short \
       -k "Protocol or Generation or Export"
```

### Pre-Publish (Full Coverage)

```bash
# Complete test suite with coverage
pytest test_mcp_comprehensive.py test_operation_matrix.py \
       -v --cov --cov-report=html --cov-report=term

# Generate HTML report
open htmlcov/index.html
```

### Continuous Integration

```bash
# Parallel execution with detailed reporting
pytest test_mcp_comprehensive.py test_operation_matrix.py \
       -n auto --tb=short --junitxml=test-results.xml \
       --cov --cov-report=xml
```

---

## Test Result Interpretation

### Success Criteria

✅ **Protocol Compliance**: 100% pass rate required
✅ **Tool Discovery**: All tools properly listed and documented
✅ **Space Groups**: >95% success rate (some edge cases may fail)
✅ **Core Operations**: >98% success rate
✅ **Scientific Correctness**: 100% pass rate
✅ **Export Formats**: 100% pass rate for implemented formats

### Expected Failures

Some tests may legitimately fail if features are not yet implemented:

- Advanced 2D materials (graphene nanoribbons, etc.)
- Twistronics features (twisted bilayers)
- Some exotic defect types
- Certain export formats

**These should fail gracefully** with clear error messages, not crash.

---

## Performance Benchmarks

### Target Performance

| Operation | Target Time | Acceptable |
|-----------|-------------|------------|
| Simple structure generation | < 1s | < 3s |
| Space group validation | < 0.5s | < 2s |
| 2x2x2 Supercell | < 2s | < 5s |
| 3x3x3 Supercell | < 5s | < 10s |
| Surface slab (10Å) | < 3s | < 8s |
| Symmetry analysis | < 1s | < 3s |
| Export to VASP | < 0.5s | < 1s |

### Measuring Performance

```bash
# Show 10 slowest tests
pytest test_mcp_comprehensive.py --durations=10

# Detailed timing
pytest test_mcp_comprehensive.py -v --durations=0
```

---

## Coverage Goals

### Minimum Acceptable Coverage

- **Protocol Layer**: 100%
- **Tool Handlers**: >95%
- **Core Generation**: >90%
- **Transformations**: >85%
- **Export Functions**: >90%
- **Error Handling**: 100%

### Generating Coverage Reports

```bash
# Terminal report
pytest test_mcp_comprehensive.py --cov --cov-report=term-missing

# HTML report (detailed)
pytest test_mcp_comprehensive.py --cov --cov-report=html
open htmlcov/index.html

# XML report (for CI)
pytest test_mcp_comprehensive.py --cov --cov-report=xml
```

---

## Debugging Failed Tests

### Verbose Output

```bash
# Maximum verbosity
pytest test_mcp_comprehensive.py -vv

# Show stdout/stderr
pytest test_mcp_comprehensive.py -v -s

# Full tracebacks
pytest test_mcp_comprehensive.py -v --tb=long
```

### Debugging Specific Test

```bash
# Run single test with debugging
pytest test_mcp_comprehensive.py::TestCrystalGeneration::test_generate_simple_cubic_si -vv -s

# Drop into debugger on failure
pytest test_mcp_comprehensive.py::TestCrystalGeneration -v --pdb
```

### Server Logs

The test client captures stderr. Look for:

```
[SERVER STDERR]: ...error details...
```

Common issues:
- Python path not found
- Missing dependencies (pyxtal, pymatgen, etc.)
- TypeScript compilation errors
- Port binding issues

---

## Pre-Publishing Checklist

### Before publishing your MCP server, ensure:

#### 1. Protocol Compliance ✓
```bash
pytest test_mcp_comprehensive.py::TestProtocolCompliance -v
# All tests must PASS
```

#### 2. Tool Schema Validation ✓
```bash
pytest test_mcp_comprehensive.py::TestToolDiscovery -v
# All tools properly documented
```

#### 3. Core Functionality ✓
```bash
pytest test_mcp_comprehensive.py::TestCrystalGeneration -v
pytest test_mcp_comprehensive.py::TestTransformations -v
# >95% pass rate required
```

#### 4. Space Group Coverage ✓
```bash
pytest test_operation_matrix.py::TestBulkStructures -v
# Test representative samples from all crystal systems
```

#### 5. Error Handling ✓
```bash
pytest test_mcp_comprehensive.py::TestErrorHandling -v
# All tests must PASS - graceful failures required
```

#### 6. Export Formats ✓
```bash
pytest test_mcp_comprehensive.py::TestExportFormats -v
# All advertised formats must work
```

#### 7. Performance Acceptable ✓
```bash
pytest test_mcp_comprehensive.py::TestPerformance -v
# No operation should timeout
```

#### 8. Scientific Validation ✓
```bash
pytest test_mcp_comprehensive.py::TestScientificCorrectness -v
# All tests must PASS
```

#### 9. Documentation Complete ✓
- [ ] README.md with usage examples
- [ ] All tools have clear descriptions
- [ ] Error messages are actionable
- [ ] Example workflows documented

#### 10. Final Comprehensive Run ✓
```bash
pytest test_mcp_comprehensive.py test_operation_matrix.py \
       -v --cov --cov-report=html \
       --junitxml=test-results.xml

# Review coverage report
open htmlcov/index.html

# Coverage should be >90% overall
```

---

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Test MCP Server

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v2
    
    - name: Setup Node
      uses: actions/setup-node@v2
      with:
        node-version: '18'
    
    - name: Setup Python
      uses: actions/setup-python@v2
      with:
        python-version: '3.10'
    
    - name: Install dependencies
      run: |
        npm install
        pip install pytest pytest-cov pytest-xdist
        pip install pyxtal pymatgen ase chgnet
    
    - name: Build server
      run: npm run build
    
    - name: Run tests
      run: |
        pytest test_mcp_comprehensive.py test_operation_matrix.py \
               -n auto -v --cov --cov-report=xml \
               --junitxml=test-results.xml
    
    - name: Upload coverage
      uses: codecov/codecov-action@v2
      with:
        file: ./coverage.xml
```

---

## Test Maintenance

### Adding New Operations

When adding a new operation to your server:

1. **Add to `test_operation_matrix.py`**:
   ```python
   def test_new_operation(self, initialized_client):
       """Test the new operation"""
       res = initialized_client.send_request("tools/call", {
           "name": "comprehensive_generate",
           "arguments": {
               "operation": "new_operation_name",
               # ... parameters
           }
       })
       
       assert "result" in res
       data = json.loads(res["result"]["content"][-1]["text"])
       assert data.get("success") is True
   ```

2. **Add integration test** if it's part of a workflow

3. **Update documentation**

### Updating Test Expectations

If you change server behavior:

1. Update corresponding test assertions
2. Document breaking changes
3. Ensure backward compatibility where possible

---

## Common Test Failures and Solutions

### Issue: "Server not built"
```
Solution: Run `npm run build` before testing
```

### Issue: "No response from server"
```
Solution: Check server stderr logs, verify Python dependencies installed
```

### Issue: "Space group X failed"
```
Solution: Some high-symmetry groups are challenging. Adjust lattice parameters
or mark as xfail if systematic PyXtal limitation
```

### Issue: "Timeout in performance test"
```
Solution: Increase timeout for complex operations, or optimize Python backend
```

### Issue: "Scientific correctness test failed"
```
Solution: Critical! Debug immediately - indicates physical error
```

---

## Advanced Testing

### Testing with MCP Inspector

```bash
# Install MCP Inspector
npx @modelcontextprotocol/inspector

# Launch inspector pointing to your server
npx @modelcontextprotocol/inspector dist/index.js

# Manually test operations through UI
```

### Load Testing

```bash
# Install locust (HTTP load testing)
pip install locust

# Create locustfile.py for your server
# Run load test
locust -f locustfile.py
```

### Profiling Performance

```bash
# Profile Python backend
python -m cProfile -o profile.stats src/python/crystal_generator.py

# Analyze results
python -m pstats profile.stats
```

---

## Getting Help

If tests fail unexpectedly:

1. **Check server logs**: Look for stderr output in test results
2. **Run single test**: Isolate the failing test
3. **Enable verbose mode**: `pytest -vv -s`
4. **Check dependencies**: Verify pyxtal, pymatgen versions
5. **Review project knowledge**: Check architecture documentation

---

## Summary

This comprehensive test suite ensures your Crystal MCP server is:

✅ **MCP Compliant**: Follows all protocol specifications
✅ **Scientifically Accurate**: Produces correct crystal structures
✅ **Robust**: Handles errors gracefully
✅ **Performant**: Completes operations in reasonable time
✅ **Complete**: Covers all 228 advertised operations
✅ **Production-Ready**: Ready for publishing and real-world use

**Recommended workflow**:
1. Development: Run quick tests frequently
2. Pre-commit: Run representative sample
3. Pre-publish: Run full suite with coverage
4. Post-deploy: Monitor with integration tests

**Test early, test often, publish confidently!** 🚀
