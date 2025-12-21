# Quick Reference - MCP Server Testing

## 🚀 Setup (One-time)

```bash
# 1. Build server
npm run build

# 2. Install test dependencies
pip install pytest pytest-cov pytest-timeout pytest-xdist

# 3. Make test runner executable
chmod +x run_tests.py
```

---

## ⚡ Quick Commands

### Development (30 seconds)
```bash
./run_tests.py --profile dev
```

### Pre-Commit (2 minutes)
```bash
./run_tests.py --profile pre-commit
```

### Full Suite (10-15 minutes)
```bash
./run_tests.py --profile full
```

### Pre-Publish (15-20 minutes)
```bash
./run_tests.py --profile publish
open htmlcov/index.html  # View coverage
```

---

## 🎯 Test Categories

### Protocol Compliance (CRITICAL)
```bash
./run_tests.py --category protocol
```

### Crystal Generation
```bash
./run_tests.py --category generation
```

### Space Groups
```bash
./run_tests.py --category spacegroup
```

### Error Handling (CRITICAL)
```bash
./run_tests.py --category error
```

### Scientific Accuracy (CRITICAL)
```bash
./run_tests.py --category scientific
```

### Performance
```bash
./run_tests.py --category performance
```

---

## 🔍 Specific Tests

### Single Test
```bash
pytest test_mcp_comprehensive.py::TestCrystalGeneration::test_generate_simple_cubic_si -v
```

### Test Class
```bash
pytest test_mcp_comprehensive.py::TestProtocolCompliance -v
```

### By Keyword
```bash
pytest test_mcp_comprehensive.py -v -k "generate"
```

---

## 📊 Coverage & Reporting

### With Coverage Report
```bash
pytest test_mcp_comprehensive.py --cov --cov-report=html
open htmlcov/index.html
```

### Show Slowest Tests
```bash
pytest test_mcp_comprehensive.py --durations=10
```

### Stop on First Failure
```bash
pytest test_mcp_comprehensive.py -v -x
```

---

## 🐛 Debugging

### Verbose Output
```bash
pytest test_mcp_comprehensive.py -vv -s
```

### With Debugger
```bash
pytest test_mcp_comprehensive.py::TestCrystalGeneration --pdb
```

### Show Full Tracebacks
```bash
pytest test_mcp_comprehensive.py -v --tb=long
```

---

## ✅ Pre-Publish Checklist

```bash
# 1. Protocol compliance
./run_tests.py --category protocol

# 2. Error handling
./run_tests.py --category error

# 3. Scientific correctness
./run_tests.py --category scientific

# 4. Full suite with coverage
./run_tests.py --profile publish

# 5. Check coverage >90%
open htmlcov/index.html
```

---

## 📋 Test Files

- `test_mcp_comprehensive.py` - Main test suite (60+ tests)
- `test_operation_matrix.py` - All 228 operations (80+ tests)
- `run_tests.py` - Unified test runner
- `pytest.ini` - Configuration
- `TESTING_GUIDE.md` - Complete documentation
- `TEST_SUITE_OVERVIEW.md` - This summary

---

## 💡 Pro Tips

1. Run `--profile dev` frequently during development
2. Always run `--profile publish` before releasing
3. Monitor coverage trends over time
4. Fix critical test failures immediately
5. Add tests when adding new features

---

## 🎯 Success Metrics

- **Protocol**: 100% pass rate
- **Core Tests**: >98% pass rate
- **Space Groups**: >95% success rate
- **Coverage**: >90% overall
- **Performance**: No timeouts

---

## 📚 Help

```bash
# List all profiles
./run_tests.py --list-profiles

# List all categories
./run_tests.py --list-categories

# Pytest help
pytest --help

# View full documentation
cat TESTING_GUIDE.md
```

---

*Keep this card handy for quick reference!*
