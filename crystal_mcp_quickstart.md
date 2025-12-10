# Crystal Structure MCP Server - Quick Start Guide

## TL;DR: First Steps to Implementation

### Week 1: Foundation Setup

#### Day 1-2: Environment Setup
```bash
# 1. Create project
mkdir crystal-mcp-server
cd crystal-mcp-server

# 2. Initialize TypeScript project
npm init -y
npm install @modelcontextprotocol/sdk zod python-shell
npm install --save-dev @types/node typescript jest @types/jest

# 3. Setup Python environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install pyxtal pymatgen spglib ase chgnet

# 4. Create basic structure
mkdir -p src/{tools,python,types,utils}
touch src/index.ts src/server.ts
```

#### Day 3: Test PyXtal Installation
```python
# test_pyxtal.py
from pyxtal import pyxtal

# Test simple structure generation
crystal = pyxtal()
crystal.from_random(3, 225, ['Na', 'Cl'], [4, 4])  # Rock salt
print(f"Generated {crystal.formula}")
print(f"Space group: {crystal.group.symbol}")
crystal.to_file("NaCl.cif")
print("✓ PyXtal working!")
```

Run it: `python test_pyxtal.py`

#### Day 4-5: First MCP Tool
Create your first tool - basic structure generation:

**src/tools/generate-crystal.ts**
```typescript
import { z } from "zod";
import { PythonShell } from "python-shell";

export const GenerateCrystalSchema = z.object({
  composition: z.array(z.string()),
  space_group: z.number().int().min(1).max(230)
});

export async function generateCrystal(input: any) {
  // Minimal implementation to test MCP flow
  const result = await PythonShell.run('crystal_generator.py', {
    args: [JSON.stringify(input)]
  });
  
  return {
    content: [{
      type: "text",
      text: `Generated structure: ${result}`
    }]
  };
}
```

**src/python/crystal_generator.py**
```python
#!/usr/bin/env python3
import sys
import json
from pyxtal import pyxtal

def main():
    params = json.loads(sys.argv[1])
    
    crystal = pyxtal()
    crystal.from_random(
        dim=3,
        group=params['space_group'],
        species=params['composition'],
        numIons=len(params['composition'])
    )
    
    result = {
        "formula": crystal.formula,
        "space_group": crystal.group.number,
        "success": True
    }
    
    print(json.dumps(result))

if __name__ == "__main__":
    main()
```

Make it executable: `chmod +x src/python/crystal_generator.py`

### Week 2-3: Core Tools Implementation Priority

#### Implementation Order (MVP)
1. **generate_crystal** - Basic structure generation ✓ (Day 1-5)
2. **analyze_symmetry** - Validate what you create (Day 6-7)
3. **make_supercell** - Essential for DFT (Day 8-9)
4. **export_structure** - Multiple format support (Day 10-11)
5. **validate_structure** - Quality assurance (Day 12-14)

#### Sample Validation Test
```typescript
// Quick test to verify your tool works
const testInput = {
  composition: ["Si", "Si"],
  space_group: 227,  // Diamond structure
  seed: 42
};

const result = await generateCrystal(testInput);
// Should produce cubic diamond Si structure
```

### Week 4-5: MLFF Integration

#### Minimal MLFF Integration
```python
# src/python/mlff_calculator.py
from ase.optimize import BFGS
from chgnet.model.model import CHGNet

def optimize_with_chgnet(structure_dict):
    # Convert to ASE atoms
    atoms = dict_to_ase(structure_dict)
    
    # Load CHGNet
    chgnet = CHGNet.load()
    atoms.calc = chgnet.get_ase_calculator()
    
    # Optimize
    opt = BFGS(atoms)
    opt.run(fmax=0.01, steps=500)
    
    energy = atoms.get_potential_energy()
    
    return {
        "energy": float(energy),
        "converged": True,
        "structure": ase_to_dict(atoms)
    }
```

### Critical Decision Points

#### When to Use Each Tool
```
User Goal                          → Tool to Use
───────────────────────────────────────────────────
Generate one structure             → generate_crystal
Test all space groups              → generate_space_group_scan  
Prepare DFT calculation            → generate_crystal + export_structure
Find ground state                  → ground_state_search
Create surface                     → generate_slab
Create larger cell                 → make_supercell
Check structure validity           → validate_structure + analyze_symmetry
```

## Testing Strategy

### Unit Tests (Essential)
```typescript
// tests/basic.test.ts
test('generates diamond Si', async () => {
  const result = await generateCrystal({
    composition: ['Si', 'Si'],
    space_group: 227
  });
  
  expect(result.structure.space_group.number).toBe(227);
  expect(result.structure.metadata.formula).toBe('Si2');
});
```

### Integration Tests (Week 5)
Test complete workflows:
1. Generate → Validate → Export
2. Generate → Optimize → Calculate Energy
3. Generate Scan → Find Lowest Energy

## Common Pitfalls & Solutions

### Pitfall 1: PyXtal Generation Failures
**Problem:** Structure generation fails silently
**Solution:** Always set max_attempts and handle errors
```python
try:
    crystal.from_random(...)
except Exception as e:
    return {
        "error": str(e),
        "suggestion": "Try increasing volume_factor or relaxing constraints"
    }
```

### Pitfall 2: Space Group Incompatibility
**Problem:** Composition doesn't fit in space group
**Solution:** Pre-validate with Wyckoff multiplicities
```python
def check_compatibility(composition, space_group):
    group = Group(space_group)
    multiplicities = [wp.multiplicity for wp in group.Wyckoff_positions]
    n_atoms = len(composition)
    return any(n_atoms % m == 0 for m in multiplicities)
```

### Pitfall 3: MLFF Memory Issues
**Problem:** Running out of memory with large systems
**Solution:** Implement batching and cleanup
```python
import gc
import torch

def optimize_large_system(structure):
    # Process in chunks
    result = optimize(structure)
    
    # Clean up
    gc.collect()
    torch.cuda.empty_cache()  # If using GPU
    
    return result
```

## MVP Feature Set (Weeks 1-6)

### Must-Have Tools
1. ✓ generate_crystal
2. ✓ analyze_symmetry  
3. ✓ validate_structure
4. ✓ make_supercell
5. ✓ export_structure
6. ✓ optimize_structure_mlff

### Nice-to-Have (Weeks 7-10)
- generate_space_group_scan
- generate_slab
- ground_state_search
- create_defect
- generate_visualization

## Debugging Checklist

When something doesn't work:

1. **Python Environment**
   ```bash
   python -c "import pyxtal; print(pyxtal.__version__)"
   python -c "import chgnet; print('CHGNet OK')"
   ```

2. **MCP Connection**
   ```bash
   npx @modelcontextprotocol/inspector node dist/index.js
   ```

3. **Tool Registration**
   - Check tool appears in ListTools response
   - Verify schema is valid JSON Schema
   - Test with simple inputs first

4. **Python Bridge**
   ```typescript
   // Add logging
   const results = await PythonShell.run(script, {
     pythonPath: 'python3',
     mode: 'text',
     scriptPath: pythonDir
   });
   console.log('Python output:', results);
   ```

## Production Readiness Checklist

### Before Deploying
- [ ] All core tools tested with real-world examples
- [ ] Error messages are helpful and actionable  
- [ ] Documentation covers common use cases
- [ ] Performance benchmarks met (see main plan)
- [ ] Memory leaks checked and fixed
- [ ] Concurrent request handling tested
- [ ] Rate limiting implemented if needed
- [ ] Security review completed

### Deployment Options

#### Option 1: Local (Development)
```bash
node dist/index.js
```

#### Option 2: Docker (Production)
```bash
docker build -t crystal-mcp .
docker run -p 3000:3000 crystal-mcp
```

#### Option 3: Cloud (Scale)
- Deploy to AWS Lambda / Google Cloud Functions
- Use containerized deployment
- Add job queue for long-running tasks

## Resource Requirements

### Minimal Setup
- CPU: 4 cores
- RAM: 8 GB
- Storage: 10 GB
- Python 3.8+
- Node.js 18+

### Recommended Setup
- CPU: 8+ cores (for parallel generation)
- RAM: 16+ GB (for MLFF calculations)
- Storage: 50 GB (for structure database)
- GPU: Optional (speeds up MLFF)

## Getting Help

### Key Resources
1. PyXtal docs: https://pyxtal.readthedocs.io
2. MCP SDK: https://github.com/modelcontextprotocol
3. Spglib docs: https://spglib.github.io/spglib/
4. ASE docs: https://wiki.fysik.dtu.dk/ase/

### Common Issues Database
```
Issue: "Cannot import pyxtal"
Solution: Check virtual environment is activated

Issue: "Structure generation timeout"  
Solution: Increase max_attempts, adjust constraints

Issue: "MLFF model not found"
Solution: Download model: chgnet.model.download_model()

Issue: "MCP server not responding"
Solution: Check server.ts error logs, verify tool schemas
```

## Next Steps After MVP

Once you have the 6 core tools working:

1. **Performance Optimization**
   - Profile and optimize hot paths
   - Add caching for expensive operations
   - Implement parallel structure generation

2. **Enhanced Features**
   - Add more MLFF models (M3GNet, MACE)
   - Implement symmetry-constrained optimization
   - Add structure similarity detection

3. **User Experience**
   - Create web interface for visualization
   - Add progress tracking for long operations
   - Implement structure database

4. **Scientific Validation**
   - Compare generated structures with experimental data
   - Validate ground state predictions
   - Benchmark MLFF accuracy

---

## Rapid Prototyping Path (3 Days)

If you want to test the concept quickly:

### Day 1: Minimal Server
```typescript
// Absolute minimum MCP server
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { generateCrystal } from "./tools/generate-crystal.js";

const server = new Server({ name: "crystal-mcp", version: "0.1.0" });

server.setRequestHandler(ListToolsRequestSchema, async () => ({
  tools: [{
    name: "generate_crystal",
    description: "Generate crystal structure",
    inputSchema: { /* minimal schema */ }
  }]
}));

server.setRequestHandler(CallToolRequestSchema, async (req) => {
  if (req.params.name === "generate_crystal") {
    return await generateCrystal(req.params.arguments);
  }
});

// Start
const transport = new StdioServerTransport();
await server.connect(transport);
```

### Day 2: Test with Real Structures
Generate these to verify:
- Si (diamond): space group 227
- NaCl (rock salt): space group 225  
- TiO2 (rutile): space group 136

### Day 3: Add MLFF Optimization
Optimize one structure and calculate energy.

If these work, you've proven the concept and can build out the full server!

---

**Remember:** Start simple, test frequently, and build incrementally. The full plan is comprehensive, but you can achieve a working prototype in 1-2 weeks by focusing on the MVP tools first.
