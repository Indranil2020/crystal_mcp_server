# Molecular Arrangement Engine

## A Production-Ready Generic Molecular Clustering System for MCP

This system synthesizes **three complementary approaches** into a unified, production-ready backend for the Molecule Generator MCP with PubChem's 130M molecule database.

---

## 🎯 Design Philosophy

The system addresses the fundamental limitation of **hardcoded stacking patterns** by providing:

1. **Infinite Flexibility**: Any arrangement conceivable by a chemist can be expressed
2. **Chemical Intelligence**: Understands H-bonds, π-stacking, ring centers, etc.
3. **Natural Language Support**: Parse requests like "Stack benzenes with π-stacking"
4. **Scientific Rigor**: Based on proper rigid body mechanics (6 DOF per molecule)
5. **Deterministic Output**: Always produces valid coordinates (no stochastic search)

---

## 🏗️ Architecture: Four-Layer Design

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Layer 4: Natural Language Interface (MCP Entry Point)                   │
│  ───────────────────────────────────────────────────────────────────────│
│  parse_arrangement_request("Stack 3 benzenes with π-stacking at 3.4Å")   │
│  → {'pattern': 'pi_stacking', 'n_molecules': 3, 'distance': 3.4}        │
├─────────────────────────────────────────────────────────────────────────┤
│  Layer 3: Chemical Intelligence (K2_plan)                                │
│  ───────────────────────────────────────────────────────────────────────│
│  AtomSelector DSL:    "0:ring_center(0)", "1:donor_h(0)"                │
│  ChemicalPatterns:    π-stacking, H-bonding, T-shaped, herringbone      │
│  Constraints:         Distance, Angle, PlaneAlignment, HBond, Clash     │
│  ConstraintSolver:    Gradient descent optimization                      │
├─────────────────────────────────────────────────────────────────────────┤
│  Layer 2: Mathematical Generators (Original Plan)                        │
│  ───────────────────────────────────────────────────────────────────────│
│  Position Generators: Linear, Circular, Helical, Grid, Spherical        │
│  Formula DSL:         "x = 3*sqrt(i)*cos(2.4*i)", "z = 3.4*i"           │
│  Orientation Gen.:    Fixed, Alternating, Incremental, FaceCenter       │
├─────────────────────────────────────────────────────────────────────────┤
│  Layer 1: Core Rigid Body Engine (z2_plan)                               │
│  ───────────────────────────────────────────────────────────────────────│
│  Position:            3D vector (x, y, z)                                │
│  Orientation:         Euler angles or quaternion                         │
│  MolecularFrame:      Local coordinate system (origin + 3 axes)          │
│  MoleculePose:        6-DOF specification (absolute or relative)         │
│  Relative Placement:  Parent-child hierarchical positioning              │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 📚 Synthesis of Three Plans

### From K2_plan (Constraint-Based)

**Key Contribution**: Chemical intelligence and declarative constraints

```python
# AtomSelector DSL - Chemistry-aware targeting
"0:centroid()"           # Molecular center of mass
"0:ring_center(0)"       # Center of first aromatic ring
"0:plane_normal()"       # Normal to molecular plane
"0:donor_h(0)"           # First H-bond donor hydrogen
"0:func_group(carbonyl)" # Carbonyl group center

# Chemical Constraints
DistanceConstraint("0:centroid()", "1:centroid()", target=3.4)
PlaneAlignmentConstraint(0, 1, mode="parallel")
HBondConstraint(donor_mol=0, acceptor_mol=1)
ClashConstraint(0, 1)  # Prevent atomic overlaps

# Constraint Solver (gradient descent)
solver = ConstraintSolver(molecules, poses, constraints)
optimized_poses = solver.solve(max_iterations=1000)
```

### From z2_plan (Pose-Based)

**Key Contribution**: Relative placement and molecular frames

```python
# MolecularFrame - Local coordinate systems
frame = MolecularFrame.from_coords(coords, method="svd")
# Returns origin + x_axis + y_axis + z_axis (plane normal)

# MoleculePose - 6 degrees of freedom
pose = MoleculePose(
    position=Position(0, 0, 3.4),
    orientation=Orientation(roll=0, pitch=0, yaw=180),
    molecule_idx=1
)

# Relative Placement - Parent-child relationships
relative_pose = MoleculePose(
    parent_idx=0,                    # Relative to molecule 0
    local_offset=Position(0, 0, 3.4),  # In parent's frame
    local_orientation=Orientation(yaw=180)
)
```

### From Original Plan (Parametric)

**Key Contribution**: Mathematical formula generators

```python
# Position generators
positions = generate_linear_positions(n=5, spacing=3.4)
positions = generate_circular_positions(n=6, radius=5.0)
positions = generate_helical_positions(n=10, radius=5.0, pitch=3.4, turns=2)
positions = generate_spherical_positions(n=20, radius=8.0)

# Formula-based positions (infinite flexibility)
positions = generate_formula_positions(
    n=10,
    x_formula="3 * sqrt(i) * cos(2.4 * i)",  # Golden spiral
    y_formula="3 * sqrt(i) * sin(2.4 * i)",
    z_formula="0.5 * i"
)

# Orientation generators
orientations = generate_fixed_orientations(n, roll=0, pitch=0, yaw=0)
orientations = generate_alternating_orientations(n, (0,0,0), (0,0,180))
orientations = generate_face_center_orientations(positions, facing="inward")
```

---

## 🚀 Usage Examples

### 1. Named Patterns (Simplest)

```python
from molecular_arrangement_engine import generate_molecular_cluster

# π-π stacking
result = generate_molecular_cluster(
    molecules=[benzene, benzene, benzene],
    arrangement="pi_pi_parallel",
    distance=3.4
)

# Circular H-bonded arrangement
result = generate_molecular_cluster(
    molecules=[water] * 6,
    arrangement="h_bonded_circular",
    distance=2.8
)

# T-shaped dimer
result = generate_molecular_cluster(
    molecules=[naphthalene, naphthalene],
    arrangement="t_shaped",
    distance=5.0
)
```

### 2. Natural Language (MCP Integration)

```python
from molecular_arrangement_engine import parse_arrangement_request, generate_molecular_cluster

# Parse user's natural language request
request = "Stack 4 benzene molecules with π-stacking at 3.5 angstroms"
spec = parse_arrangement_request(request)
# → {'pattern': 'pi_stacking', 'n_molecules': 4, 'distance': 3.5}

# Generate the cluster
result = generate_molecular_cluster(
    molecules=[benzene] * spec['n_molecules'],
    arrangement=spec['pattern'],
    distance=spec['distance']
)
```

### 3. Custom Formulas (Advanced Users)

```python
from molecular_arrangement_engine import generate_formula_positions, generate_molecular_cluster

# Archimedean spiral
positions = generate_formula_positions(
    n=12,
    x_formula="(2 + 0.5*i) * cos(i * 0.8)",
    y_formula="(2 + 0.5*i) * sin(i * 0.8)",
    z_formula="0"
)

# Double helix (DNA-like)
positions = generate_formula_positions(
    n=20,
    x_formula="10 * cos(i * 0.6)",
    y_formula="10 * sin(i * 0.6)",
    z_formula="3.4 * i"
)
```

### 4. Chemical Constraints (Research)

```python
from molecular_arrangement_engine import (
    generate_molecular_cluster,
    DistanceConstraint,
    HBondConstraint,
    PlaneAlignmentConstraint
)

# Define custom constraints
constraints = [
    DistanceConstraint("0:ring_center(0)", "1:ring_center(0)", target=3.5),
    PlaneAlignmentConstraint(0, 1, mode="parallel", tolerance=5.0),
    HBondConstraint(donor_mol=0, acceptor_mol=1)
]

result = generate_molecular_cluster(
    molecules=[mol_a, mol_b],
    arrangement="pi_pi_parallel",
    constraints=constraints,
    optimize=True  # Use solver to satisfy constraints
)
```

### 5. Constraint DSL (String-Based)

```python
result = generate_molecular_cluster(
    molecules=[...],
    arrangement="custom",
    constraints=[
        "distance(sel1=0:centroid(), sel2=1:centroid(), target=3.4)",
        "h_bond(donor=0, acceptor=1)",
        "plane_parallel(mol1=0, mol2=1)"
    ],
    optimize=True
)
```

---

## 📋 Available Patterns

| Pattern | Description | Typical Distance |
|---------|-------------|------------------|
| `pi_pi_parallel` | Face-to-face π-stacking | 3.3-3.8 Å |
| `pi_pi_antiparallel` | 180° rotated stacking | 3.3-3.8 Å |
| `pi_pi_offset` | Slip-stacked arrangement | 3.3-4.0 Å |
| `t_shaped` | Edge-to-face arrangement | 4.5-5.5 Å |
| `h_bonded_circular` | Circular H-bonded ring | 2.5-3.2 Å |
| `herringbone` | Alternating tilt pattern | 5.0-6.0 Å |
| `helical` | Helical arrangement | Variable |
| `sandwich` | A-B-A intercalated | 3.3-3.8 Å |
| `spherical` | Fibonacci lattice sphere | Variable |
| `grid` | 2D planar grid | Variable |
| `linear` | Simple linear stack | Variable |

---

## 🔬 Chemical Intelligence Features

### AtomSelector Targets

| Target | Example | Returns |
|--------|---------|---------|
| `centroid()` | `0:centroid()` | Molecular center |
| `atom(elem)` | `0:atom(O)` | First oxygen coordinates |
| `atom(elem, n)` | `0:atom(C, 2)` | Third carbon coordinates |
| `ring_center(n)` | `0:ring_center(0)` | Center of aromatic ring |
| `plane_normal()` | `0:plane_normal()` | Normal vector to plane |
| `func_group(name)` | `0:func_group(carbonyl)` | Functional group center |
| `donor_h(n)` | `0:donor_h(0)` | H-bond donor hydrogen |
| `acceptor(n)` | `0:acceptor(0)` | H-bond acceptor atom |

### Constraint Types

| Constraint | Purpose | Parameters |
|------------|---------|------------|
| `DistanceConstraint` | Control point-to-point distance | selector1, selector2, target/min/max |
| `AngleConstraint` | Control angle between 3 points | selector1, selector2, selector3, target |
| `PlaneAlignmentConstraint` | Parallel/perpendicular planes | mol1, mol2, mode, tolerance |
| `HBondConstraint` | Enforce H-bond geometry | donor_mol, acceptor_mol |
| `ClashConstraint` | Prevent atomic overlaps | mol1, mol2, scale_factor |

---

## 🔧 Formula Variables and Functions

For `generate_formula_positions()`:

| Variable | Meaning |
|----------|---------|
| `i` | Current molecule index (0, 1, 2, ...) |
| `n` | Total number of molecules |
| `t` | Normalized index = i/(n-1), range [0, 1] |

| Function | Description |
|----------|-------------|
| `sin`, `cos`, `tan` | Trigonometric |
| `sqrt`, `exp`, `log` | Mathematical |
| `abs` | Absolute value |
| `atan2` | Two-argument arctangent |
| `pi`, `e` | Constants |

---

## ✅ Validation Features

```python
result = generate_molecular_cluster(molecules, arrangement, validate=True)

# Validation results
result['validation'] = {
    'valid': True/False,
    'errors': ['Clash between mol 0 and mol 1: 1.5 Å'],
    'warnings': ['Minimum distance 2.1 Å is very small'],
    'statistics': {
        'min_intermolecular_distance': 2.1,
        'n_molecules': 5
    }
}
```

---

## 📊 Output Format

```python
result = generate_molecular_cluster(...)

# Result structure:
{
    'molecules': [
        {
            'atoms': ['C', 'C', ...],
            'coords': [[x, y, z], ...],  # Transformed coordinates
            'identifier': 'benzene'
        },
        ...
    ],
    'poses': [
        {
            'position': [0.0, 0.0, 0.0],
            'orientation': [roll, pitch, yaw],  # degrees
            'molecule_idx': 0
        },
        ...
    ],
    'validation': {...},
    'metadata': {
        'pattern': 'pi_pi_parallel',
        'n_molecules': 3,
        'distance': 3.4,
        'optimized': False
    }
}
```

---

## 🔗 MCP Integration

```python
# In your MCP tool handler:

from molecular_arrangement_engine import (
    generate_molecular_cluster,
    parse_arrangement_request
)

def handle_cluster_request(user_message: str, molecules: list):
    """
    Handle natural language clustering request from MCP.
    """
    # Parse NL to spec
    spec = parse_arrangement_request(user_message)
    
    # Generate cluster
    result = generate_molecular_cluster(
        molecules=molecules,
        arrangement=spec['pattern'],
        distance=spec.get('distance'),
        validate=True,
        optimize=True  # For best results
    )
    
    return result
```

---

## 🎓 Why This Design?

### The Core Insight

> "Stacking is fundamentally a **constraint satisfaction problem**, not a pattern enumeration problem."

Once you have:
1. A way to **select any atom/group** → `AtomSelector`
2. A way to **express any geometric relationship** → `Constraint` classes
3. A way to **place molecules anywhere** → `MoleculePose` with formulas
4. A way to **solve constraints** → `ConstraintSolver`

...you no longer need to predict every pattern. The system becomes **antifragile**: the more unusual the request, the more valuable the generic engine becomes.

### Mathematical Foundation

Every molecular arrangement is simply a point in **6N-dimensional parameter space**:
- N molecules × 6 degrees of freedom each
- 3 translational (x, y, z)
- 3 rotational (roll, pitch, yaw)

The generators provide **parameterized curves through this space**, while constraints define **feasible regions**.

---

## 📁 Files

| File | Purpose |
|------|---------|
| `molecular_arrangement_engine.py` | Complete production engine |
| `unified_molecular_arrangement.py` | Core system with unified API |
| `README_MOLECULAR_ARRANGEMENT.md` | This documentation |

---

## 🔮 Future Extensions

1. **Crystal Symmetry**: Space group operations (P21/c, C2/c, etc.)
2. **SMARTS Patterns**: RDKit integration for precise functional group detection
3. **Energy Minimization**: UFF/PM7 optimization post-arrangement
4. **Visual Editor**: Interactive constraint specification
5. **Pattern Database**: 100+ characterized motifs from CCDC

---

## 📜 License

MIT License - Free for academic and commercial use.

---

*Synthesized from K2_plan (constraints), z2_plan (poses), and original plan (generators) to create the ultimate generic molecular arrangement system.*
