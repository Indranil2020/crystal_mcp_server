
import sys
from pathlib import Path
sys.path.insert(0, str(Path("src/python").absolute()))
from crystal_generator import generate_crystal

print("Debugging generation...")
result = generate_crystal(
    composition=["Si", "Si"],
    space_group=227,
    num_atoms=8,
    seed=42,
    max_attempts=10
)

if result["success"]:
    s = result["structure"]
    print(f"Space Group: {s['space_group']}")
    print(f"Num Atoms: {s['metadata']['natoms']}")
else:
    print(f"Failed: {result['error']}")
