
import sys
from pathlib import Path

# Add src/python to path
sys.path.insert(0, str(Path("src/python").absolute()))

from crystal_generator import generate_crystal

print("Running reproduction...")
generate_crystal(
    composition=["Si", "Si"],
    space_group=227,
    seed=42,
    max_attempts=10
)
print("Success!")
