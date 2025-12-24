
from pyxtal import pyxtal
import numpy as np

print("Inspecting PyXtal Lattice...")
c = pyxtal()
c.from_random(3, 225, ["Na", "Cl"], [4, 4])
lattice = c.lattice
print(f"Lattice type: {type(lattice)}")
print(f"Dir of lattice: {dir(lattice)}")

print("Testing get_cartesian_coords...")
if hasattr(lattice, "get_cartesian_coords"):
    print(lattice.get_cartesian_coords(np.array([0.5, 0.5, 0.5])))
else:
    print("get_cartesian_coords not available")

print("Testing para (matrix)...")
if hasattr(lattice, "matrix"):
    print(lattice.matrix)
else:
    print("matrix not available")
