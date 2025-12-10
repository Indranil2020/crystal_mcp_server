
from pyxtal import pyxtal
import numpy as np

print("Inspecting PyXtal Lattice...")
c = pyxtal()
c.from_random(3, 225, ["Na", "Cl"], [4, 4])
lattice = c.lattice
print(f"Lattice type: {type(lattice)}")
print(f"Dir of lattice: {dir(lattice)}")

try:
    print("Testing get_cartesian_coords...")
    print(lattice.get_cartesian_coords(np.array([0.5, 0.5, 0.5])))
except Exception as e:
    print(f"get_cartesian_coords failed: {e}")

try:
    print("Testing para (matrix)...")
    print(lattice.matrix)
except Exception as e:
    print(f"matrix access failed: {e}")
