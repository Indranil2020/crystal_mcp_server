
from pyxtal import pyxtal
print("Inspecting Group...")
c = pyxtal()
c.from_random(3, 227, ["Si"], [8])
print(f"Group dir: {dir(c.group)}")
print(f"Lattice Type: '{c.group.lattice_type}'")
# Check PyXtal version if needed
# Check if we can get it from number
