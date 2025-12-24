
from pyxtal import pyxtal
from pyxtal.symmetry import Group

g = Group(227)
print(f"Group 227 loaded. {g.symbol}")
# Inspect if there's a compatible check function
print(dir(g))

# Try check compatibility manually if easy, otherwise we just fix the test input
from pyxtal.crystal import random_crystal
# Simulate what happens inside
