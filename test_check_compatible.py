
from pyxtal.symmetry import Group
g = Group(227)
print(f"Checking compatibility for [2] in SG 227")
try:
    result = g.check_compatible([2])
    print(f"Result for [2]: {result}")
except Exception as e:
    print(f"Error for [2]: {e}")

print(f"Checking compatibility for [8] in SG 227")
try:
    result = g.check_compatible([8])
    print(f"Result for [8]: {result}")
except Exception as e:
    print(f"Error for [8]: {e}")
