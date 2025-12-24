
from pyxtal.symmetry import Group

g = Group(227)
print("Checking compatibility for [2] in SG 227")
if hasattr(g, "check_compatible"):
    result = g.check_compatible([2])
    print(f"Result for [2]: {result}")
else:
    print("check_compatible not available on Group")

print("Checking compatibility for [8] in SG 227")
if hasattr(g, "check_compatible"):
    result = g.check_compatible([8])
    print(f"Result for [8]: {result}")
else:
    print("check_compatible not available on Group")
