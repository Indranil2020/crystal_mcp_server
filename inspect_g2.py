
from ase.collections import g2
print(f"Type of g2: {type(g2)}")
print(f"Has 'names' attr: {hasattr(g2, 'names')}")
if hasattr(g2, 'names'):
    print(f"First 10 names: {g2.names[:10]}")
    print(f"'C60' in names: {'C60' in g2.names}")
print(f"'C60' in g2: {'C60' in g2}")
