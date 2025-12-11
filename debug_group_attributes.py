
from pyxtal.symmetry import Group
g = Group(227)
print(f"Attributes of Group(227): {dir(g)}")
print(f"Symbol: {getattr(g, 'symbol', 'N/A')}")
print(f"Hall: {getattr(g, 'hall_symbol', 'N/A')}") # Check if it exists or maybe 'hall_number'
