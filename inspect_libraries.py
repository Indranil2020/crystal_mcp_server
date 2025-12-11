
import sys
try:
    from pymatgen.core import Structure, Lattice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # Check for specific transformation classes
    from pymatgen.transformations.standard_transformations import SupercellTransformation
    print("Pymatgen transformation found.")
    
    # Check if there's specific group-subgroup functionality
    # Usually this involves identifying relations or mapping.
    # We can check specific modules.
    import pymatgen.symmetry.groups as groups
    print("pymatgen.symmetry.groups found")
except ImportError as e:
    print(e)

try:
    import pyxtal
    print(f"PyXtal version: {pyxtal.__version__}")
    # accessing group-subgroup in PyXtal
    from pyxtal.symmetry import Group
    g = Group(227)
    print(f"Group 227 subgroups: {g.get_max_subgroup_numbers()}")
except ImportError as e:
    print(e)
