
import importlib.util
import sys

if (
    importlib.util.find_spec("pymatgen.core") is not None
    and importlib.util.find_spec("pymatgen.symmetry.analyzer") is not None
):
    from pymatgen.core import Structure, Lattice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # Check for specific transformation classes
    if importlib.util.find_spec("pymatgen.transformations.standard_transformations") is not None:
        from pymatgen.transformations.standard_transformations import SupercellTransformation
        print("Pymatgen transformation found.")
    else:
        print("Pymatgen transformations not available.")

    # Check if there's specific group-subgroup functionality
    # Usually this involves identifying relations or mapping.
    # We can check specific modules.
    import pymatgen.symmetry.groups as groups
    print("pymatgen.symmetry.groups found")
else:
    print("pymatgen not available")

if importlib.util.find_spec("pyxtal") is not None:
    import pyxtal
    print(f"PyXtal version: {pyxtal.__version__}")
    # accessing group-subgroup in PyXtal
    if importlib.util.find_spec("pyxtal.symmetry") is not None:
        from pyxtal.symmetry import Group
        g = Group(227)
        print(f"Group 227 subgroups: {g.get_max_subgroup_numbers()}")
    else:
        print("pyxtal.symmetry not available")
else:
    print("pyxtal not available")
