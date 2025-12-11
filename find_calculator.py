
import importlib
import pkgutil
import chgnet

print("Searching for CHGNetCalculator...")

possible_locations = [
    "chgnet.model",
    "chgnet.model.calculator",
    "chgnet.calculator",
    "chgnet.ase.calculator",
    "chgnet.interface.ase"
]

for loc in possible_locations:
    try:
        module = importlib.import_module(loc)
        if hasattr(module, "CHGNetCalculator"):
            print(f"FOUND in {loc}")
        else:
            print(f"Not in {loc}")
    except ImportError:
        print(f"Could not import {loc}")

# brute force check dir of chgnet
print(f"chgnet contents: {dir(chgnet)}")
