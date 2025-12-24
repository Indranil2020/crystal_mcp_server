
import importlib
import importlib.util
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
    if importlib.util.find_spec(loc) is None:
        print(f"Could not import {loc}")
        continue

    module = importlib.import_module(loc)
    if hasattr(module, "CHGNetCalculator"):
        print(f"FOUND in {loc}")
    else:
        print(f"Not in {loc}")

# brute force check dir of chgnet
print(f"chgnet contents: {dir(chgnet)}")
