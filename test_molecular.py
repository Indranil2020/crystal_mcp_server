from pyxtal import pyxtal

print("Testing molecular crystal generation...")
crystal = pyxtal(molecular=True)
species = ["H2O"]
numIons = [4]
molecules = ["H2O"]

try:
    print("Attempt 1: Passing molecules list to from_random")
    crystal.from_random(
        dim=3,
        group=19, # P212121
        species=species,
        numIons=numIons,
        factor=1.0 # Use smaller factor for faster generation
    )
    if crystal.valid:
        print("Success! Structure generated.")
        print(f"Formula: {crystal.formula}")
    else:
        print("Failed to generate structure.")
except Exception as e:
    print(f"Error: {e}")
