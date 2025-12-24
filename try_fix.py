
from pyxtal import pyxtal

print("Running fix attempts...")
crystal = pyxtal()

attempts = [
    ("Attempt 1: species=['Si'], numIons=[8]", [8]),
    ("Attempt 2: species=['Si'], numIons=None", None),
    ("Attempt 3: species=['Si'], numIons=8", 8),
]

for label, num_ions in attempts:
    print(f"\n{label}")
    if not isinstance(num_ions, (list, tuple)):
        print("Skipped: numIons should be a list or tuple")
        continue

    crystal.from_random(
        dim=3,
        group=227,
        species=['Si'],
        numIons=num_ions,
        factor=1.0
    )
    print("Success!")
