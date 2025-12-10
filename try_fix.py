
import sys
import traceback
from pyxtal import pyxtal

print("Running fix attempts...")
crystal = pyxtal()

try:
    print("\nAttempt 1: species=['Si'], numIons=[8]")
    crystal.from_random(
        dim=3,
        group=227,
        species=['Si'],
        numIons=[8],
        factor=1.0
    )
    print("Success 1!")
except Exception:
    print("Fail 1")
    traceback.print_exc()

try:
    print("\nAttempt 2: species=['Si'], numIons=None")
    crystal.from_random(
        dim=3,
        group=227,
        species=['Si'],
        numIons=None,
        factor=1.0
    )
    print("Success 2!")
except Exception:
    print("Fail 2")
    traceback.print_exc()

try:
    print("\nAttempt 3: species=['Si'], numIons=8")
    crystal.from_random(
        dim=3,
        group=227,
        species=['Si'],
        numIons=8,
        factor=1.0
    )
    print("Success 3!")
except Exception:
    print("Fail 3")
    traceback.print_exc()
