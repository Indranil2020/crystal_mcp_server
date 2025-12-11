
from chgnet.model import CHGNetCalculator
import inspect

print("Inspecting CHGNetCalculator...")
try:
    print(inspect.signature(CHGNetCalculator.__init__))
except Exception as e:
    print(f"Failed to inspect signature: {e}")

print(f"Doc: {CHGNetCalculator.__doc__}")
