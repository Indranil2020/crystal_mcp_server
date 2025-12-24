
from chgnet.model import CHGNetCalculator
import inspect

print("Inspecting CHGNetCalculator...")
print(inspect.signature(CHGNetCalculator.__init__))

print(f"Doc: {CHGNetCalculator.__doc__}")
