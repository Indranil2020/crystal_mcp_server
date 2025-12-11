
from chgnet.model import CHGNet
try:
    from chgnet.model.calculator import CHGNetCalculator
    print("CHGNetCalculator found in chgnet.model.calculator")
except ImportError:
    print("CHGNetCalculator NOT found in chgnet.model.calculator")
    
model = CHGNet.load()
print(f"Model type: {type(model)}")
print(f"Model dir: {dir(model)}")
