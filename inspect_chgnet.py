
import importlib.util
from chgnet.model import CHGNet

if importlib.util.find_spec("chgnet.model.calculator") is not None:
    from chgnet.model.calculator import CHGNetCalculator
    print("CHGNetCalculator found in chgnet.model.calculator")
else:
    print("CHGNetCalculator NOT found in chgnet.model.calculator")

model = CHGNet.load()
print(f"Model type: {type(model)}")
print(f"Model dir: {dir(model)}")
