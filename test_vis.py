
from ase.build import bulk
from ase.io import write
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

# Create a sample structure (Silicon)
si = bulk('Si', 'diamond', a=5.43)

# Test 1: PNG Generation via Matplotlib
fig, ax = plt.subplots()
plot_atoms(si, ax, radii=0.3, rotation=('10x,10y,10z'))
fig.savefig("test_viz.png")
print("PNG generated successfully.")

# Test 2: HTML/X3D Generation
write("test_viz.html", si, format='html')
print("HTML generated successfully.")

# Test 3: X3D specific
write("test_viz.x3d", si, format='x3d')
print("X3D generated successfully.")
