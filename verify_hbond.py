import json
import subprocess
import sys
import numpy as np
import math

# Define the complex request payload directly
payload = {
    "molecules": [
        {"identifier": "acetone", "count": 1},
        {"identifier": "formamide", "count": 1}
    ],
    "intermolecular_distance": 2.8,
    "stacking": "h_bonded",
    "constraints": [
        # Explicitly requesting the H-bond geometry
        "distance(0:carbonyl, 1:amide, 2.8)",
        "angle(0:atom(C,1), 0:carbonyl, 1:amide, 170)" 
        # Note: 0:atom(C,1) refers to carbonyl carbon of acetone (C=O)
        # 0:carbonyl is the Oxygen
        # 1:amide is the Hydrogen (via our smart selector)
    ],
    "use_new_engine": True,
    "optimize": True
}

# Write payload to tmp file
with open("/tmp/h_bond_test.json", "w") as f:
    json.dump(payload, f)

print("Running molecular_cluster_generator.py with payload...")
result = subprocess.run(
    ["python3", "src/python/molecular_cluster_generator.py", "/tmp/h_bond_test.json"],
    capture_output=True,
    text=True,
    cwd="/home/niel/git/crystal-mcp-server"
)

if result.returncode != 0:
    print("Error running generator:")
    print(result.stderr)
    sys.exit(1)

# Debug prints are in stderr, JSON output in stdout
print("Generator stderr debug output:")
for line in result.stderr.splitlines():
    if "DEBUG" in line:
        print(line)

try:
    data = json.loads(result.stdout)
except json.JSONDecodeError:
    print("Failed to decode JSON output")
    print("Stdout:", result.stdout)
    sys.exit(1)

if not data.get("success"):
    print("Generator reported failure:", data.get("error"))
    sys.exit(1)

# ANALYZE GEOMETRY
print("\n--- GEOMETRY ANALYSIS ---")
atoms = data["structure"]["atoms"]
coords = [np.array(a["coords"]) for a in atoms]
elements = [a["element"] for a in atoms]
# Convert fractional to cartesian if needed, but output usually has cartesian in species or sidecar
# The generator output structure: 'atoms' list has 'cartesian' key?
# Let's inspect the first atom to see format
first = atoms[0]
if "cartesian" in first:
    coords = [np.array(a["cartesian"]) for a in atoms]
else:
    # Need lattice to convert? Typically 'cartesian' is provided in final output
    print("Warning: No cartesian coords found, using fractional directly (might be wrong if not 1x1x1 box)")

# Identify molecules
# Acetone (C3H6O) should be first ~10 atoms
# Formamide (CH3NO) should be next ~6 atoms

# Find Carbonyl Oxygen in Acetone (Mol 0)
# Acetone: C-C(=O)-C
# We need to find the O in the first molecule
mol0_indices = range(0, 10) # Ends at 9
mol1_indices = range(10, 16) # Formamide

# Find O in Mol 0
o0_idx = next(i for i in mol0_indices if elements[i] == "O")
o0_pos = coords[o0_idx]
print(f"Acetone Carbonyl Oxygen found at index {o0_idx}: {o0_pos}")

# Find Amide Hydrogen in Formamide (Mol 1)
# Formamide: H-C(=O)-N-H2
# We need an H attached to N
n1_idx = next(i for i in mol1_indices if elements[i] == "N")
n1_pos = coords[n1_idx]

h_amide_idx = -1
for i in mol1_indices:
    if elements[i] == "H":
        dist = np.linalg.norm(coords[i] - n1_pos)
        if dist < 1.2: # N-H bond
            h_amide_idx = i
            break

if h_amide_idx == -1:
    print("Error: Could not find Amide H in Formamide")
    # Just grab first H in mol 1 as fallback for calc
    h_amide_idx = next(i for i in mol1_indices if elements[i] == "H")

h1_pos = coords[h_amide_idx]
print(f"Formamide Amide Hydrogen found at index {h_amide_idx}: {h1_pos}")

# Calculate Distance
dist = np.linalg.norm(o0_pos - h1_pos)
print(f"Measured H-Bond Distance: {dist:.4f} Å (Target: 2.8 Å)")

# Check Carbonyl Carbon for Angle
# Find C double bonded to O0
c_carbonyl_idx = -1
for i in mol0_indices:
    if elements[i] == "C":
        d = np.linalg.norm(coords[i] - o0_pos)
        if 1.1 < d < 1.3:
            c_carbonyl_idx = i
            break

if c_carbonyl_idx != -1:
    c_pos = coords[c_carbonyl_idx]
    # Vector C->O
    v1 = o0_pos - c_pos
    # Vector O->H
    v2 = h1_pos - o0_pos
    
    # Angle calculation
    # We want angle C-O...H to be 170 deg
    # dot product
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    angle_rad = np.arccos(np.clip(np.dot(-v1_u, v2_u), -1.0, 1.0)) # Angle at O
    # Wait, the request was "Place carbonyl oxygen ... at 170 angle"
    # Usually implies C=O...H angle.
    # If vectors are tail-to-tail...
    # Let's calculate typical bond angle
    
    # Vector O-C
    vec_oc = c_pos - o0_pos
    # Vector O-H
    vec_oh = h1_pos - o0_pos
    
    cos_theta = np.dot(vec_oc, vec_oh) / (np.linalg.norm(vec_oc) * np.linalg.norm(vec_oh))
    angle_deg = np.degrees(np.arccos(cos_theta))
    
    print(f"Measured Angle C-O...H: {angle_deg:.2f}° (Target: 170°)")
    
    if abs(dist - 2.8) < 0.2 and abs(angle_deg - 170) < 10:
         print("SUCCESS: Geometry matches constraints!")
    else:
         print("FAILURE: Geometry deviations too large.")
else:
    print("Could not find Carbonyl Carbon")

