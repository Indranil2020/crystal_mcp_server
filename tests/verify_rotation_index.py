
import sys
import os

# Robustness: Ensure src/python is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(current_dir, "../src/python")
sys.path.insert(0, src_path)

from generators.molecule.molecular_arrangement import ChemicalPatterns

def test_linear_rotation():
    print("Testing linear arrangement with per-molecule rotations...")
    
    # Setup: 3 molecules
    n = 3
    distance = 10.0
    
    # Test 1: Rotate middle molecule (index 1)
    rotations_middle = [
        {'x': 0, 'y': 0, 'z': 0},
        {'x': 0, 'y': 0, 'z': 90}, 
        {'x': 0, 'y': 0, 'z': 0}
    ]
    
    pos, orient = ChemicalPatterns.linear(n, distance, rotations=rotations_middle)
    
    print(f"Result orientations for middle rotation:")
    for i, o in enumerate(orient):
        print(f"  Mol {i}: {o}")

    # Verify middle is rotated
    assert orient[1].yaw == 90, f"Middle molecule should have yaw=90, got {orient[1].yaw}"
    assert orient[0].yaw == 0, "First molecule should have yaw=0"
    assert orient[2].yaw == 0, "Last molecule should have yaw=0"
    print("✅ Test 1 (Middle Rotation) PASSED")

    # Test 2: Rotate last molecule (index 2) - checking 'one rotated' logic
    rotations_last = [
        {'x': 0, 'y': 0, 'z': 0},
        {'x': 0, 'y': 0, 'z': 0},
        {'x': 0, 'y': 0, 'z': 90}
    ]
    
    pos, orient = ChemicalPatterns.linear(n, distance, rotations=rotations_last)
    
    print(f"Result orientations for last rotation:")
    for i, o in enumerate(orient):
        print(f"  Mol {i}: {o}")
        
    assert orient[2].yaw == 90, "Last molecule should have yaw=90"
    print("✅ Test 2 (Last Rotation) PASSED")
    
if __name__ == "__main__":
    try:
        test_linear_rotation()
        print("\nALL CHECKS PASSED")
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        sys.exit(1)
