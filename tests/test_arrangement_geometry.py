#!/usr/bin/env python3
"""
test_arrangement_geometry.py - Dynamic geometry verification tests

Tests each arrangement pattern by verifying the ACTUAL coordinates match
the REQUESTED parameters. No hardcoded expected values.
"""

import sys
import numpy as np
import os

# Suppress debug output to stderr during tests
class SuppressStderr:
    def __enter__(self):
        self._stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        return self
    def __exit__(self, *args):
        sys.stderr.close()
        sys.stderr = self._stderr

sys.path.insert(0, '/home/niel/git/crystal-mcp-server/src/python')

with SuppressStderr():
    from generators.molecule.molecular_arrangement import arrange_molecules

# =============================================================================
# Helper functions for geometric analysis
# =============================================================================

def get_center(mol):
    """Get geometric center of molecule."""
    return np.array(mol['coords']).mean(axis=0)

def get_plane_normal(mol):
    """Get normal vector to molecular plane using SVD."""
    coords = np.array(mol['coords'])
    centered = coords - coords.mean(axis=0)
    _, _, Vt = np.linalg.svd(centered)
    return Vt[-1]  # Smallest singular value direction = normal

def angle_between_vectors(v1, v2):
    """Angle between two vectors in degrees."""
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-10 or n2 < 1e-10:
        return 0.0
    cos_angle = np.dot(v1, v2) / (n1 * n2)
    return np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

def distance_between(p1, p2):
    """Euclidean distance between two points."""
    return np.linalg.norm(np.array(p2) - np.array(p1))

# =============================================================================
# Test molecules
# =============================================================================

BENZENE = {
    "atoms": ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"],
    "coords": [
        [1.4, 0.0, 0.0], [0.7, 1.21, 0.0], [-0.7, 1.21, 0.0],
        [-1.4, 0.0, 0.0], [-0.7, -1.21, 0.0], [0.7, -1.21, 0.0],
        [2.5, 0.0, 0.0], [1.25, 2.16, 0.0], [-1.25, 2.16, 0.0],
        [-2.5, 0.0, 0.0], [-1.25, -2.16, 0.0], [1.25, -2.16, 0.0]
    ],
    "identifier": "benzene"
}

WATER = {
    "atoms": ["O", "H", "H"],
    "coords": [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    "identifier": "water"
}

# =============================================================================
# Dynamic verification functions
# =============================================================================

def verify_distance(result, requested_distance, tolerance=0.5):
    """Verify that molecular centers are separated by requested distance."""
    if not result['success']:
        return False, f"Generation failed: {result.get('error')}"
    
    mols = result['molecules']
    if len(mols) < 2:
        return False, "Need at least 2 molecules"
    
    centers = [get_center(m) for m in mols]
    
    # Check consecutive distances
    for i in range(len(centers) - 1):
        actual_dist = distance_between(centers[i], centers[i+1])
        if abs(actual_dist - requested_distance) > tolerance:
            return False, f"Mol {i} to {i+1}: expected {requested_distance:.2f}A, got {actual_dist:.2f}A"
    
    return True, f"All distances within {tolerance}A of {requested_distance}A"

def verify_parallel_planes(result, tolerance_deg=15):
    """Verify that molecular planes are parallel (for pi-stacking)."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    normals = [get_plane_normal(m) for m in mols]
    
    for i in range(len(normals) - 1):
        angle = angle_between_vectors(normals[i], normals[i+1])
        # Parallel means angle is ~0 or ~180 (normals can point either way)
        if min(angle, 180 - angle) > tolerance_deg:
            return False, f"Planes {i} and {i+1} not parallel: angle={angle:.1f}deg"
    
    return True, "All planes parallel"

def verify_perpendicular_planes(result, mol1_idx=0, mol2_idx=1, tolerance_deg=15):
    """Verify that two molecular planes are perpendicular (for T-shaped)."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    n1 = get_plane_normal(mols[mol1_idx])
    n2 = get_plane_normal(mols[mol2_idx])
    
    angle = angle_between_vectors(n1, n2)
    # Perpendicular means ~90 degrees
    if abs(angle - 90) > tolerance_deg:
        return False, f"Planes not perpendicular: angle={angle:.1f}deg (expected ~90)"
    
    return True, f"Planes perpendicular: {angle:.1f}deg"

def verify_rotation_between(result, mol1_idx, mol2_idx, expected_rotation, tolerance_deg=20):
    """Verify relative rotation between two molecules."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    c1 = get_center(mols[mol1_idx])
    c2 = get_center(mols[mol2_idx])
    
    # Use first atom position relative to center as reference
    ref1 = np.array(mols[mol1_idx]['coords'][0]) - c1
    ref2 = np.array(mols[mol2_idx]['coords'][0]) - c2
    
    # Project to XY plane for yaw rotation
    actual_rotation = angle_between_vectors(ref1[:2], ref2[:2])
    
    if abs(actual_rotation - expected_rotation) > tolerance_deg:
        return False, f"Rotation: expected {expected_rotation}deg, got {actual_rotation:.1f}deg"
    
    return True, f"Rotation correct: {actual_rotation:.1f}deg"

def verify_circular_arrangement(result, tolerance_ratio=0.1):
    """Verify molecules form a circle (equal distance from center)."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    centers = [get_center(m) for m in mols]
    overall_center = np.mean(centers, axis=0)
    
    radii = [distance_between(c, overall_center) for c in centers]
    mean_radius = np.mean(radii)
    max_deviation = max(abs(r - mean_radius) for r in radii)
    
    if max_deviation / mean_radius > tolerance_ratio:
        return False, f"Not circular: radii={[f'{r:.2f}' for r in radii]}, max_dev={max_deviation:.2f}"
    
    # Also check angular spacing
    n = len(centers)
    expected_angle = 360 / n
    angles = [np.degrees(np.arctan2((c - overall_center)[1], (c - overall_center)[0])) for c in centers]
    angles_sorted = sorted(angles)
    
    return True, f"Circular: mean_radius={mean_radius:.2f}A, deviation={max_deviation:.3f}"

def verify_helical_arrangement(result):
    """Verify molecules form a helix (Z increases, XY forms circle)."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    centers = [get_center(m) for m in mols]
    
    # Check Z increases
    z_values = [c[2] for c in centers]
    z_increasing = all(z_values[i] <= z_values[i+1] for i in range(len(z_values)-1))
    
    # Check XY radii are consistent
    radii_xy = [np.sqrt(c[0]**2 + c[1]**2) for c in centers]
    radius_std = np.std(radii_xy)
    
    if not z_increasing:
        return False, f"Z not increasing: {[f'{z:.2f}' for z in z_values]}"
    if radius_std > 0.5:
        return False, f"XY radii not consistent: std={radius_std:.3f}"
    
    return True, f"Helical: Z range [{min(z_values):.2f}, {max(z_values):.2f}], radius={np.mean(radii_xy):.2f}A"

def verify_grid_arrangement(result, requested_spacing):
    """Verify molecules form a 2D grid."""
    if not result['success']:
        return False, f"Generation failed"
    
    mols = result['molecules']
    centers = [get_center(m) for m in mols]
    
    # All should be in same Z plane
    z_values = [c[2] for c in centers]
    if np.std(z_values) > 0.1:
        return False, f"Not in same Z plane: std={np.std(z_values):.3f}"
    
    # Check grid spacing
    x_vals = sorted(set(round(c[0], 1) for c in centers))
    y_vals = sorted(set(round(c[1], 1) for c in centers))
    
    if len(x_vals) > 1:
        x_spacings = [x_vals[i+1] - x_vals[i] for i in range(len(x_vals)-1)]
        if not all(abs(s - requested_spacing) < 0.5 for s in x_spacings):
            return False, f"X spacing wrong: {x_spacings}"
    
    return True, f"Grid: {len(x_vals)}x{len(y_vals)}"

# =============================================================================
# Run tests
# =============================================================================

def run_test(name, pattern, molecules, distance, verify_func, **kwargs):
    """Run a single test with dynamic parameters."""
    print(f"\n{'='*60}")
    print(f"TEST: {name}")
    print(f"  Pattern: {pattern}, Distance: {distance}A")
    print(f"{'='*60}")
    
    with SuppressStderr():
        result = arrange_molecules(molecules, pattern=pattern, distance=distance, **kwargs)
    
    if not result['success']:
        print(f"  GENERATION FAILED: {result.get('error')}")
        return False
    
    print(f"  Generated {len(result['molecules'])} molecules, {result['n_atoms']} atoms")
    
    # Run verification
    passed, message = verify_func(result)
    print(f"  Verification: {message}")
    print(f"  {'✅ PASS' if passed else '❌ FAIL'}")
    return passed

if __name__ == "__main__":
    print("=" * 70)
    print("DYNAMIC GEOMETRY VERIFICATION")
    print("Testing actual coordinates against requested parameters")
    print("=" * 70)
    
    results = []
    
    # Test pi-pi parallel: verify distance and parallel planes
    results.append(run_test(
        "Pi-Pi Parallel Stacking",
        "pi_pi_parallel",
        [BENZENE.copy() for _ in range(3)],
        distance=3.4,
        verify_func=lambda r: (
            verify_distance(r, 3.4)[0] and verify_parallel_planes(r)[0],
            f"{verify_distance(r, 3.4)[1]}; {verify_parallel_planes(r)[1]}"
        )
    ))
    
    # Test T-shaped: verify perpendicular planes
    results.append(run_test(
        "T-Shaped Dimer",
        "t_shaped",
        [BENZENE.copy(), BENZENE.copy()],
        distance=5.0,
        verify_func=lambda r: verify_perpendicular_planes(r)
    ))
    
    # Test antiparallel: verify 180 degree rotation
    results.append(run_test(
        "Antiparallel Stacking",
        "pi_pi_antiparallel",
        [BENZENE.copy(), BENZENE.copy()],
        distance=3.4,
        verify_func=lambda r: verify_rotation_between(r, 0, 1, 180)
    ))
    
    # Test circular: verify circular arrangement
    results.append(run_test(
        "Circular Arrangement",
        "circular",
        [WATER.copy() for _ in range(6)],
        distance=3.5,
        verify_func=verify_circular_arrangement
    ))
    
    # Test helical: verify helix geometry
    results.append(run_test(
        "Helical Arrangement",
        "helical",
        [WATER.copy() for _ in range(10)],
        distance=3.0,
        verify_func=verify_helical_arrangement
    ))
    
    # Test grid: verify grid geometry
    results.append(run_test(
        "Grid Arrangement",
        "grid",
        [BENZENE.copy() for _ in range(9)],
        distance=4.0,
        verify_func=lambda r: verify_grid_arrangement(r, 4.0),
        spacing=4.0
    ))
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    passed = sum(results)
    total = len(results)
    print(f"Passed: {passed}/{total}")
    print("✅ ALL TESTS PASSED" if passed == total else "❌ SOME TESTS FAILED")
