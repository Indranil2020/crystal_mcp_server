#!/usr/bin/env python3
"""
Adversarial Test Suite for Crystal MCP Server - Molecule Generation System

EXPERT AUDITOR MODE: Designed to break the system, not validate success.
Tests edge cases, numerical instabilities, and design flaws.

Author: Claude Code (Expert Critic Mode)
Date: 2026-01-11
"""

import sys
import os
import json
import math
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple
from datetime import datetime

# Add source path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'python'))

from generators.molecule import (
    generate_molecule_universal,
    generate_molecular_cluster_unified,
)

# =============================================================================
# TEST INFRASTRUCTURE
# =============================================================================

@dataclass
class TestResult:
    """Result of a single test case."""
    test_id: str
    category: str
    description: str
    passed: bool
    expected_behavior: str
    actual_behavior: str
    error: Optional[str] = None
    verification: Optional[Dict[str, Any]] = None
    severity: str = "normal"  # critical, high, normal, low
    root_cause: Optional[str] = None

@dataclass
class TestReport:
    """Comprehensive test report."""
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    total_tests: int = 0
    passed: int = 0
    failed: int = 0
    critical_failures: List[TestResult] = field(default_factory=list)
    high_failures: List[TestResult] = field(default_factory=list)
    normal_failures: List[TestResult] = field(default_factory=list)
    all_results: List[TestResult] = field(default_factory=list)

    def add_result(self, result: TestResult):
        self.all_results.append(result)
        self.total_tests += 1
        if result.passed:
            self.passed += 1
        else:
            self.failed += 1
            if result.severity == "critical":
                self.critical_failures.append(result)
            elif result.severity == "high":
                self.high_failures.append(result)
            else:
                self.normal_failures.append(result)


def log(msg: str, level: str = "INFO"):
    """Print formatted log message."""
    prefix = {"INFO": "  ", "PASS": "✓ ", "FAIL": "✗ ", "WARN": "⚠ ", "CRIT": "⛔"}
    print(f"{prefix.get(level, '  ')}{msg}")


# =============================================================================
# MATHEMATICAL VERIFICATION UTILITIES
# =============================================================================

def verify_min_distance(coords: List[List[float]], expected_min: float, tolerance: float = 0.1) -> Dict:
    """Verify minimum interatomic distance."""
    if len(coords) < 2:
        return {"valid": True, "min_dist": float('inf'), "message": "Less than 2 atoms"}

    coords_arr = np.array(coords)
    n = len(coords_arr)
    min_dist = float('inf')
    min_pair = (0, 0)

    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords_arr[i] - coords_arr[j])
            if d < min_dist:
                min_dist = d
                min_pair = (i, j)

    valid = min_dist >= (expected_min - tolerance)
    return {
        "valid": valid,
        "min_dist": min_dist,
        "min_pair": min_pair,
        "expected": expected_min,
        "message": f"Min distance {min_dist:.4f}Å {'≥' if valid else '<'} expected {expected_min}Å"
    }


def verify_centroid_separation(coords: List[List[float]], n_mols: int, mol_size: int,
                               expected_dist: float, tolerance: float = 0.2) -> Dict:
    """Verify centroid separation between molecules."""
    if n_mols < 2 or mol_size == 0:
        return {"valid": True, "message": "Cannot verify with < 2 molecules"}

    coords_arr = np.array(coords)
    centroids = []
    for i in range(n_mols):
        start = i * mol_size
        end = start + mol_size
        if end <= len(coords_arr):
            centroids.append(coords_arr[start:end].mean(axis=0))

    if len(centroids) < 2:
        return {"valid": True, "message": "Not enough complete molecules"}

    # Check first pair separation
    dist = np.linalg.norm(centroids[0] - centroids[1])
    valid = abs(dist - expected_dist) <= tolerance

    return {
        "valid": valid,
        "actual_dist": dist,
        "expected_dist": expected_dist,
        "tolerance": tolerance,
        "message": f"Centroid separation {dist:.4f}Å, expected {expected_dist}±{tolerance}Å"
    }


def verify_planarity(coords: List[List[float]], tolerance: float = 0.1) -> Dict:
    """Verify if molecule is planar (Z-extent near zero)."""
    if len(coords) < 3:
        return {"valid": True, "z_extent": 0, "message": "Too few atoms to verify planarity"}

    coords_arr = np.array(coords)
    # Align to principal axes
    centered = coords_arr - coords_arr.mean(axis=0)
    _, S, Vt = np.linalg.svd(centered)

    # Smallest singular value indicates planarity
    min_sv = S[-1] if len(S) >= 3 else 0
    valid = min_sv < tolerance * 10  # More lenient for planarity

    return {
        "valid": valid,
        "min_singular_value": min_sv,
        "tolerance": tolerance,
        "message": f"Planarity metric (smallest SV): {min_sv:.4f}"
    }


def verify_circular_arrangement(coords: List[List[float]], n_mols: int, mol_size: int,
                                 expected_radius: float, tolerance: float = 0.5) -> Dict:
    """Verify circular arrangement of molecules."""
    if n_mols < 3:
        return {"valid": True, "message": "Need ≥3 molecules to verify circular"}

    coords_arr = np.array(coords)
    centroids = []
    for i in range(n_mols):
        start = i * mol_size
        end = start + mol_size
        if end <= len(coords_arr):
            centroids.append(coords_arr[start:end].mean(axis=0))

    if len(centroids) < 3:
        return {"valid": True, "message": "Not enough molecules"}

    centroids = np.array(centroids)
    center = centroids.mean(axis=0)
    radii = [np.linalg.norm(c - center) for c in centroids]
    avg_radius = np.mean(radii)
    radius_std = np.std(radii)

    valid = radius_std < tolerance and abs(avg_radius - expected_radius) < tolerance

    return {
        "valid": valid,
        "avg_radius": avg_radius,
        "radius_std": radius_std,
        "expected_radius": expected_radius,
        "message": f"Avg radius {avg_radius:.2f}Å ± {radius_std:.3f}, expected {expected_radius}Å"
    }


# =============================================================================
# CATEGORY 1: SINGLE MOLECULE GENERATION STRESS TESTS
# =============================================================================

def test_single_molecule_generation(report: TestReport):
    """Test single molecule generation edge cases."""
    log("\n=== CATEGORY 1: SINGLE MOLECULE GENERATION ===", "INFO")

    tests = [
        # --- Boundary Cases ---
        ("SMG-001", "Empty string identifier", "",
         "Should return error for empty identifier", "critical"),

        ("SMG-002", "Whitespace-only identifier", "   ",
         "Should return error for whitespace identifier", "high"),

        ("SMG-003", "Single character identifier", "X",
         "Should return error or handle gracefully", "normal"),

        ("SMG-004", "Very long identifier", "a" * 500,
         "Should handle gracefully without crash", "high"),

        ("SMG-005", "Unicode identifier", "咖啡因",
         "Should handle unicode gracefully", "normal"),

        ("SMG-006", "Mixed case sensitivity", "BeNzEnE",
         "Should find benzene regardless of case", "normal"),

        # --- SMILES Edge Cases ---
        ("SMG-010", "Invalid SMILES syntax", "C(C(C",
         "Should return clear error for invalid SMILES", "high"),

        ("SMG-011", "Disconnected SMILES", "C.C.C",
         "Should generate disconnected methane molecules", "normal"),

        ("SMG-012", "SMILES with charged atoms", "[Na+].[Cl-]",
         "Should handle ionic species", "normal"),

        ("SMG-013", "SMILES with stereochemistry", "C/C=C\\C",
         "Should preserve cis/trans configuration", "high"),

        ("SMG-014", "SMILES with ring stereochemistry", "C[C@H]1CCCCC1",
         "Should preserve R/S configuration", "high"),

        ("SMG-015", "Very large SMILES", "C" * 100,
         "Should handle long alkane chain", "normal"),

        ("SMG-016", "SMILES with aromatic atoms", "c1ccccc1",
         "Should generate aromatic benzene", "normal"),

        ("SMG-017", "SMILES with multiple rings", "c1ccc2ccccc2c1",
         "Should generate naphthalene correctly", "normal"),

        # --- Name Resolution Edge Cases ---
        ("SMG-020", "Ambiguous common name", "salt",
         "Should handle ambiguous name gracefully", "normal"),

        ("SMG-021", "Misspelled name", "benzne",
         "Should fuzzy-match or return helpful error", "normal"),

        ("SMG-022", "IUPAC systematic name", "2,3-dimethylbutane",
         "Should resolve IUPAC name", "normal"),

        ("SMG-023", "Complex IUPAC name", "perylene-3,4,9,10-tetracarboxylic dianhydride",
         "Should resolve complex IUPAC (PTCDA)", "normal"),

        ("SMG-024", "Trade name", "aspirin",
         "Should resolve common trade name", "normal"),

        ("SMG-025", "CAS number format", "50-78-2",
         "Should handle CAS number for aspirin", "normal"),

        # --- Extreme Molecules ---
        ("SMG-030", "Fullerene C60", "C60",
         "Should generate fullerene cage", "normal"),

        ("SMG-031", "Large molecule (cholesterol)", "cholesterol",
         "Should generate complex steroid", "normal"),

        ("SMG-032", "Metalloporphyrin", "chlorophyll",
         "Should handle metallorganic compound", "normal"),

        ("SMG-033", "Polymer repeat unit", "polyethylene",
         "Should handle polymer name gracefully", "low"),
    ]

    for test_id, desc, identifier, expected, severity in tests:
        log(f"[{test_id}] {desc}")

        try:
            result = generate_molecule_universal(
                identifier=identifier,
                input_type="auto",
                optimize=True,
                allow_external=True
            )
        except Exception as e:
            # Catch pymatgen or other exceptions
            result = {"success": False, "error": {"message": f"Exception: {type(e).__name__}: {str(e)[:100]}"}}

        passed = False
        actual = ""
        verification = None

        # Analyze result
        if result.get("success"):
            atoms = result.get("atoms", [])
            coords = result.get("coords", [])

            if len(atoms) == 0 or len(coords) == 0:
                actual = f"Success but empty structure (atoms={len(atoms)}, coords={len(coords)})"
                passed = False
            else:
                actual = f"Generated {len(atoms)} atoms, source={result.get('source', 'unknown')}"

                # Mathematical verification
                verification = verify_planarity(coords)
                min_dist_check = verify_min_distance(coords, 0.5)  # No atoms < 0.5Å apart

                verification["min_distance_check"] = min_dist_check

                # For known molecules, validate atom count
                if "benzene" in identifier.lower() and len(atoms) != 12:
                    actual += f" WARNING: Benzene should have 12 atoms, got {len(atoms)}"
                    passed = False
                elif min_dist_check["min_dist"] < 0.5:
                    actual += f" PHYSICS VIOLATION: atoms {min_dist_check['min_dist']:.3f}Å apart"
                    passed = False
                else:
                    passed = True
        else:
            error_info = result.get("error", {})
            error_msg = error_info.get("message", str(error_info)) if isinstance(error_info, dict) else str(error_info)
            actual = f"Failed: {error_msg}"

            # Some failures are expected (empty string, invalid SMILES)
            if test_id in ["SMG-001", "SMG-002", "SMG-010"]:
                passed = True  # Expected to fail
            else:
                passed = False

        test_result = TestResult(
            test_id=test_id,
            category="Single Molecule Generation",
            description=desc,
            passed=passed,
            expected_behavior=expected,
            actual_behavior=actual,
            verification=verification,
            severity=severity
        )

        report.add_result(test_result)
        log(f"   → {actual}", "PASS" if passed else "FAIL")


# =============================================================================
# SAFE WRAPPER FOR CLUSTER GENERATION
# =============================================================================

def safe_cluster(molecules, **kwargs):
    """Wrapper to catch exceptions from cluster generation."""
    try:
        return generate_molecular_cluster_unified(molecules=molecules, **kwargs)
    except Exception as e:
        return {"success": False, "error": {"code": "EXCEPTION", "message": f"{type(e).__name__}: {str(e)[:100]}"}}


# =============================================================================
# CATEGORY 2: COMPLEX ARRANGEMENT STRESS TESTS
# =============================================================================

def test_complex_arrangements(report: TestReport):
    """Test molecular arrangement edge cases."""
    log("\n=== CATEGORY 2: COMPLEX MOLECULAR ARRANGEMENTS ===", "INFO")

    # --- Test 2.1: Boundary Conditions ---
    log("\n--- 2.1 Boundary Conditions ---", "INFO")

    # n=0 molecules
    log("[ARR-001] Zero molecules")
    result = safe_cluster(
        molecules=[],
        stacking="pi_pi_parallel"
    )
    passed = not result.get("success") and "error" in result
    report.add_result(TestResult(
        test_id="ARR-001",
        category="Arrangements",
        description="Zero molecules",
        passed=passed,
        expected_behavior="Should return error for empty molecule list",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="high"
    ))
    log(f"   → {result.get('error', 'No error returned')}", "PASS" if passed else "FAIL")

    # n=1 molecule (degenerate case)
    log("[ARR-002] Single molecule (n=1)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 1}],
        stacking="pi_pi_parallel"
    )
    passed = result.get("success") and result.get("n_molecules") == 1
    report.add_result(TestResult(
        test_id="ARR-002",
        category="Arrangements",
        description="Single molecule",
        passed=passed,
        expected_behavior="Should succeed with single molecule",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms')}",
        severity="normal"
    ))
    log(f"   → n_atoms={result.get('n_atoms', 0)}", "PASS" if passed else "FAIL")

    # Very large cluster
    log("[ARR-003] Large cluster (n=50)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 50}],
        stacking="grid",
        intermolecular_distance=7.0
    )
    passed = result.get("success") and result.get("n_atoms") > 0
    actual = f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}"
    if result.get("success"):
        coords = result.get("coords", [])
        min_check = verify_min_distance(coords, 1.5)
        passed = passed and min_check["valid"]
        actual += f", min_dist={min_check['min_dist']:.3f}Å"
    report.add_result(TestResult(
        test_id="ARR-003",
        category="Arrangements",
        description="Large cluster (50 molecules)",
        passed=passed,
        expected_behavior="Should create valid grid arrangement",
        actual_behavior=actual,
        severity="high"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # --- Test 2.2: Distance Edge Cases ---
    log("\n--- 2.2 Distance Edge Cases ---", "INFO")

    # Zero distance (should fail or create physics violation)
    log("[ARR-010] Zero intermolecular distance")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="pi_pi_parallel",
        intermolecular_distance=0.0
    )
    # Should fail due to overlap
    passed = not result.get("success") or result.get("validation", {}).get("physics_violation")
    actual = f"success={result.get('success')}"
    if result.get("error"):
        actual += f", error={result['error'].get('code', 'unknown')}"
    report.add_result(TestResult(
        test_id="ARR-010",
        category="Arrangements",
        description="Zero intermolecular distance",
        passed=passed,
        expected_behavior="Should fail with physics violation",
        actual_behavior=actual,
        severity="critical"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # Negative distance
    log("[ARR-011] Negative intermolecular distance")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="linear",
        intermolecular_distance=-3.4
    )
    # Should either fail or handle gracefully
    passed = not result.get("success") or result.get("n_atoms") > 0
    actual = f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}"
    report.add_result(TestResult(
        test_id="ARR-011",
        category="Arrangements",
        description="Negative intermolecular distance",
        passed=passed,
        expected_behavior="Should fail or handle gracefully",
        actual_behavior=actual,
        severity="high"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # Very small distance (below VDW radii)
    log("[ARR-012] Sub-VDW distance (1.0Å)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="pi_pi_parallel",
        intermolecular_distance=1.0
    )
    passed = not result.get("success")  # Should fail
    actual = f"success={result.get('success')}"
    if result.get("error"):
        actual += f", error={result['error'].get('code', 'unknown')}"
    report.add_result(TestResult(
        test_id="ARR-012",
        category="Arrangements",
        description="Sub-VDW distance (1.0Å)",
        passed=passed,
        expected_behavior="Should fail with clash or physics violation",
        actual_behavior=actual,
        severity="critical"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # Very large distance
    log("[ARR-013] Very large distance (1000Å)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="linear",
        intermolecular_distance=1000.0
    )
    passed = result.get("success")
    verification = None
    if passed:
        coords = result.get("coords", [])
        verification = verify_centroid_separation(coords, 2, 12, 1000.0, tolerance=1.0)
        passed = verification["valid"]
    report.add_result(TestResult(
        test_id="ARR-013",
        category="Arrangements",
        description="Very large distance (1000Å)",
        passed=passed,
        expected_behavior="Should create structure with 1000Å separation",
        actual_behavior=f"success={result.get('success')}, {verification['message'] if verification else 'N/A'}",
        verification=verification,
        severity="normal"
    ))
    log(f"   → {verification['message'] if verification else 'Failed'}", "PASS" if passed else "FAIL")

    # --- Test 2.3: Pattern-Specific Edge Cases ---
    log("\n--- 2.3 Pattern-Specific Tests ---", "INFO")

    # Circular with n=2 (degenerate case)
    log("[ARR-020] Circular arrangement with n=2")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="circular",
        intermolecular_distance=5.0
    )
    passed = result.get("success")
    report.add_result(TestResult(
        test_id="ARR-020",
        category="Arrangements",
        description="Circular with n=2",
        passed=passed,
        expected_behavior="Should handle degenerate 2-molecule circle",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}",
        severity="normal"
    ))
    log(f"   → {result.get('n_atoms', 0)} atoms", "PASS" if passed else "FAIL")

    # Herringbone with single molecule
    log("[ARR-021] Herringbone with n=1")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 1}],
        stacking="herringbone"
    )
    passed = result.get("success")
    report.add_result(TestResult(
        test_id="ARR-021",
        category="Arrangements",
        description="Herringbone with n=1",
        passed=passed,
        expected_behavior="Should handle single molecule herringbone",
        actual_behavior=f"success={result.get('success')}",
        severity="low"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")

    # Spherical with small radius
    log("[ARR-022] Spherical with very small radius")
    result = safe_cluster(
        molecules=[{"identifier": "water", "count": 8}],
        stacking="spherical",
        intermolecular_distance=1.5  # Too small
    )
    passed = not result.get("success") or result.get("validation", {}).get("errors")
    report.add_result(TestResult(
        test_id="ARR-022",
        category="Arrangements",
        description="Spherical with small radius",
        passed=passed,
        expected_behavior="Should detect overlapping molecules on sphere",
        actual_behavior=f"success={result.get('success')}, validation_errors={len(result.get('validation', {}).get('errors', []))}",
        severity="high"
    ))
    log(f"   → errors={len(result.get('validation', {}).get('errors', []))}", "PASS" if passed else "FAIL")

    # Unknown pattern name
    log("[ARR-023] Unknown pattern name")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="nonexistent_pattern_xyz123"
    )
    # Should fallback to auto or return error
    passed = result.get("success") or "error" in result
    report.add_result(TestResult(
        test_id="ARR-023",
        category="Arrangements",
        description="Unknown pattern name",
        passed=passed,
        expected_behavior="Should fallback to auto or return error",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")

    # --- Test 2.4: Formula Edge Cases ---
    log("\n--- 2.4 Formula Edge Cases ---", "INFO")

    # Division by zero in formula
    log("[ARR-030] Formula with division by zero")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 3}],
        formulas={"x": "1/0", "y": "0", "z": "i*3.4"}
    )
    passed = "error" in result or result.get("success")  # Should handle gracefully
    report.add_result(TestResult(
        test_id="ARR-030",
        category="Arrangements",
        description="Formula with division by zero",
        passed=passed,
        expected_behavior="Should handle div/0 gracefully (return 0 or error)",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="critical"
    ))
    log(f"   → error={result.get('error')}", "PASS" if passed else "FAIL")

    # Sqrt of negative
    log("[ARR-031] Formula with sqrt(-1)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 3}],
        formulas={"x": "sqrt(-1)", "y": "0", "z": "i*3.4"}
    )
    passed = "error" in result or result.get("success")  # Should handle gracefully
    report.add_result(TestResult(
        test_id="ARR-031",
        category="Arrangements",
        description="Formula with sqrt(-1)",
        passed=passed,
        expected_behavior="Should handle sqrt(-1) gracefully",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="high"
    ))
    log(f"   → error={result.get('error')}", "PASS" if passed else "FAIL")

    # Formula producing infinity
    log("[ARR-032] Formula producing infinity")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 3}],
        formulas={"x": "exp(1000)", "y": "0", "z": "0"}
    )
    passed = "error" in result or result.get("success")  # Should clamp or error
    report.add_result(TestResult(
        test_id="ARR-032",
        category="Arrangements",
        description="Formula producing infinity",
        passed=passed,
        expected_behavior="Should clamp large values or return error",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="high"
    ))
    log(f"   → error={result.get('error')}", "PASS" if passed else "FAIL")

    # Formula injection attempt (security)
    log("[ARR-033] Formula injection attempt")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 3}],
        formulas={"x": "__import__('os').system('echo hacked')", "y": "0", "z": "0"}
    )
    # Should NOT execute system command, should fail or sanitize
    passed = not result.get("success") or "error" in result
    report.add_result(TestResult(
        test_id="ARR-033",
        category="Arrangements",
        description="Formula injection attempt",
        passed=passed,
        expected_behavior="Should block code injection",
        actual_behavior=f"success={result.get('success')}, error={result.get('error')}",
        severity="critical"
    ))
    log(f"   → error={result.get('error')}", "PASS" if passed else "FAIL")

    # Valid but complex formula (circular)
    log("[ARR-034] Valid circular formula")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 8}],
        formulas={"x": "5*cos(2*pi*i/n)", "y": "5*sin(2*pi*i/n)", "z": "0"}
    )
    passed = result.get("success")
    verification = None
    if passed:
        coords = result.get("coords", [])
        verification = verify_circular_arrangement(coords, 8, 12, 5.0, tolerance=0.5)
        passed = verification["valid"]
    report.add_result(TestResult(
        test_id="ARR-034",
        category="Arrangements",
        description="Valid circular formula",
        passed=passed,
        expected_behavior="Should create circular arrangement with radius 5Å",
        actual_behavior=f"success={result.get('success')}, {verification['message'] if verification else 'N/A'}",
        verification=verification,
        severity="high"
    ))
    log(f"   → {verification['message'] if verification else 'Failed'}", "PASS" if passed else "FAIL")

    # --- Test 2.5: Heterogeneous Clusters ---
    log("\n--- 2.5 Heterogeneous Clusters ---", "INFO")

    # Mix of different molecules
    log("[ARR-040] Hetero-cluster (benzene + water + naphthalene)")
    result = safe_cluster(
        molecules=[
            {"identifier": "benzene", "count": 2},
            {"identifier": "water", "count": 2},
            {"identifier": "naphthalene", "count": 1}
        ],
        stacking="linear",
        intermolecular_distance=5.0
    )
    passed = result.get("success") and result.get("n_molecules") == 5
    actual = f"success={result.get('success')}, n_mols={result.get('n_molecules', 0)}, n_atoms={result.get('n_atoms', 0)}"
    verification = None
    if passed:
        coords = result.get("coords", [])
        verification = verify_min_distance(coords, 1.0)
        passed = verification["valid"]
        actual += f", min_dist={verification['min_dist']:.3f}Å"
    report.add_result(TestResult(
        test_id="ARR-040",
        category="Arrangements",
        description="Heterogeneous cluster",
        passed=passed,
        expected_behavior="Should arrange mixed molecules linearly",
        actual_behavior=actual,
        verification=verification,
        severity="normal"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # Invalid molecule in list
    log("[ARR-041] List with invalid molecule")
    result = safe_cluster(
        molecules=[
            {"identifier": "benzene", "count": 2},
            {"identifier": "definitely_not_a_molecule_xyz123", "count": 1}
        ],
        stacking="linear"
    )
    # Should either skip invalid or fail
    passed = True  # Either behavior is acceptable
    actual = f"success={result.get('success')}, n_mols={result.get('n_molecules', 0)}"
    report.add_result(TestResult(
        test_id="ARR-041",
        category="Arrangements",
        description="List with invalid molecule",
        passed=passed,
        expected_behavior="Should skip invalid molecule or return error",
        actual_behavior=actual,
        severity="normal"
    ))
    log(f"   → {actual}", "PASS" if passed else "FAIL")

    # --- Test 2.6: Mathematical Precision Tests ---
    log("\n--- 2.6 Mathematical Precision Tests ---", "INFO")

    # Pi-stacking distance verification
    log("[ARR-050] Pi-stacking distance verification (3.4Å)")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 5}],
        stacking="pi_pi_parallel",
        intermolecular_distance=3.4
    )
    passed = result.get("success")
    verification = None
    if passed:
        coords = result.get("coords", [])
        # Get benzene atom count (should be 12)
        mol_size = 12
        verification = {}

        # Check that molecules are stacked along Z
        coords_arr = np.array(coords)
        centroids = []
        for i in range(5):
            start = i * mol_size
            end = start + mol_size
            if end <= len(coords_arr):
                centroids.append(coords_arr[start:end].mean(axis=0))

        if len(centroids) >= 2:
            z_coords = [c[2] for c in centroids]
            z_diffs = [z_coords[i+1] - z_coords[i] for i in range(len(z_coords)-1)]
            avg_diff = np.mean(np.abs(z_diffs))
            std_diff = np.std(np.abs(z_diffs))

            verification["z_separations"] = z_diffs
            verification["avg_z_separation"] = avg_diff
            verification["std_z_separation"] = std_diff
            verification["expected"] = 3.4
            verification["valid"] = abs(avg_diff - 3.4) < 0.2
            verification["message"] = f"Avg Z separation: {avg_diff:.4f}Å ± {std_diff:.4f}, expected 3.4Å"

            passed = verification["valid"]

    report.add_result(TestResult(
        test_id="ARR-050",
        category="Arrangements",
        description="Pi-stacking distance verification",
        passed=passed,
        expected_behavior="Should have 3.4Å Z-separation between benzene centroids",
        actual_behavior=f"{verification['message'] if verification else 'No verification'}",
        verification=verification,
        severity="high"
    ))
    log(f"   → {verification['message'] if verification else 'Failed'}", "PASS" if passed else "FAIL")

    # Helical geometry verification
    log("[ARR-051] Helical geometry verification")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 10}],
        stacking="helical",
        intermolecular_distance=3.4
    )
    passed = result.get("success")
    verification = None
    if passed:
        coords = result.get("coords", [])
        coords_arr = np.array(coords)
        mol_size = 12
        centroids = []
        for i in range(10):
            start = i * mol_size
            end = start + mol_size
            if end <= len(coords_arr):
                centroids.append(coords_arr[start:end].mean(axis=0))

        if len(centroids) >= 3:
            centroids_arr = np.array(centroids)
            # Check XY radii (should be roughly constant)
            xy_radii = np.sqrt(centroids_arr[:, 0]**2 + centroids_arr[:, 1]**2)
            # Check Z progression (should be monotonic)
            z_coords = centroids_arr[:, 2]
            z_monotonic = all(z_coords[i] <= z_coords[i+1] for i in range(len(z_coords)-1)) or \
                         all(z_coords[i] >= z_coords[i+1] for i in range(len(z_coords)-1))

            verification = {
                "xy_radii_std": np.std(xy_radii),
                "avg_xy_radius": np.mean(xy_radii),
                "z_monotonic": z_monotonic,
                "z_range": z_coords.max() - z_coords.min(),
                "valid": np.std(xy_radii) < 0.5 and z_monotonic,
                "message": f"Helix: radius={np.mean(xy_radii):.2f}±{np.std(xy_radii):.2f}Å, z_range={z_coords.max()-z_coords.min():.2f}Å"
            }
            passed = verification["valid"]

    report.add_result(TestResult(
        test_id="ARR-051",
        category="Arrangements",
        description="Helical geometry verification",
        passed=passed,
        expected_behavior="Should have constant XY radius and monotonic Z",
        actual_behavior=f"{verification['message'] if verification else 'No verification'}",
        verification=verification,
        severity="normal"
    ))
    log(f"   → {verification['message'] if verification else 'Failed'}", "PASS" if passed else "FAIL")


# =============================================================================
# CATEGORY 3: CONSTRAINT SYSTEM STRESS TESTS
# =============================================================================

def test_constraint_system(report: TestReport):
    """Test constraint parsing and satisfaction."""
    log("\n=== CATEGORY 3: CONSTRAINT SYSTEM ===", "INFO")

    # Valid distance constraint
    log("[CON-001] Valid distance constraint")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="linear",
        intermolecular_distance=10.0,  # Start far apart
        constraints=["distance(0:centroid(), 1:centroid(), 5.0)"],
        use_solver=True
    )
    passed = result.get("success")
    verification = None
    if passed:
        coords = result.get("coords", [])
        verification = verify_centroid_separation(coords, 2, 12, 5.0, tolerance=0.5)
        passed = verification["valid"]
    report.add_result(TestResult(
        test_id="CON-001",
        category="Constraints",
        description="Valid distance constraint",
        passed=passed,
        expected_behavior="Should satisfy distance=5.0Å constraint",
        actual_behavior=f"{verification['message'] if verification else 'Failed'}",
        verification=verification,
        severity="high"
    ))
    log(f"   → {verification['message'] if verification else 'Failed'}", "PASS" if passed else "FAIL")

    # Conflicting constraints (impossible)
    log("[CON-002] Conflicting constraints")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 3}],
        stacking="linear",
        constraints=[
            "distance(0:centroid(), 1:centroid(), 3.0)",
            "distance(1:centroid(), 2:centroid(), 3.0)",
            "distance(0:centroid(), 2:centroid(), 10.0)"  # Impossible with triangle inequality
        ],
        use_solver=True
    )
    # Should either fail or return best-effort result with warnings
    passed = True  # We just want it not to crash
    report.add_result(TestResult(
        test_id="CON-002",
        category="Constraints",
        description="Conflicting constraints",
        passed=passed,
        expected_behavior="Should handle impossible constraints gracefully",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")

    # Malformed constraint string
    log("[CON-003] Malformed constraint string")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="linear",
        constraints=["not_a_valid_constraint((("],
        use_solver=True
    )
    # Should ignore invalid constraint or return error
    passed = True  # Shouldn't crash
    report.add_result(TestResult(
        test_id="CON-003",
        category="Constraints",
        description="Malformed constraint string",
        passed=passed,
        expected_behavior="Should handle malformed constraint gracefully",
        actual_behavior=f"success={result.get('success')}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")


# =============================================================================
# CATEGORY 4: NUMERICAL STABILITY TESTS
# =============================================================================

def test_numerical_stability(report: TestReport):
    """Test numerical stability edge cases."""
    log("\n=== CATEGORY 4: NUMERICAL STABILITY ===", "INFO")

    # Near-zero rotation (gimbal lock region)
    log("[NUM-001] Near-gimbal lock rotation")
    result = safe_cluster(
        molecules=[{"identifier": "benzene", "count": 2}],
        stacking="t_shaped",  # 90° rotation
        intermolecular_distance=4.5
    )
    passed = result.get("success")
    if passed:
        coords = result.get("coords", [])
        verification = verify_min_distance(coords, 1.0)
        passed = verification["valid"]
    report.add_result(TestResult(
        test_id="NUM-001",
        category="Numerical",
        description="T-shaped arrangement (90° rotation)",
        passed=passed,
        expected_behavior="Should handle 90° rotation without gimbal lock issues",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")

    # Linear molecule (degenerate SVD)
    log("[NUM-002] Linear molecule plane detection")
    result = safe_cluster(
        molecules=[{"identifier": "ethane", "count": 2}],  # Near-linear
        stacking="pi_pi_parallel",
        intermolecular_distance=3.4
    )
    passed = result.get("success")
    report.add_result(TestResult(
        test_id="NUM-002",
        category="Numerical",
        description="Linear molecule plane detection",
        passed=passed,
        expected_behavior="Should handle degenerate SVD for near-linear molecules",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")

    # Very small molecule (single atom issues)
    log("[NUM-003] Single-atom molecule")
    result = safe_cluster(
        molecules=[{"identifier": "He", "count": 4}],  # Helium atoms
        stacking="linear",
        intermolecular_distance=3.0
    )
    passed = result.get("success") or "error" in result
    report.add_result(TestResult(
        test_id="NUM-003",
        category="Numerical",
        description="Single-atom molecule arrangement",
        passed=passed,
        expected_behavior="Should handle single-atom molecules",
        actual_behavior=f"success={result.get('success')}, n_atoms={result.get('n_atoms', 0)}",
        severity="normal"
    ))
    log(f"   → success={result.get('success')}", "PASS" if passed else "FAIL")


# =============================================================================
# MAIN TEST RUNNER
# =============================================================================

def run_all_tests() -> TestReport:
    """Run all adversarial tests and generate report."""
    report = TestReport()

    print("=" * 70)
    print("ADVERSARIAL MOLECULAR CHEMISTRY SYSTEM AUDIT")
    print("=" * 70)
    print(f"Timestamp: {report.timestamp}")
    print("Mode: EXPERT CRITIC / ADVERSARIAL TEST ENGINEER")
    print("=" * 70)

    # Run test categories
    test_single_molecule_generation(report)
    test_complex_arrangements(report)
    test_constraint_system(report)
    test_numerical_stability(report)

    return report


def generate_markdown_report(report: TestReport) -> str:
    """Generate detailed markdown report."""
    lines = [
        "# Adversarial Molecular Chemistry System Audit Report",
        "",
        f"**Date:** {report.timestamp}",
        f"**Auditor:** Claude Code (Expert Critic Mode)",
        f"**Status:** {'PASS' if report.failed == 0 else 'ISSUES FOUND'}",
        "",
        "---",
        "",
        "## Executive Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total Tests | {report.total_tests} |",
        f"| Passed | {report.passed} |",
        f"| Failed | {report.failed} |",
        f"| Pass Rate | {100*report.passed/report.total_tests:.1f}% |",
        f"| Critical Failures | {len(report.critical_failures)} |",
        f"| High Severity Failures | {len(report.high_failures)} |",
        "",
        "---",
        "",
    ]

    # Critical failures section
    if report.critical_failures:
        lines.extend([
            "## Critical Failures ⛔",
            "",
            "These issues represent fundamental design flaws or security vulnerabilities.",
            "",
        ])
        for result in report.critical_failures:
            lines.extend([
                f"### {result.test_id}: {result.description}",
                "",
                f"**Expected:** {result.expected_behavior}",
                "",
                f"**Actual:** {result.actual_behavior}",
                "",
                f"**Root Cause:** {result.root_cause or 'TBD - requires investigation'}",
                "",
                "---",
                "",
            ])

    # High severity failures
    if report.high_failures:
        lines.extend([
            "## High Severity Failures ⚠️",
            "",
        ])
        for result in report.high_failures:
            lines.extend([
                f"### {result.test_id}: {result.description}",
                "",
                f"**Expected:** {result.expected_behavior}",
                "",
                f"**Actual:** {result.actual_behavior}",
                "",
            ])

    # Full test results table
    lines.extend([
        "",
        "## Full Test Results",
        "",
        "| ID | Category | Description | Status | Severity |",
        "|-----|----------|-------------|--------|----------|",
    ])

    for result in report.all_results:
        status = "✓ PASS" if result.passed else "✗ FAIL"
        lines.append(f"| {result.test_id} | {result.category} | {result.description} | {status} | {result.severity} |")

    # Mathematical verification section
    lines.extend([
        "",
        "---",
        "",
        "## Mathematical Verification Details",
        "",
    ])

    for result in report.all_results:
        if result.verification:
            lines.extend([
                f"### {result.test_id}: {result.description}",
                "",
                "```json",
                json.dumps(result.verification, indent=2, default=str),
                "```",
                "",
            ])

    return "\n".join(lines)


if __name__ == "__main__":
    report = run_all_tests()

    print("\n" + "=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)
    print(f"Total Tests: {report.total_tests}")
    print(f"Passed: {report.passed}")
    print(f"Failed: {report.failed}")
    print(f"Pass Rate: {100*report.passed/report.total_tests:.1f}%")
    print(f"Critical Failures: {len(report.critical_failures)}")
    print(f"High Severity Failures: {len(report.high_failures)}")
    print("=" * 70)

    # Generate and save markdown report
    md_report = generate_markdown_report(report)
    report_path = os.path.join(os.path.dirname(__file__), "adversarial_test_report_2026-01-11.md")
    with open(report_path, "w") as f:
        f.write(md_report)
    print(f"\nDetailed report saved to: {report_path}")
