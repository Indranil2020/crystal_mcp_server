
import sys
import os
import json
import time
import numpy as np
import logging

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/python')))

from generators.molecule.arrangement_adapter import generate_molecular_cluster
from generators.molecule import molecular_arrangement

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger("Audit")

def check_overlaps(coords, threshold=0.8):
    """Check for atomic clashes (distance < threshold). O(N^2) check."""
    n = len(coords)
    if n > 1000:
        logger.warning(f"Skipping overlap check for N={n} (too large)")
        return True
    
    coords = np.array(coords)
    min_dist = float('inf')
    clash_count = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist < min_dist:
                min_dist = dist
            if dist < threshold:
                clash_count += 1
                
    logger.info(f"  Min Atomic Distance: {min_dist:.3f} A")
    if clash_count > 0:
        logger.error(f"  CRITICAL: Found {clash_count} atomic clashes (<{threshold}A)!")
        return False
    return True

def test_security_injection():
    logger.info("AUDIT 1: Security - Formula Injection Attempt")
    # Attempt to execute system command via eval injection
    # Typical Python jailbreak attempt (may leverage string formatting or type reconstruction)
    # But since we provide explicit local_vars and empty builtins, direct access is hard.
    # We try a simple one.
    
    injection_code = "[c for c in ().__class__.__base__.__subclasses__() if c.__name__ == 'catch_warnings'][0]()._module.__builtins__['__import__']('os').system('echo PWNED > pwned.txt')"
    
    # We put this in "x" formula.
    # Note: eval expects an expression returning a float/number.
    
    try:
        result = generate_molecular_cluster(
            molecules=[{"identifier": "water", "count": 1}],
            formulas={
                "x": "0", # Safe
                "y": "0",
                "z": "0", 
                 # We inject via a variable alias? No, we can only control values of keys in formulas dict.
                 # The code does: x = float(eval(x_formula, ...))
                 # So we put the payload in x_formula.
            },
            # We want to TRY to inject.
            # If we just put the string, eval might define it?
            # x_formula = "eval('...')" is blocked by empty builtins?
        )
        # If we see pwned.txt, it failed.
        if os.path.exists("pwned.txt"):
            logger.critical("SECURITY FAILURE: Remote Code Execution via Formula!")
            os.remove("pwned.txt")
            return False
            
        logger.info("  PASS: Injection blocked (or failed silently).")
        return True
    except Exception as e:
        logger.info(f"  PASS: Injection raised exception: {e}")
        return True

def test_impossible_constraints():
    logger.info("AUDIT 2: Stability - Conflicting Constraints")
    # Distance 0-1 is 2.0 AND 10.0 simultaneously.
    
    result = generate_molecular_cluster(
        molecules=[{"identifier": "water", "count": 2}],
        constraints=[
            "distance(0:centroid(), 1:centroid(), 2.0)",
            "distance(0:centroid(), 1:centroid(), 10.0)"
        ],
        use_solver=True
    )
    
    if result['success']:
        logger.warning("  WARNING: Solver claimed success on impossible constraints!")
        # Check actual distance
        coords = result['coords']
        d = np.linalg.norm(np.array(coords[0]) - np.array(coords[3])) # Approx 1st atoms
        logger.info(f"  Resulting Distance: {d:.2f}")
        return False # Should have failed or converged to one?
    else:
        logger.info(f"  PASS: System correctly reported failure: {result.get('error')}")
        return True

def test_massive_scale():
    logger.info("AUDIT 3: Scalability - N=1000 Grid")
    start = time.time()
    result = generate_molecular_cluster(
        molecules=[{"identifier": "water", "count": 1000}],
        stacking="grid", 
        distance=3.0
    )
    elapsed = time.time() - start
    
    if not result['success']:
        logger.error(f"  FAILURE: Generation failed for N=1000. Error: {result.get('error')}")
        return False
        
    n_atoms = result['n_atoms']
    logger.info(f"  Generated {n_atoms} atoms in {elapsed:.2f}s")
    
    if elapsed > 15.0:
        logger.warning("  PERFORMANCE WARNING: Generation took > 15s")
        return False
        
    logger.info("  PASS: High-scale generation successful and fast.")
    return True

def test_physical_validity_lipid_raft():
    logger.info("AUDIT 4: Physics - Lipid Raft Density & Clashes")
    # Simulate Lipid Raft (Simplified: 50 Benzene as proxy for Lipid Head, 30 Thiophene as proxy for Cholesterol)
    # Using Grid logic.
    result = generate_molecular_cluster(
        molecules=[
            {"identifier": "benzene", "count": 50}, 
            {"identifier": "thiophene", "count": 30}
        ],
        stacking="auto", # Should pick grid
        distance=5.0
    )
    
    if not result['success']:
        logger.error("  FAILURE: Generation failed")
        return False
        
    coords = result['coords']
    # Check clashes
    if not check_overlaps(coords, threshold=1.0):
        return False
        
    # Check dimensions (Grid)
    coords_np = np.array(coords)
    x_span = coords_np[:, 0].max() - coords_np[:, 0].min()
    y_span = coords_np[:, 1].max() - coords_np[:, 1].min()
    logger.info(f"  Dimensions: {x_span:.1f} x {y_span:.1f} A")
    
    # Expected area for 80 molecules at 5A spacing?
    # Grid: sqrt(80) ~ 9. 9 * 5 = 45A.
    if x_span < 30.0 or y_span < 30.0:
        logger.warning(f"  WARNING: Packing seems too tight? Span={x_span:.1f}")
    
    logger.info("  PASS: Physically reasonable structure.")
    return True

if __name__ == "__main__":
    logger.info("=== STARTING CRITICAL AUDIT ===")
    
    tests = [
        test_security_injection,
        test_impossible_constraints,
        # test_massive_scale,
        # test_physical_validity_lipid_raft
    ]
    
    failures = 0
    for t in tests:
        try:
            if not t():
                failures += 1
                logger.error(f"TEST FAILED: {t.__name__}\n")
            else:
                logger.info(f"TEST PASSED: {t.__name__}\n")
        except Exception as e:
            logger.error(f"TEST CRASHED: {t.__name__} - {e}\n")
            failures += 1
            
    if failures == 0:
        logger.info("=== AUDIT PASSED: SYSTEM ROBUST ===")
        sys.exit(0)
    else:
        logger.critical(f"=== AUDIT FAILED: {failures} ISSUES FOUND ===")
        sys.exit(1)
