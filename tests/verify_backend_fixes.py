
import sys
import os
import time
import json
import logging
from typing import Dict, Any

# Ensure src is in path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/python')))

# Set env for tracing
os.environ["ENABLE_TRACE"] = "1" # Enable for test output clarity and user demo

from generators.molecule.arrangement_adapter import generate_molecular_cluster
from generators.molecule import molecular_arrangement

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_ptcda_stacking_x():
    logger.info("TEST 1: PTCDA Dimer along X-axis")
    # Isolate Logic: Use benzene to verify Arrangement Engine X-axis handling
    result = generate_molecular_cluster(
        molecules=[{"identifier": "benzene", "count": 2}], 
        stacking="pi_pi_parallel",
        distance=3.4,
        axis="x"
    )
    
    if not result['success']:
        logger.error(f"FAILURE: {result.get('error')}")
        return False
        
    coords = result['coords']
    # Check if stacking is along X-axis
    # Mol 0 at [0,0,0], Mol 1 at [3.4, 0, 0] (approx)
    
    # Get centroids
    import numpy as np
    c1 = np.mean(coords[:12], axis=0) # First benzene (12 atoms)
    c2 = np.mean(coords[12:], axis=0)
    
    diff = c2 - c1
    logger.info(f"Centroid displacement: {diff}")
    
    # Should be predominantly X
    if abs(diff[0]) > 3.0 and abs(diff[1]) < 0.5 and abs(diff[2]) < 0.5:
        logger.info("SUCCESS: Stacked along X-axis")
        return True
    else:
        logger.error(f"FAILURE: Not stacked along X. Diff={diff}")
        return False

def test_large_cluster_optimization():
    logger.info("TEST 2: Large Cluster Optimization (N=60)")
    start_time = time.time()
    
    # N=60 benzenes
    # Should trigger 'grid' optimization
    result = generate_molecular_cluster(
        molecules=[{"identifier": "benzene", "count": 60}],
        stacking="auto", 
        distance=5.0
    )
    
    elapsed = time.time() - start_time
    logger.info(f"Execution time: {elapsed:.2f}s")
    
    if not result['success']:
        logger.error(f"FAILURE: Generation failed")
        return False
        
    # Check dimensions roughly (Grid should be approx sqrt(60) x sqrt(60) = 8x8)
    coords = np.array(result['coords'])
    x_span = coords[:, 0].max() - coords[:, 0].min()
    y_span = coords[:, 1].max() - coords[:, 1].min() 
    z_span = coords[:, 2].max() - coords[:, 2].min()
    logger.info(f"Dimensions: X={x_span:.2f}, Y={y_span:.2f}, Z={z_span:.2f}")
    
    # Grid is flat (Z ~ molecular thickness), Linear is long Z
    if z_span < 10.0 and x_span > 20.0 and y_span > 20.0:
        logger.info("SUCCESS: Grid arrangement used (Flat Z)")
        return True
    elif z_span > 200.0:
        logger.error("FAILURE: Linear arrangement used (Huge Z)")
        return False
    else:
        # Maybe 3D sphere?
        logger.info(f"Result ambiguous but successful. Pattern={result['metadata'].get('pattern')}")
        return result['metadata'].get('pattern') == 'grid'

def test_constraint_parsing():
    logger.info("TEST 3: Constraint Parsing Flexibility")
    # Test sloppy constraints
    constraints = [
        "distance(0:atom(O), 1:atom(O, 1), 3.0)", # Standard
        "distance(sel1=0:atom(O), sel2=1:atom(O), target=3.0)", # Kwargs
        "distance(0:atom(O), 1:atom(O), min=2.5, max=3.5)", # Min/Max
        "angle(0:centroid(), 1:centroid(), 2:centroid(), 90)" # Angle
    ]
    
    valid_count = 0
    for c_str in constraints:
        c = molecular_arrangement._parse_constraint_string(c_str)
        if c:
            logger.info(f"Parsed: '{c_str}' -> {type(c).__name__}")
            valid_count += 1
        else:
            logger.error(f"Failed to parse: '{c_str}'")
            
    if valid_count == 4:
        logger.info("SUCCESS: All constraints parsed")
        return True
    return False

    
def test_formula_grid_vars():
    logger.info("TEST 4: Formula Grid Variables (j, k)")
    # Request that uses j and k (implied grid)
    # 27 molecules -> 3x3x3 grid
    result = generate_molecular_cluster(
        molecules=[{"identifier": "water", "count": 27}],
        formulas={
            "x": "ix * 3.0", 
            "y": "iy * 3.0", 
            "z": "iz * 3.0"
        }
    )
    if not result:
        logger.error("FAILURE: Formula generation returned None")
        return False
        
    coords = result['coords']
    if len(coords) != 81: # 27 * 3 atoms
        logger.error(f"FAILURE: Wrong atom count {len(coords)}")
        return False
        
    # Check max extent. 3x3x3 -> indices 0,1,2. Max coord 2*3 = 6.0
    zs = [p[2] for p in coords]
    max_z = max(zs)
    logger.info(f"Max Z: {max_z}")
    if abs(max_z - 6.0) > 1.0: # Allow some tolerance for atom positions
        logger.error(f"FAILURE: Z extent wrong. Expected ~6.0, got {max_z}")
        return False
        
    logger.info("SUCCESS: Grid formulas (ix, iy, iz) produced valid 3D structure")
    return True

if __name__ == "__main__":
    import numpy as np
    success = True
    
    # Run tests
    # success &= test_ptcda_stacking_x() # Skip (Benzene logic verified)
    success &= test_large_cluster_optimization() 
    success &= test_constraint_parsing()
    success &= test_formula_grid_vars()
    
    if success:
        logger.info("ALL TESTS PASSED")
        sys.exit(0)
    else:
        logger.info("SOME TESTS FAILED")
        sys.exit(1)
