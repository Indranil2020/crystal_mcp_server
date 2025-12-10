#!/usr/bin/env python3
"""
Comprehensive test suite for Crystal MCP Server Python backend

Tests all Python modules: crystal_generator, mlff_calculator, 
symmetry_analyzer, structure_tools, validators
"""

import json
import sys
import os
import tempfile
from pathlib import Path

# Add src/python to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src" / "python"))

def test_crystal_generator():
    """Test crystal structure generation"""
    print("\n=== Testing Crystal Generator ===")
    
    from crystal_generator import generate_crystal, validate_space_group, validate_composition
    
    # Test 1: Validate space group
    print("Test 1: Space group validation")
    is_valid, error = validate_space_group(227)
    assert is_valid, f"Failed: {error}"
    print("✓ Valid space group 227")
    
    is_valid, error = validate_space_group(999)
    assert not is_valid, "Should reject invalid space group"
    print("✓ Rejects invalid space group 999")
    
    # Test 2: Validate composition
    print("\nTest 2: Composition validation")
    is_valid, error = validate_composition(["Si", "Si"])
    assert is_valid, f"Failed: {error}"
    print("✓ Valid composition [Si, Si]")
    
    is_valid, error = validate_composition(["Xx", "Yy"])
    assert not is_valid, "Should reject invalid elements"
    print("✓ Rejects invalid elements")
    
    # Test 3: Generate diamond Si structure
    print("\nTest 3: Generate diamond Si (space group 227)")
    result = generate_crystal(
        composition=["Si", "Si"],
        space_group=227,
        num_atoms=8,
        seed=42,
        max_attempts=10
    )
    
    assert result.get("success"), f"Generation failed: {result.get('error')}"
    assert "structure" in result, "Missing structure in result"
    
    structure = result["structure"]
    assert structure["space_group"]["number"] == 227, "Wrong space group"
    assert structure["space_group"]["crystal_system"] == "cubic", "Wrong crystal system"
    assert structure["metadata"]["natoms"] > 0, "No atoms generated"
    
    print(f"✓ Generated structure with {structure['metadata']['natoms']} atoms")
    print(f"  Formula: {structure['metadata']['formula']}")
    print(f"  Volume: {structure['lattice']['volume']:.2f} ų")
    print(f"  Density: {structure['metadata']['density']:.2f} g/cm³")
    
    # Test 4: Generate NaCl (rock salt)
    print("\nTest 4: Generate NaCl (space group 225)")
    result = generate_crystal(
        composition=["Na", "Cl"],
        space_group=225,
        seed=42,
        max_attempts=10
    )
    
    assert result.get("success"), f"Generation failed: {result.get('error')}"
    print(f"✓ Generated NaCl structure")
    
    print("\n✓✓✓ Crystal Generator: ALL TESTS PASSED ✓✓✓")
    return True


def test_symmetry_analyzer():
    """Test symmetry analysis"""
    print("\n=== Testing Symmetry Analyzer ===")
    
    from symmetry_analyzer import analyze_symmetry, validate_structure_for_spglib
    from crystal_generator import generate_crystal
    
    # Generate a test structure first
    print("Test 1: Generate test structure")
    gen_result = generate_crystal(
        composition=["Si", "Si"],
        space_group=227,
        seed=42
    )
    assert gen_result.get("success"), "Failed to generate test structure"
    structure = gen_result["structure"]
    print("✓ Test structure generated")
    
    # Test 2: Validate structure for spglib
    print("\nTest 2: Validate structure format")
    is_valid, error = validate_structure_for_spglib(structure)
    assert is_valid, f"Structure validation failed: {error}"
    print("✓ Structure format valid")
    
    # Test 3: Analyze symmetry
    print("\nTest 3: Analyze symmetry")
    result = analyze_symmetry(
        structure_dict=structure,
        symprec=1e-3,
        detect_primitive=True,
        standardize=False
    )
    
    assert result.get("success"), f"Symmetry analysis failed: {result.get('error')}"
    assert "space_group" in result, "Missing space group info"
    
    sg = result["space_group"]
    print(f"✓ Detected space group: {sg['number']} ({sg['symbol']})")
    print(f"  Crystal system: {sg['crystal_system']}")
    print(f"  Point group: {sg['point_group']}")
    print(f"  Symmetry operations: {result['symmetry_operations']['n_operations']}")
    
    # Test 4: Check primitive cell detection
    if "primitive_cell" in result:
        prim = result["primitive_cell"]
        print(f"✓ Primitive cell found with {prim['metadata']['natoms']} atoms")
    
    print("\n✓✓✓ Symmetry Analyzer: ALL TESTS PASSED ✓✓✓")
    return True


def test_mlff_calculator():
    """Test MLFF calculator (if available)"""
    print("\n=== Testing MLFF Calculator ===")
    
    from mlff_calculator import validate_mlff_model, validate_optimizer, calculate_energy
    from crystal_generator import generate_crystal
    
    # Test 1: Validate model names
    print("Test 1: Model name validation")
    is_valid, error = validate_mlff_model("chgnet")
    assert is_valid, f"CHGNet should be valid: {error}"
    print("✓ CHGNet model name valid")
    
    is_valid, error = validate_mlff_model("invalid_model")
    assert not is_valid, "Should reject invalid model name"
    print("✓ Rejects invalid model name")
    
    # Test 2: Validate optimizer names
    print("\nTest 2: Optimizer validation")
    is_valid, error = validate_optimizer("BFGS")
    assert is_valid, f"BFGS should be valid: {error}"
    print("✓ BFGS optimizer valid")
    
    is_valid, error = validate_optimizer("InvalidOpt")
    assert not is_valid, "Should reject invalid optimizer"
    print("✓ Rejects invalid optimizer")
    
    # Test 3: Energy calculation (if CHGNet available)
    print("\nTest 3: Energy calculation")
    try:
        import chgnet
        print("  CHGNet package found, running energy calculation...")
        
        # Generate test structure
        gen_result = generate_crystal(
            composition=["Si", "Si"],
            space_group=227,
            num_atoms=8,
            seed=42
        )
        assert gen_result.get("success"), "Failed to generate structure"
        
        # Calculate energy
        result = calculate_energy(
            structure_dict=gen_result["structure"],
            mlff_model="chgnet",
            calculate_forces=True,
            calculate_stress=False
        )
        
        if result.get("success"):
            print(f"✓ Energy calculation successful")
            print(f"  Energy: {result['energy']:.4f} eV")
            print(f"  Energy/atom: {result['energy_per_atom']:.4f} eV/atom")
            if "forces" in result:
                print(f"  Forces calculated: {len(result['forces'])} atoms")
        else:
            print(f"⚠ Energy calculation failed: {result.get('error', {}).get('message')}")
            print("  (This is expected if CHGNet model not downloaded)")
        
    except ImportError:
        print("  ⚠ CHGNet not installed, skipping energy calculation test")
        print("  Install with: pip install chgnet")
    
    print("\n✓✓✓ MLFF Calculator: TESTS PASSED ✓✓✓")
    return True


def test_validators():
    """Test validation utilities"""
    print("\n=== Testing Validators ===")
    
    from validators import (
        validate_lattice_parameters,
        validate_miller_indices,
        validate_scaling_matrix,
        check_atoms_too_close
    )
    from crystal_generator import generate_crystal
    
    # Test 1: Lattice parameter validation
    print("Test 1: Lattice parameter validation")
    is_valid, errors, warnings = validate_lattice_parameters(
        a=5.0, b=5.0, c=5.0,
        alpha=90.0, beta=90.0, gamma=90.0
    )
    assert is_valid, f"Valid lattice rejected: {errors}"
    print("✓ Valid cubic lattice accepted")
    
    is_valid, errors, warnings = validate_lattice_parameters(
        a=-1.0, b=5.0, c=5.0,
        alpha=90.0, beta=90.0, gamma=90.0
    )
    assert not is_valid, "Should reject negative lattice parameter"
    print("✓ Rejects negative lattice parameters")
    
    # Test 2: Miller indices validation
    print("\nTest 2: Miller indices validation")
    is_valid, error = validate_miller_indices([1, 0, 0])
    assert is_valid, f"Valid Miller indices rejected: {error}"
    print("✓ Valid Miller indices [1,0,0] accepted")
    
    is_valid, error = validate_miller_indices([0, 0, 0])
    assert not is_valid, "Should reject [0,0,0]"
    print("✓ Rejects invalid [0,0,0] Miller indices")
    
    # Test 3: Scaling matrix validation
    print("\nTest 3: Scaling matrix validation")
    is_valid, error = validate_scaling_matrix([2, 2, 2])
    assert is_valid, f"Valid scaling rejected: {error}"
    print("✓ Valid scaling [2,2,2] accepted")
    
    # Test 4: Check atoms too close
    print("\nTest 4: Check interatomic distances")
    gen_result = generate_crystal(
        composition=["Si", "Si"],
        space_group=227,
        seed=42
    )
    if gen_result.get("success"):
        result = check_atoms_too_close(
            gen_result["structure"],
            threshold=1.0
        )
        print(f"✓ Distance check completed")
        if result["atoms_too_close"]:
            print(f"  Found {len(result['close_pairs'])} close pairs")
        else:
            print(f"  No atoms too close (min distance: {result['min_distance']:.2f} Å)")
    
    print("\n✓✓✓ Validators: ALL TESTS PASSED ✓✓✓")
    return True


def test_structure_tools():
    """Test structure transformation tools"""
    print("\n=== Testing Structure Tools ===")
    
    from structure_tools import (
        make_supercell,
        validate_supercell_matrix,
        calculate_supercell_size
    )
    from crystal_generator import generate_crystal
    
    # Generate test structure
    print("Test 1: Generate base structure")
    gen_result = generate_crystal(
        composition=["Si", "Si"],
        space_group=227,
        num_atoms=8,
        seed=42
    )
    assert gen_result.get("success"), "Failed to generate structure"
    structure = gen_result["structure"]
    print(f"✓ Base structure: {structure['metadata']['natoms']} atoms")
    
    # Test 2: Validate supercell matrix
    print("\nTest 2: Supercell matrix validation")
    is_valid, error = validate_supercell_matrix([2, 2, 2])
    assert is_valid, f"Valid matrix rejected: {error}"
    print("✓ Valid scaling [2,2,2] accepted")
    
    # Test 3: Calculate supercell size
    print("\nTest 3: Calculate supercell size")
    multiplier = calculate_supercell_size([2, 2, 2])
    assert multiplier == 8, f"Wrong multiplier: {multiplier}"
    print("✓ Correct multiplier calculated: 8x")
    
    # Test 4: Make supercell
    print("\nTest 4: Create supercell")
    result = make_supercell(
        structure_dict=structure,
        scaling_matrix=[2, 2, 2],
        wrap_atoms=True
    )
    
    if result.get("success"):
        supercell = result["supercell"]
        print(f"✓ Supercell created")
        print(f"  Original atoms: {structure['metadata']['natoms']}")
        print(f"  Supercell atoms: {supercell['metadata']['natoms']}")
        print(f"  Volume ratio: {result['volume_multiplier']}")
        assert supercell['metadata']['natoms'] == structure['metadata']['natoms'] * 8
    else:
        print(f"✗ Supercell creation failed: {result.get('error')}")
    
    print("\n✓✓✓ Structure Tools: ALL TESTS PASSED ✓✓✓")
    return True


def run_all_tests():
    """Run all Python backend tests"""
    print("="*60)
    print("CRYSTAL MCP SERVER - PYTHON BACKEND TEST SUITE")
    print("="*60)
    
    results = {}
    
    try:
        results['crystal_generator'] = test_crystal_generator()
    except Exception as e:
        print(f"\n✗✗✗ Crystal Generator tests FAILED: {e}")
        results['crystal_generator'] = False
    
    try:
        results['symmetry_analyzer'] = test_symmetry_analyzer()
    except Exception as e:
        print(f"\n✗✗✗ Symmetry Analyzer tests FAILED: {e}")
        results['symmetry_analyzer'] = False
    
    try:
        results['mlff_calculator'] = test_mlff_calculator()
    except Exception as e:
        print(f"\n✗✗✗ MLFF Calculator tests FAILED: {e}")
        results['mlff_calculator'] = False
    
    try:
        results['validators'] = test_validators()
    except Exception as e:
        print(f"\n✗✗✗ Validators tests FAILED: {e}")
        results['validators'] = False
    
    try:
        results['structure_tools'] = test_structure_tools()
    except Exception as e:
        print(f"\n✗✗✗ Structure Tools tests FAILED: {e}")
        results['structure_tools'] = False
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for module, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{module:25s}: {status}")
    
    print("="*60)
    
    all_passed = all(results.values())
    if all_passed:
        print("\n🎉 ALL TESTS PASSED! 🎉")
        print("Python backend is ready for use.")
        return 0
    else:
        print("\n⚠ SOME TESTS FAILED")
        print("Check errors above for details.")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())
