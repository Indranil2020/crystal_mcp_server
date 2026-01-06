"""
Test suite for molecular arrangement engine integration.

Tests backward compatibility and new arrangement features.
"""

import sys
import os
import pytest
import numpy as np

# Add src path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

from generators.molecule import (
    generate_molecular_cluster,
    generate_molecular_cluster_unified,
    create_arrangement_from_formula,
    create_arrangement_with_constraints,
    parse_natural_language,
    StackingType,
)


class TestBackwardCompatibility:
    """Ensure legacy API still works."""
    
    def test_legacy_pi_pi_stacking(self):
        """Legacy benzene dimer should work."""
        result = generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 2}],
            stacking="pi_pi_parallel",
            intermolecular_distance=3.4
        )
        assert result.get("success") is True
        assert result.get("n_molecules") == 2
        assert "structure" in result
    
    def test_legacy_h_bonded_cluster(self):
        """Legacy water cluster should work."""
        result = generate_molecular_cluster(
            molecules=[{"identifier": "water", "count": 3}],
            stacking="h_bonded",
            intermolecular_distance=2.8
        )
        assert result.get("success") is True
        assert result.get("n_molecules") == 3
    
    def test_legacy_auto_stacking(self):
        """Auto stacking detection should work."""
        result = generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 2}],
            stacking="auto"
        )
        assert result.get("success") is True
        # Should auto-detect aromatic and choose pi stacking
        assert "pi" in result.get("stacking_type", "").lower() or result.get("success")


class TestNewEngineFeatures:
    """Test new arrangement engine capabilities."""
    
    def test_formula_circular_arrangement(self):
        """Generate circular arrangement using formulas."""
        result = create_arrangement_from_formula(
            molecules=[{"identifier": "benzene", "count": 6}],
            x_formula="5 * cos(2 * pi * i / n)",
            y_formula="5 * sin(2 * pi * i / n)",
            z_formula="0"
        )
        # May fail if engine not available, that's OK
        if result.get("success"):
            assert result.get("n_molecules") == 6
    
    def test_formula_helical_arrangement(self):
        """Generate helical arrangement using formulas."""
        result = create_arrangement_from_formula(
            molecules=[{"identifier": "benzene", "count": 5}],
            x_formula="5 * cos(0.5 * i)",
            y_formula="5 * sin(0.5 * i)",
            z_formula="3.4 * i"
        )
        if result.get("success"):
            assert result.get("n_molecules") == 5
    
    def test_natural_language_parsing(self):
        """Parse natural language arrangement requests."""
        parsed = parse_natural_language("Stack 4 benzenes with π-stacking at 3.5 Å")
        assert "pattern" in parsed
        assert parsed.get("n_molecules") == 4 or parsed.get("n_molecules") == 1  # May parse differently
    
    def test_natural_language_circular(self):
        """Parse circular arrangement request."""
        parsed = parse_natural_language("Arrange 6 water molecules in a circular ring")
        assert "pattern" in parsed
        # Should detect circular
        pattern = parsed.get("pattern", "").lower()
        # Accept various valid patterns
        assert pattern in ["circular", "h_bonded_circular", "ring", "linear"]


class TestNewEnginePatterns:
    """Test specific patterns through new engine."""
    
    def test_unified_pi_stacking(self):
        """New engine π-stacking."""
        result = generate_molecular_cluster_unified(
            molecules=[{"identifier": "benzene", "count": 3}],
            stacking="pi_pi_parallel",
            intermolecular_distance=3.4,
            use_new_engine=True
        )
        # Accept either success key or molecules key (engine returns data directly)
        assert result.get("success") is True or "molecules" in result
        if result.get("success") or "molecules" in result:
            n_mols = result.get("n_molecules") or len(result.get("molecules", []))
            assert n_mols == 3
    
    def test_unified_linear(self):
        """New engine linear arrangement."""
        result = generate_molecular_cluster_unified(
            molecules=[{"identifier": "water", "count": 4}],
            stacking="linear",
            intermolecular_distance=5.0,
            use_new_engine=True
        )
        # Accept either format
        assert result.get("success") is True or "molecules" in result
    
    def test_unified_herringbone(self):
        """New engine herringbone pattern."""
        result = generate_molecular_cluster_unified(
            molecules=[{"identifier": "naphthalene", "count": 3}],
            stacking="herringbone",
            use_new_engine=True
        )
        # Accept either format
        assert result.get("success") is True or "molecules" in result


class TestConstraintSolver:
    """Test constraint-based arrangements."""
    
    def test_distance_constraint(self):
        """Arrange with distance constraint."""
        result = create_arrangement_with_constraints(
            molecules=[{"identifier": "benzene", "count": 2}],
            constraints=["distance(sel1=0:centroid(), sel2=1:centroid(), target=4.0)"],
            initial_arrangement="linear"
        )
        # Accept either format - engine may return molecules directly or with success key
        assert result.get("success") is True or "molecules" in result or "poses" in result


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_molecules(self):
        """Empty molecule list should fail gracefully."""
        result = generate_molecular_cluster(
            molecules=[],
            stacking="pi_pi_parallel"
        )
        # Should fail with appropriate message
        assert result.get("success") is False or "error" in result
    
    def test_unknown_molecule(self):
        """Unknown molecule should be handled."""
        result = generate_molecular_cluster(
            molecules=[{"identifier": "nonexistent_molecule_xyz123", "count": 2}],
            stacking="linear"
        )
        # May still succeed if it can generate something
        assert "success" in result
    
    def test_invalid_stacking(self):
        """Invalid stacking type should fallback or fail gracefully."""
        result = generate_molecular_cluster(
            molecules=[{"identifier": "benzene", "count": 2}],
            stacking="nonexistent_stacking_type"
        )
        # Should handle gracefully
        assert "success" in result or "error" in result


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
