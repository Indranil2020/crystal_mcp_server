"""
Test suite for universal molecule generation.
"""

import sys
import os
import pytest

# Add src path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

from generators.molecule import (
    generate_molecule,
    generate_molecule_universal,
    smiles_to_3d_structure,
    get_supported_molecules,
    MOLECULE_DATABASE,
    MOLECULE_ALIASES,
)


class TestLocalDatabase:
    """Tests for molecules in the local database."""
    
    def test_co2_from_local_db(self):
        """CO2 should be retrieved from local database."""
        result = generate_molecule('CO2')
        assert result['success'] is True
        assert result['source'] == 'local_database'
        assert result['n_atoms'] == 3
        assert 'C' in result['atoms']
        assert 'O' in result['atoms']
    
    def test_h2o_from_local_db(self):
        """H2O should be retrieved from local database."""
        result = generate_molecule('H2O')
        assert result['success'] is True
        assert result['source'] == 'local_database'
        assert result['n_atoms'] == 3
    
    def test_all_local_molecules(self):
        """All molecules in local DB should be retrievable."""
        for formula in list(MOLECULE_DATABASE.keys())[:5]:  # Test first 5
            result = generate_molecule(formula, allow_external=False)
            assert result['success'] is True
            assert 'atoms' in result
            assert 'coords' in result


class TestUniversalGeneration:
    """Tests for universal molecule generation."""
    
    def test_smiles_benzene(self):
        """Benzene from SMILES should generate 12 atoms (6C + 6H)."""
        result = generate_molecule_universal('c1ccccc1', input_type='smiles')
        assert result['success'] is True
        assert result['source'] == 'rdkit'
        assert result['n_atoms'] == 12  # 6 carbons + 6 hydrogens
        assert result['atoms'].count('C') == 6
        assert result['atoms'].count('H') == 6
    
    def test_smiles_ethanol(self):
        """Ethanol from SMILES should generate correctly."""
        result = generate_molecule_universal('CCO', input_type='smiles')
        assert result['success'] is True
        assert result['n_atoms'] == 9  # C2H6O = 2C + 6H + 1O
    
    def test_aspirin_from_alias(self):
        """Aspirin should be resolved from local alias."""
        result = generate_molecule_universal('aspirin')
        assert result['success'] is True
        assert result['n_atoms'] == 21  # C9H8O4
    
    def test_caffeine_from_alias(self):
        """Caffeine should be resolved from local alias."""
        result = generate_molecule_universal('caffeine')
        assert result['success'] is True
        assert result['n_atoms'] == 24  # C8H10N4O2


class TestPTCDAGeneration:
    """Tests for PTCDA - the user's specific request."""
    
    def test_ptcda_generation(self):
        """PTCDA should be generated via PubChem or alias."""
        result = generate_molecule_universal('PTCDA')
        assert result['success'] is True
        # PTCDA: C24H8O6 = 38 atoms
        assert result['n_atoms'] == 38
        assert result['formula'] == 'C24H8O6'


class TestFallbackBehavior:
    """Tests for fallback from generate_molecule to universal."""
    
    def test_unknown_molecule_fallback(self):
        """Unknown molecule should fall back to universal generation."""
        result = generate_molecule('caffeine')  # Not in local DB
        assert result['success'] is True
        assert 'source' in result
    
    def test_offline_mode(self):
        """Offline mode should only use local DB."""
        result = generate_molecule('caffeine', allow_external=False)
        assert result['success'] is False
        assert 'error' in result


class TestSupportedMolecules:
    """Tests for get_supported_molecules function."""
    
    def test_returns_capabilities(self):
        """Should return capability information."""
        info = get_supported_molecules()
        assert 'local_database' in info
        assert 'local_aliases' in info
        assert 'rdkit_available' in info
        assert 'capabilities' in info


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
