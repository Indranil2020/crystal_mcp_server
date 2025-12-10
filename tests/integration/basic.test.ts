/**
 * Integration Tests - Basic Functionality
 * 
 * Tests for core crystal generation and analysis tools.
 */

import { describe, test, expect, beforeAll } from '@jest/globals';
import { generateCrystal } from '../../src/tools/generation/generate-crystal.js';
import { analyzeSymmetry } from '../../src/tools/analysis/symmetry.js';
import { validateStructure } from '../../src/tools/analysis/validation.js';
import { makeSupercell } from '../../src/tools/transformation/supercell.js';

describe('Crystal Generation', () => {
  test('generates diamond Si structure', async () => {
    const result = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(result.success).toBe(true);
    if (result.success) {
      expect(result.data.structure.space_group.number).toBe(227);
      expect(result.data.structure.metadata.formula).toContain('Si');
      expect(result.data.validation.valid).toBe(true);
    }
  }, 30000);

  test('generates rock salt NaCl structure', async () => {
    const result = await generateCrystal({
      composition: ["Na", "Cl"],
      space_group: 225,
      seed: 42
    });

    expect(result.success).toBe(true);
    if (result.success) {
      expect(result.data.structure.space_group.number).toBe(225);
      expect(result.data.structure.metadata.natoms).toBeGreaterThan(0);
    }
  }, 30000);

  test('handles invalid space group', async () => {
    const result = await generateCrystal({
      composition: ["Si"],
      space_group: 999
    });

    expect(result.success).toBe(false);
    if (!result.success) {
      expect(result.error.code).toBeDefined();
    }
  });

  test('respects minimum distance constraints', async () => {
    const result = await generateCrystal({
      composition: ["Na", "Cl"],
      space_group: 225,
      min_distance: { "Na-Cl": 2.5 },
      seed: 42
    });

    expect(result.success).toBe(true);
    if (result.success) {
      expect(result.data.validation.valid).toBe(true);
    }
  }, 30000);
});

describe('Symmetry Analysis', () => {
  test('detects space group correctly', async () => {
    // First generate a structure
    const genResult = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    // Analyze symmetry
    const symResult = await analyzeSymmetry({
      structure: genResult.data.structure,
      symprec: 0.001
    });

    expect(symResult.success).toBe(true);
    if (symResult.success) {
      expect(symResult.data.space_group.number).toBe(227);
      expect(symResult.data.symmetry_operations.n_operations).toBeGreaterThan(0);
    }
  }, 60000);

  test('finds primitive cell', async () => {
    const genResult = await generateCrystal({
      composition: ["Na", "Cl"],
      space_group: 225,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    const symResult = await analyzeSymmetry({
      structure: genResult.data.structure,
      detect_primitive: true
    });

    expect(symResult.success).toBe(true);
    if (symResult.success && symResult.data.primitive_cell) {
      expect(symResult.data.primitive_cell.metadata.natoms).toBeLessThanOrEqual(
        genResult.data.structure.metadata.natoms
      );
    }
  }, 60000);
});

describe('Structure Validation', () => {
  test('validates correct structure', async () => {
    const genResult = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    const valResult = await validateStructure({
      structure: genResult.data.structure,
      checks: ["distances", "lattice", "stoichiometry"]
    });

    expect(valResult.success).toBe(true);
    if (valResult.success) {
      expect(valResult.data.valid).toBe(true);
      expect(valResult.data.metrics).toBeDefined();
      expect(valResult.data.metrics.density).toBeGreaterThan(0);
    }
  }, 60000);

  test('detects distance violations', async () => {
    const genResult = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    const valResult = await validateStructure({
      structure: genResult.data.structure,
      checks: ["distances"],
      min_distance: 10.0  // Unrealistic constraint
    });

    expect(valResult.success).toBe(true);
    if (valResult.success) {
      expect(valResult.data.valid).toBe(false);
      expect(valResult.data.errors.length).toBeGreaterThan(0);
    }
  }, 60000);
});

describe('Structure Transformations', () => {
  test('creates supercell', async () => {
    const genResult = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    const superResult = await makeSupercell({
      structure: genResult.data.structure,
      scaling_matrix: [2, 2, 2]
    });

    expect(superResult.success).toBe(true);
    if (superResult.success) {
      expect(superResult.data.volume_multiplier).toBe(8);
      expect(superResult.data.supercell.metadata.natoms).toBe(
        genResult.data.structure.metadata.natoms * 8
      );
    }
  }, 60000);

  test('creates non-cubic supercell', async () => {
    const genResult = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(genResult.success).toBe(true);
    if (!genResult.success) return;

    const superResult = await makeSupercell({
      structure: genResult.data.structure,
      scaling_matrix: [2, 3, 4]
    });

    expect(superResult.success).toBe(true);
    if (superResult.success) {
      expect(superResult.data.volume_multiplier).toBe(24);
    }
  }, 60000);
});

describe('Error Handling', () => {
  test('handles missing required parameters', async () => {
    const result = await generateCrystal({
      composition: ["Si"]
      // Missing space_group
    });

    expect(result.success).toBe(false);
    if (!result.success) {
      expect(result.error.suggestions.length).toBeGreaterThan(0);
    }
  });

  test('provides helpful error messages', async () => {
    const result = await generateCrystal({
      composition: ["InvalidElement"],
      space_group: 1
    });

    expect(result.success).toBe(false);
    if (!result.success) {
      expect(result.error.message).toBeDefined();
      expect(typeof result.error.message).toBe('string');
    }
  });
});

describe('File Format Generation', () => {
  test('generates all file formats', async () => {
    const result = await generateCrystal({
      composition: ["Si", "Si"],
      space_group: 227,
      seed: 42
    });

    expect(result.success).toBe(true);
    if (result.success) {
      expect(result.data.files.cif).toBeDefined();
      expect(result.data.files.poscar).toBeDefined();
      expect(result.data.files.xyz).toBeDefined();
      expect(result.data.files.json).toBeDefined();

      // Check CIF format
      expect(result.data.files.cif).toContain('data_crystal');
      expect(result.data.files.cif).toContain('_symmetry_space_group');

      // Check POSCAR format
      expect(result.data.files.poscar).toContain('Direct');
      expect(result.data.files.poscar).toContain('Si');

      // Check XYZ format
      const xyzLines = result.data.files.xyz.split('\n');
      expect(xyzLines.length).toBeGreaterThan(2);
    }
  }, 30000);
});
