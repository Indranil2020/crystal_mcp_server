/**
 * Structure Format Utilities
 *
 * Shared utilities for generating crystal structure file formats (CIF, POSCAR, XYZ, JSON).
 * Consolidated from duplicate implementations to ensure consistency.
 */

import { CrystalStructure } from "../types/crystal.js";

/**
 * Generate CIF file content from structure.
 */
export function generateCIF(structure: CrystalStructure | any): string {
  const lattice = structure?.lattice || {};
  const space_group = structure?.space_group || { symbol: 'P1', number: 1 };
  const atoms = structure?.atoms || [];

  let cif = `data_crystal\n`;
  cif += `_symmetry_space_group_name_H-M '${space_group.symbol || 'P1'}'\n`;
  cif += `_space_group_IT_number ${space_group.number || 1}\n`;

  if (lattice.a !== undefined) {
    cif += `_cell_length_a ${lattice.a.toFixed(6)}\n`;
    cif += `_cell_length_b ${(lattice.b || lattice.a).toFixed(6)}\n`;
    cif += `_cell_length_c ${(lattice.c || lattice.a).toFixed(6)}\n`;
    cif += `_cell_angle_alpha ${(lattice.alpha || 90).toFixed(4)}\n`;
    cif += `_cell_angle_beta ${(lattice.beta || 90).toFixed(4)}\n`;
    cif += `_cell_angle_gamma ${(lattice.gamma || 90).toFixed(4)}\n`;
  }

  if (lattice.volume !== undefined) {
    cif += `_cell_volume ${lattice.volume.toFixed(4)}\n`;
  }

  cif += `\nloop_\n`;
  cif += `_atom_site_label\n`;
  cif += `_atom_site_type_symbol\n`;
  cif += `_atom_site_fract_x\n`;
  cif += `_atom_site_fract_y\n`;
  cif += `_atom_site_fract_z\n`;

  atoms.forEach((atom: any, i: number) => {
    const element = atom.element || 'X';
    const coords = atom.coords || [0, 0, 0];
    cif += `${element}${i+1} ${element} ${coords[0].toFixed(6)} ${coords[1].toFixed(6)} ${coords[2].toFixed(6)}\n`;
  });

  return cif;
}

/**
 * Generate POSCAR file content from structure.
 */
export function generatePOSCAR(structure: CrystalStructure | any): string {
  const lattice = structure?.lattice || {};
  const atoms = structure?.atoms || [];
  const metadata = structure?.metadata || {};

  let poscar = `${metadata.formula || 'Structure'}\n`;
  poscar += `1.0\n`;

  // Lattice vectors
  if (lattice.matrix && Array.isArray(lattice.matrix)) {
    lattice.matrix.forEach((vec: number[]) => {
      if (vec && vec.length >= 3) {
        poscar += `  ${(vec[0] || 0).toFixed(12)}  ${(vec[1] || 0).toFixed(12)}  ${(vec[2] || 0).toFixed(12)}\n`;
      }
    });
  }

  // Count atoms by element
  const elementCounts: Record<string, number> = {};
  const elementOrder: string[] = [];

  atoms.forEach((atom: any) => {
    const element = atom.element || 'X';
    if (!elementCounts[element]) {
      elementCounts[element] = 0;
      elementOrder.push(element);
    }
    elementCounts[element]!++;
  });

  // Element names
  poscar += elementOrder.join(' ') + '\n';

  // Element counts
  poscar += elementOrder.map(el => elementCounts[el]).join(' ') + '\n';

  poscar += `Direct\n`;

  // Coordinates grouped by element
  elementOrder.forEach(element => {
    atoms.filter((a: any) => a.element === element).forEach((atom: any) => {
      const coords = atom.coords || [0, 0, 0];
      poscar += `  ${coords[0].toFixed(12)}  ${coords[1].toFixed(12)}  ${coords[2].toFixed(12)}\n`;
    });
  });

  return poscar;
}

/**
 * Generate XYZ file content from structure.
 */
export function generateXYZ(structure: CrystalStructure | any): string {
  const atoms = structure?.atoms || [];
  const metadata = structure?.metadata || {};

  let xyz = `${atoms.length}\n`;
  xyz += `${metadata.formula || 'Structure'}\n`;

  atoms.forEach((atom: any) => {
    const element = atom.element || 'X';
    const cartesian = atom.cartesian || atom.coords || [0, 0, 0];
    xyz += `${element}  ${cartesian[0].toFixed(8)}  ${cartesian[1].toFixed(8)}  ${cartesian[2].toFixed(8)}\n`;
  });

  return xyz;
}

/**
 * Generate JSON file content from structure.
 */
export function generateJSON(structure: CrystalStructure | any): string {
  return JSON.stringify(structure, null, 2);
}
