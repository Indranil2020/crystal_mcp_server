/**
 * Structure Format Utilities
 *
 * Shared utilities for generating crystal structure file formats (CIF, POSCAR, XYZ, JSON).
 * Consolidated from duplicate implementations to ensure consistency.
 */

import { CrystalStructure } from "../types/crystal.js";

/**
 * Generate CIF file content from structure.
 *
 * By default, uses P1 space group to prevent CIF readers from incorrectly
 * expanding symmetry (since we output all atoms, not asymmetric unit).
 * Set useOriginalSpaceGroup=true only if outputting asymmetric unit positions.
 *
 * @param structure - The crystal structure
 * @param useOriginalSpaceGroup - If true, use the structure's space group instead of P1
 */
export function generateCIF(structure: CrystalStructure | any, useOriginalSpaceGroup: boolean = false): string {
  const lattice = structure?.lattice || {};
  // Default to P1 to prevent symmetry expansion issues when outputting all atoms
  const original_space_group = structure?.space_group || { symbol: 'P1', number: 1 };
  const space_group = useOriginalSpaceGroup ? original_space_group : { symbol: 'P1', number: 1 };
  const atoms = structure?.atoms || [];
  const matrix = lattice.matrix;

  let aVal = lattice.a;
  let bVal = lattice.b;
  let cVal = lattice.c;
  let alphaVal = lattice.alpha;
  let betaVal = lattice.beta;
  let gammaVal = lattice.gamma;
  let volumeVal = lattice.volume;

  const hasMatrix =
    Array.isArray(matrix) &&
    matrix.length === 3 &&
    matrix.every((row: number[]) => Array.isArray(row) && row.length >= 3);

  if (hasMatrix) {
    const v1: [number, number, number] = [matrix[0][0], matrix[0][1], matrix[0][2]];
    const v2: [number, number, number] = [matrix[1][0], matrix[1][1], matrix[1][2]];
    const v3: [number, number, number] = [matrix[2][0], matrix[2][1], matrix[2][2]];
    const norm = (v: [number, number, number]) => Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2);
    const dot = (vA: [number, number, number], vB: [number, number, number]) => vA[0] * vB[0] + vA[1] * vB[1] + vA[2] * vB[2];

    const aCalc = norm(v1);
    const bCalc = norm(v2);
    const cCalc = norm(v3);

    const angle = (vA: [number, number, number], vB: [number, number, number]) => {
      const denom = norm(vA) * norm(vB);
      if (denom === 0) {
        return undefined;
      }
      const cosang = Math.min(1, Math.max(-1, dot(vA, vB) / denom));
      return (Math.acos(cosang) * 180) / Math.PI;
    };

    aVal = aVal ?? aCalc;
    bVal = bVal ?? bCalc;
    cVal = cVal ?? cCalc;
    alphaVal = alphaVal ?? angle(v2, v3);
    betaVal = betaVal ?? angle(v1, v3);
    gammaVal = gammaVal ?? angle(v1, v2);

    if (volumeVal === undefined) {
      const det =
        v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) -
        v1[1] * (v2[0] * v3[2] - v2[2] * v3[0]) +
        v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
      volumeVal = Math.abs(det);
    }
  }

  let cif = `data_crystal\n`;
  cif += `_symmetry_space_group_name_H-M '${space_group.symbol || 'P1'}'\n`;
  cif += `_space_group_IT_number ${space_group.number || 1}\n`;

  if (aVal !== undefined) {
    cif += `_cell_length_a ${aVal.toFixed(6)}\n`;
    cif += `_cell_length_b ${(bVal ?? aVal).toFixed(6)}\n`;
    cif += `_cell_length_c ${(cVal ?? aVal).toFixed(6)}\n`;
    cif += `_cell_angle_alpha ${((alphaVal ?? 90)).toFixed(4)}\n`;
    cif += `_cell_angle_beta ${((betaVal ?? 90)).toFixed(4)}\n`;
    cif += `_cell_angle_gamma ${((gammaVal ?? 90)).toFixed(4)}\n`;
  }

  if (volumeVal !== undefined) {
    cif += `_cell_volume ${volumeVal.toFixed(4)}\n`;
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
