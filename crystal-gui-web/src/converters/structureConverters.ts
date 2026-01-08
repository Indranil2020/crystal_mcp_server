/**
 * Structure Converters - Convert MCP output to viewer formats
 * 
 * Converts structure data from Python backend JSON to formats 
 * suitable for MolStar (CIF, PDB) and Kekule.js (Mol2, SMILES)
 * 
 * DEBUG: Comprehensive validation and logging enabled
 */

import type { StructureData, AtomData, LatticeData } from '../types';
import { debug, debugError } from '../debug';

/**
 * Validate structure data before conversion
 */
function validateStructureData(structure: StructureData): { valid: boolean; errors: string[] } {
    const errors: string[] = [];

    // Check lattice
    if (!structure.lattice) {
        errors.push('Missing lattice data');
    } else {
        if (structure.lattice.a <= 0 || structure.lattice.b <= 0 || structure.lattice.c <= 0) {
            errors.push(`Invalid lattice parameters: a=${structure.lattice.a}, b=${structure.lattice.b}, c=${structure.lattice.c}`);
        }
        if (structure.lattice.alpha <= 0 || structure.lattice.beta <= 0 || structure.lattice.gamma <= 0) {
            errors.push(`Invalid lattice angles: α=${structure.lattice.alpha}, β=${structure.lattice.beta}, γ=${structure.lattice.gamma}`);
        }
    }

    // Check atoms
    if (!structure.atoms || structure.atoms.length === 0) {
        errors.push('No atoms in structure');
    } else {
        structure.atoms.forEach((atom, i) => {
            if (!atom.element) {
                errors.push(`Atom ${i}: missing element`);
            }
            if (!atom.coords || atom.coords.length !== 3) {
                errors.push(`Atom ${i}: invalid coords`);
            }
        });
    }

    const valid = errors.length === 0;
    if (!valid) {
        debugError('CONVERTERS', errors.join('; '), 'Structure validation failed');
    }

    return { valid, errors };
}

/**
 * Convert cartesian coordinates to fractional using lattice matrix
 */
function cartesianToFractional(
    cartesian: [number, number, number],
    lattice: LatticeData
): [number, number, number] {
    // lattice.matrix is the cell matrix where columns are lattice vectors
    // To convert cartesian to fractional: fractional = inverse(matrix) * cartesian

    const matrix = lattice.matrix;
    if (!matrix || matrix.length !== 3 || matrix[0].length !== 3) {
        debug('CONVERTERS', '⚠️ Invalid lattice matrix, using coords as-is');
        return cartesian;
    }

    // Compute inverse of 3x3 matrix (for small matrices, direct formula is fine)
    const det =
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    if (Math.abs(det) < 1e-10) {
        debug('CONVERTERS', '⚠️ Singular lattice matrix, using coords as-is');
        return cartesian;
    }

    const invDet = 1.0 / det;
    const inv: number[][] = [
        [
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * invDet,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invDet,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invDet,
        ],
        [
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invDet,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invDet,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * invDet,
        ],
        [
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * invDet,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * invDet,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * invDet,
        ],
    ];

    const fractional: [number, number, number] = [
        inv[0][0] * cartesian[0] + inv[0][1] * cartesian[1] + inv[0][2] * cartesian[2],
        inv[1][0] * cartesian[0] + inv[1][1] * cartesian[1] + inv[1][2] * cartesian[2],
        inv[2][0] * cartesian[0] + inv[2][1] * cartesian[1] + inv[2][2] * cartesian[2],
    ];

    return fractional;
}

/**
 * Ensure atom has fractional coordinates, converting from cartesian if needed
 */
function ensureFractionalCoords(atom: AtomData, lattice: LatticeData): [number, number, number] {
    // If coords are already fractional (typically in [0,1) range), use them
    // Check if all coords are in reasonable fractional range
    const coordsAreFractional = atom.coords.every(c => c >= -0.5 && c <= 1.5);

    if (coordsAreFractional) {
        // Wrap to [0, 1) range
        const wrapped: [number, number, number] = [
            ((atom.coords[0] % 1) + 1) % 1,
            ((atom.coords[1] % 1) + 1) % 1,
            ((atom.coords[2] % 1) + 1) % 1,
        ];
        return wrapped;
    }

    // If we have explicit cartesian coords, convert them
    if (atom.cartesian) {
        debug('CONVERTERS', `  Converting cartesian to fractional for ${atom.element}`);
        return cartesianToFractional(atom.cartesian, lattice);
    }

    // Otherwise, assume coords are cartesian and convert
    debug('CONVERTERS', `  Assuming coords are cartesian for ${atom.element}, converting...`);
    return cartesianToFractional(atom.coords, lattice);
}

/**
 * Convert MCP structure to CIF format for MolStar
 */
export function structureToCif(structure: StructureData): string {
    debug('CONVERTERS', '='.repeat(60));
    debug('CONVERTERS', '🔄 CONVERTING STRUCTURE TO CIF FORMAT');
    debug('CONVERTERS', `Formula: ${structure.metadata?.formula || 'Unknown'}`);
    debug('CONVERTERS', `Atom count: ${structure.atoms.length}`);
    debug('CONVERTERS', `Lattice: a=${structure.lattice.a.toFixed(3)}, b=${structure.lattice.b.toFixed(3)}, c=${structure.lattice.c.toFixed(3)}`);
    debug('CONVERTERS', `Angles: α=${structure.lattice.alpha.toFixed(2)}°, β=${structure.lattice.beta.toFixed(2)}°, γ=${structure.lattice.gamma.toFixed(2)}°`);

    // Validate structure data
    const validation = validateStructureData(structure);
    if (!validation.valid) {
        debugError('CONVERTERS', validation.errors, '❌ Structure validation failed');
        debug('CONVERTERS', 'Attempting conversion anyway...');
    } else {
        debug('CONVERTERS', '✅ Structure validation passed');
    }

    // Log first few atoms for debugging
    debug('CONVERTERS', 'First 3 atoms:');
    structure.atoms.slice(0, 3).forEach((atom, i) => {
        debug('CONVERTERS', `  [${i}] ${atom.element}: coords=${JSON.stringify(atom.coords)}, cartesian=${JSON.stringify(atom.cartesian)}`);
    });

    const lines: string[] = [];
    const formula = structure.metadata?.formula || 'Unknown';

    // CIF header
    lines.push(`data_${formula.replace(/\s/g, '_')}`);
    lines.push(`_chemical_formula_structural '${formula}'`);
    lines.push(`_cell_length_a ${structure.lattice.a.toFixed(4)}`);
    lines.push(`_cell_length_b ${structure.lattice.b.toFixed(4)}`);
    lines.push(`_cell_length_c ${structure.lattice.c.toFixed(4)}`);
    lines.push(`_cell_angle_alpha ${structure.lattice.alpha.toFixed(2)}`);
    lines.push(`_cell_angle_beta ${structure.lattice.beta.toFixed(2)}`);
    lines.push(`_cell_angle_gamma ${structure.lattice.gamma.toFixed(2)}`);

    // Space group
    if (structure.space_group) {
        lines.push(`_symmetry_space_group_name_H-M '${structure.space_group.symbol}'`);
        lines.push(`_symmetry_Int_Tables_number ${structure.space_group.number}`);
        debug('CONVERTERS', `Space group: ${structure.space_group.symbol} (#${structure.space_group.number})`);
    } else {
        lines.push(`_symmetry_space_group_name_H-M 'P1'`);
        lines.push(`_symmetry_Int_Tables_number 1`);
        debug('CONVERTERS', 'Space group: P1 (default)');
    }

    // Atom loop - CRITICAL: Use fractional coordinates for CIF
    lines.push('');
    lines.push('loop_');
    lines.push('_atom_site_label');
    lines.push('_atom_site_type_symbol');
    lines.push('_atom_site_fract_x');
    lines.push('_atom_site_fract_y');
    lines.push('_atom_site_fract_z');
    lines.push('_atom_site_id');

    debug('CONVERTERS', 'Converting atoms to CIF format with fractional coordinates...');
    structure.atoms.forEach((atom, i) => {
        const label = `${atom.element}${i + 1}`;
        const fractCoords = ensureFractionalCoords(atom, structure.lattice);
        const [x, y, z] = fractCoords;

        lines.push(`${label} ${atom.element} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)} ${i + 1}`);

        // Log first few conversions
        if (i < 3) {
            debug('CONVERTERS', `  ${label}: fract=(${x.toFixed(4)}, ${y.toFixed(4)}, ${z.toFixed(4)})`);
        }
    });

    const cifData = lines.join('\n');
    debug('CONVERTERS', `✅ CIF generated: ${cifData.length} bytes, ${lines.length} lines`);
    debug('CONVERTERS', 'CIF preview (first 500 chars):');
    debug('CONVERTERS', cifData.substring(0, 500));
    debug('CONVERTERS', '='.repeat(60));

    return cifData;
}

/**
 * Convert MCP structure to PDB format
 */
export function structureToPdb(structure: StructureData): string {
    const lines: string[] = [];
    const formula = structure.metadata?.formula || 'MOL';

    // Header
    lines.push(`HEADER    ${formula.padEnd(40)} ${new Date().toISOString().slice(0, 10)}`);
    lines.push(`TITLE     Generated by Crystal GUI`);

    // Crystal unit cell (if applicable)
    const lat = structure.lattice;
    lines.push(
        `CRYST1${lat.a.toFixed(3).padStart(9)}${lat.b.toFixed(3).padStart(9)}${lat.c.toFixed(3).padStart(9)}` +
        `${lat.alpha.toFixed(2).padStart(7)}${lat.beta.toFixed(2).padStart(7)}${lat.gamma.toFixed(2).padStart(7)} P 1           1`
    );

    // Atoms
    structure.atoms.forEach((atom, i) => {
        const serial = (i + 1).toString().padStart(5);
        const name = atom.element.padEnd(4);
        const resName = 'MOL';
        const chainId = 'A';
        const resSeq = '1'.padStart(4);

        // Use cartesian coordinates for PDB
        const coords = atom.cartesian || atom.coords;
        const x = coords[0].toFixed(3).padStart(8);
        const y = coords[1].toFixed(3).padStart(8);
        const z = coords[2].toFixed(3).padStart(8);

        const occupancy = '1.00';
        const tempFactor = '0.00';
        const element = atom.element.padStart(2);

        lines.push(
            `ATOM  ${serial} ${name} ${resName} ${chainId}${resSeq}    ${x}${y}${z}  ${occupancy}  ${tempFactor}          ${element}`
        );
    });

    lines.push('END');
    return lines.join('\n');
}

/**
 * Convert MCP structure to XYZ format (simple Cartesian coordinates)
 */
export function structureToXyz(structure: StructureData): string {
    const lines: string[] = [];
    const natoms = structure.atoms.length;
    const formula = structure.metadata?.formula || 'Structure';

    lines.push(natoms.toString());
    lines.push(formula);

    structure.atoms.forEach(atom => {
        const coords = atom.cartesian || atom.coords;
        lines.push(`${atom.element} ${coords[0].toFixed(6)} ${coords[1].toFixed(6)} ${coords[2].toFixed(6)}`);
    });

    return lines.join('\n');
}

/**
 * Generate SMILES from structure (placeholder - would need RDKit integration)
 * For now, returns empty string as this requires actual chemistry library
 */
export function structureToSmiles(_structure: StructureData): string {
    // TODO: Integrate with RDKit or Open Babel for proper conversion
    console.warn('SMILES conversion requires RDKit - returning empty string');
    return '';
}

/**
 * Detect bonds from atom positions based on covalent radii
 */
export function detectBonds(atoms: AtomData[]): Array<[number, number]> {
    const bonds: Array<[number, number]> = [];

    // Covalent radii in Angstroms (approximate)
    const covalentRadii: Record<string, number> = {
        H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, P: 1.07, S: 1.05,
        Cl: 1.02, Br: 1.20, I: 1.39, Si: 1.11, Se: 1.20, Te: 1.38,
        B: 0.84, Al: 1.21, Ga: 1.22, In: 1.42, Tl: 1.45,
        Ge: 1.20, Sn: 1.39, Pb: 1.46, As: 1.19, Sb: 1.39, Bi: 1.48,
        Fe: 1.32, Co: 1.26, Ni: 1.24, Cu: 1.32, Zn: 1.22,
        Ru: 1.46, Rh: 1.42, Pd: 1.39, Ag: 1.45, Pt: 1.36, Au: 1.36,
    };

    const tolerance = 0.4; // Bond tolerance

    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            const a1 = atoms[i];
            const a2 = atoms[j];

            const c1 = a1.cartesian || a1.coords;
            const c2 = a2.cartesian || a2.coords;

            const dx = c1[0] - c2[0];
            const dy = c1[1] - c2[1];
            const dz = c1[2] - c2[2];
            const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

            const r1 = covalentRadii[a1.element] || 1.5;
            const r2 = covalentRadii[a2.element] || 1.5;
            const maxBondLength = r1 + r2 + tolerance;

            if (distance < maxBondLength && distance > 0.4) {
                bonds.push([i, j]);
            }
        }
    }

    return bonds;
}

/**
 * Calculate center of mass for a structure
 */
export function calculateCenterOfMass(atoms: AtomData[]): [number, number, number] {
    if (atoms.length === 0) return [0, 0, 0];

    let totalX = 0, totalY = 0, totalZ = 0;

    atoms.forEach(atom => {
        const coords = atom.cartesian || atom.coords;
        totalX += coords[0];
        totalY += coords[1];
        totalZ += coords[2];
    });

    const n = atoms.length;
    return [totalX / n, totalY / n, totalZ / n];
}

/**
 * Calculate bounding box for structure
 */
export function calculateBoundingBox(atoms: AtomData[]): { min: [number, number, number]; max: [number, number, number] } {
    if (atoms.length === 0) {
        return { min: [0, 0, 0], max: [0, 0, 0] };
    }

    const coords = atoms.map(a => a.cartesian || a.coords);

    const min: [number, number, number] = [
        Math.min(...coords.map(c => c[0])),
        Math.min(...coords.map(c => c[1])),
        Math.min(...coords.map(c => c[2])),
    ];

    const max: [number, number, number] = [
        Math.max(...coords.map(c => c[0])),
        Math.max(...coords.map(c => c[1])),
        Math.max(...coords.map(c => c[2])),
    ];

    return { min, max };
}
