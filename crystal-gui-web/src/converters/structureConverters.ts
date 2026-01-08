/**
 * Structure Converters - Convert MCP output to viewer formats
 *
 * Converts structure data from Python backend JSON to formats
 * suitable for MolStar (CIF, PDB) and Kekule.js (Mol2, SMILES)
 *
 * DEBUG: Comprehensive validation and logging enabled
 *
 * IMPORTANT: Backend always provides:
 *   - atom.coords = FRACTIONAL coordinates (already in [0,1] range)
 *   - atom.cartesian = CARTESIAN coordinates (in Angstroms)
 * No heuristic detection needed - trust the backend data!
 */

import type { StructureData, AtomData, LatticeData } from '../types';
import { debug, debugError } from '../debug';

/**
 * Validate structure data before conversion
 */
function validateStructureData(structure: StructureData): { valid: boolean; errors: string[] } {
    const errors: string[] = [];

    debug('CONVERTERS', '📋 VALIDATING STRUCTURE DATA');

    // Check lattice
    if (!structure.lattice) {
        errors.push('Missing lattice data');
        debug('CONVERTERS', '  ❌ Missing lattice');
    } else {
        debug('CONVERTERS', `  Lattice: a=${structure.lattice.a}, b=${structure.lattice.b}, c=${structure.lattice.c}`);
        debug('CONVERTERS', `  Angles: α=${structure.lattice.alpha}, β=${structure.lattice.beta}, γ=${structure.lattice.gamma}`);
        debug('CONVERTERS', `  Matrix present: ${!!structure.lattice.matrix}`);

        if (structure.lattice.a <= 0 || structure.lattice.b <= 0 || structure.lattice.c <= 0) {
            errors.push(`Invalid lattice parameters: a=${structure.lattice.a}, b=${structure.lattice.b}, c=${structure.lattice.c}`);
            debug('CONVERTERS', '  ❌ Invalid lattice parameters (<=0)');
        }
        if (structure.lattice.alpha <= 0 || structure.lattice.beta <= 0 || structure.lattice.gamma <= 0) {
            errors.push(`Invalid lattice angles: α=${structure.lattice.alpha}, β=${structure.lattice.beta}, γ=${structure.lattice.gamma}`);
            debug('CONVERTERS', '  ❌ Invalid lattice angles (<=0)');
        }
    }

    // Check atoms
    if (!structure.atoms || structure.atoms.length === 0) {
        errors.push('No atoms in structure');
        debug('CONVERTERS', '  ❌ No atoms');
    } else {
        debug('CONVERTERS', `  Atom count: ${structure.atoms.length}`);
        structure.atoms.forEach((atom, i) => {
            if (!atom.element) {
                errors.push(`Atom ${i}: missing element`);
                debug('CONVERTERS', `  ❌ Atom ${i}: missing element`);
            }
            if (!atom.coords || atom.coords.length !== 3) {
                errors.push(`Atom ${i}: invalid coords`);
                debug('CONVERTERS', `  ❌ Atom ${i}: invalid coords (${JSON.stringify(atom.coords)})`);
            }
            // Log first 3 atoms for debugging
            if (i < 3) {
                debug('CONVERTERS', `  Atom[${i}]: ${atom.element} coords=${JSON.stringify(atom.coords)} cartesian=${JSON.stringify(atom.cartesian)}`);
            }
        });
        if (structure.atoms.length > 3) {
            debug('CONVERTERS', `  ... and ${structure.atoms.length - 3} more atoms`);
        }
    }

    const valid = errors.length === 0;
    if (!valid) {
        debugError('CONVERTERS', errors.join('; '), 'Structure validation failed');
    } else {
        debug('CONVERTERS', '  ✅ Validation PASSED');
    }

    return { valid, errors };
}

/**
 * Get Cartesian coordinates from atom
 * Backend ALWAYS provides cartesian coords - use them directly
 * Falls back to fractional->cartesian conversion only if cartesian missing
 */
function getCartesianCoords(atom: AtomData, lattice: LatticeData): [number, number, number] {
    // Backend provides cartesian coordinates directly - use them!
    if (atom.cartesian && atom.cartesian.length === 3) {
        return atom.cartesian;
    }

    // Fallback: convert fractional to cartesian using lattice matrix
    debug('CONVERTERS', `  ⚠️ Missing cartesian for ${atom.element}, computing from fractional`);
    return fractionalToCartesian(atom.coords, lattice);
}

/**
 * Convert fractional coordinates to Cartesian using lattice matrix
 * cartesian = matrix * fractional (matrix rows are lattice vectors)
 */
function fractionalToCartesian(
    fractional: [number, number, number],
    lattice: LatticeData
): [number, number, number] {
    const matrix = lattice.matrix;

    // If no matrix, use simple scaling by a,b,c (assumes orthorhombic)
    if (!matrix || matrix.length !== 3 || matrix[0].length !== 3) {
        debug('CONVERTERS', '  ⚠️ No lattice matrix, using simple a*x, b*y, c*z');
        return [
            fractional[0] * lattice.a,
            fractional[1] * lattice.b,
            fractional[2] * lattice.c,
        ];
    }

    // Proper transformation: cartesian = matrix^T * fractional
    // (matrix is stored row-wise, so we need matrix[col][row] for the math)
    const cartesian: [number, number, number] = [
        matrix[0][0] * fractional[0] + matrix[1][0] * fractional[1] + matrix[2][0] * fractional[2],
        matrix[0][1] * fractional[0] + matrix[1][1] * fractional[1] + matrix[2][1] * fractional[2],
        matrix[0][2] * fractional[0] + matrix[1][2] * fractional[1] + matrix[2][2] * fractional[2],
    ];

    return cartesian;
}

/**
 * Convert cartesian coordinates to fractional using lattice matrix
 * (Exported for potential future use in structure editing)
 */
export function cartesianToFractional(
    cartesian: [number, number, number],
    lattice: LatticeData
): [number, number, number] {
    // lattice.matrix is the cell matrix where rows are lattice vectors
    // To convert cartesian to fractional: fractional = inverse(matrix) * cartesian

    const matrix = lattice.matrix;
    if (!matrix || matrix.length !== 3 || matrix[0].length !== 3) {
        debug('CONVERTERS', '⚠️ Invalid lattice matrix, using simple division');
        return [
            cartesian[0] / lattice.a,
            cartesian[1] / lattice.b,
            cartesian[2] / lattice.c,
        ];
    }

    // Compute inverse of 3x3 matrix (for small matrices, direct formula is fine)
    const det =
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    if (Math.abs(det) < 1e-10) {
        debug('CONVERTERS', '⚠️ Singular lattice matrix, using simple division');
        return [
            cartesian[0] / lattice.a,
            cartesian[1] / lattice.b,
            cartesian[2] / lattice.c,
        ];
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
 * Get fractional coordinates from atom
 * Backend ALWAYS provides fractional coords - use them directly (no heuristics!)
 */
function getFractionalCoords(atom: AtomData): [number, number, number] {
    // Backend provides coords as fractional - use them directly!
    // No heuristic detection - trust the data source
    return atom.coords;
}

/**
 * Convert MCP structure to CIF format for MolStar
 *
 * Generates mmCIF-compatible format with BOTH fractional AND Cartesian coordinates
 * This ensures Mol* can properly render the structure regardless of which it prefers.
 */
export function structureToCif(structure: StructureData): string {
    debug('CONVERTERS', '═'.repeat(60));
    debug('CONVERTERS', '🔄 CONVERTING STRUCTURE TO CIF FORMAT');
    debug('CONVERTERS', `  Formula: ${structure.metadata?.formula || 'Unknown'}`);
    debug('CONVERTERS', `  Atom count: ${structure.atoms.length}`);
    debug('CONVERTERS', `  Lattice: a=${structure.lattice.a.toFixed(3)}, b=${structure.lattice.b.toFixed(3)}, c=${structure.lattice.c.toFixed(3)}`);
    debug('CONVERTERS', `  Angles: α=${structure.lattice.alpha.toFixed(2)}°, β=${structure.lattice.beta.toFixed(2)}°, γ=${structure.lattice.gamma.toFixed(2)}°`);

    // Validate structure data
    const validation = validateStructureData(structure);
    if (!validation.valid) {
        debugError('CONVERTERS', validation.errors, '❌ Structure validation failed');
        debug('CONVERTERS', '  ⚠️ Attempting conversion anyway...');
    }

    const lines: string[] = [];
    const formula = structure.metadata?.formula || 'Unknown';
    // Sanitize formula for CIF data block name (no spaces, special chars)
    const safeFormula = formula.replace(/[^a-zA-Z0-9_]/g, '_');

    // CIF header with data block
    lines.push(`data_${safeFormula}`);
    lines.push('#');

    // Chemical formula
    lines.push(`_chemical_formula_sum '${formula}'`);
    lines.push('#');

    // Cell parameters
    lines.push('# Unit cell parameters');
    lines.push(`_cell_length_a    ${structure.lattice.a.toFixed(4)}`);
    lines.push(`_cell_length_b    ${structure.lattice.b.toFixed(4)}`);
    lines.push(`_cell_length_c    ${structure.lattice.c.toFixed(4)}`);
    lines.push(`_cell_angle_alpha ${structure.lattice.alpha.toFixed(2)}`);
    lines.push(`_cell_angle_beta  ${structure.lattice.beta.toFixed(2)}`);
    lines.push(`_cell_angle_gamma ${structure.lattice.gamma.toFixed(2)}`);
    lines.push('#');

    // Space group
    const spaceGroupSymbol = structure.space_group?.symbol || 'P 1';
    const spaceGroupNumber = structure.space_group?.number || 1;
    lines.push('# Space group');
    lines.push(`_symmetry_space_group_name_H-M '${spaceGroupSymbol}'`);
    lines.push(`_symmetry_Int_Tables_number ${spaceGroupNumber}`);
    debug('CONVERTERS', `  Space group: ${spaceGroupSymbol} (#${spaceGroupNumber})`);
    lines.push('#');

    // Atom site loop with BOTH fractional AND Cartesian coordinates
    // This is key for Mol* compatibility - it can use whichever it prefers
    lines.push('# Atom coordinates');
    lines.push('loop_');
    lines.push('_atom_site_id');
    lines.push('_atom_site_type_symbol');
    lines.push('_atom_site_label');
    lines.push('_atom_site_fract_x');
    lines.push('_atom_site_fract_y');
    lines.push('_atom_site_fract_z');
    lines.push('_atom_site_Cartn_x');
    lines.push('_atom_site_Cartn_y');
    lines.push('_atom_site_Cartn_z');
    lines.push('_atom_site_occupancy');

    debug('CONVERTERS', '  Converting atoms to CIF format...');
    structure.atoms.forEach((atom, i) => {
        const id = i + 1;
        const label = `${atom.element}${id}`;

        // Get fractional coords (from backend, already fractional)
        const fractCoords = getFractionalCoords(atom);
        const [fx, fy, fz] = fractCoords;

        // Get Cartesian coords (from backend, already cartesian)
        const cartCoords = getCartesianCoords(atom, structure.lattice);
        const [cx, cy, cz] = cartCoords;

        // Format: id type_symbol label fract_x fract_y fract_z Cartn_x Cartn_y Cartn_z occupancy
        lines.push(
            `${id} ${atom.element} ${label} ` +
            `${fx.toFixed(6)} ${fy.toFixed(6)} ${fz.toFixed(6)} ` +
            `${cx.toFixed(4)} ${cy.toFixed(4)} ${cz.toFixed(4)} ` +
            `1.0`
        );

        // Log first few atoms for debugging
        if (i < 3) {
            debug('CONVERTERS', `  Atom[${i}]: ${atom.element} fract=(${fx.toFixed(4)}, ${fy.toFixed(4)}, ${fz.toFixed(4)}) cart=(${cx.toFixed(2)}, ${cy.toFixed(2)}, ${cz.toFixed(2)})`);
        }
    });

    if (structure.atoms.length > 3) {
        debug('CONVERTERS', `  ... and ${structure.atoms.length - 3} more atoms`);
    }

    lines.push('#');
    lines.push('# End of CIF');

    const cifData = lines.join('\n');
    debug('CONVERTERS', `  ✅ CIF generated: ${cifData.length} bytes, ${lines.length} lines`);
    debug('CONVERTERS', '  CIF content:');
    debug('CONVERTERS', '  ' + cifData.split('\n').slice(0, 25).join('\n  '));
    if (lines.length > 25) {
        debug('CONVERTERS', `  ... (${lines.length - 25} more lines)`);
    }
    debug('CONVERTERS', '═'.repeat(60));

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
