import { test, expect } from '@playwright/test';
import fs from 'fs';

const FRONTEND_URL = 'http://localhost:5173';
const REPORT_FILE = 'tests/molecule_test/failure_investigation_report.json';

/**
 * FAILURE INVESTIGATION TEST SUITE
 * 
 * 100 tests specifically designed to probe the 7 root cause categories:
 * 1. Molecule lookup failures (formulas, exotic names, SMILES)
 * 2. Count handling (LLM passing wrong count)
 * 3. Invalid parameter handling (0, negative, extreme values)
 * 4. Pattern implementation gaps
 * 5. LLM tool selection issues
 * 6. Direction vector/axis handling
 * 7. Complex heterogeneous systems
 */

interface TestCase {
    id: number;
    prompt: string;
    category: string;
    root_cause: string;
    validation: { type: string; expected: Record<string, any> };
}

const FAILURE_INVESTIGATION_TESTS: TestCase[] = [
    // =========================================================================
    // GROUP 1: MOLECULE LOOKUP - Formula Resolution (1-20)
    // Root cause: LLM sends formula like C18H12 instead of name
    // =========================================================================
    { id: 1, prompt: "generate 2 C6H6 molecules 3.4 angstrom apart", category: "Formula-Lookup", root_cause: "Backend should resolve C6H6 to benzene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 2, prompt: "stack 2 C10H8 molecules along z", category: "Formula-Lookup", root_cause: "Backend should resolve C10H8 to naphthalene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 3, prompt: "create 2 C14H10 molecules 3.5 angstrom apart", category: "Formula-Lookup", root_cause: "C14H10 = anthracene or phenanthrene - ambiguous", validation: { type: 'count', expected: { count: 2 } } },
    { id: 4, prompt: "generate 2 C16H10 molecules stacked", category: "Formula-Lookup", root_cause: "C16H10 = pyrene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 5, prompt: "stack 2 C18H12 molecules 3.4 angstrom apart", category: "Formula-Lookup", root_cause: "C18H12 = tetracene/chrysene - LLM uses this", validation: { type: 'count', expected: { count: 2 } } },
    { id: 6, prompt: "create 2 C20H12 molecules along z", category: "Formula-Lookup", root_cause: "C20H12 = perylene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 7, prompt: "generate 2 C24H12 molecules 3.5 angstrom apart", category: "Formula-Lookup", root_cause: "C24H12 = coronene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 8, prompt: "stack 2 C8H8 molecules 7 angstrom apart", category: "Formula-Lookup", root_cause: "C8H8 = cubane or styrene - need disambiguation", validation: { type: 'count', expected: { count: 2 } } },
    { id: 9, prompt: "create 2 C20H20 molecules 12 angstrom apart", category: "Formula-Lookup", root_cause: "C20H20 = dodecahedrane", validation: { type: 'count', expected: { count: 2 } } },
    { id: 10, prompt: "generate 2 H2O molecules 3 angstrom apart", category: "Formula-Lookup", root_cause: "H2O should resolve to water", validation: { type: 'count', expected: { count: 2 } } },

    // Direct name tests - these names often fail
    { id: 11, prompt: "stack 2 tetracene molecules 3.4 angstrom apart", category: "Name-Lookup", root_cause: "Tetracene not in local DB, should use PubChem", validation: { type: 'count', expected: { count: 2 } } },
    { id: 12, prompt: "create 2 perylene molecules 3.5 angstrom apart", category: "Name-Lookup", root_cause: "Perylene lookup failure", validation: { type: 'count', expected: { count: 2 } } },
    { id: 13, prompt: "generate 2 triphenylene molecules stacked", category: "Name-Lookup", root_cause: "Triphenylene not resolving", validation: { type: 'count', expected: { count: 2 } } },
    { id: 14, prompt: "stack 2 chrysene molecules 3.4 angstrom apart", category: "Name-Lookup", root_cause: "Chrysene lookup failure", validation: { type: 'count', expected: { count: 2 } } },
    { id: 15, prompt: "create 2 pentacene molecules 3.5 angstrom apart", category: "Name-Lookup", root_cause: "Large PAH lookup", validation: { type: 'count', expected: { count: 2 } } },
    { id: 16, prompt: "generate 2 corannulene molecules 4 angstrom apart", category: "Name-Lookup", root_cause: "Bowl PAH lookup", validation: { type: 'count', expected: { count: 2 } } },
    { id: 17, prompt: "stack 2 TTF molecules 4 angstrom apart", category: "Name-Lookup", root_cause: "TTF abbreviation not recognized", validation: { type: 'count', expected: { count: 2 } } },
    { id: 18, prompt: "create 2 TCNQ molecules 4 angstrom apart", category: "Name-Lookup", root_cause: "TCNQ abbreviation lookup", validation: { type: 'count', expected: { count: 2 } } },
    { id: 19, prompt: "generate 2 fullerene molecules 15 angstrom apart", category: "Name-Lookup", root_cause: "Fullerene should map to C60", validation: { type: 'count', expected: { count: 2 } } },
    { id: 20, prompt: "stack 2 buckyball molecules 15 angstrom apart", category: "Name-Lookup", root_cause: "Buckyball = C60 alias", validation: { type: 'count', expected: { count: 2 } } },

    // =========================================================================
    // GROUP 2: COUNT HANDLING - LLM Sending Wrong Count (21-40)
    // Root cause: LLM sends count=1 or omits count
    // =========================================================================
    { id: 21, prompt: "generate exactly 2 benzene molecules stacked", category: "Count-Explicit", root_cause: "Explicit 'exactly 2' should enforce count=2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 22, prompt: "create a dimer of pyrene", category: "Count-Dimer", root_cause: "'Dimer' means exactly 2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 23, prompt: "generate a trimer of benzene 3.4 angstrom apart", category: "Count-Trimer", root_cause: "'Trimer' means exactly 3", validation: { type: 'count', expected: { count: 3 } } },
    { id: 24, prompt: "create a tetramer of naphthalene", category: "Count-Tetramer", root_cause: "'Tetramer' means exactly 4", validation: { type: 'count', expected: { count: 4 } } },
    { id: 25, prompt: "stack 2 molecules of adamantane 8 angstrom apart", category: "Count-Explicit", root_cause: "'2 molecules of X' should be count=2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 26, prompt: "generate three benzene molecules stacked", category: "Count-Word", root_cause: "'three' = 3, not 1", validation: { type: 'count', expected: { count: 3 } } },
    { id: 27, prompt: "create five water molecules along z", category: "Count-Word", root_cause: "'five' = 5", validation: { type: 'count', expected: { count: 5 } } },
    { id: 28, prompt: "stack ten methane molecules 4 angstrom apart", category: "Count-Word", root_cause: "'ten' = 10", validation: { type: 'count', expected: { count: 10 } } },
    { id: 29, prompt: "generate pair of pyrene molecules", category: "Count-Pair", root_cause: "'pair' = exactly 2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 30, prompt: "create couple of benzene molecules stacked", category: "Count-Colloquial", root_cause: "'couple' = ~2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 31, prompt: "stack 2 cubane at 7 angstrom", category: "Count-NoMolWord", root_cause: "'2 cubane' without 'molecules' should work", validation: { type: 'count', expected: { count: 2 } } },
    { id: 32, prompt: "generate 2 corannulene stacked", category: "Count-NoMolWord", root_cause: "'2 X' pattern", validation: { type: 'count', expected: { count: 2 } } },
    { id: 33, prompt: "create 2 of dodecahedrane 12 angstrom apart", category: "Count-Of", root_cause: "'2 of X' pattern", validation: { type: 'count', expected: { count: 2 } } },
    { id: 34, prompt: "stack both benzene and naphthalene molecules", category: "Count-Both", root_cause: "'both' implies 2 different types", validation: { type: 'count', expected: { count: 2 } } },
    { id: 35, prompt: "generate two identical benzene molecules stacked", category: "Count-Identical", root_cause: "'two identical X' = 2", validation: { type: 'count', expected: { count: 2 } } },
    { id: 36, prompt: "create double benzene stack", category: "Count-Double", root_cause: "'double X stack' = 2?", validation: { type: 'count', expected: { count: 2 } } },
    { id: 37, prompt: "stack benzene times 2", category: "Count-Times", root_cause: "'X times 2' pattern", validation: { type: 'count', expected: { count: 2 } } },
    { id: 38, prompt: "generate 2x pyrene stack", category: "Count-2x", root_cause: "'2x X' pattern", validation: { type: 'count', expected: { count: 2 } } },
    { id: 39, prompt: "create benzene x2 stacked", category: "Count-x2", root_cause: "'X x2' pattern", validation: { type: 'count', expected: { count: 2 } } },
    { id: 40, prompt: "stack benzene, benzene along z 3.4 angstrom apart", category: "Count-Repeat", root_cause: "'X, X' = 2 molecules", validation: { type: 'count', expected: { count: 2 } } },

    // =========================================================================
    // GROUP 3: INVALID PARAMETER HANDLING (41-55)
    // Root cause: System crashes on 0, negative, extreme values
    // =========================================================================
    { id: 41, prompt: "generate 2 benzene 0 angstrom apart", category: "Invalid-Zero", root_cause: "0Å should use safe default, not crash", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 42, prompt: "stack 2 pyrene at -5 angstrom", category: "Invalid-Negative", root_cause: "Negative distance should be rejected/abs'd", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 43, prompt: "create 2 benzene 0.1 angstrom apart", category: "Invalid-TooSmall", root_cause: "0.1Å < VdW contact, must adjust", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 44, prompt: "generate 2 adamantane 1 angstrom apart", category: "Invalid-TooSmall", root_cause: "1Å for 5Å cage = overlap", validation: { type: 'no_clash', expected: { min: 5.0 } } },
    { id: 45, prompt: "stack 2 coronene 0.5 angstrom apart", category: "Invalid-TooSmall", root_cause: "0.5Å for large PAH = impossible", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 46, prompt: "create 2 C60 1 angstrom apart", category: "Invalid-TooSmall", root_cause: "C60 7Å diameter, 1Å = inside", validation: { type: 'no_clash', expected: { min: 7.0 } } },
    { id: 47, prompt: "generate 1000 benzene molecules stacked", category: "Invalid-TooMany", root_cause: "1000 molecules should work or error gracefully", validation: { type: 'count', expected: { min: 2 } } },
    { id: 48, prompt: "stack 0 benzene molecules", category: "Invalid-Zero", root_cause: "0 count should error gracefully", validation: { type: 'count', expected: { min: 1 } } },
    { id: 49, prompt: "generate -2 benzene molecules", category: "Invalid-Negative", root_cause: "Negative count should error", validation: { type: 'count', expected: { min: 1 } } },
    { id: 50, prompt: "create 2 benzene 10000 angstrom apart", category: "Invalid-TooLarge", root_cause: "10000Å is valid but unusual", validation: { type: 'distance', expected: { min: 9000, max: 11000 } } },
    { id: 51, prompt: "stack 2 benzene along q axis", category: "Invalid-Axis", root_cause: "'q axis' doesn't exist, should default to z", validation: { type: 'count', expected: { count: 2 } } },
    { id: 52, prompt: "generate 2 benzene along axis 4", category: "Invalid-Axis", root_cause: "'axis 4' not valid, should default", validation: { type: 'count', expected: { count: 2 } } },
    { id: 53, prompt: "create 2 benzene along [0,0,0] 5 angstrom apart", category: "Invalid-Vector", root_cause: "[0,0,0] is null vector, should default to z", validation: { type: 'count', expected: { count: 2 } } },
    { id: 54, prompt: "stack benzene at NaN angstrom", category: "Invalid-NaN", root_cause: "NaN should be detected and defaulted", validation: { type: 'count', expected: { min: 1 } } },
    { id: 55, prompt: "generate 2 benzene at infinity angstrom", category: "Invalid-Inf", root_cause: "'infinity' should error gracefully", validation: { type: 'count', expected: { min: 1 } } },

    // =========================================================================
    // GROUP 4: PATTERN IMPLEMENTATION (56-70)
    // Root cause: Patterns mentioned but not implemented
    // =========================================================================
    { id: 56, prompt: "create helical stack of 5 benzene with 3.4 pitch", category: "Pattern-Helix", root_cause: "Helix pattern implementation", validation: { type: 'count', expected: { count: 5 } } },
    { id: 57, prompt: "generate zigzag chain of 4 benzene", category: "Pattern-Zigzag", root_cause: "Zigzag pattern implementation", validation: { type: 'count', expected: { count: 4 } } },
    { id: 58, prompt: "create tetrahedral arrangement of 4 benzene", category: "Pattern-Tetrahedral", root_cause: "3D polyhedra pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 59, prompt: "generate octahedral arrangement of 6 benzene", category: "Pattern-Octahedral", root_cause: "6-vertex polyhedra", validation: { type: 'count', expected: { count: 6 } } },
    { id: 60, prompt: "create T-shaped benzene dimer", category: "Pattern-TShape", root_cause: "T-shaped edge-face interaction", validation: { type: 'count', expected: { count: 2 } } },
    { id: 61, prompt: "generate herringbone pattern of 4 benzene", category: "Pattern-Herringbone", root_cause: "45° alternating tilt", validation: { type: 'count', expected: { count: 4 } } },
    { id: 62, prompt: "create sandwich structure with benzene around ferrocene", category: "Pattern-Sandwich", root_cause: "3-layer A-B-A structure", validation: { type: 'count', expected: { count: 3 } } },
    { id: 63, prompt: "generate slipped-parallel benzene dimer", category: "Pattern-SlippedParallel", root_cause: "Offset parallel stacking", validation: { type: 'count', expected: { count: 2 } } },
    { id: 64, prompt: "create face-to-edge benzene arrangement", category: "Pattern-FaceEdge", root_cause: "Mixed orientation", validation: { type: 'count', expected: { count: 2 } } },
    { id: 65, prompt: "generate antiparallel benzene dimer", category: "Pattern-Antiparallel", root_cause: "180° rotated stacking", validation: { type: 'count', expected: { count: 2 } } },
    { id: 66, prompt: "create circular ring of 6 benzene", category: "Pattern-Ring", root_cause: "Molecules in a circle", validation: { type: 'count', expected: { count: 6 } } },
    { id: 67, prompt: "generate 3x3 grid of benzene", category: "Pattern-Grid", root_cause: "2D grid arrangement", validation: { type: 'count', expected: { count: 9 } } },
    { id: 68, prompt: "create cubic arrangement of 8 benzene", category: "Pattern-Cubic", root_cause: "Molecules at cube vertices", validation: { type: 'count', expected: { count: 8 } } },
    { id: 69, prompt: "generate alternating up-down benzene chain", category: "Pattern-Alternating", root_cause: "Alternating orientations", validation: { type: 'count', expected: { min: 2 } } },
    { id: 70, prompt: "create star arrangement of 5 benzene around 1 center", category: "Pattern-Star", root_cause: "Central + surrounding pattern", validation: { type: 'count', expected: { count: 6 } } },

    // =========================================================================
    // GROUP 5: DIRECTION VECTOR HANDLING (71-85)
    // Root cause: Direction vectors not parsed correctly
    // =========================================================================
    { id: 71, prompt: "stack 2 benzene along direction [1,0,0]", category: "Direction-UnitX", root_cause: "[1,0,0] should work like axis='x'", validation: { type: 'count', expected: { count: 2 } } },
    { id: 72, prompt: "generate 2 benzene along direction [0,1,0]", category: "Direction-UnitY", root_cause: "[0,1,0] should work like axis='y'", validation: { type: 'count', expected: { count: 2 } } },
    { id: 73, prompt: "create 2 benzene along direction [0,0,1]", category: "Direction-UnitZ", root_cause: "[0,0,1] should work like axis='z'", validation: { type: 'count', expected: { count: 2 } } },
    { id: 74, prompt: "stack 2 benzene along [1,1,0] direction 5 angstrom apart", category: "Direction-Diagonal", root_cause: "Diagonal in XY plane", validation: { type: 'count', expected: { count: 2 } } },
    { id: 75, prompt: "generate 2 benzene along [1,1,1] direction 5 angstrom apart", category: "Direction-BodyDiag", root_cause: "Body diagonal direction", validation: { type: 'count', expected: { count: 2 } } },
    { id: 76, prompt: "create 2 benzene along [1,0,1] direction 5 angstrom apart", category: "Direction-FaceDiag", root_cause: "Face diagonal XZ", validation: { type: 'count', expected: { count: 2 } } },
    { id: 77, prompt: "stack 2 benzene along [0,1,1] direction 5 angstrom apart", category: "Direction-FaceDiag", root_cause: "Face diagonal YZ", validation: { type: 'count', expected: { count: 2 } } },
    { id: 78, prompt: "generate 2 benzene at 45 degrees in xy plane", category: "Direction-AngleXY", root_cause: "Angle-based direction", validation: { type: 'count', expected: { count: 2 } } },
    { id: 79, prompt: "create 2 benzene at 30 degrees from z axis", category: "Direction-AngleZ", root_cause: "Elevation angle from z", validation: { type: 'count', expected: { count: 2 } } },
    { id: 80, prompt: "stack 2 benzene perpendicular to x axis", category: "Direction-Perp", root_cause: "Perpendicular = yz plane", validation: { type: 'count', expected: { count: 2 } } },
    { id: 81, prompt: "generate 2 benzene perpendicular to y axis", category: "Direction-Perp", root_cause: "Perpendicular = xz plane", validation: { type: 'count', expected: { count: 2 } } },
    { id: 82, prompt: "create 2 benzene perpendicular to z axis", category: "Direction-Perp", root_cause: "Perpendicular = xy plane", validation: { type: 'count', expected: { count: 2 } } },
    { id: 83, prompt: "stack 2 benzene horizontally 5 angstrom apart", category: "Direction-Natural", root_cause: "'Horizontally' = x or y", validation: { type: 'count', expected: { count: 2 } } },
    { id: 84, prompt: "generate 2 benzene vertically 5 angstrom apart", category: "Direction-Natural", root_cause: "'Vertically' = z", validation: { type: 'count', expected: { count: 2 } } },
    { id: 85, prompt: "create 2 benzene sideways 5 angstrom apart", category: "Direction-Natural", root_cause: "'Sideways' = x or y", validation: { type: 'count', expected: { count: 2 } } },

    // =========================================================================
    // GROUP 6: SMILES/CID INPUT (86-95)
    // Root cause: SMILES and CID not resolving correctly
    // =========================================================================
    { id: 86, prompt: "generate 2 c1ccccc1 molecules stacked", category: "SMILES-Benzene", root_cause: "Benzene SMILES should work", validation: { type: 'count', expected: { count: 2 } } },
    { id: 87, prompt: "create 2 c1ccc2ccccc2c1 molecules 3.5 angstrom apart", category: "SMILES-Naphthalene", root_cause: "Naphthalene SMILES", validation: { type: 'count', expected: { count: 2 } } },
    { id: 88, prompt: "stack 2 CC(=O)O molecules 4 angstrom apart", category: "SMILES-AceticAcid", root_cause: "Acetic acid SMILES", validation: { type: 'count', expected: { count: 2 } } },
    { id: 89, prompt: "generate 2 CCO molecules 4 angstrom apart", category: "SMILES-Ethanol", root_cause: "Ethanol SMILES", validation: { type: 'count', expected: { count: 2 } } },
    { id: 90, prompt: "create 2 c1ccncc1 molecules stacked", category: "SMILES-Pyridine", root_cause: "Pyridine SMILES", validation: { type: 'count', expected: { count: 2 } } },
    { id: 91, prompt: "generate 2 CID:241 molecules 3.4 angstrom apart", category: "CID-Benzene", root_cause: "CID 241 = benzene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 92, prompt: "stack 2 PubChem:241 molecules", category: "CID-PubChem", root_cause: "PubChem:CID format", validation: { type: 'count', expected: { count: 2 } } },
    { id: 93, prompt: "create 2 [C]1[C][C][C][C][C]1 molecules stacked", category: "SMILES-Explicit", root_cause: "Explicit carbon SMILES", validation: { type: 'count', expected: { count: 2 } } },
    { id: 94, prompt: "generate 2 InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H molecules", category: "InChI-Benzene", root_cause: "InChI input format", validation: { type: 'count', expected: { count: 2 } } },
    { id: 95, prompt: "stack 2 O=C=O molecules 4 angstrom apart", category: "SMILES-CO2", root_cause: "CO2 SMILES", validation: { type: 'count', expected: { count: 2 } } },

    // =========================================================================
    // GROUP 7: LLM TOOL SELECTION (96-100)
    // Root cause: LLM picks wrong tool or no tool
    // =========================================================================
    { id: 96, prompt: "generate 2 H2 molecules 3 angstrom apart", category: "LLM-H2", root_cause: "H2 is dimer of H, needs cluster tool", validation: { type: 'count', expected: { count: 2 } } },
    { id: 97, prompt: "create 2 N2 molecules 4 angstrom apart", category: "LLM-N2", root_cause: "N2 diatomic dimer", validation: { type: 'count', expected: { count: 2 } } },
    { id: 98, prompt: "stack 2 O2 molecules 4 angstrom apart", category: "LLM-O2", root_cause: "O2 diatomic dimer", validation: { type: 'count', expected: { count: 2 } } },
    { id: 99, prompt: "generate molecular dimer of CO", category: "LLM-CO", root_cause: "CO diatomic 'dimer' = 2 CO molecules", validation: { type: 'count', expected: { count: 2 } } },
    { id: 100, prompt: "create stacked pair of chlorine molecules", category: "LLM-Cl2", root_cause: "Cl2 diatomic pair", validation: { type: 'count', expected: { count: 2 } } },
];

// Helper functions
function calculateCentroid(coords: number[][]): number[] {
    if (coords.length === 0) return [0, 0, 0];
    const sum = [0, 0, 0];
    for (const c of coords) { sum[0] += c[0]; sum[1] += c[1]; sum[2] += c[2]; }
    return [sum[0] / coords.length, sum[1] / coords.length, sum[2] / coords.length];
}

function calculateDistance(c1: number[], c2: number[]): number {
    return Math.sqrt(Math.pow(c1[0] - c2[0], 2) + Math.pow(c1[1] - c2[1], 2) + Math.pow(c1[2] - c2[2], 2));
}

function validateTest(testCase: TestCase, atoms: any[], metadata: any): { valid: boolean; details: string; measurements: any } {
    const validation = testCase.validation;
    const measurements: any = {};
    const moleculeCount = metadata?.n_molecules || 1;

    measurements.molecule_count = moleculeCount;
    measurements.total_atoms = atoms.length;

    if (moleculeCount >= 2) {
        const atomsPerMol = Math.floor(atoms.length / moleculeCount);
        const centroids: number[][] = [];
        for (let i = 0; i < Math.min(moleculeCount, 2); i++) {
            const molCoords = atoms.slice(i * atomsPerMol, (i + 1) * atomsPerMol).map((a: any) => a.cartesian || a.coords || [0, 0, 0]);
            if (molCoords.length > 0) centroids.push(calculateCentroid(molCoords));
        }
        if (centroids.length >= 2) measurements.actual_distance = calculateDistance(centroids[0], centroids[1]);
    }

    switch (validation.type) {
        case 'count':
            const minCount = validation.expected.min ?? validation.expected.count;
            const exactCount = validation.expected.count;
            if (exactCount && moleculeCount !== exactCount) {
                return { valid: false, details: `Expected ${exactCount} molecules, got ${moleculeCount}`, measurements };
            }
            if (minCount && moleculeCount < minCount) {
                return { valid: false, details: `Expected ≥${minCount} molecules, got ${moleculeCount}`, measurements };
            }
            return { valid: true, details: `✓ ${moleculeCount} molecules`, measurements };

        case 'distance':
            if (!measurements.actual_distance) return { valid: false, details: 'Need 2+ molecules for distance', measurements };
            const dist = measurements.actual_distance;
            if (dist < validation.expected.min || dist > validation.expected.max) {
                return { valid: false, details: `Distance ${dist.toFixed(2)}Å outside [${validation.expected.min}-${validation.expected.max}]`, measurements };
            }
            return { valid: true, details: `✓ Distance ${dist.toFixed(2)}Å`, measurements };

        case 'no_clash':
            if (!measurements.actual_distance) return { valid: true, details: '✓ Single mol', measurements };
            if (measurements.actual_distance < validation.expected.min) {
                return { valid: false, details: `Clash: ${measurements.actual_distance.toFixed(2)}Å < ${validation.expected.min}Å`, measurements };
            }
            return { valid: true, details: `✓ No clash: ${measurements.actual_distance.toFixed(2)}Å`, measurements };

        default:
            return { valid: atoms.length > 0, details: `${atoms.length} atoms`, measurements };
    }
}

test('Failure Investigation - 100 Root Cause Tests', async ({ page }) => {
    let results: any[] = [];
    if (fs.existsSync(REPORT_FILE)) {
        try { results = JSON.parse(fs.readFileSync(REPORT_FILE, 'utf-8')); } catch { }
    }

    test.setTimeout(14400000);
    await page.goto(FRONTEND_URL);
    await page.waitForLoadState('networkidle');

    const inputSelector = 'textarea, input[type="text"]';
    const tutorial = page.getByText("Get Started");
    if (await tutorial.isVisible({ timeout: 2000 }).catch(() => false)) {
        await tutorial.click();
        await page.waitForTimeout(300);
    }
    await page.waitForSelector(inputSelector);

    for (const tc of FAILURE_INVESTIGATION_TESTS) {
        if (results.find(r => r.id === tc.id)) continue;

        console.log(`\n[${tc.id}] ${tc.category}: ${tc.root_cause.slice(0, 50)}`);

        const timing: any = { start: Date.now() };
        let reqPayload: any = null, resPayload: any = null, atoms: any[] = [], metadata: any = null;

        const reqP = page.waitForRequest(req => req.url().includes('/mcp/call') && req.method() === 'POST', { timeout: 120000 }).catch(() => null);
        const resP = page.waitForResponse(res => res.url().includes('/mcp/call'), { timeout: 120000 }).catch(() => null);

        await page.fill(inputSelector, '');
        await page.fill(inputSelector, tc.prompt);
        await page.press(inputSelector, 'Enter');
        timing.sent = Date.now();

        try {
            const req = await reqP;
            timing.llm_done = Date.now();
            if (req) reqPayload = req.postDataJSON();

            const res = await resP;
            timing.mcp_done = Date.now();
            if (res) try { resPayload = JSON.parse(await res.text()); } catch { }

            await page.waitForTimeout(500);

            let mcpSuccess = false;
            if (resPayload?.content) {
                for (const item of resPayload.content) {
                    if (item.text?.includes('<json-data>')) {
                        const m = item.text.match(/<json-data>([\s\S]*?)<\/json-data>/);
                        if (m) {
                            try {
                                const d = JSON.parse(m[1]);
                                mcpSuccess = d.success === true;
                                atoms = d.structure?.atoms || d.atoms || [];
                                metadata = d.structure?.metadata || d.metadata || {};
                            } catch { }
                        }
                    }
                }
            }

            let testResult = { valid: false, details: 'No atoms', measurements: {} };
            if (atoms.length > 0) testResult = validateTest(tc, atoms, metadata);

            const status = mcpSuccess && testResult.valid ? 'PASSED' : (mcpSuccess ? 'SEMANTIC_FAIL' : 'MCP_FAIL');
            timing.total_ms = Date.now() - timing.start;
            timing.llm_ms = (timing.llm_done || timing.start) - timing.sent;
            timing.mcp_ms = (timing.mcp_done || timing.start) - (timing.llm_done || timing.start);

            results.push({
                id: tc.id, category: tc.category, root_cause: tc.root_cause, prompt: tc.prompt,
                llm_tool: reqPayload?.name || null, llm_args: reqPayload?.arguments || null,
                mcp_success: mcpSuccess, atom_count: atoms.length, metadata,
                test_valid: testResult.valid, test_details: testResult.details,
                measurements: testResult.measurements, position_info: metadata?.position_info || null,
                timing, status, timestamp: new Date().toISOString(),
            });

            console.log(`  ${status}: ${testResult.details} [LLM:${timing.llm_ms}ms MCP:${timing.mcp_ms}ms]`);
        } catch (e) {
            results.push({
                id: tc.id, category: tc.category, root_cause: tc.root_cause, prompt: tc.prompt,
                error: (e as Error).message, timing: { total_ms: Date.now() - timing.start }, status: 'ERROR', timestamp: new Date().toISOString()
            });
        }

        fs.writeFileSync(REPORT_FILE, JSON.stringify(results, null, 2));
    }

    // Summary by root cause
    console.log('\n' + '='.repeat(70));
    console.log('FAILURE INVESTIGATION COMPLETE');
    console.log('='.repeat(70));

    const categories: Record<string, { passed: number; total: number }> = {};
    for (const r of results) {
        const cat = r.category.split('-')[0];
        if (!categories[cat]) categories[cat] = { passed: 0, total: 0 };
        categories[cat].total++;
        if (r.status === 'PASSED') categories[cat].passed++;
    }

    console.log('\nResults by Category:');
    for (const [cat, data] of Object.entries(categories)) {
        console.log(`  ${cat}: ${data.passed}/${data.total} (${(data.passed / data.total * 100).toFixed(0)}%)`);
    }

    const p = results.filter(r => r.status === 'PASSED').length;
    console.log(`\nTotal: ${p}/${results.length} PASSED (${(p / results.length * 100).toFixed(1)}%)`);
});
