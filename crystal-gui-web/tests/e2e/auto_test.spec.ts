import { test, expect } from '@playwright/test';
import fs from 'fs';

const FRONTEND_URL = 'http://localhost:5173';
const REPORT_FILE = 'tests/molecule_test/failure_investigation_report.json';

/**
 * FAILURE INVESTIGATION TEST SUITE - v2 HARDER
 * 
 * 27 FAILED tests KEPT from previous run
 * 73 PASSED tests REPLACED with HARDER versions
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
    // GROUP 1: FORMULA LOOKUP - HARDER VERSIONS (1-10)
    // Original: 1,2,3,5,7,8,9,10 passed → now harder
    // =========================================================================
    { id: 1, prompt: "generate 3 C6H6 molecules in zigzag pattern 3.4 angstrom apart", category: "Formula-Hard", root_cause: "Formula + pattern combo", validation: { type: 'count', expected: { count: 3 } } },
    { id: 2, prompt: "stack 4 C10H8 molecules along [1,1,0] at 3.5 angstrom", category: "Formula-Hard", root_cause: "Formula + vector + count", validation: { type: 'count', expected: { count: 4 } } },
    { id: 3, prompt: "create 5 C14H10 molecules in herringbone pattern 3.5 angstrom apart", category: "Formula-Hard", root_cause: "C14H10 + pattern", validation: { type: 'count', expected: { count: 5 } } },
    // KEPT FAILED: id 4 - C16H10
    { id: 4, prompt: "generate 2 C16H10 molecules stacked", category: "Formula-Lookup", root_cause: "C16H10 = pyrene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 5, prompt: "create tetramer of C18H12 molecules antiparallel at 3.5 angstrom", category: "Formula-Hard", root_cause: "C18H12 + tetramer + pattern", validation: { type: 'count', expected: { count: 4 } } },
    // KEPT FAILED: id 6 - C20H12  
    { id: 6, prompt: "create 2 C20H12 molecules along z", category: "Formula-Lookup", root_cause: "C20H12 = perylene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 7, prompt: "generate trimer of C24H12 at [1,0,1] direction 4 angstrom", category: "Formula-Hard", root_cause: "C24H12 + trimer + vector", validation: { type: 'count', expected: { count: 3 } } },
    { id: 8, prompt: "stack 3 C8H8 molecules in T-shaped arrangement 7 angstrom", category: "Formula-Hard", root_cause: "C8H8 + count + pattern", validation: { type: 'count', expected: { count: 3 } } },
    { id: 9, prompt: "create 4 C20H20 molecules in tetrahedral arrangement 12 angstrom", category: "Formula-Hard", root_cause: "C20H20 + tetrahedral", validation: { type: 'count', expected: { count: 4 } } },
    { id: 10, prompt: "generate 5 H2O molecules in circular ring 3 angstrom apart", category: "Formula-Hard", root_cause: "H2O + ring pattern", validation: { type: 'count', expected: { count: 5 } } },

    // =========================================================================
    // GROUP 2: NAME LOOKUP - HARDER VERSIONS (11-20)
    // =========================================================================
    // KEPT FAILED: id 11 - tetracene
    { id: 11, prompt: "stack 2 tetracene molecules 3.4 angstrom apart", category: "Name-Lookup", root_cause: "Tetracene not in local DB", validation: { type: 'count', expected: { count: 2 } } },
    { id: 12, prompt: "create trimer of perylene in slipped-parallel at 3.5 angstrom", category: "Name-Hard", root_cause: "perylene + trimer + pattern", validation: { type: 'count', expected: { count: 3 } } },
    { id: 13, prompt: "generate 4 triphenylene molecules along [1,1,1] at 4 angstrom", category: "Name-Hard", root_cause: "triphenylene + vector", validation: { type: 'count', expected: { count: 4 } } },
    { id: 14, prompt: "stack tetramer of chrysene in herringbone pattern", category: "Name-Hard", root_cause: "chrysene + tetramer + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 15, prompt: "create 5 pentacene molecules in zigzag chain 3.5 angstrom", category: "Name-Hard", root_cause: "pentacene + count + pattern", validation: { type: 'count', expected: { count: 5 } } },
    { id: 16, prompt: "generate 3 corannulene molecules along [0,1,1] at 4 angstrom", category: "Name-Hard", root_cause: "corannulene + vector", validation: { type: 'count', expected: { count: 3 } } },
    { id: 17, prompt: "stack 4 TTF molecules antiparallel at 4 angstrom spacing", category: "Name-Hard", root_cause: "TTF + count + pattern", validation: { type: 'count', expected: { count: 4 } } },
    // KEPT FAILED: id 18 - TCNQ
    { id: 18, prompt: "create 2 TCNQ molecules 4 angstrom apart", category: "Name-Lookup", root_cause: "TCNQ abbreviation", validation: { type: 'count', expected: { count: 2 } } },
    { id: 19, prompt: "generate trimer of fullerene in triangular arrangement 15 angstrom", category: "Name-Hard", root_cause: "fullerene + trimer + geometry", validation: { type: 'count', expected: { count: 3 } } },
    { id: 20, prompt: "stack 4 buckyball molecules along [1,1,0] at 15 angstrom", category: "Name-Hard", root_cause: "buckyball + count + vector", validation: { type: 'count', expected: { count: 4 } } },

    // =========================================================================
    // GROUP 3: COUNT HANDLING - HARDER VERSIONS (21-40)
    // =========================================================================
    { id: 21, prompt: "generate exactly 6 benzene molecules in herringbone at 3.4 angstrom", category: "Count-Hard", root_cause: "exact count + pattern", validation: { type: 'count', expected: { count: 6 } } },
    { id: 22, prompt: "create a pentamer dimer hybrid of pyrene at 3.5 angstrom", category: "Count-Hard", root_cause: "pentamer = 5", validation: { type: 'count', expected: { count: 5 } } },
    { id: 23, prompt: "generate a hexamer of benzene in zigzag pattern 3.4 angstrom", category: "Count-Hard", root_cause: "hexamer = 6 + pattern", validation: { type: 'count', expected: { count: 6 } } },
    { id: 24, prompt: "create a pentamer of naphthalene antiparallel", category: "Count-Hard", root_cause: "pentamer = 5 + pattern", validation: { type: 'count', expected: { count: 5 } } },
    { id: 25, prompt: "stack 7 molecules of adamantane along [1,0,1] at 8 angstrom", category: "Count-Hard", root_cause: "7 + vector", validation: { type: 'count', expected: { count: 7 } } },
    { id: 26, prompt: "generate eight benzene molecules in double stack 3.4 angstrom", category: "Count-Hard", root_cause: "'eight' = 8", validation: { type: 'count', expected: { count: 8 } } },
    { id: 27, prompt: "create fifteen water molecules in 3x5 grid along z", category: "Count-Hard", root_cause: "fifteen = 15 + grid", validation: { type: 'count', expected: { count: 15 } } },
    { id: 28, prompt: "stack twenty methane molecules in helix at 4 angstrom pitch", category: "Count-Hard", root_cause: "'twenty' = 20 + helix", validation: { type: 'count', expected: { count: 20 } } },
    { id: 29, prompt: "generate half-dozen pyrene molecules in ring arrangement", category: "Count-Hard", root_cause: "'half-dozen' = 6", validation: { type: 'count', expected: { count: 6 } } },
    { id: 30, prompt: "create quartet of benzene molecules in T-shaped pattern", category: "Count-Hard", root_cause: "'quartet' = 4", validation: { type: 'count', expected: { count: 4 } } },
    // KEPT FAILED: id 31 - cubane
    { id: 31, prompt: "stack 2 cubane at 7 angstrom", category: "Count-NoMolWord", root_cause: "'2 cubane' without 'molecules'", validation: { type: 'count', expected: { count: 2 } } },
    { id: 32, prompt: "generate 6 corannulene in alternating orientations 4 angstrom", category: "Count-Hard", root_cause: "6 + pattern", validation: { type: 'count', expected: { count: 6 } } },
    { id: 33, prompt: "create 8 of dodecahedrane in cubic arrangement 12 angstrom", category: "Count-Hard", root_cause: "8 + cubic", validation: { type: 'count', expected: { count: 8 } } },
    { id: 34, prompt: "stack benzene, naphthalene, anthracene, pyrene in a line", category: "Count-Hard", root_cause: "4 different molecules", validation: { type: 'count', expected: { count: 4 } } },
    { id: 35, prompt: "generate three identical pairs of benzene at 3.4 angstrom", category: "Count-Hard", root_cause: "3 pairs = 6", validation: { type: 'count', expected: { count: 6 } } },
    { id: 36, prompt: "create double-decker benzene sandwich 3 layers thick", category: "Count-Hard", root_cause: "double-decker = 3", validation: { type: 'count', expected: { count: 3 } } },
    { id: 37, prompt: "stack benzene times 4 in slipped parallel 3.4 angstrom", category: "Count-Hard", root_cause: "'times 4' = 4 + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 38, prompt: "generate 3x pyrene stack with alternating rotation", category: "Count-Hard", root_cause: "'3x' = 3 + rotation", validation: { type: 'count', expected: { count: 3 } } },
    { id: 39, prompt: "create benzene x5 in zigzag chain 3.4 angstrom", category: "Count-Hard", root_cause: "'x5' = 5 + pattern", validation: { type: 'count', expected: { count: 5 } } },
    { id: 40, prompt: "stack benzene, benzene, benzene, benzene along [1,1,0]", category: "Count-Hard", root_cause: "4 repeated = 4", validation: { type: 'count', expected: { count: 4 } } },

    // =========================================================================
    // GROUP 4: INVALID PARAMETER HANDLING (41-55)
    // KEPT FAILED: 41,42,46,47,48,49,52,54,55
    // =========================================================================
    { id: 41, prompt: "generate 2 benzene 0 angstrom apart", category: "Invalid-Zero", root_cause: "0Å should default", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 42, prompt: "stack 2 pyrene at -5 angstrom", category: "Invalid-Negative", root_cause: "Negative distance", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 43, prompt: "create 3 benzene at 0.1 angstrom along [1,1,1] direction", category: "Invalid-Hard", root_cause: "0.1Å + vector + count", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 44, prompt: "generate 4 adamantane at 1 angstrom in herringbone pattern", category: "Invalid-Hard", root_cause: "1Å + 4 + pattern", validation: { type: 'no_clash', expected: { min: 5.0 } } },
    { id: 45, prompt: "stack 5 coronene at 0.5 angstrom in ring arrangement", category: "Invalid-Hard", root_cause: "0.5Å + 5 + ring", validation: { type: 'no_clash', expected: { min: 3.0 } } },
    { id: 46, prompt: "create 2 C60 1 angstrom apart", category: "Invalid-TooSmall", root_cause: "C60 7Å, 1Å = clash", validation: { type: 'no_clash', expected: { min: 7.0 } } },
    { id: 47, prompt: "generate 1000 benzene molecules stacked", category: "Invalid-TooMany", root_cause: "1000 molecules", validation: { type: 'count', expected: { min: 2 } } },
    { id: 48, prompt: "stack 0 benzene molecules", category: "Invalid-Zero", root_cause: "0 count", validation: { type: 'count', expected: { min: 1 } } },
    { id: 49, prompt: "generate -2 benzene molecules", category: "Invalid-Negative", root_cause: "Negative count", validation: { type: 'count', expected: { min: 1 } } },
    { id: 50, prompt: "create 3 benzene at 10000 angstrom along [0,1,1] in T-shape", category: "Invalid-Hard", root_cause: "10000Å + vector + pattern", validation: { type: 'distance', expected: { min: 9000, max: 11000 } } },
    { id: 51, prompt: "stack 4 benzene along imaginary axis at 5 angstrom", category: "Invalid-Hard", root_cause: "'imaginary axis' invalid", validation: { type: 'count', expected: { count: 4 } } },
    { id: 52, prompt: "generate 2 benzene along axis 4", category: "Invalid-Axis", root_cause: "'axis 4' invalid", validation: { type: 'count', expected: { count: 2 } } },
    { id: 53, prompt: "create 4 benzene along [0,0,0] in zigzag at 5 angstrom", category: "Invalid-Hard", root_cause: "[0,0,0] + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 54, prompt: "stack benzene at NaN angstrom", category: "Invalid-NaN", root_cause: "NaN distance", validation: { type: 'count', expected: { min: 1 } } },
    { id: 55, prompt: "generate 2 benzene at infinity angstrom", category: "Invalid-Inf", root_cause: "'infinity' distance", validation: { type: 'count', expected: { min: 1 } } },

    // =========================================================================
    // GROUP 5: PATTERN IMPLEMENTATION (56-70)
    // KEPT FAILED: 56,58,62,66,67,68,70
    // =========================================================================
    { id: 56, prompt: "create helical stack of 5 benzene with 3.4 pitch", category: "Pattern-Helix", root_cause: "Helix", validation: { type: 'count', expected: { count: 5 } } },
    { id: 57, prompt: "generate 6 naphthalene in zigzag chain along [1,0,1] at 3.5 angstrom", category: "Pattern-Hard", root_cause: "zigzag + vector", validation: { type: 'count', expected: { count: 6 } } },
    { id: 58, prompt: "create tetrahedral arrangement of 4 benzene", category: "Pattern-Tetrahedral", root_cause: "3D polyhedra", validation: { type: 'count', expected: { count: 4 } } },
    { id: 59, prompt: "generate dodecahedral arrangement of 12 water molecules 3 angstrom", category: "Pattern-Hard", root_cause: "dodecahedra = 12 vertices", validation: { type: 'count', expected: { count: 12 } } },
    { id: 60, prompt: "create T-shaped trimer of pyrene along [1,1,0]", category: "Pattern-Hard", root_cause: "T-shaped + trimer + vector", validation: { type: 'count', expected: { count: 3 } } },
    { id: 61, prompt: "generate herringbone pattern of 8 naphthalene at 3.5 angstrom along z", category: "Pattern-Hard", root_cause: "herringbone + 8 + axis", validation: { type: 'count', expected: { count: 8 } } },
    { id: 62, prompt: "create sandwich structure with benzene around ferrocene", category: "Pattern-Sandwich", root_cause: "3-layer A-B-A", validation: { type: 'count', expected: { count: 3 } } },
    { id: 63, prompt: "generate 4 slipped-parallel anthracene along [1,0,1] at 3.6 angstrom", category: "Pattern-Hard", root_cause: "slipped + 4 + vector", validation: { type: 'count', expected: { count: 4 } } },
    { id: 64, prompt: "create 6 face-to-edge pyrene molecules in alternating pattern", category: "Pattern-Hard", root_cause: "face-edge + 6 + alternating", validation: { type: 'count', expected: { count: 6 } } },
    { id: 65, prompt: "generate antiparallel hexamer of benzene along [0,1,1]", category: "Pattern-Hard", root_cause: "antiparallel + 6 + vector", validation: { type: 'count', expected: { count: 6 } } },
    { id: 66, prompt: "create circular ring of 6 benzene", category: "Pattern-Ring", root_cause: "Ring", validation: { type: 'count', expected: { count: 6 } } },
    { id: 67, prompt: "generate 3x3 grid of benzene", category: "Pattern-Grid", root_cause: "2D grid", validation: { type: 'count', expected: { count: 9 } } },
    { id: 68, prompt: "create cubic arrangement of 8 benzene", category: "Pattern-Cubic", root_cause: "Cube vertices", validation: { type: 'count', expected: { count: 8 } } },
    { id: 69, prompt: "generate alternating tetramer of naphthalene at 3.5 angstrom along [1,1,1]", category: "Pattern-Hard", root_cause: "alternating + 4 + vector", validation: { type: 'count', expected: { count: 4 } } },
    { id: 70, prompt: "create star arrangement of 5 benzene around 1 center", category: "Pattern-Star", root_cause: "Star", validation: { type: 'count', expected: { count: 6 } } },

    // =========================================================================
    // GROUP 6: DIRECTION VECTOR HANDLING (71-85)
    // KEPT FAILED: 80,81,82
    // =========================================================================
    { id: 71, prompt: "stack 4 pyrene along direction [1,0,0] in herringbone at 3.5 angstrom", category: "Direction-Hard", root_cause: "[1,0,0] + 4 + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 72, prompt: "generate 5 naphthalene along direction [0,1,0] in zigzag at 3.5 angstrom", category: "Direction-Hard", root_cause: "[0,1,0] + 5 + zigzag", validation: { type: 'count', expected: { count: 5 } } },
    { id: 73, prompt: "create 6 anthracene along direction [0,0,1] antiparallel at 3.6 angstrom", category: "Direction-Hard", root_cause: "[0,0,1] + 6 + antiparallel", validation: { type: 'count', expected: { count: 6 } } },
    { id: 74, prompt: "stack 4 coronene along [1,1,0] direction in slipped-parallel 4 angstrom", category: "Direction-Hard", root_cause: "[1,1,0] + 4 + slipped", validation: { type: 'count', expected: { count: 4 } } },
    { id: 75, prompt: "generate 5 perylene along [1,1,1] direction in T-shape 4 angstrom", category: "Direction-Hard", root_cause: "[1,1,1] + 5 + T-shape", validation: { type: 'count', expected: { count: 5 } } },
    { id: 76, prompt: "create 3 benzene along [1,0,1] with alternating rotation 5 angstrom", category: "Direction-Hard", root_cause: "[1,0,1] + 3 + rotation", validation: { type: 'count', expected: { count: 3 } } },
    { id: 77, prompt: "stack 4 pyrene along [0,1,1] in face-to-edge pattern 4 angstrom", category: "Direction-Hard", root_cause: "[0,1,1] + 4 + face-edge", validation: { type: 'count', expected: { count: 4 } } },
    { id: 78, prompt: "generate 5 molecules at 45 degrees with herringbone pattern 3.5 angstrom", category: "Direction-Hard", root_cause: "45° + 5 + herringbone", validation: { type: 'count', expected: { count: 5 } } },
    { id: 79, prompt: "create 4 benzene at 30 degrees from z in zigzag at 3.5 angstrom", category: "Direction-Hard", root_cause: "30° + 4 + zigzag", validation: { type: 'count', expected: { count: 4 } } },
    { id: 80, prompt: "stack 2 benzene perpendicular to x axis", category: "Direction-Perp", root_cause: "Perpendicular = yz", validation: { type: 'count', expected: { count: 2 } } },
    { id: 81, prompt: "generate 2 benzene perpendicular to y axis", category: "Direction-Perp", root_cause: "Perpendicular = xz", validation: { type: 'count', expected: { count: 2 } } },
    { id: 82, prompt: "create 2 benzene perpendicular to z axis", category: "Direction-Perp", root_cause: "Perpendicular = xy", validation: { type: 'count', expected: { count: 2 } } },
    { id: 83, prompt: "stack 5 naphthalene horizontally with alternating tilt 3.5 angstrom", category: "Direction-Hard", root_cause: "horizontal + 5 + tilt", validation: { type: 'count', expected: { count: 5 } } },
    { id: 84, prompt: "generate 4 pyrene vertically in slipped-parallel 4 angstrom", category: "Direction-Hard", root_cause: "vertical + 4 + slipped", validation: { type: 'count', expected: { count: 4 } } },
    { id: 85, prompt: "create 3 benzene sideways with T-shaped arrangement 5 angstrom", category: "Direction-Hard", root_cause: "sideways + 3 + T-shape", validation: { type: 'count', expected: { count: 3 } } },

    // =========================================================================
    // GROUP 7: SMILES/CID INPUT (86-95)
    // KEPT FAILED: 91,92
    // =========================================================================
    { id: 86, prompt: "generate 4 c1ccccc1 molecules in herringbone pattern 3.4 angstrom", category: "SMILES-Hard", root_cause: "benzene SMILES + 4 + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 87, prompt: "create 5 c1ccc2ccccc2c1 molecules in zigzag at 3.5 angstrom", category: "SMILES-Hard", root_cause: "naphthalene SMILES + 5 + zigzag", validation: { type: 'count', expected: { count: 5 } } },
    { id: 88, prompt: "stack 3 CC(=O)O molecules along [1,1,0] at 4 angstrom", category: "SMILES-Hard", root_cause: "acetic acid + 3 + vector", validation: { type: 'count', expected: { count: 3 } } },
    { id: 89, prompt: "generate 4 CCO molecules in alternating pattern 4 angstrom", category: "SMILES-Hard", root_cause: "ethanol + 4 + alternating", validation: { type: 'count', expected: { count: 4 } } },
    { id: 90, prompt: "create 6 c1ccncc1 molecules along [0,1,1] in T-shape 3.5 angstrom", category: "SMILES-Hard", root_cause: "pyridine + 6 + vector + T", validation: { type: 'count', expected: { count: 6 } } },
    { id: 91, prompt: "generate 2 CID:241 molecules 3.4 angstrom apart", category: "CID-Benzene", root_cause: "CID 241 = benzene", validation: { type: 'count', expected: { count: 2 } } },
    { id: 92, prompt: "stack 2 PubChem:241 molecules", category: "CID-PubChem", root_cause: "PubChem:CID", validation: { type: 'count', expected: { count: 2 } } },
    { id: 93, prompt: "create 4 [C]1[C][C][C][C][C]1 along [1,0,1] in herringbone 4 angstrom", category: "SMILES-Hard", root_cause: "explicit C + 4 + vector + pattern", validation: { type: 'count', expected: { count: 4 } } },
    { id: 94, prompt: "generate 3 InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H molecules in zigzag 3.5 angstrom", category: "SMILES-Hard", root_cause: "InChI + 3 + zigzag", validation: { type: 'count', expected: { count: 3 } } },
    { id: 95, prompt: "stack 5 O=C=O molecules along [1,1,1] in slipped-parallel 4 angstrom", category: "SMILES-Hard", root_cause: "CO2 + 5 + vector + slipped", validation: { type: 'count', expected: { count: 5 } } },

    // =========================================================================
    // GROUP 8: LLM TOOL SELECTION (96-100)
    // KEPT FAILED: 96
    // =========================================================================
    { id: 96, prompt: "generate 2 H2 molecules 3 angstrom apart", category: "LLM-H2", root_cause: "H2 dimer", validation: { type: 'count', expected: { count: 2 } } },
    { id: 97, prompt: "create 5 N2 molecules in alternating pattern 4 angstrom", category: "LLM-Hard", root_cause: "N2 + 5 + alternating", validation: { type: 'count', expected: { count: 5 } } },
    { id: 98, prompt: "stack 4 O2 molecules along [1,0,1] at 4 angstrom", category: "LLM-Hard", root_cause: "O2 + 4 + vector", validation: { type: 'count', expected: { count: 4 } } },
    { id: 99, prompt: "generate 3 CO molecules in T-shaped arrangement 4 angstrom", category: "LLM-Hard", root_cause: "CO + 3 + T-shape", validation: { type: 'count', expected: { count: 3 } } },
    { id: 100, prompt: "create 6 Cl2 molecules in zigzag chain 4 angstrom", category: "LLM-Hard", root_cause: "Cl2 + 6 + zigzag", validation: { type: 'count', expected: { count: 6 } } },
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

test('Failure Investigation v2 - 100 HARDER Tests', async ({ page }) => {
    test.setTimeout(30 * 60 * 1000);  // 30 minutes for 100 tests
    let results: any[] = [];

    await page.goto(FRONTEND_URL);
    await page.waitForTimeout(3000);

    // Dismiss welcome modal if present
    const getStartedButton = page.locator('button:has-text("Get Started")');
    if (await getStartedButton.isVisible({ timeout: 2000 }).catch(() => false)) {
        await getStartedButton.click();
        await page.waitForTimeout(500);
    }

    console.log('\n' + '='.repeat(80));
    console.log('FAILURE INVESTIGATION v2 - 27 KEPT FAILURES + 73 HARDER TESTS');
    console.log('='.repeat(80) + '\n');

    for (const testCase of FAILURE_INVESTIGATION_TESTS) {
        const startTime = Date.now();
        let llmTool: string | null = null;
        let llmArgs: any = null;
        let mcpStartTime = 0;
        let mcpTime = 0;

        // Clear previous input
        const textarea = page.locator('textarea');
        await textarea.fill('');
        await textarea.fill(testCase.prompt);

        // Capture LLM response
        const llmPromise = page.waitForResponse(
            (r) => r.url().includes('/api/chat') || r.url().includes(':11434'),
            { timeout: 120000 }
        ).catch(() => null);

        // Submit
        const sendButton = page.locator('button:has-text("Send"), button[type="submit"]').first();
        await sendButton.click();

        const llmResponse = await llmPromise;
        const llmTime = Date.now() - startTime;
        mcpStartTime = Date.now();

        if (llmResponse) {
            try {
                const body = await llmResponse.json();
                if (body.message?.tool_calls?.[0]) {
                    llmTool = body.message.tool_calls[0].function?.name;
                    llmArgs = body.message.tool_calls[0].function?.arguments;
                    if (typeof llmArgs === 'string') llmArgs = JSON.parse(llmArgs);
                }
            } catch { }
        }

        // Wait for MCP response / Tool Execution
        // Increased from 2000ms to 10000ms to allow for Python backend time
        await page.waitForTimeout(10000);
        mcpTime = Date.now() - mcpStartTime;

        // Initialize LLM Client
        await page.evaluate(() => {
            // (window as any).llmClient = new LLMClient(); // Singleton already exists
        });

        // Filter for specific test if needed
        const hardTests = [
            {
                id: "DEBUG-001",
                category: "Debug",
                prompt: "Create a benzene dimer separated by 3.5 Angstroms",
                root_cause: "Debug"
            }
        ];

        // Get structure data
        const structureData = await page.evaluate(() => {
            const state = (window as any).__REDUX_STATE__;
            if (!state) return { error: "Redux state not found on window" };
            if (!state.structure) return { error: "structure slice not found" };

            console.log('Structure state:', {
                structuresCount: state.structure.structures?.length,
                activeId: state.structure.activeStructureId
            });

            if (state.structure.structures?.length > 0) {
                const latest = state.structure.structures[state.structure.structures.length - 1];
                return { atoms: latest.data?.atoms || [], metadata: latest.data?.metadata || {} };
            }
            return null;
        });

        if (structureData && (structureData as any).error) {
            console.log("DEBUG STATE ERROR:", (structureData as any).error);
        }

        // Validate
        let status = 'MCP_FAIL';
        let details = 'No atoms';

        if (structureData && structureData.atoms.length > 0) {
            const validation = validateTest(testCase, structureData.atoms, structureData.metadata);
            status = validation.valid ? 'PASSED' : 'SEMANTIC_FAIL';
            details = validation.details;
        }

        const result = {
            id: testCase.id,
            prompt: testCase.prompt,
            category: testCase.category,
            root_cause: testCase.root_cause,
            status,
            test_details: details,
            llm_tool: llmTool,
            llm_args: llmArgs,
            llm_time_ms: llmTime,
            mcp_time_ms: mcpTime,
            atom_count: structureData?.atoms?.length || 0,
        };

        results.push(result);

        const icon = status === 'PASSED' ? '' : (status === 'SEMANTIC_FAIL' ? '⚠' : '✗');
        console.log(`\n[${testCase.id}] ${testCase.category}: ${testCase.root_cause}`);
        console.log(`  ${status}: ${icon} ${details} [LLM:${llmTime}ms MCP:${mcpTime}ms]`);

        // Clear for next test
        await page.evaluate(() => {
            const state = (window as any).__REDUX_STATE__;
            if (state?.structure) state.structure.structures = [];
        });
    }

    // Write results
    fs.writeFileSync(REPORT_FILE, JSON.stringify(results, null, 2));

    // Summary
    const passed = results.filter(r => r.status === 'PASSED').length;
    const semantic = results.filter(r => r.status === 'SEMANTIC_FAIL').length;
    const mcp = results.filter(r => r.status === 'MCP_FAIL').length;

    console.log('\n' + '='.repeat(80));
    console.log('HARDER TEST SUITE COMPLETE');
    console.log('='.repeat(80));
    console.log(`\nPASSED: ${passed}/100 (${(passed).toFixed(0)}%)`);
    console.log(`SEMANTIC_FAIL: ${semantic}`);
    console.log(`MCP_FAIL: ${mcp}`);

    // Category breakdown
    const cats: Record<string, { p: number; t: number }> = {};
    for (const r of results) {
        const cat = r.category.split('-')[0];
        if (!cats[cat]) cats[cat] = { p: 0, t: 0 };
        cats[cat].t++;
        if (r.status === 'PASSED') cats[cat].p++;
    }
    console.log('\nBy Category:');
    for (const [cat, d] of Object.entries(cats).sort()) {
        console.log(`  ${cat}: ${d.p}/${d.t} (${(d.p / d.t * 100).toFixed(0)}%)`);
    }
});
