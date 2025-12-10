/**
 * Generate Crystal Tool
 *
 * Main tool for generating crystal structures with PyXtal.
 * Follows defensive programming with Result<T> pattern.
 */
import { GenerateCrystalSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
/**
 * Generate a crystal structure with specified composition and space group.
 *
 * @param input - Generation parameters
 * @returns Result containing structure data or error
 */
export async function generateCrystal(input) {
    // Validate input schema
    const parsed = GenerateCrystalSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", {
            zodErrors: parsed.error.errors,
            receivedInput: input
        }, ["Check input parameters match the required schema", "See documentation for valid parameter types"], true));
    }
    const params = parsed.data;
    // Execute Python crystal generator
    const result = await executePythonWithJSON("crystal_generator.py", params, {
        timeout: 60000 // 60 second timeout
    });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data.data;
    // Check if Python execution succeeded
    if (!pythonResult.success) {
        const error = pythonResult.error;
        return createFailure(createError(error.code || CrystalErrorCode.GENERATION_FAILED, error.message, error.details, error.suggestions || [], error.recoverable !== false));
    }
    // Extract structure data
    const structure = pythonResult.structure;
    const validation = pythonResult.validation;
    const metadata = pythonResult.metadata;
    // Create result
    const generationResult = {
        structure,
        validation: {
            valid: validation.valid,
            issues: validation.errors || [],
            warnings: validation.warnings || []
        },
        generation_metadata: {
            attempts: metadata?.attempts || 1,
            seed: params.seed,
            generation_time_ms: metadata?.generation_time_ms || 0
        },
        files: {
            cif: generateCIFContent(structure),
            poscar: generatePOSCARContent(structure),
            xyz: generateXYZContent(structure),
            json: JSON.stringify(structure, null, 2)
        }
    };
    return createSuccess(generationResult);
}
/**
 * Generate CIF file content from structure.
 */
function generateCIFContent(structure) {
    const { lattice, space_group, atoms } = structure;
    let cif = `data_crystal\n`;
    cif += `_symmetry_space_group_name_H-M '${space_group.symbol}'\n`;
    cif += `_space_group_IT_number ${space_group.number}\n`;
    cif += `_cell_length_a ${lattice.a.toFixed(6)}\n`;
    cif += `_cell_length_b ${lattice.b.toFixed(6)}\n`;
    cif += `_cell_length_c ${lattice.c.toFixed(6)}\n`;
    cif += `_cell_angle_alpha ${lattice.alpha.toFixed(4)}\n`;
    cif += `_cell_angle_beta ${lattice.beta.toFixed(4)}\n`;
    cif += `_cell_angle_gamma ${lattice.gamma.toFixed(4)}\n`;
    cif += `_cell_volume ${lattice.volume.toFixed(4)}\n`;
    cif += `\nloop_\n`;
    cif += `_atom_site_label\n`;
    cif += `_atom_site_type_symbol\n`;
    cif += `_atom_site_fract_x\n`;
    cif += `_atom_site_fract_y\n`;
    cif += `_atom_site_fract_z\n`;
    atoms.forEach((atom, i) => {
        cif += `${atom.element}${i + 1} ${atom.element} ${atom.coords[0].toFixed(6)} ${atom.coords[1].toFixed(6)} ${atom.coords[2].toFixed(6)}\n`;
    });
    return cif;
}
/**
 * Generate POSCAR file content from structure.
 */
function generatePOSCARContent(structure) {
    const { lattice, atoms, metadata } = structure;
    let poscar = `${metadata.formula}\n`;
    poscar += `1.0\n`;
    // Lattice vectors
    lattice.matrix.forEach(vec => {
        poscar += `  ${vec[0].toFixed(12)}  ${vec[1].toFixed(12)}  ${vec[2].toFixed(12)}\n`;
    });
    // Count atoms by element
    const elementCounts = {};
    const elementOrder = [];
    atoms.forEach(atom => {
        const element = atom.element;
        if (element) {
            if (!elementCounts[element]) {
                elementCounts[element] = 0;
                elementOrder.push(element);
            }
            elementCounts[element]++;
        }
    });
    // Element names
    poscar += elementOrder.join(' ') + '\n';
    // Element counts
    poscar += elementOrder.map(el => elementCounts[el]).join(' ') + '\n';
    poscar += `Direct\n`;
    // Coordinates grouped by element
    elementOrder.forEach(element => {
        atoms.filter(a => a.element === element).forEach(atom => {
            poscar += `  ${atom.coords[0].toFixed(12)}  ${atom.coords[1].toFixed(12)}  ${atom.coords[2].toFixed(12)}\n`;
        });
    });
    return poscar;
}
/**
 * Generate XYZ file content from structure.
 */
function generateXYZContent(structure) {
    const { atoms, metadata } = structure;
    let xyz = `${atoms.length}\n`;
    xyz += `${metadata.formula}\n`;
    atoms.forEach(atom => {
        xyz += `${atom.element}  ${atom.cartesian[0].toFixed(8)}  ${atom.cartesian[1].toFixed(8)}  ${atom.cartesian[2].toFixed(8)}\n`;
    });
    return xyz;
}
/**
 * MCP tool handler for generate_crystal.
 */
export async function handleGenerateCrystal(args) {
    const result = await generateCrystal(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Crystal Generation Failed**\n\n**Error:** ${result.error.message}\n\n**Code:** ${result.error.code}\n\n**Suggestions:**\n${result.error.suggestions.map(s => `- ${s}`).join('\n')}`
                }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.structure, data.validation);
    return {
        content: [{
                type: "text",
                text: outputText
            }]
    };
}
//# sourceMappingURL=generate-crystal.js.map