/**
 * Generate Slab Tool
 *
 * Tool for creating surface slabs for DFT calculations.
 */
import { GenerateSlabSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function generateSlab(input) {
    const parsed = GenerateSlabSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check miller_indices format: [h, k, l]", "Ensure thickness >= 1 and vacuum >= 0"], true));
    }
    // Rename structure -> structure_dict for Python compatibility
    const { structure, ...rest } = parsed.data;
    const params = { structure_dict: structure, ...rest, operation: "slab" };
    const result = await executePythonWithJSON("structure_tools.py", params, { timeout: 120000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error.code, pythonResult.error.message, pythonResult.error.details));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateSlab(args) {
    const result = await generateSlab(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `**Slab Generation Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const data = result.data;
    let output = `## Surface Slab Generated\n\n`;
    output += `**Miller Indices:** (${data.miller_indices.join(' ')})\n`;
    output += `**Slab Thickness:** ${data.slab_thickness.toFixed(3)} Angstroms\n`;
    output += `**Vacuum Thickness:** ${data.vacuum_thickness.toFixed(3)} Angstroms\n`;
    output += `**Surface Area:** ${data.surface_area.toFixed(3)} Angstrom^2\n`;
    if (data.fixed_atoms && data.fixed_atoms.length > 0) {
        output += `**Fixed Atoms:** ${data.fixed_atoms.length} atoms (bottom layers)\n`;
    }
    output += `\n${formatStructureOutput(data.slab, undefined)}`;
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.slab,
        slab_thickness: data.slab_thickness,
        vacuum_thickness: data.vacuum_thickness,
        surface_area: data.surface_area,
        miller_indices: data.miller_indices
    });
    return {
        content: [
            {
                type: "text",
                text: output
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=slab.js.map