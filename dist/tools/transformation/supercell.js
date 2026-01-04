/**
 * Make Supercell Tool
 *
 * Tool for creating supercells from crystal structures.
 */
import { MakeSupercellSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function makeSupercell(input) {
    const parsed = MakeSupercellSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check scaling_matrix format: [[nx,0,0],[0,ny,0],[0,0,nz]] or [nx,ny,nz]"], true));
    }
    const params = { ...parsed.data, operation: "supercell" };
    const result = await executePythonWithJSON("structure_tools.py", params, { timeout: 60000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error.code, pythonResult.error.message, pythonResult.error.details));
    }
    return createSuccess(pythonResult);
}
export async function handleMakeSupercell(args) {
    const result = await makeSupercell(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `**Supercell Generation Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.supercell, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.supercell,
        volume_multiplier: data.volume_multiplier,
        scaling_matrix: data.scaling_matrix
    });
    return {
        content: [
            {
                type: "text",
                text: `## 🔲 Supercell Generated\n\n**Volume Multiplier:** ${data.volume_multiplier}x\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=supercell.js.map