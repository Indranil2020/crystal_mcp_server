/**
 * Create Heterostructure Tool
 *
 * Tool for stacking structures.
 */
import { CreateHeterostructureSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function createHeterostructure(input) {
    const parsed = CreateHeterostructureSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Ensure both substrate and overlayer are valid structures"], true));
    }
    // Rename substrate/overlayer to _dict variants for Python compatibility
    const { substrate, overlayer, ...rest } = parsed.data;
    const params = { substrate_dict: substrate, overlayer_dict: overlayer, ...rest, operation: "heterostructure" };
    const result = await executePythonWithJSON("structure_tools.py", params, { timeout: 60000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error?.code || CrystalErrorCode.EXECUTION_ERROR, pythonResult.error?.message || "Heterostructure generation failed", pythonResult.error?.details));
    }
    return createSuccess(pythonResult);
}
export async function handleCreateHeterostructure(args) {
    const result = await createHeterostructure(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `**Heterostructure Creation Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const data = result.data;
    const warnings = data.warnings || [];
    const warningText = warnings.length > 0 ? `> ⚠️ **Warnings:**\n> ${warnings.join('\n> ')}\n\n` : "";
    const outputText = formatStructureOutput(data.heterostructure, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.heterostructure,
        warnings: data.warnings
    });
    return {
        content: [
            {
                type: "text",
                text: `## 🥞 Heterostructure Created\n\n${warningText}${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=heterostructure.js.map