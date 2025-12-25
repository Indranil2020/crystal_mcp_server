/**
 * Create Alloy Tool
 *
 * Tool for creating substitutional alloys.
 */
import { CreateAlloySchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function createAlloy(input) {
    const parsed = CreateAlloySchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check substitutions format"], true));
    }
    // Rename structure -> structure_dict for Python compatibility
    const { structure, ...rest } = parsed.data;
    const params = { structure_dict: structure, ...rest, operation: "alloy" };
    const result = await executePythonWithJSON("structure_tools.py", params, { timeout: 60000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error?.code || CrystalErrorCode.EXECUTION_ERROR, pythonResult.error?.message || "Alloy generation failed", pythonResult.error?.details));
    }
    return createSuccess(pythonResult);
}
export async function handleCreateAlloy(args) {
    const result = await createAlloy(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Alloy Creation Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.alloy_structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 🧪 Alloy Created\n\n${outputText}`
            }]
    };
}
//# sourceMappingURL=alloy.js.map