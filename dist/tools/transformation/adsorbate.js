/**
 * Add Adsorbate Tool
 *
 * Tool for adding adsorbates to surfaces.
 */
import { AddAdsorbateSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function addAdsorbate(input) {
    const parsed = AddAdsorbateSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, [], true));
    }
    // Rename structure -> structure_dict for Python compatibility
    const { structure, ...rest } = parsed.data;
    const params = { structure_dict: structure, ...rest, operation: "adsorbate" };
    const result = await executePythonWithJSON("structure_tools.py", params, { timeout: 60000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error?.code || CrystalErrorCode.EXECUTION_ERROR, pythonResult.error?.message || "Adsorbate addition failed", pythonResult.error?.details));
    }
    return createSuccess(pythonResult);
}
export async function handleAddAdsorbate(args) {
    const result = await addAdsorbate(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Adsorbate Addition Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 📎 Adsorbate Added\n\n${outputText}`
            }]
    };
}
//# sourceMappingURL=adsorbate.js.map