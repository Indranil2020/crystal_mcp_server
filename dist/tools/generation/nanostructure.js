/**
 * Nanostructure Generation Tool
 *
 * Tool for creating nanotubes, graphene sheets, molecular crystals (C60), etc.
 */
// import { z } from "zod";
import { GenerateNanostructureSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
export async function generateNanostructure(input) {
    const parsed = GenerateNanostructureSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid nanostructure parameters", { zodErrors: parsed.error.errors }));
    }
    // Map 'ribbon_type' to 'type' for Python backend if generic type is nanoribbon
    const pythonParams = { ...parsed.data.params };
    if (parsed.data.type === "nanoribbon" && parsed.data.params.ribbon_type) {
        pythonParams.type = parsed.data.params.ribbon_type;
    }
    const result = await executePythonWithJSON("nanostructure_generator.py", {
        type: parsed.data.type,
        params: pythonParams
    });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateNanostructure(args) {
    const result = await generateNanostructure(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Nanostructure Generation Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    return {
        content: [{
                type: "text",
                text: `## 🧪 Nanostructure Generated\n\n${formatStructureOutput(result.data.structure)}`
            }]
    };
}
//# sourceMappingURL=nanostructure.js.map