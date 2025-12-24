/**
 * Apply Strain Tool
 * 
 * Tool for apply elastic strain.
 */

import { ApplyStrainSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";

export async function applyStrain(input: unknown): Promise<Result<any>> {
    const parsed = ApplyStrainSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            [],
            true
        ));
    }

    const params = { ...parsed.data, operation: "strain" };

    const result = await executePythonWithJSON<typeof params, any>(
        "structure_tools.py",
        params,
        { timeout: 60000 }
    );

    if (!result.success) {
        return createFailure(result.error);
    }

    const pythonResult = result.data;

    if (!pythonResult.success) {
        return createFailure(createError(
            pythonResult.error?.code as CrystalErrorCode || CrystalErrorCode.EXECUTION_ERROR,
            pythonResult.error?.message || "Strain application failed",
            pythonResult.error?.details
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleApplyStrain(args: unknown): Promise<any> {
    const result = await applyStrain(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `❌ **Strain Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const outputText = formatStructureOutput(data.strained_structure, undefined);

    return {
        content: [{
            type: "text",
            text: `## 📏 Strain Applied\n\n**Strain Tensor:** ${JSON.stringify(data.strain_tensor)}\n\n${outputText}`
        }]
    };
}
