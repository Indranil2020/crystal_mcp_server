/**
 * Create Heterostructure Tool
 * 
 * Tool for stacking structures.
 */

import { CreateHeterostructureSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";

export async function createHeterostructure(input: unknown): Promise<Result<any>> {
    const parsed = CreateHeterostructureSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            ["Ensure both substrate and overlayer are valid structures"],
            true
        ));
    }

    const params = { ...parsed.data, operation: "heterostructure" };

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
            pythonResult.error?.message || "Heterostructure generation failed",
            pythonResult.error?.details
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleCreateHeterostructure(args: unknown): Promise<any> {
    const result = await createHeterostructure(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `❌ **Heterostructure Creation Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const warnings = data.warnings || [];
    const warningText = warnings.length > 0 ? `> ⚠️ **Warnings:**\n> ${warnings.join('\n> ')}\n\n` : "";

    const outputText = formatStructureOutput(data.heterostructure, undefined);

    return {
        content: [{
            type: "text",
            text: `## 🥞 Heterostructure Created\n\n${warningText}${outputText}`
        }]
    };
}
