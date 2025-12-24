/**
 * Create Defect Tool
 * 
 * Tool for creating point defects in crystal structures.
 */

import { CreateDefectSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";

export async function createDefect(input: unknown): Promise<Result<any>> {
    const parsed = CreateDefectSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            ["Check defect_type and defect_site", "Ensure defect_species is provided for substitution/interstitial"],
            true
        ));
    }

    const params = { ...parsed.data, operation: "defect" };

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
            pythonResult.error.code as CrystalErrorCode,
            pythonResult.error.message,
            pythonResult.error.details
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleCreateDefect(args: unknown): Promise<any> {
    const result = await createDefect(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `❌ **Defect Creation Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    let output = `## 🧩 Defect Created\n\n`;
    output += `**Type:** ${data.defect_type}\n`;
    output += `**Site Index:** ${data.defect_site}\n`;

    if (data.defect_species) {
        output += `**Species:** ${data.defect_species}\n`;
    }

    output += `\n${formatStructureOutput(data.defected_structure, undefined)}`;

    return {
        content: [{
            type: "text",
            text: output
        }]
    };
}
