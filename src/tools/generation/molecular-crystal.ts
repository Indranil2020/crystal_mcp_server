/**
 * Generate Molecular Crystal Tool
 * 
 * Tool for generating molecular crystal structures using PyXtal.
 */

import { GenerateMolecularCrystalSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";

export async function generateMolecularCrystal(input: unknown): Promise<Result<any>> {
    const parsed = GenerateMolecularCrystalSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            ["Check molecules list and space group"],
            true
        ));
    }

    const params = { ...parsed.data, operation: "generate_molecular" };

    const result = await executePythonWithJSON<typeof params, any>(
        "crystal_generator.py",
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

export async function handleGenerateMolecularCrystal(args: unknown): Promise<any> {
    const result = await generateMolecularCrystal(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `**Molecular Crystal Generation Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const outputText = formatStructureOutput(data.structure, data.validation);

    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        validation: data.validation
    });

    return {
        content: [
            {
                type: "text",
                text: `## Molecular Crystal Generated\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
