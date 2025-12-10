/**
 * Analyze Symmetry Tool
 *
 * Tool for detecting and analyzing crystal symmetry.
 */
import { AnalyzeSymmetrySchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatSymmetryOutput } from "../../utils/formatting.js";
export async function analyzeSymmetry(input) {
    const parsed = AnalyzeSymmetrySchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check symprec and angle_tolerance values"], true));
    }
    const params = parsed.data;
    const result = await executePythonWithJSON("symmetry_analyzer.py", params, { timeout: 60000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error.code, pythonResult.error.message, pythonResult.error.details));
    }
    return createSuccess(pythonResult);
}
export async function handleAnalyzeSymmetry(args) {
    const result = await analyzeSymmetry(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Symmetry Analysis Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const outputText = formatSymmetryOutput(result.data);
    return {
        content: [{
                type: "text",
                text: outputText
            }]
    };
}
//# sourceMappingURL=symmetry.js.map