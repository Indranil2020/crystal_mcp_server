import { VisualizationSchema } from "../../types/tools.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
export async function generateVisualization(input) {
    const parsed = VisualizationSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check structure object and format"], true));
    }
    const result = await executePythonWithJSON("visualization.py", parsed.data, { timeout: 30000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, pythonResult.error || "Visualization failed", pythonResult));
    }
    return createSuccess(pythonResult.result);
}
export async function handleVisualization(args) {
    const result = await generateVisualization(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `** Visualization Failed **\n\n${result.error.message} `
                }],
            isError: true
        };
    }
    const data = result.data;
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        data: data
    });
    return {
        content: [
            {
                type: "text",
                text: JSON.stringify(data, null, 2)
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=visualization.js.map