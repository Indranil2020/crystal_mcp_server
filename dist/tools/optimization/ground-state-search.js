/**
 * Ground State Search Tool
 *
 * Tool for finding ground state structures across multiple space groups.
 */
import { GroundStateSearchSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatGroundStateSearchOutput } from "../../utils/formatting.js";
export async function groundStateSearch(input) {
    const parsed = GroundStateSearchSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check composition and space_groups parameters", "Verify MLFF model is available"], true));
    }
    const params = parsed.data;
    const result = await executePythonWithJSON("ground_state_searcher.py", params, { timeout: 1800000 } // 30 minute timeout for extensive searches
    );
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error.code, pythonResult.error.message, pythonResult.error.details));
    }
    return createSuccess(pythonResult);
}
export async function handleGroundStateSearch(args) {
    const result = await groundStateSearch(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `**Ground State Search Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const outputText = formatGroundStateSearchOutput(result.data);
    // Include raw JSON data for the frontend viewer
    // Python returns: { best_structure: {...}, min_energy_per_atom: float, all_results: [...] }
    const jsonData = JSON.stringify({
        success: true,
        structure: result.data.best_structure,
        ground_state: result.data.best_structure,
        energy_ranking: result.data.all_results
    });
    return {
        content: [
            {
                type: "text",
                text: outputText
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=ground-state-search.js.map