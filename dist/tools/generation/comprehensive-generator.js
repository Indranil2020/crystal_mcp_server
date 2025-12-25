import { ComprehensiveGenerateSchema } from "../../types/tools.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
const isRecord = (value) => typeof value === "object" && value !== null;
/**
 * Handle comprehensive_generate tool execution
 *
 * Routes requests to the unified Python router (comprehensive_structures.py),
 * which then delegates to the specific generator module.
 */
export async function handleComprehensiveGenerate(args) {
    // Validate input with Zod schema
    const parsed = ComprehensiveGenerateSchema.safeParse(args);
    if (!parsed.success) {
        return {
            content: [
                {
                    type: "text",
                    text: `Invalid input parameters:\n${parsed.error.errors.map(e => `- ${e.path.join('.')}: ${e.message}`).join('\n')}`
                }
            ],
            isError: true
        };
    }
    const validatedArgs = parsed.data;
    // Execute Python script
    const result = await executePythonWithJSON("comprehensive_structures.py", validatedArgs);
    if (!result.success) {
        return {
            content: [
                {
                    type: "text",
                    text: `Python execution failed: ${typeof result.error === 'object' ? JSON.stringify(result.error) : result.error || "Unknown error"}`
                }
            ],
            isError: true
        };
    }
    // Handle errors returned from Python
    const data = result.data;
    if (!isRecord(data)) {
        return {
            content: [
                {
                    type: "text",
                    text: "Python execution returned an unexpected response shape."
                }
            ],
            isError: true
        };
    }
    if (data.success === false) {
        // Return JSON with success: false for proper API response
        // This allows clients to parse and check the error details
        return {
            content: [
                {
                    type: "text",
                    text: JSON.stringify(data, null, 2)
                }
            ],
            isError: true
        };
    }
    // Structure generation successful
    const structure = data;
    const operation = validatedArgs.operation;
    // Format success message
    let summary = `Successfully executed operation '${operation}'`;
    // Add operation-specific details if available
    if (typeof structure.formula === "string") {
        summary += `\nGenerated structure: ${structure.formula}`;
    }
    if (typeof structure.n_atoms === "number") {
        summary += `\nNumber of atoms: ${structure.n_atoms}`;
    }
    if (typeof structure.spacegroup_symbol === "string") {
        summary += `\nSpace Group: ${structure.spacegroup_symbol}`;
    }
    // Handle listing operations
    if (operation === "list_all" || operation === "list_category" || validatedArgs.list_available) {
        return {
            content: [
                {
                    type: "text",
                    text: JSON.stringify(data, null, 2)
                }
            ]
        };
    }
    // Return formatted result
    return {
        content: [
            {
                type: "text",
                text: summary
            },
            {
                type: "text",
                text: JSON.stringify(structure, null, 2)
            }
        ]
    };
}
//# sourceMappingURL=comprehensive-generator.js.map