import { CallToolResult } from "@modelcontextprotocol/sdk/types.js";
import { ComprehensiveGenerateSchema } from "../../types/tools.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

/**
 * Handle comprehensive_generate tool execution
 *
 * Routes requests to the unified Python router (comprehensive_structures.py),
 * which then delegates to the specific generator module.
 */
export async function handleComprehensiveGenerate(
    args: unknown
): Promise<CallToolResult> {
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
    const result = await executePythonWithJSON(
        "comprehensive_structures.py",
        validatedArgs
    );

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
    // Type assertion needed as result.data is generic unknown type
    const data = result.data as any;

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
    if (structure.formula) {
        summary += `\nGenerated structure: ${structure.formula}`;
    }
    if (structure.n_atoms) {
        summary += `\nNumber of atoms: ${structure.n_atoms}`;
    }
    if (structure.spacegroup_symbol) {
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
