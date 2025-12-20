import { executePythonWithJSON } from "../../utils/python-bridge.js";
/**
 * Handle comprehensive_generate tool execution
 *
 * Routes requests to the unified Python router (comprehensive_structures.py),
 * which then delegates to the specific generator module.
 */
export async function handleComprehensiveGenerate(args) {
    // Execute Python script
    const result = await executePythonWithJSON("comprehensive_structures.py", args);
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
    const data = result.data;
    if (data.success === false) {
        // Format error message from Python generator
        const errorInfo = data.error || {};
        let errorMessage = `Generator failed: ${errorInfo.message || "Unknown error"}`;
        if (errorInfo.hint) {
            errorMessage += `\nHint: ${errorInfo.hint}`;
        }
        if (errorInfo.available_categories) {
            errorMessage += `\nAvailable categories: ${errorInfo.available_categories.join(", ")}`;
        }
        if (errorInfo.available_operations) {
            errorMessage += `\nAvailable operations: ${Object.keys(errorInfo.available_operations).join(", ")}`;
        }
        return {
            content: [
                {
                    type: "text",
                    text: errorMessage
                }
            ],
            isError: true
        };
    }
    // Structure generation successful
    const structure = data;
    const operation = args.operation;
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
    if (operation === "list_all" || operation === "list_category" || args.list_available) {
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