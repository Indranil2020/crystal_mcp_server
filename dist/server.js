/**
 * Crystal MCP Server
 *
 * Main server implementation with tool registration and request handling.
 */
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import { CallToolRequestSchema, ListToolsRequestSchema } from "@modelcontextprotocol/sdk/types.js";
import { zodToJsonSchema } from "zod-to-json-schema";
// Import tool handlers
import { handleGenerateCrystal } from "./tools/generation/generate-crystal.js";
import { handleSpaceGroupScan } from "./tools/generation/space-group-scan.js";
import { handleMakeSupercell } from "./tools/transformation/supercell.js";
import { handleGenerateSlab } from "./tools/transformation/slab.js";
import { handleAnalyzeSymmetry } from "./tools/analysis/symmetry.js";
import { handleValidateStructure } from "./tools/analysis/validation.js";
import { handleOptimizeStructureMLFF } from "./tools/optimization/mlff-optimize.js";
import { handleCalculateEnergyMLFF } from "./tools/optimization/calculate-energy-mlff.js";
import { handleGroundStateSearch } from "./tools/optimization/ground-state-search.js";
import { handleExportStructure } from "./tools/export/export-structure.js";
// Import tool definitions
import { TOOL_DEFINITIONS } from "./types/tools.js";
/**
 * Create and configure the MCP server.
 */
export function createServer() {
    const server = new Server({
        name: "crystal-structure-generator",
        version: "1.0.0"
    }, {
        capabilities: {
            tools: {}
        }
    });
    // Register ListTools handler
    server.setRequestHandler(ListToolsRequestSchema, async () => {
        const tools = TOOL_DEFINITIONS.map(tool => ({
            name: tool.name,
            description: tool.description,
            inputSchema: zodToJsonSchema(tool.inputSchema)
        }));
        return {
            tools: tools
        };
    });
    // Register CallTool handler
    server.setRequestHandler(CallToolRequestSchema, async (request) => {
        const { name, arguments: args } = request.params;
        switch (name) {
            case "generate_crystal":
                return await handleGenerateCrystal(args);
            case "generate_space_group_scan":
                return await handleSpaceGroupScan(args);
            case "make_supercell":
                return await handleMakeSupercell(args);
            case "generate_slab":
                return await handleGenerateSlab(args);
            case "analyze_symmetry":
                return await handleAnalyzeSymmetry(args);
            case "validate_structure":
                return await handleValidateStructure(args);
            case "optimize_structure_mlff":
                return await handleOptimizeStructureMLFF(args);
            case "calculate_energy_mlff":
                return await handleCalculateEnergyMLFF(args);
            case "ground_state_search":
                return await handleGroundStateSearch(args);
            case "export_structure":
                return await handleExportStructure(args);
            default:
                return {
                    content: [{
                            type: "text",
                            text: `❌ Unknown tool: ${name}\n\nAvailable tools:\n${TOOL_DEFINITIONS.map(t => `- ${t.name}`).join('\n')}`
                        }],
                    isError: true
                };
        }
    });
    return server;
}
/**
 * Start the MCP server with stdio transport.
 */
export async function startServer() {
    const server = createServer();
    const transport = new StdioServerTransport();
    await server.connect(transport);
    console.error("Crystal Structure Generator MCP Server running on stdio");
    console.error("Version: 1.0.0");
    console.error("Available tools:", TOOL_DEFINITIONS.length);
}
//# sourceMappingURL=server.js.map