/**
 * Crystal MCP Server
 * 
 * Main server implementation with tool registration and request handling.
 */

import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import {
  CallToolRequestSchema,
  ListToolsRequestSchema,
  Tool
} from "@modelcontextprotocol/sdk/types.js";
import { zodToJsonSchema } from "zod-to-json-schema";
import { createRequire } from "module";
const require = createRequire(import.meta.url);
const pkg = require("../package.json") as { version: string };

// Import tool handlers
import { handleComprehensiveGenerate } from "./tools/generation/comprehensive-generator.js";
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
import { handleVisualization } from "./tools/export/visualization.js";
import { handleCreateDefect } from "./tools/transformation/defect.js";
import { handleGenerateMolecularCrystal } from "./tools/generation/molecular-crystal.js";
import { handleGenerateNanostructure } from "./tools/generation/nanostructure.js";
import { handleExploreSymmetryRelations } from "./tools/analysis/symmetry-relations.js";
import { handleBuildMolecule } from "./tools/generation/build-molecule.js";
import { handleBuildMolecularCluster } from "./tools/generation/build-molecular-cluster.js";
import { handleCreateAlloy } from "./tools/transformation/alloy.js";
import { handleCreateHeterostructure } from "./tools/transformation/heterostructure.js";
import { handleAddAdsorbate } from "./tools/transformation/adsorbate.js";
import { handleApplyStrain } from "./tools/transformation/strain.js";
import {
  handleGeneratePrototype,
  handleGenerateTwistedBilayer,
  handleGenerateHighEntropyAlloy,
  handleGenerate2DMaterial,
  handleGenerateMOF,
  handleGenerateCage
} from "./tools/generation/advanced-structures.js";

// Import tool definitions
import { TOOL_DEFINITIONS } from "./types/tools.js";


/**
 * Create and configure the MCP server.
 */
export function createServer(): Server {
  const serverVersion = pkg.version;
  const server = new Server(
    {
      name: "crystal-structure-generator",
      version: serverVersion
    },
    {
      capabilities: {
        tools: {}
      }
    }
  );

  // Register ListTools handler
  const isToolInputSchema = (value: unknown): value is Tool["inputSchema"] => {
    if (!value || typeof value !== "object") {
      return false;
    }
    if ("type" in value) {
      const typeValue = (value as { type?: unknown }).type;
      return typeValue === "object";
    }
    return "anyOf" in value || "oneOf" in value;
  };

  server.setRequestHandler(ListToolsRequestSchema, async () => {
    const tools = TOOL_DEFINITIONS.map(tool => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const inputSchema = zodToJsonSchema(tool.inputSchema as any) as unknown;
      if (!isToolInputSchema(inputSchema)) {
        throw new Error(`Invalid input schema for tool: ${tool.name}`);
      }
      return {
        name: tool.name,
        description: tool.description,
        inputSchema,
        annotations: tool.annotations
      };
    });

    return {
      tools: tools as Tool[]
    };
  });

  // Register CallTool handler
  server.setRequestHandler(CallToolRequestSchema, async (request) => {
    const { name, arguments: args } = request.params;

    // DEBUG: Log all tool calls
    console.error(`[MCP DEBUG] 📞 Request received: Tool='${name}'`);
    console.error(`[MCP DEBUG] 📦 Arguments: ${JSON.stringify(args, null, 2)}`);

    switch (name) {
      case "comprehensive_generate":
        return await handleComprehensiveGenerate(args);

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

      case "generate_visualization":
        return await handleVisualization(args);

      case "create_defect":
        return await handleCreateDefect(args);

      case "generate_molecular_crystal":
        return await handleGenerateMolecularCrystal(args);

      case "generate_nanostructure":
        return await handleGenerateNanostructure(args);

      case "explore_symmetry_relations":
        return await handleExploreSymmetryRelations(args);

      case "build_molecule":
        return await handleBuildMolecule(args);

      case "build_molecular_cluster":
        return await handleBuildMolecularCluster(args);

      case "create_alloy":
        return await handleCreateAlloy(args);

      case "create_heterostructure":
        return await handleCreateHeterostructure(args);

      case "add_adsorbate":
        return await handleAddAdsorbate(args);

      case "apply_strain":
        return await handleApplyStrain(args);

      // New advanced structure tools
      case "generate_prototype":
        return await handleGeneratePrototype(args);

      case "generate_twisted_bilayer":
        return await handleGenerateTwistedBilayer(args);

      case "generate_high_entropy_alloy":
        return await handleGenerateHighEntropyAlloy(args);

      case "generate_2d_material":
        return await handleGenerate2DMaterial(args);

      case "generate_mof":
        return await handleGenerateMOF(args);

      case "generate_cage":
        return await handleGenerateCage(args);

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
export async function startServer(): Promise<void> {
  const server = createServer();
  const transport = new StdioServerTransport();

  await server.connect(transport);
  // Server is now running - no stderr output to avoid interfering with MCP protocol
}
