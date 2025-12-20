import { CallToolResult } from "@modelcontextprotocol/sdk/types.js";
import { ComprehensiveGenerateInput } from "../../types/tools.js";
/**
 * Handle comprehensive_generate tool execution
 *
 * Routes requests to the unified Python router (comprehensive_structures.py),
 * which then delegates to the specific generator module.
 */
export declare function handleComprehensiveGenerate(args: ComprehensiveGenerateInput): Promise<CallToolResult>;
//# sourceMappingURL=comprehensive-generator.d.ts.map