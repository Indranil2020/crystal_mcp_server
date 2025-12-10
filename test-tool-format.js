// Check if the tools are in the correct MCP format
import { TOOL_DEFINITIONS } from './dist/types/tools.js';

const firstTool = TOOL_DEFINITIONS[0];
console.log('First tool structure:');
console.log(JSON.stringify(firstTool, null, 2));

// Check if inputSchema is a Zod schema object
console.log('\ninputSchema type:', typeof firstTool.inputSchema);
console.log('Has _def:', firstTool.inputSchema && typeof firstTool.inputSchema._def !== 'undefined');

// MCP expects inputSchema to be a JSON Schema, not a Zod schema
// Let's see if we need to convert it
