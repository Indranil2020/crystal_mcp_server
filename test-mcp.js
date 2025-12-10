// Quick test to see what the server returns for ListTools
import { createServer } from './dist/server.js';

const server = createServer();

// Mock the request
const mockRequest = {
  params: {},
  method: 'tools/list'
};

// Try to get the handler
console.log('Server info:', {
  hasListToolsHandler: typeof server._requestHandlers !== 'undefined'
});

// Try calling listTools directly
import { TOOL_DEFINITIONS } from './dist/types/tools.js';
console.log('\nTool count:', TOOL_DEFINITIONS.length);
console.log('\nTool names:');
TOOL_DEFINITIONS.forEach(tool => {
  console.log(`  - ${tool.name}`);
});
