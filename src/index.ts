#!/usr/bin/env node

/**
 * Crystal MCP Server Entry Point
 * 
 * Main entry point with environment verification.
 */

import { startServer } from "./server.js";

/**
 * Main entry point.
 */
async function main(): Promise<void> {
  // Start server immediately - Python validation happens on first tool call
  await startServer();
  
  // Keep the process alive - MCP server runs indefinitely
  // The transport will handle stdin/stdout and keep event loop active
}

// Run main function
main().catch(() => {
  process.exit(1);
});
