#!/usr/bin/env node
/**
 * Crystal MCP Server Entry Point
 *
 * Main entry point with environment verification.
 * CRITICAL: No output to stdout except JSON-RPC messages.
 */
import { startServer } from "./server.js";
// Suppress all uncaught errors from going to stdout/stderr
// This is critical for MCP protocol compliance
process.on('uncaughtException', () => {
    process.exit(1);
});
process.on('unhandledRejection', () => {
    process.exit(1);
});
/**
 * Main entry point.
 */
async function main() {
    // Start server immediately - Python validation happens on first tool call
    await startServer();
    // Keep the process alive - MCP server runs indefinitely
    // The transport will handle stdin/stdout and keep event loop active
}
// Run main function
main().catch(() => {
    process.exit(1);
});
//# sourceMappingURL=index.js.map