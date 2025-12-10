#!/usr/bin/env node
/**
 * Crystal MCP Server Entry Point
 *
 * Main entry point with environment verification.
 */
import { startServer } from "./server.js";
import { checkPythonAvailable } from "./utils/python-bridge.js";
/**
 * Verify Python environment before starting server.
 */
async function verifyEnvironment() {
    // Check Python availability
    const pythonCheck = await checkPythonAvailable();
    if (!pythonCheck.success) {
        console.error("❌ Python environment check failed:");
        console.error(pythonCheck.error.message);
        console.error("\nPlease ensure:");
        console.error("- Python 3.8+ is installed");
        console.error("- Required packages are installed (run: pip install -r requirements.txt)");
        return false;
    }
    console.error("✅ Python environment verified");
    return true;
}
/**
 * Main entry point.
 */
async function main() {
    console.error("🔬 Crystal Structure Generator MCP Server");
    console.error("========================================\n");
    // Verify environment
    const envOk = await verifyEnvironment();
    if (!envOk) {
        process.exit(1);
    }
    console.error("");
    // Start server
    await startServer();
}
// Run main function
main().catch((error) => {
    console.error("Fatal error:", error);
    process.exit(1);
});
//# sourceMappingURL=index.js.map