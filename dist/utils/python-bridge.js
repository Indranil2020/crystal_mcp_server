/**
 * Python Bridge Utility
 *
 * This module provides a safe interface for executing Python scripts
 * from TypeScript, with comprehensive error handling and validation.
 * Uses defensive programming without try/catch blocks.
 */
import { PythonShell } from "python-shell";
import * as fs from "fs";
import * as path from "path";
import * as os from "os";
import { fileURLToPath } from "url";
import { spawnSync } from "child_process";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../types/errors.js";
/**
 * Get the directory path in ES modules
 */
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
/**
 * Configuration for Python environment
 */
export class PythonConfig {
    static instance = null;
    pythonPath = process.env.PYTHON_PATH ??
        process.env.PYTHON ??
        (process.platform === "win32" ? "python" : "python3");
    scriptsDirectory;
    constructor() {
        // Default to src/python directory
        // Check possible locations
        const possiblePaths = [
            path.resolve(__dirname, "..", "python"), // Relative to src/utils
            path.resolve(__dirname, "..", "..", "src", "python"), // Relative to dist/utils
        ];
        this.scriptsDirectory = possiblePaths.find(p => fs.existsSync(p)) || possiblePaths[0];
    }
    static getInstance() {
        if (PythonConfig.instance === null) {
            PythonConfig.instance = new PythonConfig();
        }
        return PythonConfig.instance;
    }
    setPythonPath(pythonPath) {
        this.pythonPath = pythonPath;
    }
    getPythonPath() {
        return this.pythonPath;
    }
    setScriptsDirectory(dir) {
        this.scriptsDirectory = dir;
    }
    getScriptsDirectory() {
        return this.scriptsDirectory;
    }
    /**
     * Get full path to a Python script
     */
    getScriptPath(scriptName) {
        return path.join(this.scriptsDirectory, scriptName);
    }
}
/**
 * Verify Python script exists
 */
function verifyScriptExists(scriptPath) {
    if (!fs.existsSync(scriptPath)) {
        return createFailure(createError(CrystalErrorCode.SCRIPT_NOT_FOUND, `Python script not found: ${scriptPath}`, { scriptPath }, [
            "Check that the script file exists",
            "Verify the scripts directory is correctly configured",
            "Ensure the script name is spelled correctly"
        ]));
    }
    const stats = fs.statSync(scriptPath);
    if (!stats.isFile()) {
        return createFailure(createError(CrystalErrorCode.SCRIPT_NOT_FOUND, `Script path is not a file: ${scriptPath}`, { scriptPath }, ["Provide a valid file path"]));
    }
    return createSuccess(scriptPath);
}
/**
 * Validate JSON string
 */
function isValidJSON(str) {
    if (str.trim().length === 0) {
        return false;
    }
    const firstChar = str.trim()[0];
    if (firstChar !== "{" && firstChar !== "[") {
        return false;
    }
    return true;
}
/**
 * Parse JSON output from Python script
 */
function parseJSONOutput(output) {
    if (!isValidJSON(output)) {
        return createFailure(createError(CrystalErrorCode.INVALID_JSON_RESPONSE, "Python script did not return valid JSON", { output: output.substring(0, 200) }, [
            "Ensure Python script prints only JSON to stdout",
            "Check for print statements or warnings in the script",
            "Verify script completed successfully"
        ]));
    }
    // Parse JSON with proper error handling
    let parsed;
    try {
        parsed = JSON.parse(output);
    }
    catch (parseError) {
        const errorMessage = parseError instanceof Error ? parseError.message : 'Unknown parse error';
        return createFailure(createError(CrystalErrorCode.INVALID_JSON_RESPONSE, `Failed to parse JSON from Python output: ${errorMessage}`, { output: output.substring(0, 200) }, [
            "Check JSON syntax in Python output",
            "Ensure proper JSON formatting",
            "Look for stray print statements or warnings in Python script"
        ]));
    }
    return createSuccess(parsed);
}
function isPythonResult(value) {
    if (!value || typeof value !== "object") {
        return false;
    }
    if (!("success" in value)) {
        return false;
    }
    const successValue = value.success;
    return typeof successValue === "boolean";
}
/**
 * Execute a Python script
 */
export async function executePython(options) {
    const config = PythonConfig.getInstance();
    // Get script path
    const scriptPath = config.getScriptPath(options.scriptName);
    // Verify script exists
    const scriptVerification = verifyScriptExists(scriptPath);
    if (!scriptVerification.success) {
        return createFailure(scriptVerification.error);
    }
    // Prepare execution options
    const pythonOptions = {
        mode: "text",
        pythonPath: options.pythonPath ?? config.getPythonPath(),
        scriptPath: path.dirname(scriptPath),
        args: options.args ? [...options.args] : [],
    };
    if (options.env) {
        pythonOptions.env = { ...process.env, ...options.env };
    }
    const startTime = Date.now();
    let timedOut = false;
    let timeoutId = null;
    // Set up timeout if specified
    return new Promise((resolve) => {
        const shell = new PythonShell(path.basename(scriptPath), pythonOptions);
        const stdoutLines = [];
        const stderrLines = [];
        // Set up timeout with proper process termination
        if (options.timeout) {
            timeoutId = setTimeout(() => {
                timedOut = true;
                // Kill the Python process when timeout occurs
                shell.terminate();
            }, options.timeout);
        }
        // Collect stdout
        shell.on("message", (message) => {
            if (!timedOut) {
                stdoutLines.push(message);
            }
        });
        // Collect stderr
        shell.on("stderr", (stderr) => {
            if (!timedOut) {
                stderrLines.push(stderr);
            }
        });
        // Handle completion
        shell.end((err, code) => {
            if (timeoutId) {
                clearTimeout(timeoutId);
            }
            const executionTime = Date.now() - startTime;
            // Check for timeout
            if (timedOut) {
                resolve(createFailure(createError(CrystalErrorCode.TIMEOUT_ERROR, `Python script execution timed out after ${options.timeout}ms`, { scriptName: options.scriptName, timeout: options.timeout }, [
                    "Increase timeout value",
                    "Optimize the Python script",
                    "Check for infinite loops or blocking operations"
                ])));
                return;
            }
            // Check for execution error
            if (err) {
                resolve(createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, `Python script execution failed: ${err.message}`, {
                    scriptName: options.scriptName,
                    error: err.message,
                    stderr: stderrLines.join("\n")
                }, [
                    "Check Python installation",
                    "Verify script syntax",
                    "Check for missing dependencies",
                    "Review error messages in stderr"
                ])));
                return;
            }
            // Check exit code
            const exitCode = code ?? 0;
            if (exitCode !== 0) {
                resolve(createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, `Python script exited with code ${exitCode}`, {
                    scriptName: options.scriptName,
                    exitCode,
                    stderr: stderrLines.join("\n")
                }, [
                    "Check script for errors",
                    "Review stderr output for details",
                    "Verify input parameters are correct"
                ])));
                return;
            }
            // Parse stdout as JSON if available
            const combinedOutput = stdoutLines.join("\n").trim();
            let parsedData = undefined;
            if (combinedOutput.length > 0) {
                const parseResult = parseJSONOutput(combinedOutput);
                if (!parseResult.success) {
                    resolve(createFailure(parseResult.error));
                    return;
                }
                parsedData = parseResult.data;
            }
            // Return successful result
            resolve(createSuccess({
                stdout: stdoutLines,
                stderr: stderrLines,
                data: parsedData,
                exitCode,
                executionTime
            }));
        });
    });
}
/**
 * Execute Python script with JSON input
 */
export async function executePythonWithJSON(scriptName, input, options) {
    // Serialize input to JSON
    let jsonInput;
    // Manual JSON stringification
    const stringifyResult = (() => {
        const result = JSON.stringify(input);
        return { success: true, data: result };
    })();
    if (!stringifyResult.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_OPERATION, "Failed to serialize input to JSON", { input }, ["Ensure input is JSON-serializable"]));
    }
    jsonInput = stringifyResult.data;
    // Write input to temporary file (cross-platform)
    const tempDir = os.tmpdir();
    const tempFileName = `crystal_mcp_${Date.now()}_${Math.random().toString(36).substring(7)}.json`;
    const tempFilePath = path.join(tempDir, tempFileName);
    // Write file synchronously with error checking
    fs.writeFileSync(tempFilePath, jsonInput, "utf-8");
    const writeResult = { success: true };
    if (!writeResult.success) {
        return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, "Failed to write temporary input file", { tempFilePath }, ["Check write permissions", "Ensure /tmp directory exists"]));
    }
    // Execute Python script with temp file as argument
    const execResult = await executePython({
        scriptName,
        args: [tempFilePath],
        ...options
    });
    // Clean up temp file
    if (fs.existsSync(tempFilePath)) {
        fs.unlinkSync(tempFilePath);
    }
    if (!execResult.success) {
        return execResult;
    }
    // Return parsed data
    if (execResult.data.data === undefined) {
        return createFailure(createError(CrystalErrorCode.INVALID_JSON_RESPONSE, "Python script did not return any data", { scriptName }, ["Ensure script prints JSON output"]));
    }
    if (!isPythonResult(execResult.data.data)) {
        return createFailure(createError(CrystalErrorCode.INVALID_JSON_RESPONSE, "Python response missing required 'success' flag", { scriptName, response: execResult.data.data }, ["Ensure Python scripts return {'success': true|false, ...}"]));
    }
    return createSuccess(execResult.data.data);
}
/**
 * Check if Python is available
 */
export function checkPythonAvailable(pythonPath = PythonConfig.getInstance().getPythonPath()) {
    // Check if python executable exists
    const result = spawnSync(pythonPath, ["--version"], {
        encoding: "utf-8",
        timeout: 5000
    });
    const checkResult = {
        success: result.status === 0,
        stdout: result.stdout,
        stderr: result.stderr,
        error: result.error
    };
    if (!checkResult.success) {
        return createFailure(createError(CrystalErrorCode.PYTHON_NOT_FOUND, `Python not found at: ${pythonPath}`, { pythonPath, error: checkResult.error }, [
            "Install Python 3.8 or later",
            "Ensure Python is in system PATH",
            "Specify correct python path in configuration"
        ]));
    }
    const version = checkResult.stdout || checkResult.stderr || "";
    const requiredModules = ["numpy", "pymatgen", "ase", "pyxtal", "spglib", "matplotlib"];
    const moduleCheck = spawnSync(pythonPath, [
        "-c",
        `import importlib.util, json; mods=${JSON.stringify(requiredModules)}; missing=[m for m in mods if importlib.util.find_spec(m) is None]; print(json.dumps(missing))`
    ], { encoding: "utf-8", timeout: 5000 });
    if (moduleCheck.status !== 0) {
        return createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, "Failed to verify Python dependencies", { pythonPath, stderr: moduleCheck.stderr }, ["Ensure required Python packages are installed"]));
    }
    const missingRaw = moduleCheck.stdout.trim();
    let missing = [];
    if (missingRaw.length > 0) {
        try {
            missing = JSON.parse(missingRaw);
        }
        catch {
            return createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, "Failed to parse Python dependency check output", { output: missingRaw }, ["Re-run dependency check or verify Python output"]));
        }
    }
    if (missing.length > 0) {
        return createFailure(createError(CrystalErrorCode.PYTHON_EXECUTION_FAILED, `Missing Python dependencies: ${missing.join(", ")}`, { missing }, ["Install missing packages with pip", "Verify active Python environment"]));
    }
    return createSuccess(version.trim());
}
/**
 * Get Python configuration instance
 */
export function getPythonConfig() {
    return PythonConfig.getInstance();
}
//# sourceMappingURL=python-bridge.js.map