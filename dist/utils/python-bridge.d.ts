/**
 * Python Bridge Utility
 *
 * This module provides a safe interface for executing Python scripts
 * from TypeScript, with comprehensive error handling and validation.
 * Uses defensive programming without try/catch blocks.
 */
import type { Result } from "../types/errors.js";
/**
 * Python execution options
 */
export interface PythonExecutionOptions {
    readonly scriptName: string;
    readonly args?: readonly string[];
    readonly pythonPath?: string;
    readonly cwd?: string;
    readonly timeout?: number;
    readonly env?: Record<string, string>;
}
/**
 * Python execution result
 */
export interface PythonExecutionResult<T = unknown> {
    readonly stdout: readonly string[];
    readonly stderr: readonly string[];
    readonly data?: T;
    readonly exitCode: number;
    readonly executionTime: number;
}
/**
 * Configuration for Python environment
 */
export declare class PythonConfig {
    private static instance;
    private pythonPath;
    private scriptsDirectory;
    private constructor();
    static getInstance(): PythonConfig;
    setPythonPath(pythonPath: string): void;
    getPythonPath(): string;
    setScriptsDirectory(dir: string): void;
    getScriptsDirectory(): string;
    /**
     * Get full path to a Python script
     */
    getScriptPath(scriptName: string): string;
}
/**
 * Execute a Python script
 */
export declare function executePython<T = unknown>(options: PythonExecutionOptions): Promise<Result<PythonExecutionResult<T>>>;
/**
 * Execute Python script with JSON input
 */
export declare function executePythonWithJSON<TInput, TOutput>(scriptName: string, input: TInput, options?: Partial<PythonExecutionOptions>): Promise<Result<TOutput>>;
/**
 * Check if Python is available
 */
export declare function checkPythonAvailable(pythonPath?: string): Result<string>;
/**
 * Get Python configuration instance
 */
export declare function getPythonConfig(): PythonConfig;
//# sourceMappingURL=python-bridge.d.ts.map