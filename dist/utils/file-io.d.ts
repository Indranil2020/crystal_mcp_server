/**
 * File I/O Utilities
 *
 * This module provides safe file reading and writing operations
 * with comprehensive validation and error handling.
 */
import type { Result } from "../types/errors.js";
/**
 * Check if a file exists and is readable
 */
export declare function fileExists(filePath: string): boolean;
/**
 * Check if a directory exists
 */
export declare function directoryExists(dirPath: string): boolean;
/**
 * Read file contents as string
 */
export declare function readFileSync(filePath: string): Result<string>;
/**
 * Write string contents to file
 */
export declare function writeFileSync(filePath: string, content: string): Result<void>;
/**
 * Ensure directory exists, create if it doesn't
 */
export declare function ensureDirectory(dirPath: string): Result<void>;
/**
 * Delete file if it exists
 */
export declare function deleteFileIfExists(filePath: string): Result<boolean>;
/**
 * Read JSON file and parse
 */
export declare function readJSONFile<T>(filePath: string): Result<T>;
/**
 * Write object to JSON file
 */
export declare function writeJSONFile<T>(filePath: string, data: T, pretty?: boolean): Result<void>;
/**
 * Get file extension
 */
export declare function getFileExtension(filePath: string): string;
/**
 * Check if file has specific extension
 */
export declare function hasExtension(filePath: string, extension: string): boolean;
/**
 * Generate unique filename in directory
 */
export declare function generateUniqueFilename(directory: string, baseName: string, extension: string): string;
//# sourceMappingURL=file-io.d.ts.map