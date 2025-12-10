/**
 * File I/O Utilities
 *
 * This module provides safe file reading and writing operations
 * with comprehensive validation and error handling.
 */
import * as fs from "fs";
import * as path from "path";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../types/errors.js";
/**
 * Check if a file exists and is readable
 */
export function fileExists(filePath) {
    if (!filePath || typeof filePath !== "string") {
        return false;
    }
    if (!fs.existsSync(filePath)) {
        return false;
    }
    const stats = fs.statSync(filePath);
    return stats.isFile();
}
/**
 * Check if a directory exists
 */
export function directoryExists(dirPath) {
    if (!dirPath || typeof dirPath !== "string") {
        return false;
    }
    if (!fs.existsSync(dirPath)) {
        return false;
    }
    const stats = fs.statSync(dirPath);
    return stats.isDirectory();
}
/**
 * Read file contents as string
 */
export function readFileSync(filePath) {
    if (!filePath || typeof filePath !== "string") {
        return createFailure(createError(CrystalErrorCode.FILE_READ_ERROR, "Invalid file path", { filePath }, ["Provide a valid file path string"]));
    }
    if (!fs.existsSync(filePath)) {
        return createFailure(createError(CrystalErrorCode.FILE_NOT_FOUND, `File not found: ${filePath}`, { filePath }, [
            "Check that the file exists",
            "Verify the file path is correct",
            "Ensure you have read permissions"
        ]));
    }
    const stats = fs.statSync(filePath);
    if (!stats.isFile()) {
        return createFailure(createError(CrystalErrorCode.FILE_READ_ERROR, `Path is not a file: ${filePath}`, { filePath }, ["Provide a path to a file, not a directory"]));
    }
    // Check read permissions
    fs.accessSync(filePath, fs.constants.R_OK);
    const content = fs.readFileSync(filePath, "utf-8");
    return createSuccess(content);
}
/**
 * Write string contents to file
 */
export function writeFileSync(filePath, content) {
    if (!filePath || typeof filePath !== "string") {
        return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, "Invalid file path", { filePath }, ["Provide a valid file path string"]));
    }
    if (typeof content !== "string") {
        return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, "Content must be a string", { contentType: typeof content }, ["Provide string content to write"]));
    }
    // Create directory if it doesn't exist
    const dirPath = path.dirname(filePath);
    if (!directoryExists(dirPath)) {
        const mkdirResult = ensureDirectory(dirPath);
        if (!mkdirResult.success) {
            return mkdirResult;
        }
    }
    fs.writeFileSync(filePath, content, "utf-8");
    return createSuccess(undefined);
}
/**
 * Ensure directory exists, create if it doesn't
 */
export function ensureDirectory(dirPath) {
    if (!dirPath || typeof dirPath !== "string") {
        return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, "Invalid directory path", { dirPath }, ["Provide a valid directory path string"]));
    }
    if (fs.existsSync(dirPath)) {
        const stats = fs.statSync(dirPath);
        if (!stats.isDirectory()) {
            return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, `Path exists but is not a directory: ${dirPath}`, { dirPath }, ["Use a different path or remove the existing file"]));
        }
        return createSuccess(undefined);
    }
    fs.mkdirSync(dirPath, { recursive: true });
    return createSuccess(undefined);
}
/**
 * Delete file if it exists
 */
export function deleteFileIfExists(filePath) {
    if (!filePath || typeof filePath !== "string") {
        return createFailure(createError(CrystalErrorCode.FILE_WRITE_ERROR, "Invalid file path", { filePath }, ["Provide a valid file path string"]));
    }
    if (!fs.existsSync(filePath)) {
        return createSuccess(false);
    }
    fs.unlinkSync(filePath);
    return createSuccess(true);
}
/**
 * Read JSON file and parse
 */
export function readJSONFile(filePath) {
    const readResult = readFileSync(filePath);
    if (!readResult.success) {
        return readResult;
    }
    const content = readResult.data.trim();
    if (content.length === 0) {
        return createFailure(createError(CrystalErrorCode.CORRUPTED_DATA, `File is empty: ${filePath}`, { filePath }, ["Ensure file contains valid JSON data"]));
    }
    const parsed = JSON.parse(content);
    return createSuccess(parsed);
}
/**
 * Write object to JSON file
 */
export function writeJSONFile(filePath, data, pretty = true) {
    const jsonString = JSON.stringify(data, null, pretty ? 2 : 0);
    return writeFileSync(filePath, jsonString);
}
/**
 * Get file extension
 */
export function getFileExtension(filePath) {
    if (!filePath || typeof filePath !== "string") {
        return "";
    }
    const ext = path.extname(filePath);
    return ext.toLowerCase().replace(".", "");
}
/**
 * Check if file has specific extension
 */
export function hasExtension(filePath, extension) {
    const fileExt = getFileExtension(filePath);
    const targetExt = extension.toLowerCase().replace(".", "");
    return fileExt === targetExt;
}
/**
 * Generate unique filename in directory
 */
export function generateUniqueFilename(directory, baseName, extension) {
    const ext = extension.startsWith(".") ? extension : `.${extension}`;
    let filePath = path.join(directory, `${baseName}${ext}`);
    if (!fs.existsSync(filePath)) {
        return filePath;
    }
    let counter = 1;
    while (fs.existsSync(filePath)) {
        filePath = path.join(directory, `${baseName}_${counter}${ext}`);
        counter++;
    }
    return filePath;
}
//# sourceMappingURL=file-io.js.map