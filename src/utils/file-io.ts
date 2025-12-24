/**
 * File I/O Utilities
 * 
 * This module provides safe file reading and writing operations
 * with comprehensive validation and error handling.
 */

import * as fs from "fs";
import * as path from "path";
import type { Result } from "../types/errors.js";
import { 
  createSuccess, 
  createFailure, 
  createError, 
  CrystalErrorCode 
} from "../types/errors.js";

/**
 * Check if a file exists and is readable
 */
export function fileExists(filePath: string): boolean {
  if (!filePath || typeof filePath !== "string") {
    return false;
  }

  try {
    if (!fs.existsSync(filePath)) {
      return false;
    }
    const stats = fs.statSync(filePath);
    return stats.isFile();
  } catch {
    return false;
  }
}

/**
 * Check if a directory exists
 */
export function directoryExists(dirPath: string): boolean {
  if (!dirPath || typeof dirPath !== "string") {
    return false;
  }

  try {
    if (!fs.existsSync(dirPath)) {
      return false;
    }
    const stats = fs.statSync(dirPath);
    return stats.isDirectory();
  } catch {
    return false;
  }
}

/**
 * Read file contents as string
 */
export function readFileSync(filePath: string): Result<string> {
  if (!filePath || typeof filePath !== "string") {
    return createFailure(createError(
      CrystalErrorCode.FILE_READ_ERROR,
      "Invalid file path",
      { filePath },
      ["Provide a valid file path string"]
    ));
  }

  try {
    if (!fs.existsSync(filePath)) {
      return createFailure(createError(
        CrystalErrorCode.FILE_NOT_FOUND,
        `File not found: ${filePath}`,
        { filePath },
        [
          "Check that the file exists",
          "Verify the file path is correct",
          "Ensure you have read permissions"
        ]
      ));
    }

    const stats = fs.statSync(filePath);
    if (!stats.isFile()) {
      return createFailure(createError(
        CrystalErrorCode.FILE_READ_ERROR,
        `Path is not a file: ${filePath}`,
        { filePath },
        ["Provide a path to a file, not a directory"]
      ));
    }

    // Check read permissions
    fs.accessSync(filePath, fs.constants.R_OK);

    const content = fs.readFileSync(filePath, "utf-8");
    return createSuccess(content);
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : 'Unknown error';
    return createFailure(createError(
      CrystalErrorCode.FILE_READ_ERROR,
      `Failed to read file: ${errorMessage}`,
      { filePath },
      ["Check file permissions", "Ensure file is not locked by another process"]
    ));
  }
}

/**
 * Write string contents to file
 */
export function writeFileSync(
  filePath: string, 
  content: string
): Result<void> {
  if (!filePath || typeof filePath !== "string") {
    return createFailure(createError(
      CrystalErrorCode.FILE_WRITE_ERROR,
      "Invalid file path",
      { filePath },
      ["Provide a valid file path string"]
    ));
  }
  
  if (typeof content !== "string") {
    return createFailure(createError(
      CrystalErrorCode.FILE_WRITE_ERROR,
      "Content must be a string",
      { contentType: typeof content },
      ["Provide string content to write"]
    ));
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
export function ensureDirectory(dirPath: string): Result<void> {
  if (!dirPath || typeof dirPath !== "string") {
    return createFailure(createError(
      CrystalErrorCode.FILE_WRITE_ERROR,
      "Invalid directory path",
      { dirPath },
      ["Provide a valid directory path string"]
    ));
  }
  
  if (fs.existsSync(dirPath)) {
    const stats = fs.statSync(dirPath);
    if (!stats.isDirectory()) {
      return createFailure(createError(
        CrystalErrorCode.FILE_WRITE_ERROR,
        `Path exists but is not a directory: ${dirPath}`,
        { dirPath },
        ["Use a different path or remove the existing file"]
      ));
    }
    return createSuccess(undefined);
  }
  
  fs.mkdirSync(dirPath, { recursive: true });
  return createSuccess(undefined);
}

/**
 * Delete file if it exists
 */
export function deleteFileIfExists(filePath: string): Result<boolean> {
  if (!filePath || typeof filePath !== "string") {
    return createFailure(createError(
      CrystalErrorCode.FILE_WRITE_ERROR,
      "Invalid file path",
      { filePath },
      ["Provide a valid file path string"]
    ));
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
export function readJSONFile<T>(filePath: string): Result<T> {
  const readResult = readFileSync(filePath);
  if (!readResult.success) {
    return readResult;
  }
  
  const content = readResult.data.trim();
  if (content.length === 0) {
    return createFailure(createError(
      CrystalErrorCode.CORRUPTED_DATA,
      `File is empty: ${filePath}`,
      { filePath },
      ["Ensure file contains valid JSON data"]
    ));
  }
  
  const parsed = JSON.parse(content);
  return createSuccess(parsed as T);
}

/**
 * Write object to JSON file
 */
export function writeJSONFile<T>(
  filePath: string, 
  data: T,
  pretty: boolean = true
): Result<void> {
  const jsonString = JSON.stringify(data, null, pretty ? 2 : 0);
  return writeFileSync(filePath, jsonString);
}

/**
 * Get file extension
 */
export function getFileExtension(filePath: string): string {
  if (!filePath || typeof filePath !== "string") {
    return "";
  }
  
  const ext = path.extname(filePath);
  return ext.toLowerCase().replace(".", "");
}

/**
 * Check if file has specific extension
 */
export function hasExtension(filePath: string, extension: string): boolean {
  const fileExt = getFileExtension(filePath);
  const targetExt = extension.toLowerCase().replace(".", "");
  return fileExt === targetExt;
}

/**
 * Generate unique filename in directory
 */
export function generateUniqueFilename(
  directory: string,
  baseName: string,
  extension: string
): string {
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
