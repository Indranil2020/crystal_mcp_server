/**
 * Validation Utilities
 *
 * This module provides comprehensive validation functions for crystal structures,
 * input parameters, and data integrity checks. All validations follow defensive
 * programming principles without using try/catch blocks.
 */
import type { LatticeParameters, CrystalSystem } from "../types/crystal.js";
import type { ValidationResult } from "../types/errors.js";
/**
 * Validate space group number
 */
export declare function validateSpaceGroup(spaceGroup: number | string): ValidationResult;
/**
 * Validate chemical composition
 */
export declare function validateComposition(composition: readonly string[]): ValidationResult;
/**
 * Validate lattice parameters
 */
export declare function validateLatticeParameters(params: Partial<LatticeParameters>): ValidationResult;
/**
 * Validate volume factor
 */
export declare function validateVolumeFactor(factor: number): ValidationResult;
/**
 * Validate minimum distance constraints
 */
export declare function validateMinDistance(minDistance: Record<string, number>): ValidationResult;
/**
 * Validate number of atoms
 */
export declare function validateNumAtoms(numAtoms: number | undefined): ValidationResult;
/**
 * Validate Vector3
 */
export declare function validateVector3(vec: unknown, name?: string): ValidationResult;
/**
 * Validate Matrix3x3
 */
export declare function validateMatrix3x3(mat: unknown, name?: string): ValidationResult;
/**
 * Get expected minimum distance between two elements
 */
export declare function getExpectedMinDistance(elem1: string, elem2: string, scale?: number): number;
/**
 * Check if a value is a valid crystal system
 */
export declare function isValidCrystalSystem(system: string): system is CrystalSystem;
/**
 * Validate crystal system
 */
export declare function validateCrystalSystem(system: string): ValidationResult;
//# sourceMappingURL=validation.d.ts.map