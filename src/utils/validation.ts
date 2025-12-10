/**
 * Validation Utilities
 * 
 * This module provides comprehensive validation functions for crystal structures,
 * input parameters, and data integrity checks. All validations follow defensive
 * programming principles without using try/catch blocks.
 */

import type {
  LatticeParameters,
  CrystalSystem
} from "../types/crystal.js";
import type { 
  CrystalError, 
  ValidationResult 
} from "../types/errors.js";
import { 
  CrystalErrorCode, 
  createError,
  createValidationSuccess,
  createValidationFailure,
  ERROR_MESSAGES
} from "../types/errors.js";

/**
 * Chemical elements data for validation
 */
const VALID_ELEMENTS = new Set([
  "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
  "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
  "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]);

/**
 * Covalent radii in Angstroms for minimum distance calculations
 */
const COVALENT_RADII: Record<string, number> = {
  "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84,
  "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
  "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
  "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76,
  "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39,
  "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
  "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
  "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
  "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42,
  "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39,
  "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
};

/**
 * Validate space group number
 */
export function validateSpaceGroup(spaceGroup: number | string): ValidationResult {
  const errors: CrystalError[] = [];
  
  if (typeof spaceGroup === "number") {
    if (!Number.isInteger(spaceGroup)) {
      errors.push(createError(
        CrystalErrorCode.INVALID_SPACE_GROUP,
        "Space group number must be an integer",
        { spaceGroup },
        ["Provide an integer space group number between 1 and 230"]
      ));
    } else if (spaceGroup < 1 || spaceGroup > 230) {
      const msg = ERROR_MESSAGES.INVALID_SPACE_GROUP(spaceGroup);
      errors.push(createError(
        CrystalErrorCode.INVALID_SPACE_GROUP,
        msg.message,
        { spaceGroup },
        msg.suggestions
      ));
    }
  } else if (typeof spaceGroup === "string") {
    if (spaceGroup.trim().length === 0) {
      errors.push(createError(
        CrystalErrorCode.INVALID_SPACE_GROUP,
        "Space group symbol cannot be empty",
        { spaceGroup },
        ["Provide a valid Hermann-Mauguin symbol or space group number"]
      ));
    }
  } else {
    errors.push(createError(
      CrystalErrorCode.INVALID_SPACE_GROUP,
      "Space group must be a number or string",
      { spaceGroup, type: typeof spaceGroup },
      ["Provide a space group number (1-230) or Hermann-Mauguin symbol"]
    ));
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors);
  }
  
  return createValidationSuccess();
}

/**
 * Validate chemical composition
 */
export function validateComposition(composition: readonly string[]): ValidationResult {
  const errors: CrystalError[] = [];
  const warnings: string[] = [];
  
  if (composition.length === 0) {
    errors.push(createError(
      CrystalErrorCode.INVALID_COMPOSITION,
      "Composition cannot be empty",
      { composition },
      ["Provide at least one element symbol"]
    ));
    return createValidationFailure(errors, warnings);
  }
  
  for (let i = 0; i < composition.length; i++) {
    const element = composition[i];
    
    if (!element || typeof element !== "string") {
      errors.push(createError(
        CrystalErrorCode.INVALID_COMPOSITION,
        `Invalid element at index ${i}`,
        { element, index: i },
        ["Provide valid chemical element symbols (e.g., 'Si', 'O', 'Na')"]
      ));
      continue;
    }
    
    const trimmed = element.trim();
    if (trimmed.length === 0) {
      errors.push(createError(
        CrystalErrorCode.INVALID_COMPOSITION,
        `Empty element at index ${i}`,
        { element, index: i },
        ["Remove empty strings from composition"]
      ));
      continue;
    }
    
    if (!VALID_ELEMENTS.has(trimmed)) {
      errors.push(createError(
        CrystalErrorCode.INVALID_COMPOSITION,
        `Unknown element '${trimmed}' at index ${i}`,
        { element: trimmed, index: i },
        [
          "Check element symbol spelling",
          "Ensure proper capitalization (e.g., 'Si' not 'SI')",
          `Valid elements: ${Array.from(VALID_ELEMENTS).slice(0, 20).join(", ")}...`
        ]
      ));
    }
  }
  
  // Check for unusual element combinations
  const uniqueElements = new Set(composition);
  if (uniqueElements.size > 6) {
    warnings.push(
      `Composition contains ${uniqueElements.size} different elements. ` +
      "Complex compositions may be difficult to generate."
    );
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors, warnings);
  }
  
  return { valid: true, errors: [], warnings };
}

/**
 * Validate lattice parameters
 */
export function validateLatticeParameters(
  params: Partial<LatticeParameters>
): ValidationResult {
  const errors: CrystalError[] = [];
  const warnings: string[] = [];
  
  // Check lattice lengths
  const lengths = ["a", "b", "c"] as const;
  for (const key of lengths) {
    const value = params[key];
    if (value !== undefined) {
      if (typeof value !== "number" || !Number.isFinite(value)) {
        errors.push(createError(
          CrystalErrorCode.INVALID_LATTICE_PARAMS,
          `Lattice parameter ${key} must be a finite number`,
          { [key]: value },
          [`Provide a positive number for ${key} in Angstroms`]
        ));
      } else if (value <= 0) {
        errors.push(createError(
          CrystalErrorCode.INVALID_LATTICE_PARAMS,
          `Lattice parameter ${key} must be positive`,
          { [key]: value },
          [`${key} should be > 0 Angstroms`]
        ));
      } else if (value < 2.0) {
        warnings.push(
          `Lattice parameter ${key} = ${value.toFixed(2)} Å is very small (< 2 Å)`
        );
      } else if (value > 50.0) {
        warnings.push(
          `Lattice parameter ${key} = ${value.toFixed(2)} Å is very large (> 50 Å)`
        );
      }
    }
  }
  
  // Check angles
  const angles = ["alpha", "beta", "gamma"] as const;
  for (const key of angles) {
    const value = params[key];
    if (value !== undefined) {
      if (typeof value !== "number" || !Number.isFinite(value)) {
        errors.push(createError(
          CrystalErrorCode.INVALID_LATTICE_PARAMS,
          `Angle ${key} must be a finite number`,
          { [key]: value },
          [`Provide an angle for ${key} in degrees (0-180)`]
        ));
      } else if (value <= 0 || value >= 180) {
        errors.push(createError(
          CrystalErrorCode.INVALID_LATTICE_PARAMS,
          `Angle ${key} must be between 0 and 180 degrees`,
          { [key]: value },
          [`${key} should be in range (0, 180) degrees`]
        ));
      } else if (value < 30 || value > 150) {
        warnings.push(
          `Angle ${key} = ${value.toFixed(1)}° is unusual (< 30° or > 150°)`
        );
      }
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors, warnings);
  }
  
  return { valid: true, errors: [], warnings };
}

/**
 * Validate volume factor
 */
export function validateVolumeFactor(factor: number): ValidationResult {
  const errors: CrystalError[] = [];
  const warnings: string[] = [];
  
  if (typeof factor !== "number" || !Number.isFinite(factor)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_VOLUME_FACTOR,
      "Volume factor must be a finite number",
      { factor },
      ["Provide a positive number (typically 0.8 to 2.0)"]
    ));
  } else if (factor <= 0) {
    errors.push(createError(
      CrystalErrorCode.INVALID_VOLUME_FACTOR,
      "Volume factor must be positive",
      { factor },
      ["Use values > 0, typically around 1.0"]
    ));
  } else {
    if (factor < 0.5) {
      warnings.push(
        `Volume factor ${factor} is very small (< 0.5). ` +
        "Atoms may be too close together."
      );
    } else if (factor > 3.0) {
      warnings.push(
        `Volume factor ${factor} is very large (> 3.0). ` +
        "Structure may be very sparse."
      );
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors, warnings);
  }
  
  return { valid: true, errors: [], warnings };
}

/**
 * Validate minimum distance constraints
 */
export function validateMinDistance(
  minDistance: Record<string, number>
): ValidationResult {
  const errors: CrystalError[] = [];
  const warnings: string[] = [];
  
  for (const [pair, distance] of Object.entries(minDistance)) {
    // Validate pair format
    const elements = pair.split("-");
    if (elements.length !== 2) {
      errors.push(createError(
        CrystalErrorCode.INVALID_MIN_DISTANCE,
        `Invalid element pair format: '${pair}'`,
        { pair },
        ["Use format 'Element1-Element2', e.g., 'Si-O' or 'Na-Cl'"]
      ));
      continue;
    }
    
    const elem1 = elements[0];
    const elem2 = elements[1];

    if (!elem1 || !elem2 || !VALID_ELEMENTS.has(elem1) || !VALID_ELEMENTS.has(elem2)) {
      errors.push(createError(
        CrystalErrorCode.INVALID_MIN_DISTANCE,
        `Invalid elements in pair '${pair}'`,
        { pair, elem1, elem2 },
        ["Use valid element symbols"]
      ));
      continue;
    }

    // Validate distance value
    if (typeof distance !== "number" || !Number.isFinite(distance)) {
      errors.push(createError(
        CrystalErrorCode.INVALID_MIN_DISTANCE,
        `Invalid distance for pair '${pair}'`,
        { pair, distance },
        ["Provide a positive number in Angstroms"]
      ));
      continue;
    }

    if (distance <= 0) {
      errors.push(createError(
        CrystalErrorCode.INVALID_MIN_DISTANCE,
        `Distance for pair '${pair}' must be positive`,
        { pair, distance },
        ["Use positive values in Angstroms"]
      ));
      continue;
    }

    // Check against covalent radii
    const r1 = COVALENT_RADII[elem1] ?? 1.5;
    const r2 = COVALENT_RADII[elem2] ?? 1.5;
    const expectedMin = (r1 + r2) * 0.7;  // 70% of sum of covalent radii
    
    if (distance < expectedMin) {
      warnings.push(
        `Minimum distance for ${pair} (${distance.toFixed(2)} Å) is ` +
        `smaller than expected (${expectedMin.toFixed(2)} Å)`
      );
    }
    
    if (distance > (r1 + r2) * 2.5) {
      warnings.push(
        `Minimum distance for ${pair} (${distance.toFixed(2)} Å) is ` +
        "unusually large and may prevent structure generation"
      );
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors, warnings);
  }
  
  return { valid: true, errors: [], warnings };
}

/**
 * Validate number of atoms
 */
export function validateNumAtoms(numAtoms: number | undefined): ValidationResult {
  const errors: CrystalError[] = [];
  const warnings: string[] = [];
  
  if (numAtoms === undefined) {
    return createValidationSuccess();
  }
  
  if (typeof numAtoms !== "number" || !Number.isFinite(numAtoms)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_NUM_ATOMS,
      "Number of atoms must be a finite number",
      { numAtoms },
      ["Provide a positive integer"]
    ));
  } else if (!Number.isInteger(numAtoms)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_NUM_ATOMS,
      "Number of atoms must be an integer",
      { numAtoms },
      ["Provide a whole number"]
    ));
  } else if (numAtoms <= 0) {
    errors.push(createError(
      CrystalErrorCode.INVALID_NUM_ATOMS,
      "Number of atoms must be positive",
      { numAtoms },
      ["Use a positive integer"]
    ));
  } else {
    if (numAtoms > 1000) {
      warnings.push(
        `Large number of atoms (${numAtoms}) may be slow to generate. ` +
        "Consider generating smaller cells and using supercells."
      );
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors, warnings);
  }
  
  return { valid: true, errors: [], warnings };
}

/**
 * Validate Vector3
 */
export function validateVector3(vec: unknown, name: string = "vector"): ValidationResult {
  const errors: CrystalError[] = [];
  
  if (!Array.isArray(vec)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_OPERATION,
      `${name} must be an array`,
      { value: vec },
      ["Provide a 3-element array of numbers"]
    ));
    return createValidationFailure(errors);
  }
  
  if (vec.length !== 3) {
    errors.push(createError(
      CrystalErrorCode.INVALID_OPERATION,
      `${name} must have exactly 3 elements`,
      { value: vec, length: vec.length },
      ["Provide exactly 3 numbers"]
    ));
    return createValidationFailure(errors);
  }
  
  for (let i = 0; i < 3; i++) {
    if (typeof vec[i] !== "number" || !Number.isFinite(vec[i])) {
      errors.push(createError(
        CrystalErrorCode.INVALID_OPERATION,
        `${name}[${i}] must be a finite number`,
        { value: vec[i], index: i },
        ["All elements must be finite numbers"]
      ));
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors);
  }
  
  return createValidationSuccess();
}

/**
 * Validate Matrix3x3
 */
export function validateMatrix3x3(mat: unknown, name: string = "matrix"): ValidationResult {
  const errors: CrystalError[] = [];
  
  if (!Array.isArray(mat)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_OPERATION,
      `${name} must be an array`,
      { value: mat },
      ["Provide a 3x3 matrix as nested arrays"]
    ));
    return createValidationFailure(errors);
  }
  
  if (mat.length !== 3) {
    errors.push(createError(
      CrystalErrorCode.INVALID_OPERATION,
      `${name} must have exactly 3 rows`,
      { rows: mat.length },
      ["Provide a 3x3 matrix"]
    ));
    return createValidationFailure(errors);
  }
  
  for (let i = 0; i < 3; i++) {
    if (!Array.isArray(mat[i])) {
      errors.push(createError(
        CrystalErrorCode.INVALID_OPERATION,
        `${name}[${i}] must be an array`,
        { value: mat[i] },
        ["Each row must be an array of 3 numbers"]
      ));
      continue;
    }
    
    if (mat[i].length !== 3) {
      errors.push(createError(
        CrystalErrorCode.INVALID_OPERATION,
        `${name}[${i}] must have exactly 3 elements`,
        { length: mat[i].length },
        ["Each row must have 3 numbers"]
      ));
      continue;
    }
    
    for (let j = 0; j < 3; j++) {
      if (typeof mat[i][j] !== "number" || !Number.isFinite(mat[i][j])) {
        errors.push(createError(
          CrystalErrorCode.INVALID_OPERATION,
          `${name}[${i}][${j}] must be a finite number`,
          { value: mat[i][j] },
          ["All matrix elements must be finite numbers"]
        ));
      }
    }
  }
  
  if (errors.length > 0) {
    return createValidationFailure(errors);
  }
  
  return createValidationSuccess();
}

/**
 * Get expected minimum distance between two elements
 */
export function getExpectedMinDistance(
  elem1: string, 
  elem2: string,
  scale: number = 0.8
): number {
  const r1 = COVALENT_RADII[elem1] ?? 1.5;
  const r2 = COVALENT_RADII[elem2] ?? 1.5;
  return (r1 + r2) * scale;
}

/**
 * Check if a value is a valid crystal system
 */
export function isValidCrystalSystem(system: string): system is CrystalSystem {
  return [
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic"
  ].includes(system);
}

/**
 * Validate crystal system
 */
export function validateCrystalSystem(system: string): ValidationResult {
  const errors: CrystalError[] = [];
  
  if (!isValidCrystalSystem(system)) {
    errors.push(createError(
      CrystalErrorCode.INVALID_OPERATION,
      `Invalid crystal system: ${system}`,
      { system },
      [
        "Valid crystal systems: triclinic, monoclinic, orthorhombic, " +
        "tetragonal, trigonal, hexagonal, cubic"
      ]
    ));
    return createValidationFailure(errors);
  }
  
  return createValidationSuccess();
}
