/**
 * Structure Types - TypeScript interfaces for molecular structures
 */

import type { StructureData, LatticeData, AtomData, SpaceGroupData, StructureMetadata } from './mcp';

// Re-export from mcp.ts
export type { StructureData, LatticeData, AtomData, SpaceGroupData, StructureMetadata };

/** Structure with unique ID for state management */
export interface Structure {
    id: string;
    name: string;
    data: StructureData;
    source: 'mcp' | 'kekule' | 'file' | 'manual';
    createdAt: number;
    modifiedAt: number;
    visible: boolean;
}

/** Structure representation mode for MolStar */
export type RepresentationMode =
    | 'ball-and-stick'
    | 'spacefill'
    | 'cartoon'
    | 'licorice'
    | 'surface'
    | 'wireframe';

/** Color scheme for structure visualization */
export type ColorScheme =
    | 'element'
    | 'chain'
    | 'residue'
    | 'secondary-structure'
    | 'uniform'
    | 'property';

/** Viewer settings */
export interface ViewerSettings {
    representation: RepresentationMode;
    colorScheme: ColorScheme;
    showUnitCell: boolean;
    showAxes: boolean;
    showLabels: boolean;
    backgroundColor: string;
}

/** Measurement types */
export type MeasurementType = 'distance' | 'angle' | 'dihedral';

/** Measurement data */
export interface Measurement {
    id: string;
    type: MeasurementType;
    atomIndices: number[];
    value: number;
    unit: string;
    label?: string;
}

/** Selection state */
export interface SelectionState {
    structureId: string;
    atomIndices: number[];
    selectionMode: 'single' | 'range' | 'element' | 'residue';
}
