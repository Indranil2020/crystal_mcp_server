/**
 * Structure Slice - Redux state for molecular structures
 */

import { createSlice, type PayloadAction } from '@reduxjs/toolkit';
import type { Structure, ViewerSettings, SelectionState, Measurement } from '../types';

interface StructureState {
    structures: Structure[];
    activeStructureId: string | null;
    viewerSettings: ViewerSettings;
    selection: SelectionState | null;
    measurements: Measurement[];
    isLoading: boolean;
}

const initialState: StructureState = {
    structures: [],
    activeStructureId: null,
    viewerSettings: {
        representation: 'ball-and-stick',
        colorScheme: 'element',
        showUnitCell: true,
        showAxes: false,
        showLabels: false,
        backgroundColor: '#0f172a',
    },
    selection: null,
    measurements: [],
    isLoading: false,
};

export const structureSlice = createSlice({
    name: 'structure',
    initialState,
    reducers: {
        addStructure: (state, action: PayloadAction<Structure>) => {
            state.structures.push(action.payload);
            state.activeStructureId = action.payload.id;
        },
        removeStructure: (state, action: PayloadAction<string>) => {
            state.structures = state.structures.filter(s => s.id !== action.payload);
            if (state.activeStructureId === action.payload) {
                state.activeStructureId = state.structures[0]?.id || null;
            }
        },
        setActiveStructure: (state, action: PayloadAction<string>) => {
            state.activeStructureId = action.payload;
        },
        toggleStructureVisibility: (state, action: PayloadAction<string>) => {
            const structure = state.structures.find(s => s.id === action.payload);
            if (structure) {
                structure.visible = !structure.visible;
            }
        },
        updateViewerSettings: (state, action: PayloadAction<Partial<ViewerSettings>>) => {
            state.viewerSettings = { ...state.viewerSettings, ...action.payload };
        },
        setSelection: (state, action: PayloadAction<SelectionState | null>) => {
            state.selection = action.payload;
            // Also store logic for updating visualization if needed?
            // Pure state update here.
        },
        addMeasurement: (state, action: PayloadAction<Measurement>) => {
            state.measurements.push(action.payload);
        },
        removeMeasurement: (state, action: PayloadAction<string>) => {
            state.measurements = state.measurements.filter(m => m.id !== action.payload);
        },
        clearMeasurements: (state) => {
            state.measurements = [];
        },
        setLoading: (state, action: PayloadAction<boolean>) => {
            state.isLoading = action.payload;
        },
        clearAllStructures: (state) => {
            state.structures = [];
            state.activeStructureId = null;
            state.selection = null;
            state.measurements = [];
        },
        deleteAtoms: (state, action: PayloadAction<{ structureId: string, atomIndices: number[] }>) => {
            const { structureId, atomIndices } = action.payload;
            const structure = state.structures.find(s => s.id === structureId);

            if (structure) {
                const indicesSet = new Set(atomIndices);
                structure.data.atoms = structure.data.atoms.filter((_atom, index) => !indicesSet.has(index));
                structure.modifiedAt = Date.now();

                if (state.selection?.structureId === structureId) {
                    state.selection = null;
                }
            }
        },
    },
});

export const {
    addStructure,
    removeStructure,
    setActiveStructure,
    toggleStructureVisibility,
    updateViewerSettings,
    setSelection,
    addMeasurement,
    removeMeasurement,
    clearMeasurements,
    setLoading,
    clearAllStructures,
    deleteAtoms,
} = structureSlice.actions;

export default structureSlice.reducer;
