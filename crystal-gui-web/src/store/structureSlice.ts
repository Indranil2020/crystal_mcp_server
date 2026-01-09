/**
 * Structure Slice - Redux state for molecular structures
 * 
 * Includes undo/redo history management with debug logging.
 * Max history depth: 50 states.
 */

import { createSlice, type PayloadAction } from '@reduxjs/toolkit';
import type { Structure, ViewerSettings, SelectionState, Measurement } from '../types';

// Debug logger for undo/redo operations
function debugLog(action: string, details?: string): void {
    const timestamp = new Date().toISOString().split('T')[1].slice(0, 12);
    console.log(`%c[REDUX/UNDO] ${timestamp} ${action}${details ? ': ' + details : ''}`, 'color: #FF5722; font-weight: bold');
}

// Maximum history depth to prevent memory issues
const MAX_HISTORY_DEPTH = 50;

// Snapshot of structure state for undo/redo (excludes history stacks)
interface StructureSnapshot {
    structures: Structure[];
    activeStructureId: string | null;
    selection: SelectionState | null;
    measurements: Measurement[];
}

interface StructureState {
    structures: Structure[];
    activeStructureId: string | null;
    viewerSettings: ViewerSettings;
    selection: SelectionState | null;
    measurements: Measurement[];
    isLoading: boolean;
    loadingMessage: string | null;
    loadingProgress: number | null;
    // Undo/Redo history
    undoStack: StructureSnapshot[];
    redoStack: StructureSnapshot[];
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
    loadingMessage: null,
    loadingProgress: null,
    undoStack: [],
    redoStack: [],
};

// Helper to create a snapshot of the current state (for undo/redo)
function createSnapshot(state: StructureState): StructureSnapshot {
    return {
        structures: JSON.parse(JSON.stringify(state.structures)),
        activeStructureId: state.activeStructureId,
        selection: state.selection ? JSON.parse(JSON.stringify(state.selection)) : null,
        measurements: JSON.parse(JSON.stringify(state.measurements)),
    };
}

// Helper to restore state from a snapshot
function restoreFromSnapshot(state: StructureState, snapshot: StructureSnapshot): void {
    state.structures = snapshot.structures;
    state.activeStructureId = snapshot.activeStructureId;
    state.selection = snapshot.selection;
    state.measurements = snapshot.measurements;
}

// Helper to push current state to undo stack (called before destructive operations)
function pushToUndoStack(state: StructureState): void {
    const snapshot = createSnapshot(state);
    state.undoStack.push(snapshot);

    // Limit history depth
    if (state.undoStack.length > MAX_HISTORY_DEPTH) {
        state.undoStack.shift();
        debugLog('TRIM', `Undo stack trimmed to ${MAX_HISTORY_DEPTH} entries`);
    }

    // Clear redo stack on new action
    state.redoStack = [];

    debugLog('PUSH', `Undo stack size: ${state.undoStack.length}, Redo stack cleared`);
}

export const structureSlice = createSlice({
    name: 'structure',
    initialState,
    reducers: {
        // === Structure Management ===
        addStructure: (state, action: PayloadAction<Structure>) => {
            pushToUndoStack(state);
            debugLog('ACTION', `addStructure: ${action.payload.name}`);
            state.structures.push(action.payload);
            state.activeStructureId = action.payload.id;
        },

        removeStructure: (state, action: PayloadAction<string>) => {
            pushToUndoStack(state);
            debugLog('ACTION', `removeStructure: ${action.payload}`);
            state.structures = state.structures.filter(s => s.id !== action.payload);
            if (state.activeStructureId === action.payload) {
                state.activeStructureId = state.structures[0]?.id || null;
            }
        },

        setActiveStructure: (state, action: PayloadAction<string>) => {
            // No undo for this - it's just a view change
            state.activeStructureId = action.payload;
        },

        toggleStructureVisibility: (state, action: PayloadAction<string>) => {
            const structure = state.structures.find(s => s.id === action.payload);
            if (structure) {
                structure.visible = !structure.visible;
            }
        },

        // === Viewer Settings ===
        updateViewerSettings: (state, action: PayloadAction<Partial<ViewerSettings>>) => {
            // No undo for settings changes
            state.viewerSettings = { ...state.viewerSettings, ...action.payload };
        },

        // === Selection ===
        setSelection: (state, action: PayloadAction<SelectionState | null>) => {
            // No undo for selection changes
            state.selection = action.payload;
        },

        // === Measurements ===
        addMeasurement: (state, action: PayloadAction<Measurement>) => {
            pushToUndoStack(state);
            debugLog('ACTION', `addMeasurement: ${action.payload.type}`);
            state.measurements.push(action.payload);
        },

        removeMeasurement: (state, action: PayloadAction<string>) => {
            pushToUndoStack(state);
            debugLog('ACTION', `removeMeasurement: ${action.payload}`);
            state.measurements = state.measurements.filter(m => m.id !== action.payload);
        },

        clearMeasurements: (state) => {
            if (state.measurements.length > 0) {
                pushToUndoStack(state);
                debugLog('ACTION', 'clearMeasurements');
            }
            state.measurements = [];
        },

        // === Loading State ===
        setLoading: (state, action: PayloadAction<boolean>) => {
            state.isLoading = action.payload;
            if (!action.payload) {
                state.loadingMessage = null;
                state.loadingProgress = null;
            }
        },

        setLoadingState: (state, action: PayloadAction<{ message?: string; progress?: number | null }>) => {
            if (action.payload.message !== undefined) {
                state.loadingMessage = action.payload.message;
            }
            if (action.payload.progress !== undefined) {
                state.loadingProgress = action.payload.progress;
            }
        },

        // === Clear All ===
        clearAllStructures: (state) => {
            if (state.structures.length > 0) {
                pushToUndoStack(state);
                debugLog('ACTION', 'clearAllStructures');
            }
            state.structures = [];
            state.activeStructureId = null;
            state.selection = null;
            state.measurements = [];
        },

        // === Atom Operations ===
        deleteAtoms: (state, action: PayloadAction<{ structureId: string, atomIndices: number[] }>) => {
            const { structureId, atomIndices } = action.payload;
            const structure = state.structures.find(s => s.id === structureId);

            if (structure && atomIndices.length > 0) {
                pushToUndoStack(state);
                debugLog('ACTION', `deleteAtoms: ${atomIndices.length} atoms from ${structureId}`);

                const indicesSet = new Set(atomIndices);
                structure.data.atoms = structure.data.atoms.filter((_atom, index) => !indicesSet.has(index));
                structure.modifiedAt = Date.now();

                if (state.selection?.structureId === structureId) {
                    state.selection = null;
                }
            }
        },

        // === Undo/Redo ===
        undo: (state) => {
            if (state.undoStack.length === 0) {
                debugLog('UNDO', 'Nothing to undo');
                return;
            }

            // Save current state to redo stack
            const currentSnapshot = createSnapshot(state);
            state.redoStack.push(currentSnapshot);

            // Pop and restore from undo stack
            const previousSnapshot = state.undoStack.pop()!;
            restoreFromSnapshot(state, previousSnapshot);

            debugLog('UNDO', `Restored state. Undo: ${state.undoStack.length}, Redo: ${state.redoStack.length}`);
        },

        redo: (state) => {
            if (state.redoStack.length === 0) {
                debugLog('REDO', 'Nothing to redo');
                return;
            }

            // Save current state to undo stack
            const currentSnapshot = createSnapshot(state);
            state.undoStack.push(currentSnapshot);

            // Pop and restore from redo stack
            const nextSnapshot = state.redoStack.pop()!;
            restoreFromSnapshot(state, nextSnapshot);

            debugLog('REDO', `Restored state. Undo: ${state.undoStack.length}, Redo: ${state.redoStack.length}`);
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
    setLoadingState,
    clearAllStructures,
    deleteAtoms,
    undo,
    redo,
} = structureSlice.actions;

export default structureSlice.reducer;
