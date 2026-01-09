/**
 * useKeyboardShortcuts - Global keyboard shortcuts for Crystal GUI
 * 
 * Comprehensive keyboard shortcuts with debug logging at every action.
 * 
 * Shortcuts:
 *   Ctrl+Z / Cmd+Z         - Undo
 *   Ctrl+Shift+Z / Cmd+Shift+Z - Redo
 *   Ctrl+S / Cmd+S         - Export screenshot
 *   Delete / Backspace     - Delete selected atoms
 *   Escape                 - Clear selection
 *   1-5                    - Quick representation (1=ball-stick, 2=spacefill, etc.)
 *   C                      - Toggle unit cell
 *   M                      - Toggle measurement/selection mode
 *   R                      - Reset camera view
 */

import { useEffect, useCallback } from 'react';
import { useAppDispatch, useAppSelector } from '../store/hooks';
import {
    updateViewerSettings,
    setSelection,
    deleteAtoms,
    undo,
    redo,
} from '../store/structureSlice';
import type { RepresentationMode } from '../types';

// Debug logger
function debugLog(action: string, details?: string): void {
    const timestamp = new Date().toISOString().split('T')[1].slice(0, 12);
    console.log(`%c[KEYBOARD] ${timestamp} ${action}${details ? ': ' + details : ''}`, 'color: #9C27B0; font-weight: bold');
}

interface KeyboardShortcutsConfig {
    onExportScreenshot?: () => void;
    onResetCamera?: () => void;
    onToggleMeasurementMode?: () => void;
    enabled?: boolean;
}

export function useKeyboardShortcuts(config: KeyboardShortcutsConfig = {}): void {
    const dispatch = useAppDispatch();
    const { activeStructureId, selection, viewerSettings } = useAppSelector(state => state.structure);

    const {
        onExportScreenshot,
        onResetCamera,
        onToggleMeasurementMode,
        enabled = true,
    } = config;

    const handleKeyDown = useCallback((event: KeyboardEvent) => {
        // Skip if typing in an input field
        const target = event.target as HTMLElement;
        if (target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || target.isContentEditable) {
            debugLog('SKIP', 'Input field focused');
            return;
        }

        if (!enabled) {
            debugLog('SKIP', 'Shortcuts disabled');
            return;
        }

        const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
        const modKey = isMac ? event.metaKey : event.ctrlKey;
        const key = event.key.toLowerCase();

        debugLog('KEY_DOWN', `key="${event.key}" ctrl=${event.ctrlKey} meta=${event.metaKey} shift=${event.shiftKey}`);

        // === Ctrl/Cmd shortcuts ===
        if (modKey) {
            // Ctrl+Z / Cmd+Z - Undo
            if (key === 'z' && !event.shiftKey) {
                event.preventDefault();
                debugLog('ACTION', 'UNDO triggered');
                dispatch(undo());
                return;
            }

            // Ctrl+Shift+Z / Cmd+Shift+Z - Redo
            if (key === 'z' && event.shiftKey) {
                event.preventDefault();
                debugLog('ACTION', 'REDO triggered');
                dispatch(redo());
                return;
            }

            // Ctrl+Y / Cmd+Y - Also Redo (Windows convention)
            if (key === 'y') {
                event.preventDefault();
                debugLog('ACTION', 'REDO (Ctrl+Y) triggered');
                dispatch(redo());
                return;
            }

            // Ctrl+S / Cmd+S - Export screenshot
            if (key === 's') {
                event.preventDefault();
                debugLog('ACTION', 'EXPORT SCREENSHOT triggered');
                if (onExportScreenshot) {
                    onExportScreenshot();
                } else {
                    debugLog('WARN', 'No screenshot handler provided');
                }
                return;
            }

            // Don't intercept other Ctrl combinations
            return;
        }

        // === Single key shortcuts (no modifier) ===
        switch (key) {
            case 'delete':
            case 'backspace':
                // Delete selected atoms
                if (selection && selection.atomIndices.length > 0 && activeStructureId) {
                    event.preventDefault();
                    debugLog('ACTION', `DELETE ${selection.atomIndices.length} atoms from structure ${activeStructureId}`);
                    dispatch(deleteAtoms({
                        structureId: activeStructureId,
                        atomIndices: selection.atomIndices,
                    }));
                } else {
                    debugLog('SKIP', 'No atoms selected to delete');
                }
                break;

            case 'escape':
                // Clear selection
                debugLog('ACTION', 'CLEAR SELECTION');
                dispatch(setSelection(null));
                break;

            case '1':
                // Ball-and-stick representation
                debugLog('ACTION', 'Set representation: ball-and-stick');
                dispatch(updateViewerSettings({ representation: 'ball-and-stick' as RepresentationMode }));
                break;

            case '2':
                // Spacefill representation
                debugLog('ACTION', 'Set representation: spacefill');
                dispatch(updateViewerSettings({ representation: 'spacefill' as RepresentationMode }));
                break;

            case '3':
                // Licorice representation
                debugLog('ACTION', 'Set representation: licorice');
                dispatch(updateViewerSettings({ representation: 'licorice' as RepresentationMode }));
                break;

            case '4':
                // Wireframe representation
                debugLog('ACTION', 'Set representation: wireframe');
                dispatch(updateViewerSettings({ representation: 'wireframe' as RepresentationMode }));
                break;

            case '5':
                // Surface representation
                debugLog('ACTION', 'Set representation: surface');
                dispatch(updateViewerSettings({ representation: 'surface' as RepresentationMode }));
                break;

            case 'c':
                // Toggle unit cell
                debugLog('ACTION', `Toggle unit cell: ${!viewerSettings.showUnitCell}`);
                dispatch(updateViewerSettings({ showUnitCell: !viewerSettings.showUnitCell }));
                break;

            case 'm':
                // Toggle measurement/selection mode
                debugLog('ACTION', 'Toggle measurement mode');
                if (onToggleMeasurementMode) {
                    onToggleMeasurementMode();
                } else {
                    debugLog('WARN', 'No measurement mode handler provided');
                }
                break;

            case 'r':
                // Reset camera view
                debugLog('ACTION', 'RESET CAMERA VIEW');
                if (onResetCamera) {
                    onResetCamera();
                } else {
                    debugLog('WARN', 'No reset camera handler provided');
                }
                break;

            default:
                // Unknown key - no action
                break;
        }
    }, [dispatch, enabled, activeStructureId, selection, viewerSettings.showUnitCell, onExportScreenshot, onResetCamera, onToggleMeasurementMode]);

    // Register global keyboard listener
    useEffect(() => {
        debugLog('INIT', 'Registering global keyboard shortcuts');
        window.addEventListener('keydown', handleKeyDown);

        return () => {
            debugLog('CLEANUP', 'Removing global keyboard shortcuts');
            window.removeEventListener('keydown', handleKeyDown);
        };
    }, [handleKeyDown]);
}

export default useKeyboardShortcuts;
