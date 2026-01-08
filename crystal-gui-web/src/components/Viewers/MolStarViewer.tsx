/**
 * MolStarViewer - 3D Molecular Visualization with MolStar
 *
 * Professional-grade 3D viewer for crystal structures and molecules.
 * Features: Multiple representations, unit cell, measurements, export.
 *
 * DEBUG: Comprehensive logging enabled at every step of structure loading.
 */

import { useEffect, useRef, useCallback, useState } from 'react';
import { useAppSelector, useAppDispatch } from '../../store/hooks';
import { updateViewerSettings, setSelection, deleteAtoms } from '../../store/structureSlice';
import { structureToCif } from '../../converters';
import { debug } from '../../debug';
import type { Structure, RepresentationMode } from '../../types';

// MolStar imports
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import {
    StructureElement,
    StructureProperties,
    StructureSelection,
    QueryContext
} from 'molstar/lib/mol-model/structure';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { compile } from 'molstar/lib/mol-script/runtime/query/compiler';

// MolStar CSS
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

// Debug helper to log Mol* state
function logPluginState(plugin: PluginContext, label: string): void {
    debug('VIEWERS', `[${label}] Mol* plugin state:`);
    // State tree may have complex structure, just log structure count
    debug('VIEWERS', `  - Structures loaded: ${plugin.managers.structure.hierarchy.current.structures.length}`);
    debug('VIEWERS', `  - Canvas3D initialized: ${!!plugin.canvas3d}`);
}

interface Props {
    className?: string;
}

export default function MolStarViewer({ className = '' }: Props) {
    const dispatch = useAppDispatch();
    const { structures, activeStructureId, viewerSettings } = useAppSelector(state => state.structure);
    const activeStructure = structures.find(s => s.id === activeStructureId);

    const containerRef = useRef<HTMLDivElement>(null);
    const pluginRef = useRef<PluginContext | null>(null);
    const isInitializingRef = useRef(false); // Track async initialization in progress
    const [isInitialized, setIsInitialized] = useState(false);
    const [error, setError] = useState<string | null>(null);

    // Check WebGL support
    const checkWebGLSupport = (): boolean => {
        const canvas = document.createElement('canvas');
        const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
        return !!gl;
    };

    // Initialize MolStar plugin
    useEffect(() => {
        console.log('[MolStarViewer] Component MOUNTED');
        console.log('[MolStarViewer] Current state - pluginRef:', !!pluginRef.current, 'isInitializing:', isInitializingRef.current, 'container:', !!containerRef.current);

        // Guard: Don't initialize if already initialized or currently initializing
        if (pluginRef.current) {
            console.log('[MolStarViewer] Plugin already exists, skipping initialization');
            return;
        }

        if (isInitializingRef.current) {
            console.log('[MolStarViewer] Initialization already in progress, skipping');
            return;
        }

        if (!containerRef.current) {
            console.log('[MolStarViewer] No container ref, skipping initialization');
            return;
        }

        // Check WebGL support first
        if (!checkWebGLSupport()) {
            console.warn('[MolStarViewer] WebGL not supported, skipping initialization');
            setError('WebGL not available (headless browser or GPU disabled)');
            return;
        }

        let isMounted = true;
        let plugin: PluginContext | null = null;

        async function initPlugin() {
            const container = containerRef.current;
            if (!container || !isMounted) {
                console.log('[MolStarViewer] Container or mount check failed before init');
                return;
            }

            console.log('[MolStarViewer] ========================================');
            console.log('[MolStarViewer] STARTING plugin initialization');
            console.log('[MolStarViewer] Container element:', container.tagName, container.className);
            console.log('[MolStarViewer] Container children count:', container.childElementCount);

            isInitializingRef.current = true;

            try {
                console.log('[MolStarViewer] Calling createPluginUI...');
                plugin = await createPluginUI({
                    target: container,
                    render: renderReact18,
                    spec: {
                        ...DefaultPluginUISpec(),
                        layout: {
                            initial: {
                                isExpanded: false,
                                showControls: false,
                                controlsDisplay: 'reactive',
                            },
                        },
                    },
                });

                console.log('[MolStarViewer] createPluginUI returned:', !!plugin);

                if (!isMounted || !plugin) {
                    console.log('[MolStarViewer] Component unmounted or plugin null, disposing');
                    plugin?.dispose();
                    isInitializingRef.current = false;
                    return;
                }

                pluginRef.current = plugin;
                setIsInitialized(true);
                setError(null);
                console.log('[MolStarViewer] ========================================');
                console.log('[MolStarViewer] Plugin initialized SUCCESSFULLY');
                console.log('[MolStarViewer] Container children count after init:', container.childElementCount);
                console.log('[MolStarViewer] ========================================');
            } catch (err) {
                console.error('[MolStarViewer] ========================================');
                console.error('[MolStarViewer] FAILED to initialize plugin:', err);
                console.error('[MolStarViewer] Error type:', err instanceof Error ? err.constructor.name : typeof err);
                console.error('[MolStarViewer] ========================================');
                if (isMounted) {
                    setError('Failed to initialize 3D viewer');
                }
            } finally {
                isInitializingRef.current = false;
                console.log('[MolStarViewer] Initialization process complete, isInitializing set to false');
            }
        }

        initPlugin();

        return () => {
            console.log('[MolStarViewer] ========================================');
            console.log('[MolStarViewer] Component UNMOUNTING - cleanup starting');
            console.log('[MolStarViewer] Plugin exists:', !!pluginRef.current);

            isMounted = false;

            if (pluginRef.current) {
                console.log('[MolStarViewer] Disposing plugin...');
                try {
                    pluginRef.current.dispose();
                    console.log('[MolStarViewer] Plugin disposed successfully');
                } catch (err) {
                    console.error('[MolStarViewer] Error during plugin disposal:', err);
                }
                pluginRef.current = null;
            }

            // CRITICAL: Use innerHTML to clear, NOT removeChild
            // This prevents React from trying to remove children it doesn't own
            if (containerRef.current) {
                console.log('[MolStarViewer] Clearing container via innerHTML...');
                try {
                    containerRef.current.innerHTML = '';
                    console.log('[MolStarViewer] Container cleared successfully');
                } catch (err) {
                    console.error('[MolStarViewer] Error clearing container:', err);
                }
            }

            isInitializingRef.current = false;
            console.log('[MolStarViewer] Cleanup complete');
            console.log('[MolStarViewer] ========================================');
        };
    }, []);

    // Load structure when active structure changes
    useEffect(() => {
        if (!isInitialized || !pluginRef.current || !activeStructure) return;

        debug('VIEWERS', '🎯 Active structure changed, loading...');
        loadStructure(activeStructure);
    }, [isInitialized, activeStructure]);

    // Update representation when settings change
    useEffect(() => {
        if (!isInitialized || !pluginRef.current) return;

        updateRepresentation(viewerSettings.representation);
    }, [isInitialized, viewerSettings.representation]);

    // Load structure into MolStar
    const loadStructure = useCallback(async (structure: Structure) => {
        const plugin = pluginRef.current;
        if (!plugin) {
            debug('VIEWERS', '[ERROR] Cannot load structure - plugin not initialized');
            return;
        }

        debug('VIEWERS', '═'.repeat(60));
        debug('VIEWERS', '[MOLSTAR] LOADING STRUCTURE INTO MOL*');
        debug('VIEWERS', `  Structure ID: ${structure.id}`);
        debug('VIEWERS', `  Name: ${structure.name}`);
        debug('VIEWERS', `  Atom count: ${structure.data.atoms.length}`);
        debug('VIEWERS', `  Source: ${structure.source}`);
        debug('VIEWERS', `  Lattice: a=${structure.data.lattice.a.toFixed(2)}, b=${structure.data.lattice.b.toFixed(2)}, c=${structure.data.lattice.c.toFixed(2)}`);
        debug('VIEWERS', `  Formula: ${structure.data.metadata?.formula || 'unknown'}`);

        // Log first few atoms
        structure.data.atoms.slice(0, 3).forEach((atom, i) => {
            debug('VIEWERS', `  Atom[${i}]: ${atom.element} coords=${JSON.stringify(atom.coords)} cartesian=${JSON.stringify(atom.cartesian)}`);
        });
        if (structure.data.atoms.length > 3) {
            debug('VIEWERS', `  ... and ${structure.data.atoms.length - 3} more atoms`);
        }

        logPluginState(plugin, 'Before clear');

        // Clear existing structures
        debug('VIEWERS', '  Clearing existing structures...');
        await plugin.clear();
        logPluginState(plugin, 'After clear');

        // Convert to CIF format
        debug('VIEWERS', '  [CIF] Converting to CIF format...');
        const cifData = structureToCif(structure.data);
        debug('VIEWERS', `  [CIF] Generated: ${cifData.length} bytes`);

        // Try loading with CIF parser (not mmcif - our format is standard CIF)
        debug('VIEWERS', '  [LOAD] Loading CIF into Mol*...');
        debug('VIEWERS', '  [LOAD] Step 1: Loading raw data...');

        const dataResult = await plugin.builders.data.rawData({
            data: cifData,
            label: structure.name,
        });

        if (!dataResult.ref) {
            debug('VIEWERS', '  [ERROR] rawData returned no ref');
            setError('Failed to load raw data');
            return;
        }
        debug('VIEWERS', `  [OK] Raw data loaded, ref=${dataResult.ref}`);
        logPluginState(plugin, 'After rawData');

        // Parse as CIF trajectory
        // Mol* supports: mmcif, cifCore, pdb, xyz, mol, sdf, mol2, etc.
        // Use 'cifCore' for standard crystallographic CIF files
        debug('VIEWERS', '  [LOAD] Step 2: Parsing trajectory (format: cifCore)...');

        let trajectoryResult = await plugin.builders.structure.parseTrajectory(dataResult, 'cifCore');

        if (!trajectoryResult.ref) {
            debug('VIEWERS', '  [WARN] cifCore format failed, trying mmcif...');

            // Fallback: try mmcif format (for macromolecular-style CIF)
            trajectoryResult = await plugin.builders.structure.parseTrajectory(dataResult, 'mmcif');
            if (!trajectoryResult.ref) {
                debug('VIEWERS', '  [ERROR] mmcif format also failed');
                setError('Failed to parse CIF data');
                return;
            }
            debug('VIEWERS', `  [OK] mmcif fallback succeeded, ref=${trajectoryResult.ref}`);
        } else {
            debug('VIEWERS', `  [OK] Trajectory parsed with cifCore, ref=${trajectoryResult.ref}`);
        }

        logPluginState(plugin, 'After parseTrajectory');

        // Apply default preset for visualization
        debug('VIEWERS', '  [LOAD] Step 3: Applying visualization preset...');
        const structures = plugin.managers.structure.hierarchy.current.structures;
        debug('VIEWERS', `  Structures in hierarchy: ${structures.length}`);

        if (structures.length === 0) {
            debug('VIEWERS', '  [WARN] No structures in hierarchy after parsing');
            debug('VIEWERS', '  [LOAD] Attempting to apply preset to trajectory...');

            // Apply preset to the trajectory we already parsed
            if (trajectoryResult.ref) {
                await plugin.builders.structure.hierarchy.applyPreset(trajectoryResult, 'default');
                logPluginState(plugin, 'After applyPreset');
            }
        }

        // Check final state
        const finalStructures = plugin.managers.structure.hierarchy.current.structures;
        debug('VIEWERS', `  Final structures count: ${finalStructures.length}`);

        if (finalStructures.length > 0) {
            debug('VIEWERS', '  [SUCCESS] STRUCTURE LOADED SUCCESSFULLY');
            debug('VIEWERS', `  [SUCCESS] Structure "${structure.name}" is now visible`);
            setError(null);

            // Reset camera to show the structure
            debug('VIEWERS', '  Resetting camera view...');
            plugin.managers.camera.reset();
        } else {
            debug('VIEWERS', '  [ERROR] Structure loading failed - no structures in viewer');
            setError(`Failed to render: ${structure.name}`);
        }

        debug('VIEWERS', '═'.repeat(60));
    }, []);

    // Update visual representation
    const updateRepresentation = useCallback(async (mode: RepresentationMode) => {
        const plugin = pluginRef.current;
        if (!plugin) return;

        // Map our representation modes to MolStar presets
        const presetMap: Record<RepresentationMode, string> = {
            'ball-and-stick': 'ball-and-stick',
            'spacefill': 'spacefill',
            'cartoon': 'cartoon',
            'licorice': 'ball-and-stick', // Similar
            'surface': 'molecular-surface',
            'wireframe': 'ball-and-stick',
        };

        const structures = plugin.managers.structure.hierarchy.current.structures;
        for (const s of structures) {
            await plugin.builders.structure.representation.applyPreset(
                s.cell.transform.ref,
                presetMap[mode] || 'ball-and-stick'
            );
        }
    }, []);

    // Handle selection changes (Redux -> MolStar)
    const { selection } = useAppSelector(state => state.structure);

    useEffect(() => {
        if (!pluginRef.current || !activeStructure) return;

        const plugin = pluginRef.current;
        const manager = plugin.managers.structure.selection;

        if (!selection || selection.structureId !== activeStructure.id || selection.atomIndices.length === 0) {
            manager.clear();
            return;
        }

        // Apply selection from Redux
        const structureRef = plugin.managers.structure.hierarchy.current.structures[0];

        // Explicitly check for structure reference
        if (!structureRef) return;

        // Select by atom indices (using the IDs we put in the CIF: 1-based index)
        const ids = selection.atomIndices.map(i => i + 1); // 0-based to 1-based

        // Explicitly check for structure data
        const data = structureRef.cell.obj?.data;
        if (!data) return;

        // Construct query using MolScript Builder
        const expression = MS.struct.generator.atomGroups({
            'atom-test': MS.core.set.has([
                MS.set(...ids),
                MS.ammp('id') // Matches _atom_site.id property
            ])
        });

        // Compile and run query
        // We assume compilation is safe given valid expression construction
        const query = compile(expression);
        const result = query(new QueryContext(data));

        // Convert to Loci and apply
        const loci = StructureSelection.toLociWithSourceUnits(result);

        manager.fromLoci('set', loci);

    }, [selection, activeStructure]);

    // Initialize selection behavior (MolStar -> Redux)
    useEffect(() => {
        if (!isInitialized || !pluginRef.current) return;
        const plugin = pluginRef.current;

        const sub = plugin.behaviors.interaction.click.subscribe(({ current }) => {
            if (activeStructureId && current.loci.kind === 'element-loci') {
                const loci = current.loci;

                // Extract atom ID from Loci
                const loc = StructureElement.Location.create(loci.structure);
                let atomIndex = -1;

                // Get the first selected atom
                if (StructureElement.Loci.getFirstLocation(loci, loc)) {
                    // Try to get the ID we embedded in CIF (_atom_site_id)
                    // Usually mapped to StructureProperties.atom.id OR .auth_seq_id
                    // Default mmCIF parser might map _atom_site_id to atom.id

                    // Direct access (no try/catch)
                    const id = StructureProperties.atom.id(loc);
                    atomIndex = id - 1; // Convert back to 0-based
                }

                if (atomIndex >= 0) {
                    dispatch(setSelection({
                        structureId: activeStructureId,
                        atomIndices: [atomIndex],
                        selectionMode: 'single'
                    }));
                }
            } else {
                // Clicked empty space
                dispatch(setSelection(null));
            }
        });

        return () => sub.unsubscribe();
    }, [isInitialized, activeStructureId, dispatch]);

    // Export screenshot
    const exportScreenshot = useCallback(async () => {
        const plugin = pluginRef.current;
        if (!plugin?.canvas3d) return;

        await plugin.helpers.viewportScreenshot?.getImageDataUri()
            .then(imageData => {
                if (imageData) {
                    const link = document.createElement('a');
                    link.download = `${activeStructure?.name || 'structure'}.png`;
                    link.href = imageData;
                    link.click();
                }
            })
            .catch(err => {
                console.error('[MolStarViewer] Screenshot failed:', err);
            });
    }, [activeStructure]);

    // Representation selector
    const handleRepresentationChange = (mode: RepresentationMode) => {
        dispatch(updateViewerSettings({ representation: mode }));
    };

    if (error) {
        return (
            <div className={`flex items-center justify-center h-full bg-slate-900 ${className}`}>
                <div className="text-center text-red-400">
                    <p className="text-lg">[ERROR] Viewer Error</p>
                    <p className="text-sm mt-2">{error}</p>
                </div>
            </div>
        );
    }

    return (
        <div className={`flex flex-col h-full ${className}`}>
            {/* Viewer Controls */}
            <div className="px-3 py-2 bg-slate-800 border-b border-slate-700 flex items-center justify-between">
                <div className="flex items-center gap-2">
                    <span className="text-sm text-slate-300">[3D Viewer]</span>
                    {activeStructure && (
                        <span className="text-xs text-slate-500">
                            {activeStructure.name} ({activeStructure.data.atoms.length} atoms)
                        </span>
                    )}
                </div>

                <div className="flex items-center gap-2">
                    {/* Delete selection */}
                    <button
                        onClick={() => selection && activeStructureId && dispatch(deleteAtoms({
                            structureId: activeStructureId,
                            atomIndices: selection.atomIndices
                        }))}
                        disabled={!selection || selection.atomIndices.length === 0}
                        className="text-xs px-2 py-1 rounded bg-red-900/50 text-red-200 hover:bg-red-800 disabled:opacity-50 disabled:bg-slate-700 disabled:text-slate-500"
                        title="Delete selected atoms (Del)"
                    >
                        DEL
                    </button>

                    <div className="w-px h-4 bg-slate-600 mx-1" />

                    {/* Representation dropdown */}
                    <select
                        value={viewerSettings.representation}
                        onChange={e => handleRepresentationChange(e.target.value as RepresentationMode)}
                        className="text-xs bg-slate-700 text-slate-200 rounded px-2 py-1 border border-slate-600"
                    >
                        <option value="ball-and-stick">Ball & Stick</option>
                        <option value="spacefill">Spacefill</option>
                        <option value="licorice">Licorice</option>
                        <option value="wireframe">Wireframe</option>
                        <option value="surface">Surface</option>
                    </select>

                    {/* Unit cell toggle */}
                    <button
                        onClick={() => dispatch(updateViewerSettings({ showUnitCell: !viewerSettings.showUnitCell }))}
                        className={`text-xs px-2 py-1 rounded ${viewerSettings.showUnitCell
                            ? 'bg-blue-600 text-white'
                            : 'bg-slate-700 text-slate-300'
                            }`}
                    >
                        Cell
                    </button>

                    {/* Screenshot */}
                    <button
                        onClick={exportScreenshot}
                        disabled={!activeStructure}
                        className="text-xs px-2 py-1 rounded bg-slate-700 text-slate-300 hover:bg-slate-600 disabled:opacity-50"
                    >
                        IMG
                    </button>
                </div>
            </div>

            {/* MolStar Container - completely managed by Mol*, React must not touch its children */}
            <div className="flex-1 relative">
                <div
                    ref={containerRef}
                    className="absolute inset-0"
                    style={{ minHeight: '300px' }}
                />

                {/* Overlay messages - separate from Mol* container */}
                {!isInitialized && !error && (
                    <div className="absolute inset-0 flex items-center justify-center bg-slate-900 pointer-events-none">
                        <div className="text-center text-slate-400">
                            <div className="w-8 h-8 border-2 border-blue-500 border-t-transparent rounded-full animate-spin mx-auto mb-2" />
                            <p className="text-sm">Initializing MolStar...</p>
                        </div>
                    </div>
                )}

                {isInitialized && !activeStructure && (
                    <div className="absolute inset-0 flex items-center justify-center pointer-events-none">
                        <div className="text-center text-slate-500">
                            <p className="text-lg">No structure loaded</p>
                            <p className="text-sm mt-2">Generate a molecule using the chat</p>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
}
