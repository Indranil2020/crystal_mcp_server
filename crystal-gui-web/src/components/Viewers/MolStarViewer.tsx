/**
 * MolStarViewer - 3D Molecular Visualization with MolStar
 * 
 * Professional-grade 3D viewer for crystal structures and molecules.
 * Features: Multiple representations, unit cell, measurements, export.
 */

import { useEffect, useRef, useCallback, useState } from 'react';
import { useAppSelector, useAppDispatch } from '../../store/hooks';
import { updateViewerSettings, setSelection, deleteAtoms } from '../../store/structureSlice';
import { structureToCif } from '../../converters';
import { debug, debugError } from '../../debug';
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

interface Props {
    className?: string;
}

export default function MolStarViewer({ className = '' }: Props) {
    const dispatch = useAppDispatch();
    const { structures, activeStructureId, viewerSettings } = useAppSelector(state => state.structure);
    const activeStructure = structures.find(s => s.id === activeStructureId);

    const containerRef = useRef<HTMLDivElement>(null);
    const pluginRef = useRef<PluginContext | null>(null);
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
        if (!containerRef.current || pluginRef.current) return;

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
            if (!container || !isMounted) return;

            console.log('[MolStarViewer] Initializing plugin...');

            try {
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

                if (!isMounted || !plugin) {
                    plugin?.dispose();
                    return;
                }

                pluginRef.current = plugin;
                setIsInitialized(true);
                setError(null);
                console.log('[MolStarViewer] Plugin initialized successfully');
            } catch (err) {
                console.error('[MolStarViewer] Failed to initialize:', err);
                if (isMounted) {
                    setError('Failed to initialize 3D viewer');
                }
            }
        }

        initPlugin();

        return () => {
            isMounted = false;
            if (pluginRef.current) {
                pluginRef.current.dispose();
                pluginRef.current = null;
            }
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
        if (!plugin) return;

        debug('VIEWERS', '='.repeat(60));
        debug('VIEWERS', '🔬 LOADING STRUCTURE INTO MOL*');
        debug('VIEWERS', `Structure ID: ${structure.id}`);
        debug('VIEWERS', `Structure name: ${structure.name}`);
        debug('VIEWERS', `Atom count: ${structure.data.atoms.length}`);
        debug('VIEWERS', `Source: ${structure.source}`);

        // Log full structure data for debugging
        debug('VIEWERS', 'Full structure data:');
        debug('VIEWERS', `  Lattice: a=${structure.data.lattice.a}, b=${structure.data.lattice.b}, c=${structure.data.lattice.c}`);
        debug('VIEWERS', `  Formula: ${structure.data.metadata?.formula}`);
        debug('VIEWERS', `  First atom: ${JSON.stringify(structure.data.atoms[0])}`);

        // Clear existing structures
        debug('VIEWERS', 'Clearing existing structures...');
        await plugin.clear();

        // Convert to CIF format
        debug('VIEWERS', '📝 Converting structure to CIF format...');
        try {
            const cifData = structureToCif(structure.data);
            debug('VIEWERS', `✅ CIF conversion complete: ${cifData.length} bytes`);
            debug('VIEWERS', `CIF preview (first 300 chars): ${cifData.substring(0, 300)}`);

            // Load CIF data into Mol*
            debug('VIEWERS', '⏳ Loading CIF into Mol* plugin...');
            await plugin.builders.data.rawData({
                data: cifData,
                label: structure.name,
            })
                .then(data => {
                    debug('VIEWERS', '✅ Raw data loaded, parsing trajectory...');
                    return plugin.builders.structure.parseTrajectory(data, 'mmcif');
                })
                .then(trajectory => {
                    debug('VIEWERS', '✅ Trajectory parsed, applying preset...');
                    return plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
                })
                .then(() => {
                    debug('VIEWERS', '✅✅✅ STRUCTURE LOADED SUCCESSFULLY INTO MOL*');
                    debug('VIEWERS', `Structure "${structure.name}" is now visible in viewer`);
                    setError(null); // Clear any previous errors
                })
                .catch(err => {
                    debugError('VIEWERS', err, '❌ Failed to load structure into Mol*');
                    debug('VIEWERS', 'Error details:', err);
                    debug('VIEWERS', 'Error stack:', err.stack);
                    setError(`Failed to load: ${structure.name}`);
                });
        } catch (err) {
            debugError('VIEWERS', err, '❌ CIF conversion failed');
            setError(`Failed to convert structure: ${structure.name}`);
        }

        debug('VIEWERS', '='.repeat(60));
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
                    <p className="text-lg">⚠️ Viewer Error</p>
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
                    <span className="text-sm text-slate-300">🔬 3D Viewer</span>
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
                        🗑️
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
                        📷
                    </button>
                </div>
            </div>

            {/* MolStar Container */}
            <div
                ref={containerRef}
                className="flex-1 relative"
                style={{ minHeight: '300px' }}
                suppressHydrationWarning
            >
                {!isInitialized && !error && (
                    <div className="absolute inset-0 flex items-center justify-center bg-slate-900">
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
