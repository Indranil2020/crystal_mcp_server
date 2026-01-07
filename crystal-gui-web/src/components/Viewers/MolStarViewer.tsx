/**
 * MolStarViewer - 3D Molecular Visualization with MolStar
 * 
 * Professional-grade 3D viewer for crystal structures and molecules.
 * Features: Multiple representations, unit cell, measurements, export.
 */

import { useEffect, useRef, useCallback, useState } from 'react';
import { useAppSelector, useAppDispatch } from '../../store/hooks';
import { updateViewerSettings } from '../../store/structureSlice';
import { structureToCif } from '../../converters';
import type { Structure, RepresentationMode } from '../../types';

// MolStar imports
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { DefaultPluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { PluginContext } from 'molstar/lib/mol-plugin/context';

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
        try {
            const canvas = document.createElement('canvas');
            const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
            return !!gl;
        } catch {
            return false;
        }
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

        async function initPlugin() {
            try {
                const plugin = await createPluginUI({
                    target: containerRef.current!,
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

                pluginRef.current = plugin;
                setIsInitialized(true);
                console.log('[MolStarViewer] Plugin initialized');
            } catch (err) {
                console.error('[MolStarViewer] Failed to initialize:', err);
                setError('Failed to initialize 3D viewer (WebGL error)');
            }
        }

        initPlugin();

        return () => {
            if (pluginRef.current) {
                pluginRef.current.dispose();
                pluginRef.current = null;
            }
        };
    }, []);

    // Load structure when active structure changes
    useEffect(() => {
        if (!isInitialized || !pluginRef.current || !activeStructure) return;

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

        try {
            // Clear existing structures
            await plugin.clear();

            // Convert to CIF format
            const cifData = structureToCif(structure.data);

            // Load CIF data
            const data = await plugin.builders.data.rawData({
                data: cifData,
                label: structure.name,
            });

            const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

            // Structure loaded

            console.log('[MolStarViewer] Loaded structure:', structure.name);
        } catch (err) {
            console.error('[MolStarViewer] Failed to load structure:', err);
            setError(`Failed to load: ${structure.name}`);
        }
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

        try {
            // Apply representation preset
            const structures = plugin.managers.structure.hierarchy.current.structures;
            for (const s of structures) {
                await plugin.builders.structure.representation.applyPreset(
                    s.cell.transform.ref,
                    presetMap[mode] || 'ball-and-stick'
                );
            }
        } catch (err) {
            console.error('[MolStarViewer] Failed to update representation:', err);
        }
    }, []);

    // Export screenshot
    const exportScreenshot = useCallback(async () => {
        const plugin = pluginRef.current;
        if (!plugin?.canvas3d) return;

        try {
            const imageData = await plugin.helpers.viewportScreenshot?.getImageDataUri();
            if (imageData) {
                const link = document.createElement('a');
                link.download = `${activeStructure?.name || 'structure'}.png`;
                link.href = imageData;
                link.click();
            }
        } catch (err) {
            console.error('[MolStarViewer] Screenshot failed:', err);
        }
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
