/**
 * KekuleEditor - 2D Molecular Editor with Kekule.js
 *
 * Professional-grade 2D structure editor for drawing and editing molecules.
 * Features: Drawing tools, templates, SMILES I/O, sync with 3D viewer.
 *
 * DEBUG: Comprehensive logging for initialization and operation tracking.
 */

import { useEffect, useRef, useState, useCallback } from 'react';
import { useAppDispatch } from '../../store/hooks';
import { addStructure } from '../../store/structureSlice';
import { v4 as uuidv4 } from 'uuid';
import { debug, debugError } from '../../debug';
import type { StructureData } from '../../types';

// Kekule.js imports - using global loaded from CDN for better compatibility
declare global {
    interface Window {
        Kekule: unknown;
    }
}

interface Props {
    className?: string;
    onStructureChange?: (smiles: string) => void;
}

export default function KekuleEditor({ className = '', onStructureChange }: Props) {
    const dispatch = useAppDispatch();
    const editorRef = useRef<HTMLDivElement>(null);
    const composerRef = useRef<any>(null);
    const [isLoaded, setIsLoaded] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [currentSmiles, setCurrentSmiles] = useState<string>('');

    // Load Kekule.js from CDN
    useEffect(() => {
        debug('VIEWERS', '═'.repeat(50));
        debug('VIEWERS', '✏️ KEKULE EDITOR INITIALIZATION START');

        let isMounted = true;

        if (window.Kekule) {
            debug('VIEWERS', '  ℹ️ Kekule.js global object already exists on window');
            // Even if it exists, checking if it's fully ready might be safer, but usually it implies loaded
            // We verify the subsystems we need
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const K = window.Kekule as any;
            debug('VIEWERS', `  🔍 Global State Check: Widget=${!!K.Widget}, Editor=${!!K.Editor}, ChemWidget=${!!(K.Editor && K.Editor.Composer)}`);

            if (K.Editor && K.Editor.Composer) {
                debug('VIEWERS', '  ✓ Kekule seems fully loaded, initializing editor immediately');
                initEditor();
            } else {
                debugError('VIEWERS', '  ⚠️ Kekule global exists but Scheduler/Composer is missing. This might be a partial load.', 'KekuleEditor');
                // In a partial load scenario, we might want to try reloading or just fail gracefully. 
                // For now, let's attempt to proceed to let the specific checks fail if they must.
                initEditor();
            }
            return;
        }

        debug('VIEWERS', '  ⬇️ Loading Kekule.js from CDN (Not found on window)...');

        // Load CSS
        const css = document.createElement('link');
        css.rel = 'stylesheet';
        css.href = 'https://cdn.jsdelivr.net/npm/kekule/dist/themes/default/kekule.css';
        document.head.appendChild(css);
        debug('VIEWERS', '  ✓ CSS tag appended');

        // Load main script
        const script = document.createElement('script');
        script.src = 'https://cdn.jsdelivr.net/npm/kekule/dist/kekule.min.js';

        script.onload = () => {
            if (!isMounted) {
                debug('VIEWERS', '  ⚠️ Script loaded but component unmounted, skipping init');
                return;
            }
            debug('VIEWERS', '  ✓ Kekule.js script onload triggered');
            debug('VIEWERS', '  🕒 Waiting 100ms for script parse/execution just to be safe...');
            // Small breathing room for script execution if needed, though onload usually suffices
            setTimeout(() => {
                if (isMounted) initEditor();
            }, 100);
        };

        script.onerror = (e) => {
            if (!isMounted) return;
            debugError('VIEWERS', `  ❌ Failed to load Kekule.js script: ${e}`, 'KekuleEditor');
            setError('Failed to load Kekule.js from CDN. Check internet connection.');
        };

        document.head.appendChild(script);
        debug('VIEWERS', '  ✓ Script tag appended, waiting for browser to load...');

        return () => {
            isMounted = false;
        };
    }, []);

    // Initialize the Kekule Composer
    const initEditor = useCallback(() => {
        debug('VIEWERS', '  🛠️ initEditor() called');

        if (!editorRef.current) {
            debug('VIEWERS', '  ❌ DOM Error: editorRef.current is null. Cannot attach Composer.');
            setError('Editor DOM element not found');
            return;
        }

        if (!window.Kekule) {
            debugError('VIEWERS', '  ❌ Critical: initEditor called but window.Kekule is undefined', 'KekuleEditor');
            setError('Kekule.js library failed to initialize');
            return;
        }

        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const Kekule = window.Kekule as any;

        debug('VIEWERS', '  🔍 Verifying Kekule Subsystems...');
        const hasWidget = !!Kekule.Widget;
        const hasComposer = !!(Kekule.Editor && Kekule.Editor.Composer);

        debug('VIEWERS', `    - Kekule.Widget: ${hasWidget ? 'OK' : 'MISSING'}`);
        debug('VIEWERS', `    - Kekule.Editor.Composer: ${hasComposer ? 'OK' : 'MISSING'}`);

        if (!hasWidget || !hasComposer) {
            debugError('VIEWERS', '  ❌ Required Kekule subsystems are missing.', 'KekuleEditor');
            setError('Kekule.js loaded but required components (Widget/Composer) are missing.');
            return;
        }

        try {
            debug('VIEWERS', '  🚀 Creating Kekule.Editor.Composer instance...');
            const composer = new Kekule.Editor.Composer(editorRef.current);
            debug('VIEWERS', '  ✓ Composer instance created');

            composer.setDimension('100%', '100%');
            debug('VIEWERS', '  ✓ Dimensions set to 100%');

            // Configure toolbar
            debug('VIEWERS', '  ⚙️ Configuring Toolbars...');
            if (composer.setCommonToolButtons) {
                composer.setCommonToolButtons([
                    'newDoc', 'loadData', 'saveData',
                    'undo', 'redo',
                    'copy', 'cut', 'paste',
                    'zoomIn', 'zoomOut', 'reset',
                ]);
            } else {
                debug('VIEWERS', '  ⚠️ setCommonToolButtons missing from composer instance');
            }

            if (composer.setChemToolButtons) {
                composer.setChemToolButtons([
                    'manipulate', 'erase', 'bond', 'atom',
                    'ring', 'charge', 'glyph',
                ]);
            } else {
                debug('VIEWERS', '  ⚠️ setChemToolButtons missing from composer instance');
            }

            // Listen for changes
            composer.addEventListener('editObjsUpdated', () => {
                const mol = composer.getChemObj();
                if (mol && Kekule.IO && Kekule.IO.saveFormatData) {
                    const smiles = Kekule.IO.saveFormatData(mol, 'smi');
                    // debug('VIEWERS', `  Molecule updated: SMILES=${smiles}`); // Commented to reduce noise on every edit
                    setCurrentSmiles(smiles);
                    onStructureChange?.(smiles);
                }
            });

            composerRef.current = composer;
            setIsLoaded(true); // This removes the loading spinner
            setError(null);    // Clear any previous errors

            debug('VIEWERS', '  ✅ KEKULE EDITOR FULLY INITIALIZED AND READY');
            debug('VIEWERS', '═'.repeat(50));

        } catch (e) {
            debugError('VIEWERS', `  ❌ Exception during Composer instantiation: ${e}`, 'KekuleEditor');
            setError(`Failed to create editor: ${e instanceof Error ? e.message : String(e)}`);
        }
    }, [onStructureChange]);

    // Load molecule from SMILES
    const loadFromSmiles = useCallback((smiles: string) => {
        if (!composerRef.current || !window.Kekule) return;
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const Kekule = window.Kekule as any;

        if (Kekule.IO && Kekule.IO.loadFormatData) {
            const mol = Kekule.IO.loadFormatData(smiles, 'smi');
            if (mol) {
                composerRef.current.setChemObj(mol);
            }
        }
    }, []);

    // Push structure to 3D viewer
    const pushTo3D = useCallback(() => {
        if (!composerRef.current || !window.Kekule) return;

        const mol = composerRef.current.getChemObj();

        if (!mol) {
            console.warn('[KekuleEditor] No molecule to push');
            return;
        }

        // Get coordinates from Kekule molecule
        const atoms: StructureData['atoms'] = [];
        const nodeCount = mol.getNodeCount?.() || 0;

        for (let i = 0; i < nodeCount; i++) {
            const node = mol.getNodeAt(i);
            const coord = node.getCoord2D?.() || { x: 0, y: 0 };
            const element = node.getSymbol?.() || 'C';

            atoms.push({
                element,
                coords: [coord.x, coord.y, 0],
                cartesian: [coord.x, coord.y, 0],
            });
        }

        if (atoms.length === 0) {
            console.warn('[KekuleEditor] No atoms in molecule');
            return;
        }

        // Create structure for Redux
        const smiles = currentSmiles || 'drawn-molecule';
        const structure = {
            id: uuidv4(),
            name: `Drawn: ${smiles.slice(0, 20)}`,
            data: {
                lattice: { a: 20, b: 20, c: 20, alpha: 90, beta: 90, gamma: 90, matrix: [[20, 0, 0], [0, 20, 0], [0, 0, 20]] },
                atoms,
                metadata: { formula: smiles, natoms: atoms.length },
            },
            source: 'kekule' as const,
            createdAt: Date.now(),
            modifiedAt: Date.now(),
            visible: true,
        };

        dispatch(addStructure(structure));
        console.log('[KekuleEditor] Pushed to 3D:', structure.name);
    }, [dispatch, currentSmiles]);

    // Clear editor
    const clearEditor = useCallback(() => {
        if (composerRef.current) {
            composerRef.current.newDoc();
            setCurrentSmiles('');
        }
    }, []);

    // Template molecules
    const loadTemplate = useCallback((template: string) => {
        const templates: Record<string, string> = {
            benzene: 'c1ccccc1',
            cyclohexane: 'C1CCCCC1',
            pyridine: 'c1ccncc1',
            naphthalene: 'c1ccc2ccccc2c1',
            phenol: 'c1ccc(O)cc1',
            aniline: 'c1ccc(N)cc1',
        };

        if (templates[template]) {
            loadFromSmiles(templates[template]);
        }
    }, [loadFromSmiles]);

    if (error) {
        return (
            <div className={`flex items-center justify-center h-full bg-slate-900 ${className}`}>
                <div className="text-center text-amber-400 p-4">
                    <p className="text-lg font-bold">⚠️ 2D Editor Unavailable</p>
                    <p className="text-sm mt-2 text-slate-300">{error}</p>
                    <p className="text-xs mt-4 text-slate-500 font-mono">Check console for 'VIEWERS' logs</p>
                </div>
            </div>
        );
    }

    return (
        <div className={`flex flex-col h-full ${className}`}>
            {/* Editor Controls */}
            <div className="px-3 py-2 bg-slate-800 border-b border-slate-700 flex items-center justify-between">
                <div className="flex items-center gap-2">
                    <span className="text-sm text-slate-300">✏️ 2D Editor</span>
                    {currentSmiles && (
                        <span className="text-xs text-slate-500 truncate max-w-32">
                            {currentSmiles}
                        </span>
                    )}
                </div>

                <div className="flex items-center gap-1">
                    {/* Templates */}
                    <select
                        onChange={e => e.target.value && loadTemplate(e.target.value)}
                        className="text-xs bg-slate-700 text-slate-200 rounded px-2 py-1 border border-slate-600"
                        defaultValue=""
                    >
                        <option value="">Templates</option>
                        <option value="benzene">Benzene</option>
                        <option value="cyclohexane">Cyclohexane</option>
                        <option value="pyridine">Pyridine</option>
                        <option value="naphthalene">Naphthalene</option>
                        <option value="phenol">Phenol</option>
                        <option value="aniline">Aniline</option>
                    </select>

                    {/* Clear */}
                    <button
                        onClick={clearEditor}
                        className="text-xs px-2 py-1 rounded bg-slate-700 text-slate-300 hover:bg-slate-600"
                    >
                        Clear
                    </button>

                    {/* Push to 3D */}
                    <button
                        onClick={pushTo3D}
                        className="text-xs px-2 py-1 rounded bg-blue-600 text-white hover:bg-blue-500"
                    >
                        → 3D
                    </button>
                </div>
            </div>

            {/* Kekule Editor Container */}
            <div
                ref={editorRef}
                className="flex-1 bg-white"
                style={{ minHeight: '250px' }}
            >
                {!isLoaded && !error && (
                    <div className="flex items-center justify-center h-full">
                        <div className="text-center text-slate-600">
                            <div className="w-8 h-8 border-2 border-blue-500 border-t-transparent rounded-full animate-spin mx-auto mb-2" />
                            <p className="text-sm">Loading Kekule.js...</p>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
}
