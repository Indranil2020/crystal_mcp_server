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

    // Load molecule from SMILES (reserved for future use)
    // const _loadFromSmiles = useCallback((smiles: string) => {
    //     if (!composerRef.current || !window.Kekule) return;
    //     // eslint-disable-next-line @typescript-eslint/no-explicit-any
    //     const Kekule = window.Kekule as any;
    //
    //     if (Kekule.IO && Kekule.IO.loadFormatData) {
    //         const mol = Kekule.IO.loadFormatData(smiles, 'smi');
    //         if (mol) {
    //             composerRef.current.setChemObj(mol);
    //         }
    //     }
    // }, []);

    // Push structure to 3D viewer via MDL MOL export (robust, preserves bonds)
    const pushTo3D = useCallback(() => {
        console.log('[KekuleEditor] pushTo3D() called');

        if (!composerRef.current || !window.Kekule) {
            console.warn('[KekuleEditor] Editor not ready');
            return;
        }

        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const Kekule = window.Kekule as any;

        let mol = composerRef.current.getChemObj();
        if (!mol) {
            console.warn('[KekuleEditor] No molecule to push');
            return;
        }

        // Debug object type
        const className = mol.getClassName ? mol.getClassName() : 'Unknown';
        console.log(`[KekuleEditor] Object type: ${className}`);

        // Unwrap ChemDocument/ChemSpace wrapper to get the actual Molecule
        if (className === 'Kekule.ChemDocument' || className === 'Kekule.ChemSpace') {
            const childCount = mol.getChildCount ? mol.getChildCount() : 0;
            console.log(`[KekuleEditor] Wrapper has ${childCount} children`);

            let foundMol = null;
            for (let i = 0; i < childCount; i++) {
                const child = mol.getChildAt(i);
                const childClass = child?.getClassName?.() || '';
                console.log(`[KekuleEditor]   Child ${i}: ${childClass}`);
                if (childClass === 'Kekule.Molecule') {
                    foundMol = child;
                    break;
                }
            }

            if (!foundMol) {
                console.warn('[KekuleEditor] No Molecule found in wrapper');
                return;
            }
            mol = foundMol;
            console.log('[KekuleEditor] Extracted Molecule from wrapper');
        }

        // Validate molecule has atoms
        const nodeCount = mol.getNodeCount?.() || 0;
        console.log(`[KekuleEditor] Molecule has ${nodeCount} atoms`);
        if (nodeCount === 0) {
            console.warn('[KekuleEditor] Molecule has no atoms');
            return;
        }

        // Export as MDL MOL V2000 format (preserves bonds, stereochemistry)
        let molData: string | null = null;
        if (Kekule.IO?.saveFormatData) {
            molData = Kekule.IO.saveFormatData(mol, 'mol');
            console.log('[KekuleEditor] Exported MOL data:', molData?.substring(0, 100) + '...');
        } else {
            console.error('[KekuleEditor] Kekule.IO.saveFormatData not available');
            return;
        }

        // Validate MOL export
        if (!molData || !molData.includes('V2000')) {
            console.error('[KekuleEditor] MOL export failed or invalid format');
            return;
        }

        // Extract formula for naming
        let formula = 'Molecule';
        if (mol.calcFormula) {
            formula = mol.calcFormula() || formula;
        }
        console.log(`[KekuleEditor] Formula: ${formula}`);

        // Create structure for Redux with raw MOL data
        const structure = {
            id: uuidv4(),
            name: `Drawn: ${formula}`,
            data: {
                lattice: { a: 20, b: 20, c: 20, alpha: 90, beta: 90, gamma: 90, matrix: [[20, 0, 0], [0, 20, 0], [0, 0, 20]] },
                atoms: [],  // Empty - viewer will use molData instead
                metadata: { formula, natoms: nodeCount },
            },
            molData,         // Raw MDL MOL string
            format: 'mol' as const,  // Format identifier for viewer
            source: 'kekule' as const,
            createdAt: Date.now(),
            modifiedAt: Date.now(),
            visible: true,
        };

        dispatch(addStructure(structure));
        console.log('[KekuleEditor] ✅ Pushed to 3D:', structure.name);
    }, [dispatch]);

    // Clear editor
    const clearEditor = useCallback(() => {
        if (composerRef.current) {
            composerRef.current.newDoc();
            setCurrentSmiles('');
        }
    }, []);

    // Template molecules (using MDL Molfile V2000 format for reliability without OpenBabel)
    const loadTemplate = useCallback((template: string) => {
        // Simple Benzene MolBlock
        const BENZENE_MOL = `
  Kekule.js   

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END`;

        // Cyclohexane (Chair confirmation approximation)
        const CYCLOHEXANE_MOL = `
  Kekule.js

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  1  1  0  0  0  0
M  END`;

        const templates: Record<string, string> = {
            benzene: BENZENE_MOL,
            cyclohexane: CYCLOHEXANE_MOL,
            pyridine: '', // TODO: Add others
            naphthalene: '',
            phenol: '',
            aniline: '',
        };

        if (templates[template] && composerRef.current && window.Kekule) {
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const Kekule = window.Kekule as any;
            if (Kekule.IO && Kekule.IO.loadFormatData) {
                try {
                    const mol = Kekule.IO.loadFormatData(templates[template], 'mol');
                    if (mol) {
                        composerRef.current.setChemObj(mol);
                    }
                } catch (e) {
                    debugError('VIEWERS', `Failed to load template: ${e}`, 'KekuleEditor');
                }
            }
        }
    }, []);

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
