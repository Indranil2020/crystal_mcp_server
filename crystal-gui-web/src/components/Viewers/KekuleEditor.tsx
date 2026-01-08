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
        debug('VIEWERS', '✏️ KEKULE EDITOR INITIALIZATION');

        // Check if already loaded
        if (window.Kekule) {
            debug('VIEWERS', '  Kekule.js already loaded, initializing editor...');
            initEditor();
            return;
        }

        debug('VIEWERS', '  Loading Kekule.js from CDN...');

        // Load CSS
        const css = document.createElement('link');
        css.rel = 'stylesheet';
        css.href = 'https://cdn.jsdelivr.net/npm/kekule/dist/themes/default/kekule.css';
        document.head.appendChild(css);
        debug('VIEWERS', '  CSS loaded');

        // Load main script
        const script = document.createElement('script');
        script.src = 'https://cdn.jsdelivr.net/npm/kekule/dist/kekule.min.js';
        script.onload = () => {
            debug('VIEWERS', '  ✓ Kekule.js script loaded successfully');
            initEditor();
        };
        script.onerror = () => {
            debugError('VIEWERS', 'Failed to load Kekule.js from CDN', 'KekuleEditor');
            setError('Failed to load Kekule.js');
        };
        document.head.appendChild(script);
        debug('VIEWERS', '  Script tag added, waiting for load...');
        debug('VIEWERS', '═'.repeat(50));
    }, []);

    // Initialize the Kekule Composer
    const initEditor = useCallback(() => {
        debug('VIEWERS', '  initEditor called');

        if (!editorRef.current) {
            debug('VIEWERS', '  ❌ editorRef not ready');
            return;
        }
        if (!window.Kekule) {
            debug('VIEWERS', '  ❌ window.Kekule not available');
            return;
        }

        const checkReadiness = (attempts = 0) => {
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const Kekule = window.Kekule as any;

            // Check specific dependencies required for Composer
            const isKekuleReady = !!Kekule;
            const isWidgetReady = !!(Kekule && Kekule.Widget);
            const isClassReady = !!(Kekule && Kekule.Widget && Kekule.Widget.getPreferredWidgetClass);

            debug('VIEWERS', `  Readiness check #${attempts}: Kekule=${isKekuleReady}, Widget=${isWidgetReady}, Class=${isClassReady}`);

            // If any dependency is missing, retry or determine unsupported
            if (!isKekuleReady || !isWidgetReady || !isClassReady) {
                if (attempts > 20) {
                    debug('VIEWERS', '  ⚠️ Kekule Widget subsystem unavailable after 20 attempts');
                    debug('VIEWERS', '  This usually indicates a headless browser environment');
                    setError('Kekule Widget system not available in this environment');
                    return;
                }
                setTimeout(() => checkReadiness(attempts + 1), 100);
                return;
            }

            debug('VIEWERS', '  ✓ All Kekule systems ready, creating Composer...');

            // All systems go - initialize directly
            const composer = new Kekule.Editor.Composer(editorRef.current);
            composer.setDimension('100%', '100%');
            debug('VIEWERS', '  Composer created');

            // Configure toolbar - strict feature checks instead of try/catch
            if (typeof composer.setCommonToolButtons === 'function') {
                composer.setCommonToolButtons([
                    'newDoc', 'loadData', 'saveData',
                    'undo', 'redo',
                    'copy', 'cut', 'paste',
                    'zoomIn', 'zoomOut', 'reset',
                ]);
                debug('VIEWERS', '  Common tool buttons configured');
            }

            if (typeof composer.setChemToolButtons === 'function') {
                composer.setChemToolButtons([
                    'manipulate', 'erase', 'bond', 'atom',
                    'ring', 'charge', 'glyph',
                ]);
                debug('VIEWERS', '  Chem tool buttons configured');
            }

            // Listen for changes
            composer.addEventListener('editObjsUpdated', () => {
                const mol = composer.getChemObj();
                if (mol && Kekule.IO && Kekule.IO.saveFormatData) {
                    const smiles = Kekule.IO.saveFormatData(mol, 'smi');
                    debug('VIEWERS', `  Molecule updated: SMILES=${smiles}`);
                    setCurrentSmiles(smiles);
                    onStructureChange?.(smiles);
                }
            });

            composerRef.current = composer;
            setIsLoaded(true);
            debug('VIEWERS', '  ✅ KEKULE COMPOSER INITIALIZED SUCCESSFULLY');
        };

        checkReadiness();
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
                <div className="text-center text-amber-400">
                    <p className="text-lg">⚠️ 2D Editor Unavailable</p>
                    <p className="text-sm mt-2 text-slate-400">Kekule.js requires a graphical environment.</p>
                    <p className="text-xs mt-1 text-slate-500">(Headless browser detected)</p>
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
