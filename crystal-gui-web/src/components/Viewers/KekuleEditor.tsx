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

            // Configure Common Toolbar (file ops, edit ops, view ops, settings)
            debug('VIEWERS', '  ⚙️ Configuring Toolbars...');
            if (composer.setCommonToolButtons) {
                composer.setCommonToolButtons([
                    'newDoc', 'loadData', 'saveData',  // File operations
                    'undo', 'redo',                      // Edit history
                    'copy', 'cut', 'paste',              // Clipboard
                    'zoomIn', 'zoomOut', 'reset',        // View controls
                    'config',                            // Configuration dialog
                    'objInspector',                      // Object inspector panel
                ]);
                debug('VIEWERS', '    ✓ Common Toolbar configured (12 buttons)');
            } else {
                debug('VIEWERS', '  ⚠️ setCommonToolButtons missing from composer instance');
            }

            // Configure Chem Toolbar (all chemistry tools)
            if (composer.setChemToolButtons) {
                composer.setChemToolButtons([
                    'manipulate',        // Selection and move
                    'erase',             // Delete objects
                    'bond',              // Draw bonds (single, double, triple, etc.)
                    'atomAndFormula',    // Input atom symbol or formula
                    'ring',              // Common ring structures
                    'charge',            // Add charge (+/-)
                    'glyph',             // Reaction arrows, symbols
                    'textImage',         // Text and image annotations
                ]);
                debug('VIEWERS', '    ✓ Chem Toolbar configured (8 buttons)');
            } else {
                debug('VIEWERS', '  ⚠️ setChemToolButtons missing from composer instance');
            }

            // Configure Style Toolbar (text/color styling)
            if (composer.setStyleToolComponentNames) {
                composer.setStyleToolComponentNames([
                    'fontName',          // Font family selector
                    'fontSize',          // Font size selector
                    'color',             // Color picker
                    'textDirection',     // Text direction (LTR/RTL)
                    'textAlign',         // Text alignment
                ]);
                debug('VIEWERS', '    ✓ Style Toolbar configured (5 components)');
            } else {
                debug('VIEWERS', '  ⚠️ setStyleToolComponentNames not available');
            }

            // Enable the Object Inspector by default for advanced editing
            if (composer.setEnableObjModifier) {
                composer.setEnableObjModifier(true);
                debug('VIEWERS', '    ✓ Object Modifier enabled');
            }

            // Listen for changes and sync SMILES
            composer.addEventListener('editObjsUpdated', () => {
                const mol = composer.getChemObj();
                if (mol && Kekule.IO && Kekule.IO.saveFormatData) {
                    try {
                        const smiles = Kekule.IO.saveFormatData(mol, 'smi');
                        setCurrentSmiles(smiles || '');
                        onStructureChange?.(smiles || '');
                    } catch {
                        // SMILES export may fail for complex objects, ignore
                    }
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

    // Push structure to 3D viewer via backend RDKit (ensures chemically valid structures with H atoms)
    const pushTo3D = useCallback(async () => {
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

        // Collect all molecules for SMILES export
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const molecules: any[] = [];

        // Handle ChemDocument/ChemSpace wrapper
        if (className === 'Kekule.ChemDocument' || className === 'Kekule.ChemSpace') {
            const childCount = mol.getChildCount ? mol.getChildCount() : 0;
            console.log(`[KekuleEditor] Wrapper has ${childCount} children`);

            for (let i = 0; i < childCount; i++) {
                const child = mol.getChildAt(i);
                const childClass = child?.getClassName?.() || '';
                console.log(`[KekuleEditor]   Child ${i}: ${childClass}`);
                if (childClass === 'Kekule.Molecule') {
                    molecules.push(child);
                }
            }

            if (molecules.length === 0) {
                console.warn('[KekuleEditor] No Molecule found in wrapper');
                return;
            }

            // For single molecule, use it directly
            if (molecules.length === 1) {
                mol = molecules[0];
                console.log('[KekuleEditor] Single molecule - using directly');
            }
        } else {
            molecules.push(mol);
        }

        // Export as SMILES for each molecule
        console.log(`[KekuleEditor] Exporting ${molecules.length} molecule(s) to SMILES...`);
        const smilesStrings: string[] = [];

        for (let i = 0; i < molecules.length; i++) {
            const m = molecules[i];
            if (Kekule.IO?.saveFormatData) {
                const smiles = Kekule.IO.saveFormatData(m, 'smi');
                if (smiles && smiles.trim()) {
                    smilesStrings.push(smiles.trim());
                    console.log(`[KekuleEditor]   Molecule ${i}: ${smiles.trim()}`);
                }
            }
        }

        if (smilesStrings.length === 0) {
            console.error('[KekuleEditor] Failed to export SMILES');
            // Fallback to direct MOL export (without hydrogens)
            console.warn('[KekuleEditor] Falling back to direct MOL export (no hydrogen saturation)');
            await pushTo3DFallback(mol, Kekule);
            return;
        }

        // Combine SMILES using dot notation for disconnected fragments
        const combinedSmiles = smilesStrings.join('.');
        console.log(`[KekuleEditor] Combined SMILES: ${combinedSmiles}`);

        // Call backend to generate 3D structure with hydrogens
        console.log('[KekuleEditor] Calling backend /chemistry/smiles-to-3d...');

        const response = await fetch('http://localhost:8080/chemistry/smiles-to-3d', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                smiles: combinedSmiles,
                optimize: true,
                name: 'Drawn Molecule',
            }),
        });

        if (!response.ok) {
            const errorText = await response.text();
            console.error('[KekuleEditor] Backend error:', response.status, errorText);
            // Fallback to direct MOL export
            console.warn('[KekuleEditor] Falling back to direct MOL export (backend unavailable)');
            await pushTo3DFallback(mol, Kekule);
            return;
        }

        const result = await response.json();
        console.log('[KekuleEditor] Backend response:', result);
        console.log(`[KekuleEditor] Formula: ${result.formula}, Atoms: ${result.n_atoms} (${result.n_heavy_atoms} heavy)`);

        // Log validation results
        if (result.validation) {
            console.log(`[KekuleEditor] Validation: valid=${result.validation.valid}`);
            if (result.validation.warnings?.length > 0) {
                console.warn('[KekuleEditor] Warnings:', result.validation.warnings);
            }
            if (result.validation.errors?.length > 0) {
                console.error('[KekuleEditor] Errors:', result.validation.errors);
            }
        }

        // Create structure for Redux with SDF data from backend
        const structure = {
            id: uuidv4(),
            name: `Drawn: ${result.formula}`,
            data: {
                lattice: { a: 20, b: 20, c: 20, alpha: 90, beta: 90, gamma: 90, matrix: [[20, 0, 0], [0, 20, 0], [0, 0, 20]] },
                atoms: [],
                metadata: {
                    formula: result.formula,
                    natoms: result.n_atoms,
                    molecular_weight: result.molecular_weight,
                    smiles: result.canonical_smiles,
                    validation: result.validation,
                },
            },
            molData: result.sdf,  // SDF with 3D coords and explicit H
            format: 'mol' as const,
            source: 'kekule' as const,
            createdAt: Date.now(),
            modifiedAt: Date.now(),
            visible: true,
        };

        console.log(`[KekuleEditor] Structure: ${structure.name} (${result.n_atoms} atoms, ${result.sdf.length} bytes)`);
        dispatch(addStructure(structure));
        console.log('[KekuleEditor] ✅ Pushed to 3D with hydrogens:', structure.name);
    }, [dispatch]);

    // Fallback: Direct MOL export without backend (no hydrogen saturation)
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const pushTo3DFallback = async (mol: any, Kekule: any) => {
        console.log('[KekuleEditor] pushTo3DFallback() - exporting without hydrogen saturation');

        let molData: string | null = null;
        if (Kekule.IO?.saveFormatData) {
            molData = Kekule.IO.saveFormatData(mol, 'mol');
        }

        if (!molData || !molData.includes('V2000')) {
            console.error('[KekuleEditor] Fallback MOL export failed');
            return;
        }

        // Extract formula
        let formula = 'Molecule';
        if (mol.calcFormula) {
            const formulaResult = mol.calcFormula();
            if (formulaResult) {
                formula = typeof formulaResult === 'string'
                    ? formulaResult
                    : (formulaResult.getText?.() || formulaResult.toString?.() || 'Molecule');
            }
        }

        const structure = {
            id: uuidv4(),
            name: `Drawn: ${formula} (no H)`,
            data: {
                lattice: { a: 20, b: 20, c: 20, alpha: 90, beta: 90, gamma: 90, matrix: [[20, 0, 0], [0, 20, 0], [0, 0, 20]] },
                atoms: [],
                metadata: { formula, natoms: mol.getNodeCount?.() || 0 },
            },
            molData,
            format: 'mol' as const,
            source: 'kekule' as const,
            createdAt: Date.now(),
            modifiedAt: Date.now(),
            visible: true,
        };

        dispatch(addStructure(structure));
        console.log('[KekuleEditor] ✅ Pushed to 3D (fallback, no H):', structure.name);
    };

    // Clear editor
    const clearEditor = useCallback(() => {
        if (composerRef.current) {
            composerRef.current.newDoc();
            setCurrentSmiles('');
        }
    }, []);

    // Template molecules (using MDL Molfile V2000 format for reliability)
    const loadTemplate = useCallback((template: string) => {
        // Benzene (C6H6) - aromatic ring
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

        // Cyclohexane (C6H12) - saturated ring
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

        // Pyridine (C5H5N) - 6-membered ring with N
        const PYRIDINE_MOL = `
  Kekule.js

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
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

        // Naphthalene (C10H8) - fused benzene rings
        const NAPHTHALENE_MOL = `
  Kekule.js

 10 11  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7320    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5980    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5980   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7320   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10  5  1  0  0  0  0
M  END`;

        // Phenol (C6H5OH) - benzene with OH group
        const PHENOL_MOL = `
  Kekule.js

  7  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
M  END`;

        // Aniline (C6H5NH2) - benzene with NH2 group
        const ANILINE_MOL = `
  Kekule.js

  7  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
M  END`;

        const templates: Record<string, string> = {
            benzene: BENZENE_MOL,
            cyclohexane: CYCLOHEXANE_MOL,
            pyridine: PYRIDINE_MOL,
            naphthalene: NAPHTHALENE_MOL,
            phenol: PHENOL_MOL,
            aniline: ANILINE_MOL,
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

                    <div className="w-px h-4 bg-slate-600 mx-1" />

                    {/* Export MOL */}
                    <button
                        onClick={() => {
                            if (!composerRef.current || !window.Kekule) return;
                            // eslint-disable-next-line @typescript-eslint/no-explicit-any
                            const Kekule = window.Kekule as any;
                            const mol = composerRef.current.getChemObj();
                            if (mol && Kekule.IO?.saveFormatData) {
                                const molData = Kekule.IO.saveFormatData(mol, 'mol');
                                const blob = new Blob([molData], { type: 'chemical/x-mdl-molfile' });
                                const url = URL.createObjectURL(blob);
                                const a = document.createElement('a');
                                a.href = url;
                                a.download = 'molecule.mol';
                                a.click();
                                URL.revokeObjectURL(url);
                                console.log('[KekuleEditor] Exported MOL file');
                            }
                        }}
                        disabled={!isLoaded}
                        className="text-xs px-2 py-1 rounded bg-slate-700 text-slate-300 hover:bg-slate-600 disabled:opacity-50"
                        title="Export as MDL MOL file"
                    >
                        Export
                    </button>

                    {/* Import MOL */}
                    <label className="text-xs px-2 py-1 rounded bg-slate-700 text-slate-300 hover:bg-slate-600 cursor-pointer">
                        Import
                        <input
                            type="file"
                            accept=".mol,.sdf"
                            className="hidden"
                            onChange={(e) => {
                                const file = e.target.files?.[0];
                                if (!file || !composerRef.current || !window.Kekule) return;
                                // eslint-disable-next-line @typescript-eslint/no-explicit-any
                                const Kekule = window.Kekule as any;
                                const reader = new FileReader();
                                reader.onload = (ev) => {
                                    const data = ev.target?.result as string;
                                    if (data && Kekule.IO?.loadFormatData) {
                                        const mol = Kekule.IO.loadFormatData(data, 'mol');
                                        if (mol) {
                                            composerRef.current?.setChemObj(mol);
                                            console.log('[KekuleEditor] Imported MOL file');
                                        }
                                    }
                                };
                                reader.readAsText(file);
                                e.target.value = ''; // Reset for re-import
                            }}
                        />
                    </label>

                    <div className="w-px h-4 bg-slate-600 mx-1" />

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
