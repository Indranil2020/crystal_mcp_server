/**
 * ExportMenu - Structure and image export functionality
 */

import { useState, useRef, useEffect } from 'react';
import { useAppSelector } from '../../store/hooks';
import { structureToCif, structureToPdb, structureToXyz } from '../../converters';

export default function ExportMenu() {
    const [isOpen, setIsOpen] = useState(false);
    const menuRef = useRef<HTMLDivElement>(null);
    const { structures, activeStructureId } = useAppSelector(state => state.structure);
    const activeStructure = structures.find(s => s.id === activeStructureId);

    // Close menu on outside click
    useEffect(() => {
        function handleClickOutside(event: MouseEvent) {
            if (menuRef.current && !menuRef.current.contains(event.target as Node)) {
                setIsOpen(false);
            }
        }
        document.addEventListener('mousedown', handleClickOutside);
        return () => document.removeEventListener('mousedown', handleClickOutside);
    }, []);

    const downloadFile = (content: string, filename: string, mimeType: string) => {
        const blob = new Blob([content], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        URL.revokeObjectURL(url);
        setIsOpen(false);
    };

    const exportCif = () => {
        if (!activeStructure) return;
        const cif = structureToCif(activeStructure.data);
        const filename = `${activeStructure.name.replace(/\s+/g, '_')}.cif`;
        downloadFile(cif, filename, 'chemical/x-cif');
    };

    const exportPdb = () => {
        if (!activeStructure) return;
        const pdb = structureToPdb(activeStructure.data);
        const filename = `${activeStructure.name.replace(/\s+/g, '_')}.pdb`;
        downloadFile(pdb, filename, 'chemical/x-pdb');
    };

    const exportXyz = () => {
        if (!activeStructure) return;
        const xyz = structureToXyz(activeStructure.data);
        const filename = `${activeStructure.name.replace(/\s+/g, '_')}.xyz`;
        downloadFile(xyz, filename, 'chemical/x-xyz');
    };

    const exportJson = () => {
        if (!activeStructure) return;
        const json = JSON.stringify(activeStructure.data, null, 2);
        const filename = `${activeStructure.name.replace(/\s+/g, '_')}.json`;
        downloadFile(json, filename, 'application/json');
    };

    return (
        <div ref={menuRef} className="relative">
            <button
                onClick={() => setIsOpen(!isOpen)}
                disabled={!activeStructure}
                className="px-3 py-1 rounded text-sm bg-green-600 text-white 
                   hover:bg-green-500 disabled:opacity-50 disabled:cursor-not-allowed
                   flex items-center gap-1"
            >
                <span>📥</span>
                <span>Export</span>
            </button>

            {isOpen && activeStructure && (
                <div className="absolute right-0 mt-1 w-48 bg-slate-800 border border-slate-600 
                        rounded-lg shadow-lg overflow-hidden z-20">
                    <div className="px-3 py-2 text-xs text-slate-400 border-b border-slate-700">
                        Structure Formats
                    </div>
                    <button
                        onClick={exportCif}
                        className="w-full px-3 py-2 text-left text-sm text-slate-200 hover:bg-slate-700 flex items-center gap-2"
                    >
                        <span>📄</span> CIF (.cif)
                    </button>
                    <button
                        onClick={exportPdb}
                        className="w-full px-3 py-2 text-left text-sm text-slate-200 hover:bg-slate-700 flex items-center gap-2"
                    >
                        <span>📄</span> PDB (.pdb)
                    </button>
                    <button
                        onClick={exportXyz}
                        className="w-full px-3 py-2 text-left text-sm text-slate-200 hover:bg-slate-700 flex items-center gap-2"
                    >
                        <span>📄</span> XYZ (.xyz)
                    </button>
                    <button
                        onClick={exportJson}
                        className="w-full px-3 py-2 text-left text-sm text-slate-200 hover:bg-slate-700 flex items-center gap-2"
                    >
                        <span>📋</span> JSON (.json)
                    </button>

                    <div className="px-3 py-2 text-xs text-slate-400 border-t border-slate-700">
                        Image Formats
                    </div>
                    <button
                        onClick={() => {
                            // Trigger screenshot from MolStar viewer
                            const event = new CustomEvent('export-screenshot');
                            window.dispatchEvent(event);
                            setIsOpen(false);
                        }}
                        className="w-full px-3 py-2 text-left text-sm text-slate-200 hover:bg-slate-700 flex items-center gap-2"
                    >
                        <span>🖼️</span> PNG Screenshot
                    </button>
                </div>
            )}
        </div>
    );
}
