/**
 * Toolbar - Main application toolbar with actions
 */

import { useAppSelector, useAppDispatch } from '../../store/hooks';
import { updateViewerSettings } from '../../store/structureSlice';
import { clearMessages } from '../../store/chatSlice';
import ExportMenu from './ExportMenu';

export default function Toolbar() {
    const dispatch = useAppDispatch();
    const { structures, activeStructureId, viewerSettings } = useAppSelector(state => state.structure);
    const activeStructure = structures.find(s => s.id === activeStructureId);

    return (
        <div className="bg-slate-800 border-b border-slate-700 px-4 py-2 flex items-center justify-between">
            {/* Left: Logo & Title */}
            <div className="flex items-center gap-3">
                <div className="text-xl">🔮</div>
                <div>
                    <h1 className="text-lg font-bold text-slate-100">Crystal GUI</h1>
                    <p className="text-xs text-slate-400">Molecular Structure Generator</p>
                </div>
            </div>

            {/* Center: Structure Info */}
            <div className="flex items-center gap-4">
                {activeStructure && (
                    <div className="text-sm text-slate-300">
                        <span className="text-slate-500">Active: </span>
                        <span className="font-medium">{activeStructure.name}</span>
                        <span className="text-slate-500 ml-2">
                            ({activeStructure.data.atoms?.length || 0} atoms)
                        </span>
                    </div>
                )}
            </div>

            {/* Right: Actions */}
            <div className="flex items-center gap-2">
                {/* Representation dropdown */}
                <select
                    value={viewerSettings.representation}
                    onChange={e => dispatch(updateViewerSettings({
                        representation: e.target.value as typeof viewerSettings.representation
                    }))}
                    className="bg-slate-700 text-slate-200 text-sm rounded px-2 py-1 
                     border border-slate-600 focus:outline-none focus:ring-1 focus:ring-blue-500"
                >
                    <option value="ball-and-stick">Ball & Stick</option>
                    <option value="spacefill">Spacefill</option>
                    <option value="licorice">Licorice</option>
                    <option value="wireframe">Wireframe</option>
                </select>

                {/* Unit cell toggle */}
                <button
                    onClick={() => dispatch(updateViewerSettings({ showUnitCell: !viewerSettings.showUnitCell }))}
                    className={`px-3 py-1 rounded text-sm ${viewerSettings.showUnitCell
                        ? 'bg-blue-600 text-white'
                        : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
                        }`}
                >
                    Unit Cell
                </button>

                {/* Clear chat */}
                <button
                    onClick={() => dispatch(clearMessages())}
                    className="px-3 py-1 rounded text-sm bg-slate-700 text-slate-300 
                     hover:bg-slate-600 hover:text-white transition-colors"
                >
                    Clear
                </button>

                {/* Export menu */}
                <ExportMenu />
            </div>
        </div>
    );
}

