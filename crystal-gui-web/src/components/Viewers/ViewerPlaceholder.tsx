/**
 * ViewerPlaceholder - Placeholder for MolStar/Kekule viewers
 * 
 * To be replaced with actual viewer implementations in Phase 2
 */

import { useAppSelector } from '../../store/hooks';

interface Props {
    type: '3D' | '2D';
}

export default function ViewerPlaceholder({ type }: Props) {
    const { structures, activeStructureId, viewerSettings } = useAppSelector(state => state.structure);
    const activeStructure = structures.find(s => s.id === activeStructureId);

    return (
        <div className="flex-1 flex flex-col bg-slate-900">
            {/* Header */}
            <div className="px-4 py-2 border-b border-slate-700 flex items-center justify-between">
                <h3 className="text-sm font-medium text-slate-300">
                    {type === '3D' ? '🔬 3D Viewer (MolStar)' : '✏️ 2D Editor (Kekule.js)'}
                </h3>
                <span className="text-xs text-slate-500">
                    {viewerSettings.representation}
                </span>
            </div>

            {/* Viewer area */}
            <div className="flex-1 flex items-center justify-center p-4">
                {activeStructure ? (
                    <div className="text-center">
                        {/* Placeholder visualization */}
                        <div className="w-64 h-64 mx-auto mb-4 rounded-lg bg-slate-800 border border-slate-700 
                          flex items-center justify-center">
                            <div className="text-6xl opacity-50">
                                {type === '3D' ? '🧬' : '🔗'}
                            </div>
                        </div>

                        {/* Structure info */}
                        <div className="text-slate-300 font-medium">
                            {activeStructure.name}
                        </div>
                        <div className="text-sm text-slate-500 mt-1">
                            {activeStructure.data.atoms?.length || 0} atoms
                            {activeStructure.data.metadata?.formula && (
                                <span className="ml-2">• {activeStructure.data.metadata.formula}</span>
                            )}
                        </div>
                        <div className="text-xs text-slate-600 mt-2">
                            {type === '3D'
                                ? 'MolStar viewer will be integrated in Phase 2'
                                : 'Kekule.js editor will be integrated in Phase 2'
                            }
                        </div>

                        {/* Lattice info for crystals */}
                        {activeStructure.data.lattice && (
                            <div className="mt-4 text-xs text-slate-500">
                                <div>
                                    a={activeStructure.data.lattice.a.toFixed(2)}Å,
                                    b={activeStructure.data.lattice.b.toFixed(2)}Å,
                                    c={activeStructure.data.lattice.c.toFixed(2)}Å
                                </div>
                                {activeStructure.data.space_group && (
                                    <div className="mt-1">
                                        Space group: {activeStructure.data.space_group.symbol}
                                        (#{activeStructure.data.space_group.number})
                                    </div>
                                )}
                            </div>
                        )}
                    </div>
                ) : (
                    <div className="text-center text-slate-500">
                        <div className="text-4xl mb-4 opacity-50">
                            {type === '3D' ? '🔬' : '✏️'}
                        </div>
                        <p className="text-lg">No structure loaded</p>
                        <p className="text-sm mt-2">
                            Use the chat to generate a molecule or crystal structure
                        </p>
                    </div>
                )}
            </div>
        </div>
    );
}
