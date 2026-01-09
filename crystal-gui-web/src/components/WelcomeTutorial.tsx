/**
 * WelcomeTutorial.tsx
 * 
 * First-run tutorial overlay showing key features and shortcuts.
 */

import { useState, useEffect } from 'react';
import { X, Keyboard, MousePointer, Info } from 'lucide-react';

export function WelcomeTutorial() {
    const [isVisible, setIsVisible] = useState(false);

    useEffect(() => {
        // Check if tutorial has been shown before
        const hasSeenTutorial = localStorage.getItem('crystal-gui-seen-tutorial');
        if (!hasSeenTutorial) {
            setIsVisible(true);
        }
    }, []);

    const handleDismiss = () => {
        setIsVisible(false);
        localStorage.setItem('crystal-gui-seen-tutorial', 'true');
    };

    if (!isVisible) return null;

    return (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/80 backdrop-blur-sm p-4">
            <div className="bg-slate-800 rounded-xl max-w-2xl w-full shadow-2xl border border-slate-700 overflow-hidden flex flex-col max-h-[90vh]">

                {/* Header */}
                <div className="bg-gradient-to-r from-blue-600 to-indigo-600 p-6 flex justify-between items-start">
                    <div>
                        <h2 className="text-2xl font-bold text-white mb-2">Welcome to Crystal GUI</h2>
                        <p className="text-blue-100">Your advanced molecular visualization workspace</p>
                    </div>
                    <button
                        onClick={handleDismiss}
                        className="text-white/70 hover:text-white transition-colors"
                    >
                        <X size={24} />
                    </button>
                </div>

                {/* Content */}
                <div className="p-6 overflow-y-auto space-y-8 flex-1">

                    {/* Workflow Section */}
                    <div className="flex gap-4">
                        <div className="bg-blue-500/10 p-3 rounded-lg h-fit">
                            <Info className="text-blue-400" size={24} />
                        </div>
                        <div>
                            <h3 className="text-lg font-semibold text-slate-200 mb-2">How it works</h3>
                            <p className="text-slate-400 leading-relaxed">
                                Use the <strong>Chat Panel</strong> on the left to request molecules ("generate aspirin")
                                or complex crystals ("create gold nanocluster"). The AI will use tools to build,
                                optimize, and display them in the center viewer.
                            </p>
                        </div>
                    </div>

                    {/* Controls Grid */}
                    <div className="grid md:grid-cols-2 gap-6">
                        {/* Mouse Controls */}
                        <div className="bg-slate-900/50 p-4 rounded-lg border border-slate-700/50">
                            <div className="flex items-center gap-2 mb-4 text-emerald-400">
                                <MousePointer size={20} />
                                <h4 className="font-medium">Mouse Controls</h4>
                            </div>
                            <ul className="space-y-3 text-sm text-slate-300">
                                <li className="flex justify-between">
                                    <span>Rotate</span>
                                    <span className="text-slate-500">Left Click + Drag</span>
                                </li>
                                <li className="flex justify-between">
                                    <span>Zoom</span>
                                    <span className="text-slate-500">Scroll Wheel</span>
                                </li>
                                <li className="flex justify-between">
                                    <span>Pan</span>
                                    <span className="text-slate-500">Right Click + Drag</span>
                                </li>
                                <li className="flex justify-between">
                                    <span>Select</span>
                                    <span className="text-slate-500">Click Atom</span>
                                </li>
                            </ul>
                        </div>

                        {/* Keyboard Shortcuts */}
                        <div className="bg-slate-900/50 p-4 rounded-lg border border-slate-700/50">
                            <div className="flex items-center gap-2 mb-4 text-purple-400">
                                <Keyboard size={20} />
                                <h4 className="font-medium">Keyboard Shortcuts</h4>
                            </div>
                            <ul className="space-y-3 text-sm text-slate-300">
                                <li className="flex justify-between">
                                    <span>Undo / Redo</span>
                                    <kbd className="bg-slate-700 px-2 rounded text-xs py-0.5">Ctrl + Z / Y</kbd>
                                </li>
                                <li className="flex justify-between">
                                    <span>Delete Selection</span>
                                    <kbd className="bg-slate-700 px-2 rounded text-xs py-0.5">Delete</kbd>
                                </li>
                                <li className="flex justify-between">
                                    <span>Toggle Unit Cell</span>
                                    <kbd className="bg-slate-700 px-2 rounded text-xs py-0.5">C</kbd>
                                </li>
                                <li className="flex justify-between">
                                    <span>Reset View</span>
                                    <kbd className="bg-slate-700 px-2 rounded text-xs py-0.5">R</kbd>
                                </li>
                                <li className="flex justify-between">
                                    <span>Representations</span>
                                    <kbd className="bg-slate-700 px-2 rounded text-xs py-0.5">1 - 5</kbd>
                                </li>
                            </ul>
                        </div>
                    </div>

                    {/* Pro Tip */}
                    <div className="bg-amber-900/20 border border-amber-900/50 rounded-lg p-4">
                        <p className="text-sm text-amber-200/80">
                            <strong>Pro Tip:</strong> Toggle the 2D Editor (Right Panel) to draw custom chemical structures
                            and send them to the 3D viewer.
                        </p>
                    </div>
                </div>

                {/* Footer */}
                <div className="p-6 bg-slate-900/50 border-t border-slate-700 flex justify-end">
                    <button
                        onClick={handleDismiss}
                        className="px-6 py-2.5 bg-blue-600 hover:bg-blue-500 text-white rounded-lg font-medium transition-colors shadow-lg shadow-blue-900/20"
                    >
                        Get Started
                    </button>
                </div>
            </div>
        </div>
    );
}

export default WelcomeTutorial;
