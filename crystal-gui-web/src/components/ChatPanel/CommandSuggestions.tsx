/**
 * CommandSuggestions - Autocomplete suggestions for chat input
 */

import { useState, useEffect } from 'react';

interface Props {
    input: string;
    onSelect: (suggestion: string) => void;
    isVisible: boolean;
}

const SUGGESTIONS = [
    { text: 'Generate benzene', description: 'Create a benzene ring molecule' },
    { text: 'Generate aspirin', description: 'Create acetylsalicylic acid' },
    { text: 'Generate caffeine', description: 'Create caffeine molecule' },
    { text: 'Create benzene dimer', description: 'Pi-pi stacked benzene pair' },
    { text: 'Create water tetramer', description: 'Hydrogen-bonded water cluster' },
    { text: 'Generate silicon crystal', description: 'Silicon diamond structure' },
    { text: 'Generate NaCl', description: 'Sodium chloride rock salt' },
    { text: 'Make supercell 2x2x2', description: 'Expand current structure' },
    { text: 'Generate slab (111)', description: 'Create surface slab' },
    { text: 'List available tools', description: 'Show all MCP tools' },
];

export default function CommandSuggestions({ input, onSelect, isVisible }: Props) {
    const [filteredSuggestions, setFilteredSuggestions] = useState(SUGGESTIONS);
    const [selectedIndex, setSelectedIndex] = useState(0);

    useEffect(() => {
        if (!input.trim()) {
            setFilteredSuggestions(SUGGESTIONS.slice(0, 5));
            return;
        }

        const query = input.toLowerCase();
        const filtered = SUGGESTIONS.filter(s =>
            s.text.toLowerCase().includes(query) ||
            s.description.toLowerCase().includes(query)
        ).slice(0, 5);

        setFilteredSuggestions(filtered);
        setSelectedIndex(0);
    }, [input]);

    if (!isVisible || filteredSuggestions.length === 0) return null;

    return (
        <div className="absolute bottom-full left-0 right-0 mb-1 bg-slate-800 border border-slate-600 rounded-lg shadow-lg overflow-hidden z-10">
            {filteredSuggestions.map((suggestion, index) => (
                <button
                    key={suggestion.text}
                    onClick={() => onSelect(suggestion.text)}
                    className={`w-full px-4 py-2 text-left flex items-center justify-between
            ${index === selectedIndex ? 'bg-slate-700' : 'hover:bg-slate-700/50'}
            transition-colors`}
                >
                    <span className="text-sm text-slate-200">{suggestion.text}</span>
                    <span className="text-xs text-slate-500">{suggestion.description}</span>
                </button>
            ))}
        </div>
    );
}
