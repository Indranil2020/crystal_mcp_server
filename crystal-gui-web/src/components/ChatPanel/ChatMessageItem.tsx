/**
 * ChatMessageItem - Enhanced message display with markdown and structure previews
 */

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import type { ChatMessage } from '../../types';
import { useAppSelector } from '../../store/hooks';

interface Props {
    message: ChatMessage;
}

export default function ChatMessageItem({ message }: Props) {
    const isUser = message.role === 'user';
    const isAssistant = message.role === 'assistant';
    const structures = useAppSelector(state => state.structure.structures);

    // Find linked structure if this message has a tool result
    const linkedStructure = message.toolResult?.structureId
        ? structures.find(s => s.id === message.toolResult?.structureId)
        : null;

    return (
        <div className={`flex ${isUser ? 'justify-end' : 'justify-start'}`}>
            <div
                className={`max-w-[90%] rounded-lg px-4 py-3 ${isUser
                        ? 'bg-blue-600 text-white'
                        : isAssistant
                            ? 'bg-slate-700 text-slate-100'
                            : 'bg-slate-600 text-slate-200'
                    }`}
            >
                {/* Role indicator for non-user messages */}
                {!isUser && (
                    <div className="text-xs text-slate-400 mb-1 capitalize flex items-center gap-2">
                        <span>{message.role === 'assistant' ? '🤖' : '⚙️'}</span>
                        <span>{message.role}</span>
                    </div>
                )}

                {/* Message content with markdown */}
                <div className="prose prose-sm prose-invert max-w-none">
                    <ReactMarkdown
                        remarkPlugins={[remarkGfm]}
                        components={{
                            // Style code blocks
                            code: ({ className, children, ...props }) => {
                                const isInline = !className;
                                return isInline ? (
                                    <code className="bg-slate-800 px-1 py-0.5 rounded text-blue-300" {...props}>
                                        {children}
                                    </code>
                                ) : (
                                    <code className="block bg-slate-900 p-2 rounded text-xs overflow-x-auto" {...props}>
                                        {children}
                                    </code>
                                );
                            },
                            // Style links
                            a: ({ children, ...props }) => (
                                <a className="text-blue-400 hover:underline" {...props}>{children}</a>
                            ),
                            // Style lists
                            ul: ({ children }) => (
                                <ul className="list-disc list-inside space-y-1">{children}</ul>
                            ),
                            ol: ({ children }) => (
                                <ol className="list-decimal list-inside space-y-1">{children}</ol>
                            ),
                        }}
                    >
                        {message.content}
                    </ReactMarkdown>
                </div>

                {/* Structure preview card */}
                {linkedStructure && (
                    <div className="mt-3 p-2 bg-slate-800/50 rounded-lg border border-slate-600">
                        <div className="flex items-center gap-2">
                            <span className="text-xl">🧬</span>
                            <div>
                                <div className="text-sm font-medium text-slate-200">
                                    {linkedStructure.name}
                                </div>
                                <div className="text-xs text-slate-400">
                                    {linkedStructure.data.atoms.length} atoms
                                    {linkedStructure.data.metadata?.formula && (
                                        <span className="ml-2">• {linkedStructure.data.metadata.formula}</span>
                                    )}
                                </div>
                            </div>
                        </div>
                    </div>
                )}

                {/* Tool result indicator */}
                {message.toolResult && !linkedStructure && (
                    <div
                        className={`mt-2 pt-2 border-t ${message.toolResult.success ? 'border-green-500/30' : 'border-red-500/30'
                            }`}
                    >
                        <div className="flex items-center gap-2 text-xs">
                            <span className={message.toolResult.success ? 'text-green-400' : 'text-red-400'}>
                                {message.toolResult.success ? '✓' : '✗'}
                            </span>
                            <span className="text-slate-400 font-mono">{message.toolResult.toolName}</span>
                        </div>
                    </div>
                )}

                {/* Timestamp */}
                <div className="text-xs text-slate-500 mt-2 text-right">
                    {new Date(message.timestamp).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                </div>
            </div>
        </div>
    );
}
