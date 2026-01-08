# Services Documentation

> Detailed documentation for the service layer in `crystal-gui-web`.

---

## Overview

```mermaid
graph LR
    App["App.tsx"]
    
    subgraph Services["services/"]
        mcp["mcp.ts"]
        llm["llm.ts"]
        tools["tools.ts"]
        index["index.ts"]
    end
    
    subgraph External["External Services"]
        Bridge["Python Bridge :8080"]
        Ollama["Ollama :11434"]
    end
    
    App --> index
    index --> mcp
    index --> llm
    index --> tools
    
    mcp --> Bridge
    llm --> Ollama
    tools --> mcp
```

---

## MCP Client (`mcp.ts`)

Handles communication with the Crystal MCP Server via the Python bridge.

### API

| Method | Purpose |
|--------|---------|
| `initialize()` | Establish connection to MCP |
| `listTools()` | Fetch available tools |
| `callTool(name, params)` | Execute a tool |
| `checkConnection()` | Verify bridge is accessible |

### Communication Flow

```mermaid
sequenceDiagram
    participant GUI
    participant MCP as mcp.ts
    participant Bridge as :8080
    participant Server as MCP Server
    
    GUI->>MCP: callTool("build_molecule", {...})
    MCP->>Bridge: POST /call/build_molecule
    Bridge->>Server: JSON-RPC call
    Server-->>Bridge: Result
    Bridge-->>MCP: JSON response
    MCP-->>GUI: Structure data
```

---

## LLM Client (`llm.ts`)

Handles communication with Ollama for natural language processing.

### API

| Method | Purpose |
|--------|---------|
| `checkConnection()` | Verify Ollama is running |
| `chat(messages, tools)` | Send chat with tool definitions |
| `streamChat(messages, callback)` | Streaming responses |

### Message Flow

```mermaid
sequenceDiagram
    participant User
    participant Chat
    participant LLM as llm.ts
    participant Ollama
    
    User->>Chat: "Build benzene"
    Chat->>LLM: chat([message], tools)
    LLM->>Ollama: POST /api/chat
    Ollama-->>LLM: Stream response
    LLM-->>Chat: Tool call detected
    Chat->>Chat: Execute tool
```

---

## Tool Orchestrator (`tools.ts`)

Bridges LLM tool calls with MCP execution.

### API

| Method | Purpose |
|--------|---------|
| `setTools(tools)` | Register available tools |
| `executeTool(name, params)` | Execute via MCP |
| `formatToolsForLLM()` | Format for Ollama |

### Orchestration Flow

```mermaid
graph LR
    LLM["LLM Response"]
    Parse["Parse Tool Call"]
    Map["Map to MCP Tool"]
    Execute["Execute via MCP"]
    Result["Return Result"]
    
    LLM --> Parse --> Map --> Execute --> Result
```

---

## Usage Example

```typescript
import { mcpClient, llmClient, toolOrchestrator } from './services';

// Initialize
await mcpClient.initialize();
const tools = await mcpClient.listTools();
toolOrchestrator.setTools(tools);

// Execute tool from LLM response
const result = await toolOrchestrator.executeTool('build_molecule', { 
  formula: 'H2O' 
});
```
