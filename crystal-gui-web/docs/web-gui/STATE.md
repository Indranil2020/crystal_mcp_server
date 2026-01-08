# State Management Documentation

> Redux state architecture for `crystal-gui-web`.

---

## Store Architecture

```mermaid
graph TB
    subgraph Store["Redux Store"]
        subgraph Slices
            mcpSlice["mcpSlice"]
            structureSlice["structureSlice"]
            chatSlice["chatSlice"]
        end
    end
    
    subgraph Hooks
        useAppSelector
        useAppDispatch
    end
    
    subgraph Components
        App
        ChatPanel
        MolStarViewer
        StatusBar
    end
    
    Slices --> Hooks
    Hooks --> Components
```

---

## Slice Details

### mcpSlice

**File:** `store/mcpSlice.ts`

```mermaid
classDiagram
    class MCPState {
        status: 'disconnected' | 'connecting' | 'connected' | 'error'
        tools: Tool[]
        error: string | null
    }
```

| Action | Payload | Effect |
|--------|---------|--------|
| `setConnectionStatus` | status string | Update connection state |
| `setTools` | Tool[] | Store available tools |
| `setError` | string | Set error message |

---

### structureSlice

**File:** `store/structureSlice.ts`

```mermaid
classDiagram
    class StructureState {
        structures: Structure[]
        activeStructureId: string | null
        viewerSettings: ViewerSettings
        selection: Selection | null
        measurements: Measurement[]
    }
    
    class ViewerSettings {
        representation: RepresentationMode
        colorScheme: string
        showUnitCell: boolean
        measurementMode: string | null
    }
    
    class Selection {
        structureId: string
        atomIndices: number[]
        selectionMode: string
    }
    
    StructureState --> ViewerSettings
    StructureState --> Selection
```

| Action | Purpose |
|--------|---------|
| `addStructure` | Add new structure |
| `setActiveStructure` | Select for viewing |
| `updateViewerSettings` | Change representation/coloring |
| `setSelection` | Update selected atoms |
| `deleteAtoms` | Remove atoms from structure |
| `addMeasurement` | Add distance/angle measurement |

---

### chatSlice

**File:** `store/chatSlice.ts`

```mermaid
classDiagram
    class ChatState {
        messages: Message[]
        isStreaming: boolean
        currentToolCall: ToolCall | null
        selectedModel: string
    }
    
    class Message {
        id: string
        role: 'user' | 'assistant' | 'system'
        content: string
        timestamp: number
        toolCalls: ToolCall[]
    }
    
    class ToolCall {
        id: string
        name: string
        arguments: object
        status: 'pending' | 'running' | 'success' | 'error'
        result: any
    }
    
    ChatState --> Message
    Message --> ToolCall
```

---

## Data Flow

```mermaid
sequenceDiagram
    participant Component
    participant Dispatch
    participant Reducer
    participant State
    participant Selector
    
    Component->>Dispatch: dispatch(action)
    Dispatch->>Reducer: action
    Reducer->>State: new state
    State->>Selector: useAppSelector
    Selector->>Component: selected data
    Component->>Component: re-render
```

---

## Usage Examples

### Reading State

```typescript
import { useAppSelector } from '../store/hooks';

function MyComponent() {
  const { structures, activeStructureId } = useAppSelector(
    state => state.structure
  );
  const activeStructure = structures.find(s => s.id === activeStructureId);
  // ...
}
```

### Dispatching Actions

```typescript
import { useAppDispatch } from '../store/hooks';
import { addStructure, setActiveStructure } from '../store/structureSlice';

function MyComponent() {
  const dispatch = useAppDispatch();
  
  const handleNewStructure = (structure: Structure) => {
    dispatch(addStructure(structure));
    dispatch(setActiveStructure(structure.id));
  };
}
```
