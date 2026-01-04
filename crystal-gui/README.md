# Crystal Structure Generator GUI

A Rust-based desktop application for crystal structure generation with:
- **Chat Interface**: Natural language interaction with LLM (via Ollama)
- **MCP Client**: Direct communication with crystal-mcp-server
- **3D Visualization**: Interactive crystal structure viewer (Materials Project style)
- **Structure Editing**: Modify and refine generated structures

## Features

### Two Operating Modes

1. **Standalone MCP Server**: Connect the crystal-mcp-server to any LLM or chat interface
2. **Complete Package**: This GUI with integrated LLM, visualization, and editing

### Capabilities

- Generate crystal structures with 27 specialized tools
- Interactive 3D visualization with rotation, zoom, and atom scaling
- Quick actions for common structures (Si Diamond, NaCl, Perovskite)
- Real-time chat with AI assistant for structure generation guidance

## Prerequisites

1. **Rust** (1.70+): https://rustup.rs/
2. **Node.js** (18+): For the MCP server
3. **Ollama** (optional): For LLM chat - https://ollama.ai/

## Installation

```bash
# Build the MCP server first
cd /path/to/crystal-mcp-server
npm install
npm run build

# Build the GUI
cd crystal-gui
cargo build --release
```

## Running

```bash
# Start Ollama (optional, for chat)
ollama serve

# Pull a small model (optional)
ollama pull llama3.2:1b

# Run the GUI
./target/release/crystal-gui
```

## Usage

### Connecting

1. Click **Connection → Connect MCP** to start the crystal-mcp-server
2. Click **Connection → Connect Ollama** to enable chat (requires Ollama running)

### Generating Structures

**Via Quick Actions:**
- Click "Generate Si (Diamond)" for silicon with space group 227
- Click "Generate NaCl (Rocksalt)" for sodium chloride
- Click "Generate Perovskite" for CaTiO3

**Via Tool Panel:**
1. Select a tool from the dropdown
2. Enter parameters as JSON
3. Click "Execute Tool"

**Via Chat:**
- Type natural language requests like "Create a silicon crystal with diamond structure"
- The AI will guide you through the process

### Visualization Controls

- **Drag**: Rotate structure
- **Scroll**: Zoom in/out
- **Atom Scale slider**: Adjust atom sizes
- **View menu**: Toggle unit cell and bonds

## Architecture

```
crystal-gui/
├── src/
│   ├── main.rs           # Entry point
│   ├── app.rs            # Main application state and UI
│   ├── mcp_client.rs     # MCP JSON-RPC client
│   ├── llm_client.rs     # Ollama API client
│   └── crystal_viewer.rs # 3D visualization
└── Cargo.toml
```

## Configuration

Settings are accessible via **File → Settings**:

- **MCP Server Path**: Path to `dist/index.js`
- **Ollama URL**: Default `http://localhost:11434`
- **Ollama Model**: Default `llama3.2:1b`

## Troubleshooting

### MCP Connection Failed
- Ensure the MCP server is built: `npm run build` in the parent directory
- Check the server path in Settings

### Ollama Connection Failed
- Ensure Ollama is running: `ollama serve`
- Check if the model is pulled: `ollama list`

### No Structure Displayed
- Check the tool output for errors
- Ensure the structure JSON is valid

## License

MIT
