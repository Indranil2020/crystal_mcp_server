# Connect Crystal MCP Server to LLM Agents (incl. Windsurf)

Step-by-step instructions for first-time users who want to run this server locally and plug it into any MCP-capable agent (e.g., Windsurf, Claude desktop, custom SDK clients).

## Quick Start: Direct Git Connection (Recommended)

**No manual installation required!** Use npx to run directly from the git repository:

### For Claude Desktop
Edit your Claude Desktop config file (usually at `~/Library/Application Support/Claude/claude_desktop_config.json` on Mac):

```json
{
  "mcpServers": {
    "crystal-structure-generator": {
      "command": "npx",
      "args": [
        "-y",
        "git+https://github.com/Indranil2020/crystal_mcp_server.git"
      ],
      "env": {
        "PYTHON_PATH": "python"
      }
    }
  }
}
```

### For Windsurf
1. Open Windsurf settings > MCP (or "Connections")
2. Add a new server named `crystal-structure-generator`
3. Command: `npx`
4. Args: `-y git+https://github.com/Indranil2020/crystal_mcp_server.git`
5. Environment variables: `PYTHON_PATH=python`
6. Save and restart Windsurf

### First Run
The first time the server starts, it will automatically:
- Download the repository
- Install Node.js dependencies
- Install Python dependencies
- Build the TypeScript code

This may take 1-2 minutes. Subsequent runs will be instant.

---

## Alternative: Manual Installation

If you prefer to install manually or need more control:

## 1) Prerequisites
- Node.js 18+ and npm (`node -v`, `npm -v`)
- Python 3.9+ with `pip`
- Git (to download the repo)

## 2) Download the server
```bash
git clone https://github.com/Indranil2020/crystal_mcp_server.git
cd crystal-mcp-server
```

## 3) Install dependencies
```bash
# Node/TypeScript pieces
npm install

# Python side (use a venv if you like)
pip install -r requirements.txt
```

## 4) Build the TypeScript
```bash
npm run build    # creates dist/index.js (the MCP entrypoint)
```

## 5) Quick smoke test (optional)
In one terminal:
```bash
node dist/index.js
```
You should see a message like "Crystal Structure Generator MCP Server running on stdio". Leave this running while you test.

## 6) Wire it into an MCP-aware client
MCP clients spawn the server via a command and talk over stdio. Point your client at `node dist/index.js` in the repo directory.

### Generic MCP config shape
Most clients use a JSON/YAML entry similar to:
```json
{
  "mcpServers": {
    "crystal-mcp": {
      "command": "node",
      "args": ["dist/index.js"],
      "cwd": "/absolute/path/to/crystal-mcp-server"
    }
  }
}
```
- `command`: the executable (Node)
- `args`: the script to run
- `cwd`: repo root so the server can find its Python scripts

### Windsurf
1. Open Windsurf settings > MCP (or “Connections” depending on your version).
2. Add a new server named `crystal-mcp`.
3. Command: `node`
4. Args: `dist/index.js`
5. Working directory: absolute path to your `crystal-mcp-server` folder.
6. Save and reconnect; Windsurf should list the tools exposed by this server.

### Claude desktop (example)
If your Claude desktop uses an `mcpServers` config file, add:
```json
"crystal-mcp": {
  "command": "node",
  "args": ["dist/index.js"],
  "cwd": "/absolute/path/to/crystal-mcp-server"
}
```
Restart Claude to pick it up.

## 7) Troubleshooting tips
- If the server logs “Python environment check failed”, ensure Python 3.9+ is on your PATH and rerun `pip install -r requirements.txt`.
- If a client cannot find tools, confirm you built with `npm run build` and that the client’s `cwd` points to this repo.
- Keep only one active instance; MCP clients will start/stop their own copy when configured correctly.
