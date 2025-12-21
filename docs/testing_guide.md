
# Testing Guide: Crystal MCP Server

This guide provides step-by-step instructions on how to verify and test the Crystal MCP Server as a standalone service.

## 1. Prerequisites
Ensure you have the following installed:
- **Node.js**: v18 or higher.
- **Python**: 3.9+ with `pymatgen`, `pyxtal`, and `numpy`.
- **NPM Dependencies**: Run `npm install` in the root directory.

---

## 2. Step-by-Step Verification

### Step 1: Build the Project
First, compile the TypeScript source into executable JavaScript in the `dist` folder.
```bash
npm run build
```

### Step 2: Automated E2E Test
I have provided a specialized test harness that simulates a real MCP client connecting via stdio. This is the most "proper" way to test the full stack (TS Bridge -> Python Router -> Generator).
```bash
python3 test_mcp_e2e.py
```
**Expectation**: You should see "SUCCESS: End-to-End Test Passed."

### Step 3: Using the MCP Inspector (Recommended)
The MCP Inspector is a visual tool provided by the Model Context Protocol team to debug servers.
```bash
npx @modelcontextprotocol/inspector node dist/index.js
```
1. Open the URL provided in the terminal (usually `http://localhost:5173`).
2. Click **"List Tools"**.
3. You should see 21 tools, including **`comprehensive_generate`**.
4. Select `comprehensive_generate` and provide arguments:
   ```json
   {
     "operation": "generate_from_spacegroup",
     "spacegroup": 227,
     "elements": ["Si"],
     "composition": [8]
   }
   ```
5. Click **"Run Tool"**.

### Step 4: Connecting to Claude Desktop
To use the server live in your Claude Desktop app:
1. Open your `claude_desktop_config.json` (Mac: `~/Library/Application Support/Claude/claude_desktop_config.json`, Windows: `%APPDATA%\Claude\claude_desktop_config.json`).
2. Add the following entry:
   ```json
   {
     "mcpServers": {
       "crystal-gen": {
         "command": "node",
         "args": ["/path/to/your/git/crystal-mcp-server/dist/index.js"]
       }
     }
   }
   ```
3. Restart Claude Desktop. You should see a hammer icon indicating the server is active.

---

## 3. Advanced Operational Testing
To verify the "super-tool" architecture (21 tools vs 228 operations):

### List All Categories
Run `comprehensive_generate` with:
```json
{ "operation": "list_all" }
```

### Explore a Category
Explore the `twist` category:
```json
{ "operation": "list_category", "category": "twist" }
```

### Chain Operations
1. Generate a structure.
2. Copy the resulting `structure` JSON.
3. Pass it to `export_vasp` or `generate_vacancy`.
