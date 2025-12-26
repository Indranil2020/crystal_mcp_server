# Direct Git Repository Connection

This MCP server supports direct connection from git repositories without manual installation. This guide shows how to configure various MCP clients to use the git-direct method.

## Benefits of Git-Direct Connection

- ✅ **No manual installation** - Just configure and run
- ✅ **Auto-updates** - Pull latest changes automatically
- ✅ **Zero maintenance** - Dependencies managed automatically
- ✅ **Cross-platform** - Works on Windows, Mac, Linux
- ✅ **Quick setup** - One configuration, ready to use

## How It Works

When you use `npx -y git+https://github.com/Indranil2020/crystal_mcp_server.git`, the following happens automatically:

1. **Download**: Repository is cloned to npx cache
2. **Install**: Node.js dependencies are installed
3. **Python Setup**: Python packages are installed via pip/pip3
4. **Build**: TypeScript is compiled to JavaScript
5. **Run**: MCP server starts on stdio transport

First run takes 1-2 minutes. Subsequent runs are instant (uses cached installation).

## Configuration Examples

### Claude Desktop

**macOS/Linux**
Edit: `~/.config/Claude/claude_desktop_config.json`

**Windows**
Edit: `%APPDATA%\Claude\claude_desktop_config.json`

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

After editing, restart Claude Desktop.

### Windsurf IDE

1. Open **Settings** > **MCP** (or **Connections**)
2. Click **Add Server**
3. Configure:
   - **Name**: `crystal-structure-generator`
   - **Command**: `npx`
   - **Arguments**: `-y git+https://github.com/Indranil2020/crystal_mcp_server.git`
   - **Environment**: `PYTHON_PATH=python`
4. Save and restart Windsurf

### Cline (VSCode Extension)

Add to your VSCode settings (`settings.json`):

```json
{
  "cline.mcpServers": {
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

### Continue.dev

Add to `.continue/config.json`:

```json
{
  "mcpServers": [
    {
      "name": "crystal-structure-generator",
      "command": "npx",
      "args": [
        "-y",
        "git+https://github.com/Indranil2020/crystal_mcp_server.git"
      ],
      "env": {
        "PYTHON_PATH": "python"
      }
    }
  ]
}
```

### Generic MCP Client

Most MCP clients support stdio transport. Use this general pattern:

```json
{
  "command": "npx",
  "args": ["-y", "git+https://github.com/Indranil2020/crystal_mcp_server.git"],
  "env": {"PYTHON_PATH": "python"}
}
```

## Using Specific Git Branches

To use a specific branch, tag, or commit:

### Specific Branch
```json
{
  "command": "npx",
  "args": [
    "-y",
    "git+https://github.com/Indranil2020/crystal_mcp_server.git#dev-branch"
  ]
}
```

### Specific Tag/Release
```json
{
  "command": "npx",
  "args": [
    "-y",
    "git+https://github.com/Indranil2020/crystal_mcp_server.git#v1.2.3"
  ]
}
```

### Specific Commit
```json
{
  "command": "npx",
  "args": [
    "-y",
    "git+https://github.com/Indranil2020/crystal_mcp_server.git#abc123f"
  ]
}
```

## Environment Variables

### PYTHON_PATH
Specifies which Python interpreter to use:

```json
"env": {
  "PYTHON_PATH": "python"     // Use default python
}
```

or

```json
"env": {
  "PYTHON_PATH": "python3"    // Explicitly use python3
}
```

or

```json
"env": {
  "PYTHON_PATH": "/path/to/venv/bin/python"  // Use virtual environment
}
```

### NODE_ENV (optional)
Set to production to suppress development warnings:

```json
"env": {
  "PYTHON_PATH": "python",
  "NODE_ENV": "production"
}
```

## Troubleshooting

### First Run Takes Long Time
This is normal. The first run downloads and builds everything. Subsequent runs use the cached version and start instantly.

### Python Dependencies Fail to Install
Make sure you have Python 3.8+ and pip installed:

```bash
python --version   # or python3 --version
pip --version      # or pip3 --version
```

If using a virtual environment, set `PYTHON_PATH` to point to it.

### npx Command Not Found
Make sure Node.js 18+ and npm are installed:

```bash
node --version
npm --version
```

### Clear npx Cache (Force Reinstall)
To force a fresh installation:

```bash
# Clear npx cache
npx clear-npx-cache

# Or remove specific package
rm -rf ~/.npm/_npx
```

### Permission Errors on Windows
Run your terminal/IDE as administrator, or install Node.js with "for all users" option.

### Python Packages Not Found at Runtime
The `prepare` script should install Python dependencies automatically. If it fails:

1. Check Python/pip are in PATH
2. Manually install: `pip install -r requirements.txt`
3. Check the server logs for error messages

## Performance Notes

### Cache Location
npx caches packages at:
- **Linux/Mac**: `~/.npm/_npx/`
- **Windows**: `%LOCALAPPDATA%\npm-cache\_npx\`

### First Run Timing
- Download repository: ~10 seconds
- Install Node.js deps: ~30 seconds
- Install Python deps: ~60 seconds
- Build TypeScript: ~10 seconds
- **Total**: ~1-2 minutes

### Subsequent Runs
- Start time: <1 second (uses cached build)

## Alternative: npm Package (Future)

Once published to npm registry, you can use:

```json
{
  "command": "npx",
  "args": ["-y", "crystal-mcp-server"]
}
```

This will be faster as pre-built packages will be available.

## Comparison: Git vs Manual Installation

| Feature | Git-Direct | Manual Install |
|---------|-----------|----------------|
| Setup time | 1-2 min (once) | 5 min |
| Updates | Change git hash | `git pull && npm install` |
| Disk space | Minimal (npx cache) | Full clone |
| Best for | End users | Developers |
| Customization | Limited | Full control |

## Next Steps

After connecting the server:
- See [API Reference](./api_reference.md) for available tools
- Check [Workflows](./workflows.md) for usage examples
- Read [Troubleshooting Guide](./troubleshooting_guide.md) if issues arise

## Support

If you encounter issues with git-direct connection:
- Check MCP client logs for error messages
- Verify Python and Node.js are installed
- Open an issue on GitHub with error details
