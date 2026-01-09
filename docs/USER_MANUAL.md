# Crystal GUI Web - User Manual

## Introduction
Crystal GUI Web is a modern, advanced molecular visualization workspace powered by AI. It combines the power of large language models with rigorous scientific tools (RDKit, ASE, Pymatgen) to assist researchers in designing and analyzing chemical structures.

## Getting Started

### Accessing the Application
Open the application in your web browser (default: `http://localhost:5173`).
On first launch, you will be greeted by the **Welcome Tutorial**, which provides a quick overview of controls.

### The Interface
The workspace is divided into three main areas:
1. **Left Panel (Chat & History):** Interact with the AI agent, see your conversation history, and view active tool executions.
2. **Center Panel (3D Viewer):** The main visualization area using Mol* to display molecules and crystals in 3D.
3. **Right Panel (2D Editor):** A Kekule.js-powered editor for drawing chemical structures manually.

## Features & Workflows

### 1. Generating Molecules & Crystals
Use the chat panel to request structures using natural language.
- **Molecules:** "Generate aspirin", "Show me a caffeine molecule", "Create a water molecule"
- **Crystals:** "Create a gold crystal structure", "Show me silicon in diamond structure"
- **Clusters:** "Generate a 13-atom gold nanocluster"

The AI will parse your request, choose the appropriate tool (internal database, PubChem, or algorithmic generation), and display the result.

### 2. 3D Visualization Controls (Mol*)
Interact with the 3D structure using mouse and keyboard.

**Mouse Controls:**
- **Rotate:** Left Click + Drag
- **Zoom:** Scroll Wheel
- **Pan:** Right Click + Drag
- **Select:** Click on an atom (when Selection Mode is active)

**Keyboard Shortcuts:**
- `1` - `5`: Change representation style (Ball & Stick, Spacefill, Licorice, Wireframe, Surface)
- `C`: Toggle Crystal Unit Cell box
- `Deleting`: Press `Delete` or `Backspace` to remove selected atoms
- `R`: Reset Camera
- `M`: Toggle Measurement Controls panel
- `Ctrl + S`: Save Screenshot
- `Ctrl + Z`: Undo last action
- `Ctrl + Y` (or `Ctrl + Shift + Z`): Redo

### 3. Measurements
To measure distances, angles, or dihedrals:
1. Click the **Measure** button in the toolbar (or press `M`).
2. This opens the Mol* controls panel on the right side of the viewer.
3. Use the **Selection Mode** button (or `Select` in toolbar) to pick atoms.
4. Click two atoms to measure distance.
5. Click three atoms for an angle.
6. Click four for a dihedral.
7. Use the **CSV** button in the toolbar to export your measurements.

### 4. 2D to 3D Workflow
1. Open the **2D Editor** (toggle via the Layout button in status bar or resize the right panel).
2. Draw your molecule using the tool palette.
3. Click **"Push to 3D"**.
4. The backend will optimize the geometry (add hydrogens, minimize energy) and load it into the 3D viewer.

### 5. Undo/Redo System
The application maintains a history of your actions (50 steps).
- Accidentally deleted an atom? Press `Ctrl + Z`.
- Want to go back? Press `Ctrl + Y`.
History tracks: structure additions, deletions, measurements, and clearing.

### 6. Exporting Data
- **Screenshots:** Click the `IMG` button or press `Ctrl + S`.
- **Measurements:** Click `CSV` to download, or `Right Click` -> `CSV` (or `Shift + Click`) to copy to clipboard.
- **Session:** Click `Save` to store your current view state in the browser (local storage). Click `Load` to restore it.

## Troubleshooting

### "Optimization Failed"
If the 2D-to-3D conversion fails, try simplifying the structure or checking for valency violations in your 2D drawing.

### "WebGL Context Lost" or Viewer Issues
If the 3D viewer behaves strangely:
- Click the **Retry** button in the error overlay if visible.
- Refresh the page (Session state is not automatically saved on refresh unless you clicked "Save").
- Ensure your browser supports WebGL 2.0 (Chrome/Edge/Firefox suggested).

### Missing Dependencies
The backend uses `dependency_check.py` to manage tools. If a specific feature (e.g., crystal generation via `pyxtal`) is unavailable, it means the Python package is missing on the server. The basic features (RDKit, internal database) will still work.
