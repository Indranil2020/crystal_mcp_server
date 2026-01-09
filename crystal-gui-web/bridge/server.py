"""
MCP Bridge Server - FastAPI HTTP/WebSocket bridge to MCP Server

Bridges browser frontend to the MCP server which uses stdio communication.
Provides both HTTP endpoints (fallback) and WebSocket (streaming) support.
Also provides direct chemistry endpoints using RDKit for 3D structure generation.
"""

import asyncio
import json
import subprocess
import sys
import logging
import importlib.util
from pathlib import Path
from typing import Optional, List

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Check RDKit availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None
logger.info(f"[Bridge] RDKit available: {RDKIT_AVAILABLE}")

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# Path to the MCP server
MCP_SERVER_PATH = Path(__file__).parent.parent.parent / "dist" / "index.js"

app = FastAPI(
    title="MCP Bridge Server",
    description="HTTP/WebSocket bridge for MCP Server communication",
    version="1.0.0",
)

# CORS for local development - allow all localhost ports (Vite may use different ports)
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"^http://(localhost|127\.0\.0\.1):\d+$",  # Any localhost port
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class McpProcess:
    """Manages a single MCP server subprocess."""
    
    def __init__(self):
        self.process: Optional[subprocess.Popen] = None
        self.request_id = 0
        self.initialized = False
        self._lock = asyncio.Lock()
    
    async def start(self):
        """Start the MCP server process."""
        if self.process is not None:
            return
        
        self.process = subprocess.Popen(
            ["node", str(MCP_SERVER_PATH)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    
    async def stop(self):
        """Stop the MCP server process."""
        if self.process:
            self.process.terminate()
            self.process = None
            self.initialized = False
    
    def _next_id(self) -> int:
        self.request_id += 1
        return self.request_id
    
    async def send_request(self, method: str, params: dict = None) -> dict:
        """Send a JSON-RPC request and wait for response."""
        async with self._lock:
            if not self.process:
                await self.start()
            
            request = {
                "jsonrpc": "2.0",
                "id": self._next_id(),
                "method": method,
                "params": params or {},
            }
            
            request_str = json.dumps(request) + "\n"
            self.process.stdin.write(request_str)
            self.process.stdin.flush()
            
            # Read response
            response_line = self.process.stdout.readline()
            if not response_line:
                raise HTTPException(status_code=500, detail="Empty response from MCP server")
            
            response = json.loads(response_line)
            
            if "error" in response:
                raise HTTPException(
                    status_code=400,
                    detail=response["error"].get("message", "MCP error")
                )
            
            return response.get("result", {})
    
    async def initialize(self):
        """Initialize the MCP connection."""
        if self.initialized:
            return
        
        await self.send_request("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {
                "name": "crystal-gui-bridge",
                "version": "1.0.0",
            },
        })
        
        # Send initialized notification
        if self.process:
            notification = json.dumps({
                "jsonrpc": "2.0",
                "method": "notifications/initialized",
            }) + "\n"
            self.process.stdin.write(notification)
            self.process.stdin.flush()
        
        self.initialized = True


# Global MCP process
mcp = McpProcess()


# Pydantic models
class InitializeRequest(BaseModel):
    protocolVersion: str = "2024-11-05"
    capabilities: dict = {}
    clientInfo: dict = {"name": "crystal-gui-web", "version": "1.0.0"}


class ToolCallRequest(BaseModel):
    name: str
    arguments: dict = {}


class SmilesTo3DRequest(BaseModel):
    """Request to convert SMILES to 3D structure with hydrogens."""
    smiles: str
    optimize: bool = True
    name: str = "Molecule"


# ============================================================
# Chemistry Endpoints (Direct RDKit, no MCP)
# ============================================================

def validate_chemistry(mol) -> dict:
    """
    Validate chemical correctness of a molecule.
    
    Checks:
    - All atoms have valid valency
    - No disconnected fragments (optional warning)
    - Kekulization is possible for aromatic systems
    
    Returns dict with 'valid', 'warnings', 'errors' keys.
    """
    result = {"valid": True, "warnings": [], "errors": []}
    
    if mol is None:
        result["valid"] = False
        result["errors"].append("Invalid molecule (could not parse)")
        return result
    
    # Check for multiple fragments (disconnected structures)
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        result["warnings"].append(f"Molecule has {len(frags)} disconnected fragments")
    
    # Try Kekulization (validates aromaticity)
    mol_copy = Chem.RWMol(mol)
    kekulize_ok = True
    kekulize_error = None
    # Kekulize can raise, so we check via sanitization flags instead
    san_result = Chem.SanitizeMol(mol_copy, catchErrors=True)
    if san_result != 0:
        kekulize_ok = False
        kekulize_error = str(san_result)
    
    if not kekulize_ok:
        result["warnings"].append(f"Kekulization issue: {kekulize_error}")
    
    # Check valency per atom
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        valence = atom.GetTotalValence()
        default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
        
        # Default valence can be a tuple for multi-valent atoms (e.g., S can be 2,4,6)
        if isinstance(default_valence, tuple):
            valid_valences = default_valence
        else:
            valid_valences = (default_valence,)
        
        # Allow some flexibility for charged atoms
        formal_charge = atom.GetFormalCharge()
        adjusted_valences = [v - formal_charge for v in valid_valences]
        
        if valence not in valid_valences and valence not in adjusted_valences:
            # This might be a radical or unusual bonding
            if atom.GetNumRadicalElectrons() > 0:
                result["warnings"].append(f"Atom {symbol}{atom.GetIdx()+1} appears to be a radical")
            else:
                result["warnings"].append(
                    f"Atom {symbol}{atom.GetIdx()+1} has unusual valence {valence} "
                    f"(expected one of {valid_valences})"
                )
    
    return result


def smiles_to_3d_structure(smiles: str, optimize: bool = True) -> dict:
    """
    Convert SMILES to 3D structure with explicit hydrogens.
    
    Uses RDKit to:
    1. Parse SMILES
    2. Add explicit hydrogens
    3. Generate 3D conformer (ETKDG)
    4. Optimize geometry (MMFF94/UFF)
    5. Validate chemical correctness
    
    Returns dict with atoms, coords, formula, validation results.
    """
    logger.info(f"[Chemistry] Converting SMILES to 3D: {smiles}")
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error(f"[Chemistry] Failed to parse SMILES: {smiles}")
        return {"success": False, "error": f"Invalid SMILES: {smiles}"}
    
    logger.info(f"[Chemistry] Parsed SMILES, heavy atoms: {mol.GetNumAtoms()}")
    
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    logger.info(f"[Chemistry] Added hydrogens, total atoms: {mol.GetNumAtoms()}")
    
    # Generate 3D conformer
    embed_params = AllChem.ETKDGv3()
    embed_params.randomSeed = 42
    embed_params.useRandomCoords = True
    
    conf_id = AllChem.EmbedMolecule(mol, embed_params)
    
    if conf_id == -1:
        # Fallback: try with random coordinates
        logger.warning("[Chemistry] ETKDG failed, trying with random coords")
        conf_id = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    
    if conf_id == -1:
        # Last resort: 2D coordinates
        logger.warning("[Chemistry] 3D embedding failed, using 2D coordinates")
        AllChem.Compute2DCoords(mol)
        conf_id = 0
    
    logger.info(f"[Chemistry] Conformer generated, conf_id: {conf_id}")
    
    # Optimize geometry
    if optimize and conf_id != -1:
        logger.info("[Chemistry] Optimizing geometry with MMFF94...")
        mmff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        if mmff_result != 0:
            logger.warning("[Chemistry] MMFF94 failed, trying UFF")
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    
    # Extract coordinates
    conformer = mol.GetConformer()
    atoms = []
    coords = []
    
    for i, atom in enumerate(mol.GetAtoms()):
        atoms.append(atom.GetSymbol())
        pos = conformer.GetAtomPosition(i)
        coords.append([round(pos.x, 4), round(pos.y, 4), round(pos.z, 4)])
    
    # Calculate properties
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mol_weight = Descriptors.MolWt(mol)
    canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    
    logger.info(f"[Chemistry] Formula: {formula}, Atoms: {len(atoms)}, MW: {mol_weight:.2f}")
    
    # Validate chemistry
    validation = validate_chemistry(mol)
    logger.info(f"[Chemistry] Validation: valid={validation['valid']}, warnings={len(validation['warnings'])}")
    
    # Generate SDF for Mol* (includes bonds)
    sdf_data = Chem.MolToMolBlock(mol)
    logger.info(f"[Chemistry] Generated SDF: {len(sdf_data)} bytes")
    
    return {
        "success": True,
        "smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "formula": formula,
        "molecular_weight": round(mol_weight, 4),
        "n_atoms": len(atoms),
        "n_heavy_atoms": mol.GetNumHeavyAtoms(),
        "atoms": atoms,
        "coords": coords,
        "sdf": sdf_data,  # SDF/MOL block with 3D coords and bonds
        "optimized": optimize,
        "validation": validation,
    }


@app.post("/chemistry/smiles-to-3d")
async def convert_smiles_to_3d(request: SmilesTo3DRequest):
    """
    Convert SMILES to 3D structure with explicit hydrogens.
    
    Uses RDKit to generate chemically valid 3D structures with:
    - Explicit hydrogen atoms
    - 3D conformer generation (ETKDG algorithm)
    - Geometry optimization (MMFF94 force field)
    - Chemical validity checking
    
    Returns SDF data that can be loaded directly into Mol*.
    """
    if not RDKIT_AVAILABLE:
        raise HTTPException(
            status_code=503,
            detail="RDKit not available on server. Install rdkit-pypi."
        )
    
    logger.info(f"[API] /chemistry/smiles-to-3d called with: {request.smiles[:50]}...")
    
    result = smiles_to_3d_structure(request.smiles, optimize=request.optimize)
    
    if not result.get("success"):
        raise HTTPException(status_code=400, detail=result.get("error", "Conversion failed"))
    
    result["name"] = request.name
    return result


# Health check
@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "ok", 
        "mcp_initialized": mcp.initialized,
        "rdkit_available": RDKIT_AVAILABLE,
    }


# MCP endpoints
@app.post("/mcp/initialize")
async def mcp_initialize(request: InitializeRequest):
    """Initialize MCP connection."""
    await mcp.initialize()
    return {"success": True, "initialized": True}


@app.get("/mcp/tools")
async def mcp_list_tools():
    """List available MCP tools."""
    if not mcp.initialized:
        await mcp.initialize()
    
    result = await mcp.send_request("tools/list", {})
    return {"tools": result.get("tools", [])}


@app.post("/mcp/call")
async def mcp_call_tool(request: ToolCallRequest):
    """Call an MCP tool."""
    if not mcp.initialized:
        await mcp.initialize()
    
    result = await mcp.send_request("tools/call", {
        "name": request.name,
        "arguments": request.arguments,
    })
    
    return result


# WebSocket for streaming (future enhancement)
@app.websocket("/mcp/ws")
async def mcp_websocket(websocket: WebSocket):
    """WebSocket endpoint for streaming MCP communication."""
    await websocket.accept()
    
    try:
        while True:
            data = await websocket.receive_json()
            
            method = data.get("method")
            params = data.get("params", {})
            
            if method == "tools/list":
                if not mcp.initialized:
                    await mcp.initialize()
                result = await mcp.send_request("tools/list", {})
                await websocket.send_json({"type": "tools", "data": result})
            
            elif method == "tools/call":
                if not mcp.initialized:
                    await mcp.initialize()
                result = await mcp.send_request("tools/call", params)
                await websocket.send_json({"type": "result", "data": result})
            
            else:
                await websocket.send_json({"type": "error", "message": f"Unknown method: {method}"})
    
    except WebSocketDisconnect:
        pass


# Startup/shutdown events
@app.on_event("startup")
async def startup():
    """Start MCP process on server startup."""
    await mcp.start()


@app.on_event("shutdown")
async def shutdown():
    """Stop MCP process on server shutdown."""
    await mcp.stop()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)
