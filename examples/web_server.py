#!/usr/bin/env python3
"""
Crystal MCP Server - Web Interface (No Dependencies)

A pure Python web interface using only standard library.
Connects to the Crystal MCP Server via stdio using REAL MCP protocol.
"""

import subprocess
import json
import sys
import os
import threading
import queue
import time
import re
from pathlib import Path
from http.server import HTTPServer, BaseHTTPRequestHandler

PROJECT_ROOT = Path(__file__).parent.parent

mcp_client = None


class MCPWebClient:
    def __init__(self):
        self.process = None
        self.request_id = 0
        self.response_queue = queue.Queue()
        self.running = False
        self.initialized = False

    def start(self):
        server_path = PROJECT_ROOT / "dist" / "index.js"
        if not server_path.exists():
            raise FileNotFoundError(f"Server not built: {server_path}")

        env = os.environ.copy()
        env["PYTHON_PATH"] = sys.executable

        self.process = subprocess.Popen(
            ["node", str(server_path)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=str(PROJECT_ROOT), env=env, text=True, bufsize=1
        )
        self.running = True
        threading.Thread(target=self._read_stderr, daemon=True).start()
        threading.Thread(target=self._read_stdout, daemon=True).start()
        time.sleep(2)
        self._initialize()
        return True

    def _read_stderr(self):
        while self.running and self.process:
            line = self.process.stderr.readline()
            if line:
                print(f"[MCP] {line.rstrip()}")
            elif self.process.poll() is not None:
                break

    def _read_stdout(self):
        while self.running and self.process:
            line = self.process.stdout.readline()
            if line:
                try:
                    self.response_queue.put(json.loads(line.strip()))
                except json.JSONDecodeError:
                    pass
            elif self.process.poll() is not None:
                break

    def _send(self, method, params=None):
        self.request_id += 1
        request = {"jsonrpc": "2.0", "id": self.request_id, "method": method}
        if params:
            request["params"] = params
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        try:
            return self.response_queue.get(timeout=120)
        except queue.Empty:
            return {"error": "Timeout"}

    def _initialize(self):
        response = self._send("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {"roots": {"listChanged": True}},
            "clientInfo": {"name": "crystal-web-client", "version": "1.0.0"}
        })
        if "result" in response:
            self.process.stdin.write(json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}) + "\n")
            self.process.stdin.flush()
            self.initialized = True
            print("[MCP] Server initialized")
        return response

    def call_tool(self, tool_name, arguments):
        print(f"[MCP] Calling tool: {tool_name}")
        return self._send("tools/call", {"name": tool_name, "arguments": arguments})

    def stop(self):
        self.running = False
        if self.process:
            self.process.terminate()


# Common elements for parsing
ELEMENTS = {
    'h': 'H', 'he': 'He', 'li': 'Li', 'be': 'Be', 'b': 'B', 'c': 'C', 'n': 'N', 'o': 'O', 'f': 'F', 'ne': 'Ne',
    'na': 'Na', 'mg': 'Mg', 'al': 'Al', 'si': 'Si', 'p': 'P', 's': 'S', 'cl': 'Cl', 'ar': 'Ar',
    'k': 'K', 'ca': 'Ca', 'sc': 'Sc', 'ti': 'Ti', 'v': 'V', 'cr': 'Cr', 'mn': 'Mn', 'fe': 'Fe', 'co': 'Co', 'ni': 'Ni',
    'cu': 'Cu', 'zn': 'Zn', 'ga': 'Ga', 'ge': 'Ge', 'as': 'As', 'se': 'Se', 'br': 'Br', 'kr': 'Kr',
    'rb': 'Rb', 'sr': 'Sr', 'y': 'Y', 'zr': 'Zr', 'nb': 'Nb', 'mo': 'Mo', 'ru': 'Ru', 'rh': 'Rh', 'pd': 'Pd', 'ag': 'Ag',
    'cd': 'Cd', 'in': 'In', 'sn': 'Sn', 'sb': 'Sb', 'te': 'Te', 'i': 'I', 'xe': 'Xe',
    'cs': 'Cs', 'ba': 'Ba', 'la': 'La', 'ce': 'Ce', 'pt': 'Pt', 'au': 'Au', 'pb': 'Pb', 'bi': 'Bi',
    'silicon': 'Si', 'carbon': 'C', 'oxygen': 'O', 'nitrogen': 'N', 'iron': 'Fe', 'gold': 'Au',
    'silver': 'Ag', 'copper': 'Cu', 'zinc': 'Zn', 'titanium': 'Ti', 'barium': 'Ba', 'sodium': 'Na', 'chlorine': 'Cl'
}

# Space group multiplicities (minimum atoms needed for common space groups)
SPACE_GROUP_MULTIPLICITIES = {
    1: 1, 2: 2, 14: 4, 62: 4, 136: 2, 139: 2, 141: 4, 166: 3, 186: 2, 194: 2,
    216: 4, 221: 1, 225: 4, 227: 8, 229: 2, 230: 8
}


def parse_query(query):
    """Parse natural language query to extract element and space group."""
    q = query.lower().strip()

    # Predefined shortcuts for common structures
    if "silicon" in q and "diamond" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["Si"] * 8, "space_group": 227, "seed": 42}}
    elif "nacl" in q or "rocksalt" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["Na"] * 4 + ["Cl"] * 4, "space_group": 225, "seed": 42}}
    elif "perovskite" in q or "batio3" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["Ba", "Ti", "O", "O", "O"], "space_group": 221, "seed": 42}}
    elif "graphite" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["C"] * 4, "space_group": 194, "seed": 42}}
    elif "hbn" in q or "boron nitride" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["B", "B", "N", "N"], "space_group": 194, "seed": 42}}
    elif "iron" in q and "bcc" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["Fe", "Fe"], "space_group": 229, "seed": 42}}
    elif "gold" in q and not any(c.isdigit() for c in q):
        return {"tool": "generate_crystal", "params": {"composition": ["Au"] * 4, "space_group": 225, "seed": 42}}
    elif ("zno" in q or "zinc oxide" in q) and not any(c.isdigit() for c in q):
        return {"tool": "generate_crystal", "params": {"composition": ["Zn", "Zn", "O", "O"], "space_group": 186, "seed": 42}}
    elif "tio2" in q or "rutile" in q:
        return {"tool": "generate_crystal", "params": {"composition": ["Ti", "Ti", "O", "O", "O", "O"], "space_group": 136, "seed": 42}}

    # Flexible parsing: Extract element and space group number
    # Look for space group number (1-230)
    space_group = None
    sg_match = re.search(r'(?:space\s*group|sg|spg|#)?\s*(\d{1,3})\b', q)
    if sg_match:
        sg_num = int(sg_match.group(1))
        if 1 <= sg_num <= 230:
            space_group = sg_num

    # Look for element symbols or names
    element = None
    words = re.findall(r'[a-zA-Z]+', q)
    for word in words:
        word_lower = word.lower()
        if word_lower in ELEMENTS:
            element = ELEMENTS[word_lower]
            break

    # If we found both element and space group, create the structure
    if element and space_group:
        # Determine multiplicity based on space group
        multiplicity = SPACE_GROUP_MULTIPLICITIES.get(space_group, 4)
        composition = [element] * multiplicity
        return {"tool": "generate_crystal", "params": {"composition": composition, "space_group": space_group, "seed": 42}}

    # If only element found, use a default space group
    if element:
        return {"tool": "generate_crystal", "params": {"composition": [element] * 4, "space_group": 225, "seed": 42}}

    return None


def make_xyz(structure):
    atoms = structure.get("atoms", [])
    lines = [str(len(atoms)), "Crystal"]
    for a in atoms:
        x, y, z = a.get("cartesian", [0, 0, 0])
        lines.append(f"{a['element']}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(lines)


HTML_PAGE = '''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Crystal MCP Server</title>
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:system-ui,sans-serif;background:#1a1a2e;color:#fff;min-height:100vh}
.hdr{background:linear-gradient(135deg,#667eea,#764ba2);padding:20px;text-align:center}
.hdr h1{font-size:24px}.hdr p{opacity:.9;font-size:14px;margin-top:5px}
.main{display:flex;height:calc(100vh - 80px);padding:20px;gap:20px}
.pnl{background:rgba(255,255,255,.05);border-radius:12px;padding:20px}
.left{flex:1;display:flex;flex-direction:column}.right{width:450px}
.ex{margin-bottom:15px}
.ex button{background:rgba(255,255,255,.1);border:none;color:#fff;padding:8px 12px;margin:3px;border-radius:15px;cursor:pointer}
.ex button:hover{background:rgba(255,255,255,.2)}
#msgs{flex:1;overflow-y:auto;padding:10px 0}
.m{margin:10px 0;padding:12px 16px;border-radius:12px;max-width:85%}
.m.u{background:linear-gradient(135deg,#667eea,#764ba2);margin-left:auto;text-align:right}
.m.b{background:rgba(255,255,255,.1)}
.ir{display:flex;gap:10px;margin-top:15px}
.ir input{flex:1;padding:12px 18px;border:none;border-radius:25px;font-size:15px;background:rgba(255,255,255,.1);color:#fff}
.ir input::placeholder{color:rgba(255,255,255,.5)}
.ir button{padding:12px 25px;border:none;border-radius:25px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;cursor:pointer;font-size:15px}
#vw{height:350px;background:#fff;border-radius:8px}
.ps{margin-top:15px}.ps h3{color:#667eea;margin-bottom:10px}
.p{display:flex;justify-content:space-between;padding:6px 0;border-bottom:1px solid rgba(255,255,255,.1)}
.pl{opacity:.7}
</style>
</head>
<body>
<div class="hdr"><h1>Crystal MCP Server</h1><p>Real MCP Protocol - Natural Language to Crystal Structure</p></div>
<div class="main">
<div class="pnl left">
<div class="ex">
<button id="b1">Silicon Diamond</button>
<button id="b2">NaCl Rocksalt</button>
<button id="b3">BaTiO3</button>
<button id="b4">Graphite</button>
<button id="b5">Iron BCC</button>
<button id="b6">ZnO</button>
</div>
<div id="msgs"></div>
<div class="ir"><input id="q" placeholder="Describe the crystal structure..."><button id="go">Generate</button></div>
</div>
<div class="pnl right"><div id="vw"></div><div class="ps" id="ps"><h3>Properties</h3><p style="opacity:.5">Generate a structure</p></div></div>
</div>
<script>
var viewer;
document.addEventListener('DOMContentLoaded',function(){
viewer=$3Dmol.createViewer(document.getElementById('vw'),{backgroundColor:'white'});
viewer.render();
addM('Welcome! Type a request like "Generate silicon diamond structure" to create crystals via MCP.',0);

document.getElementById('b1').onclick=function(){document.getElementById('q').value='Generate silicon diamond structure'};
document.getElementById('b2').onclick=function(){document.getElementById('q').value='Create NaCl rocksalt'};
document.getElementById('b3').onclick=function(){document.getElementById('q').value='BaTiO3 perovskite'};
document.getElementById('b4').onclick=function(){document.getElementById('q').value='Graphite carbon'};
document.getElementById('b5').onclick=function(){document.getElementById('q').value='Iron BCC'};
document.getElementById('b6').onclick=function(){document.getElementById('q').value='ZnO wurtzite'};
document.getElementById('go').onclick=send;
document.getElementById('q').onkeypress=function(e){if(e.key==='Enter')send()};
});

function addM(t,u){
var d=document.createElement('div');
d.className='m '+(u?'u':'b');
d.appendChild(document.createTextNode(t));
document.getElementById('msgs').appendChild(d);
document.getElementById('msgs').scrollTop=99999;
return d;
}

function updV(xyz,s){
viewer.removeAllModels();
viewer.addModel(xyz,'xyz');
viewer.setStyle({},{stick:{radius:.15},sphere:{scale:.25}});
viewer.addUnitCell(viewer.getModel());
viewer.zoomTo();
viewer.render();

var ps=document.getElementById('ps');
while(ps.firstChild)ps.removeChild(ps.firstChild);
var h=document.createElement('h3');
h.appendChild(document.createTextNode(s.metadata.formula));
ps.appendChild(h);

var props=[['Space Group',s.space_group.symbol+' (#'+s.space_group.number+')'],
['System',s.space_group.crystal_system],['Atoms',s.metadata.natoms],
['Volume',s.lattice.volume.toFixed(2)+' A^3'],
['a',s.lattice.a.toFixed(3)+' A'],['b',s.lattice.b.toFixed(3)+' A'],['c',s.lattice.c.toFixed(3)+' A']];

for(var i=0;i<props.length;i++){
var row=document.createElement('div');row.className='p';
var lbl=document.createElement('span');lbl.className='pl';lbl.appendChild(document.createTextNode(props[i][0]));
var val=document.createElement('span');val.appendChild(document.createTextNode(props[i][1]));
row.appendChild(lbl);row.appendChild(val);ps.appendChild(row);
}
}

function send(){
var q=document.getElementById('q').value.trim();
if(!q)return;
addM(q,1);
document.getElementById('q').value='';
var ld=addM('Sending MCP request...',0);

fetch('/api/generate',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({query:q})})
.then(function(r){return r.json()})
.then(function(d){
ld.parentNode.removeChild(ld);
if(d.success){
addM('Generated '+d.structure.metadata.formula+' via MCP tool: '+d.tool,0);
updV(d.xyz,d.structure);
}else{addM('Error: '+d.error,0);}
})
.catch(function(e){ld.parentNode.removeChild(ld);addM('Error: '+e.message,0);});
}
</script>
</body>
</html>'''


class Handler(BaseHTTPRequestHandler):
    def log_message(self, f, *a):
        pass

    def do_GET(self):
        if self.path == '/' or self.path == '/index.html':
            self.send_response(200)
            self.send_header('Content-Type', 'text/html')
            self.end_headers()
            self.wfile.write(HTML_PAGE.encode())
        else:
            self.send_error(404)

    def do_POST(self):
        if self.path == '/api/generate':
            length = int(self.headers['Content-Length'])
            data = json.loads(self.rfile.read(length).decode())
            query = data.get('query', '')
            parsed = parse_query(query)

            if not parsed:
                self.resp({"success": False, "error": f"Could not understand: '{query}'"})
                return

            if not mcp_client or not mcp_client.initialized:
                self.resp({"success": False, "error": "MCP Server not connected"})
                return

            try:
                response = mcp_client.call_tool(parsed["tool"], parsed["params"])
                if "result" in response and "content" in response["result"]:
                    for item in response["result"]["content"]:
                        if item.get("type") == "text":
                            text = item["text"]
                            # Parse JSON from <json-data> tags (new format)
                            if "<json-data>" in text:
                                start = text.find("<json-data>") + 11
                                end = text.find("</json-data>")
                                result = json.loads(text[start:end].strip())
                                if result.get("success"):
                                    self.resp({"success": True, "structure": result["structure"], "xyz": make_xyz(result["structure"]), "tool": parsed["tool"]})
                                    return
                            # Fallback: parse JSON from ```json blocks (legacy format)
                            elif "```json" in text:
                                start = text.find("```json") + 7
                                end = text.find("```", start)
                                result = json.loads(text[start:end].strip())
                                if result.get("success"):
                                    self.resp({"success": True, "structure": result["structure"], "xyz": make_xyz(result["structure"]), "tool": parsed["tool"]})
                                    return
                self.resp({"success": False, "error": "Failed to parse response"})
            except Exception as e:
                self.resp({"success": False, "error": str(e)})
        else:
            self.send_error(404)

    def resp(self, data):
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(data).encode())


def main():
    global mcp_client
    print("=" * 60)
    print("   Crystal MCP Server - Web Interface")
    print("   Using REAL MCP Protocol")
    print("=" * 60)

    print("\nStarting MCP Server...")
    mcp_client = MCPWebClient()

    try:
        mcp_client.start()
        print("[OK] MCP Server connected")
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        print("Run 'npm run build' first!")
        return

    print("\n" + "=" * 60)
    print("   Open http://localhost:8080 in your browser")
    print("=" * 60 + "\n")

    server = HTTPServer(('0.0.0.0', 8080), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down...")
    finally:
        mcp_client.stop()
        server.shutdown()


if __name__ == '__main__':
    main()
