import os
import sys
import time
import json
import subprocess
import requests
import re
from playwright.sync_api import sync_playwright

# Configuration
SERVER_PORT = 8080
FRONTEND_URL = "http://localhost:5173"
# Path relative to project root
SERVER_SCRIPT = "crystal-gui-web/bridge/server.py"
LOG_FILE = "tests/molecule_test/server_output.log"

TEST_CASES = [
    # --- Group 1: Basic Linear Translations ---
    {"id": 101, "prompt": "generate 2 planar ptcda separated by 5 angstrom along x axis", "category": "Linear-X"},
    {"id": 102, "prompt": "generate 2 planar ptcda separated by 10 angstrom along x axis", "category": "Linear-X"},
    {"id": 103, "prompt": "generate 2 planar ptcda separated by 5 angstrom along y axis", "category": "Linear-Y"},
    {"id": 104, "prompt": "generate 3 planar ptcda separated by 4 angstrom along z axis", "category": "Linear-Z"},
    {"id": 105, "prompt": "generate 2 benzene molecules separated by 8 angstrom along x axis", "category": "Linear-X-Small"},

    # --- Group 2: Large Counts & Distances ---
    {"id": 201, "prompt": "generate 5 benzene molecules separated by 6 angstrom along x axis", "category": "Scale-Count"},
    {"id": 202, "prompt": "generate 10 water molecules separated by 3 angstrom along y axis", "category": "Scale-Count"},
    {"id": 203, "prompt": "generate 2 planar ptcda separated by 50 angstrom along x axis", "category": "Scale-Dist"},
    {"id": 204, "prompt": "generate 2 planar ptcda separated by 1.5 angstrom along x axis", "category": "Scale-Dist-Tiny"},

    # --- Group 3: Angular Translations (XY Plane) ---
    {"id": 301, "prompt": "generate 2 planar ptcda separated by 10 angstrom along 30 degree direction in xy plane", "category": "Angle-XY"},
    {"id": 302, "prompt": "generate 2 planar ptcda separated by 10 angstrom along 45 degree direction in xy plane", "category": "Angle-XY"},
    {"id": 303, "prompt": "generate 2 planar ptcda separated by 10 angstrom along 60 degree direction in xy plane", "category": "Angle-XY"},
    {"id": 304, "prompt": "generate 2 planar ptcda separated by 10 angstrom along 135 degree direction in xy plane", "category": "Angle-XY-Obtuse"},

    # --- Group 4: Angular Translations (YZ & XZ Planes) ---
    {"id": 401, "prompt": "generate 2 molecules separated by 15 angstrom along 45 degree direction in xz plane", "category": "Angle-XZ"},
    {"id": 402, "prompt": "generate 2 molecules separated by 15 angstrom along 30 degree direction in xz plane", "category": "Angle-XZ"},
    {"id": 403, "prompt": "generate 2 molecules separated by 15 angstrom along 60 degree direction in yz plane", "category": "Angle-YZ"},

    # --- Group 5: Rotations & Orientations ---
    {"id": 501, "prompt": "generate 2 ptcda separated by 8 angstrom along x, rotated 90 degrees around z", "category": "Rotation-Z"},
    {"id": 502, "prompt": "generate 2 ptcda separated by 8 angstrom along y, rotated 45 degrees around x", "category": "Rotation-X"},
    {"id": 503, "prompt": "generate 2 ptcda separated by 4 angstrom along z, rotated 180 degrees around y", "category": "Rotation-Y"},

    # --- Group 6: Hetero-Molecular Systems ---
    {"id": 601, "prompt": "generate a dimer of ptcda and ntcda separated by 7 angstrom along x axis", "category": "Hetero-Dimer"},
    {"id": 602, "prompt": "generate a trimer of ptcda, ntcda, ptcda separated by 5 angstrom along y axis", "category": "Hetero-Trimer"},
    {"id": 603, "prompt": "generate alternating chain of 4 molecules: benzene, thiophene, benzene, thiophene separated by 6 angstrom along x", "category": "Hetero-Alternating"},

    # --- Group 7: Complex/3D ---
    {"id": 701, "prompt": "generate 2 molecules separated by 5 angstrom along x and 2 angstrom along z", "category": "Complex-Slip"},
    {"id": 702, "prompt": "generate 2 molecules separated by 10 angstrom along the vector [1, 1, 1]", "category": "Complex-Vector"},
    {"id": 703, "prompt": "generate 2 planar ptcda separated by 15 angstrom along 60 degree direction in xz plane", "category": "Failure-Reproduction"},
]

def kill_server_on_port(port):
    """Kill any process running on the specified port."""
    print(f"Checking port {port}...")
    try:
        # Fuser is reliable on Linux
        subprocess.run(f"fuser -k {port}/tcp", shell=True, check=False)
        time.sleep(2)
    except Exception as e:
        print(f"Warning: Failed to kill process on port {port}: {e}")

def start_server():
    """Start the bridge server in background."""
    print("Starting server...")
    kill_server_on_port(SERVER_PORT)
    
    log_fp = open(LOG_FILE, "w")
    process = subprocess.Popen(
        ["/usr/bin/python3", SERVER_SCRIPT],
        stdout=log_fp,
        stderr=subprocess.STDOUT,
        cwd=os.getcwd()
    )
    
    # Wait for health check
    print("Waiting for server health check...")
    for i in range(20):
        try:
            resp = requests.get(f"http://localhost:{SERVER_PORT}/health")
            if resp.status_code == 200:
                print("Server started successfully.")
                return process, log_fp
        except:
            pass
        time.sleep(0.5)
    
    print("Failed to start server.")
    process.terminate()
    return None, log_fp

def get_last_mcp_call(log_path, start_offset):
    """Read the log file from start_offset and extract the last [MCP DEBUG] call."""
    try:
        with open(log_path, "r") as f:
            f.seek(start_offset)
            content = f.read()
            
        # Regex to capture JSON block after "Arguments: "
        # It handles multiline JSON
        matches = list(re.finditer(r"\[MCP DEBUG\] 📦 Arguments: ({[\s\S]*?})\n", content))
        if matches:
            last_match = matches[-1]
            try:
                return json.loads(last_match.group(1))
            except:
                return {"raw": last_match.group(1), "error": "JSON parse failed"}
    except Exception as e:
        return {"error": str(e)}
        
    return None

def run_tests():
    server_process, log_fp = start_server()
    if not server_process:
        sys.exit(1)
        
    results = []
    
    try:
        with sync_playwright() as p:
            print("Launching browser...")
            # Headless=True for CI/Background
            browser = p.chromium.launch(headless=True)
            page = browser.new_page()
            
            print(f"Navigating to {FRONTEND_URL}...")
            # Increase timeout for initial load if needed
            try:
                page.goto(FRONTEND_URL, timeout=10000)
                page.wait_for_load_state("networkidle", timeout=5000)
            except Exception as e:
                print(f"Navigation warning: {e}")
            
            # Locate chat input
            # Try multiple selectors
            input_selector = None
            for sel in ["textarea[placeholder*='Ask']", "textarea", "input[type='text']"]:
                if page.is_visible(sel):
                    input_selector = sel
                    break
            
            if not input_selector:
                print("Could not find chat input. Dumping page content snippet...")
                print(page.content()[:500])
                return

            print(f"Found input selector: {input_selector}")

            for test in TEST_CASES:
                print(f"Running Test {test['id']}: {test['prompt']}")
                
                # Mark log position
                log_fp.flush()
                log_offset = os.stat(LOG_FILE).st_size
                
                # Clear input first if needed (usually sending clears it)
                page.fill(input_selector, "")
                page.fill(input_selector, test["prompt"])
                page.press(input_selector, "Enter")
                
                # Wait for response
                # We expect an interaction. 
                # Wait for 4 seconds to be safe given local processing
                time.sleep(4)
                
                # Capture backend logs
                log_fp.flush()
                mcp_args = get_last_mcp_call(LOG_FILE, log_offset)
                
                status = "Captured" if mcp_args else "No Call Detected"
                print(f"  -> Status: {status}")
                if mcp_args and "error" in mcp_args:
                    print(f"  -> JSON Error: {mcp_args['error']}")
                
                result_entry = {
                    "id": test["id"],
                    "prompt": test["prompt"],
                    "llm_args": mcp_args,
                    "status": status
                }
                results.append(result_entry)
                
                # Wait a bit before next
                time.sleep(1)

    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()
    finally:
        print("Stopping server...")
        server_process.terminate()
        log_fp.close()
        
        output_path = "tests/molecule_test/test_report_raw.json"
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Done. Results saved to {output_path}")

if __name__ == "__main__":
    run_tests()
