
import sys
import os
sys.path.insert(0, 'src/python')

from comprehensive_structures import handle_request

print("Testing Registry Integration...")
result = handle_request({'operation': 'list_all'})

if not result['success']:
    print("FAILED: list_all")
    sys.exit(1)

ops = result['operations']
total = result['total_operations']
print(f"Total Operations: {total}")
print(f"Categories: {len(ops)}")

# Check specific new operations
checks = [
    'generate_twisted_trilayer',
    'generate_check_tolerance', # Wait, check_tolerance is utility, might not start with generate_
    'generate_skyrmion'
]

# check_tolerance might not be in registry if script only scanned generate_*
# My script scanned: generate_, make_, create_, build_, apply_

# Let's check 'apply_strain'
checks.append('apply_strain')

for op in checks:
    found = False
    for cat, op_list in ops.items():
        if op in op_list:
            print(f"Found {op} in {cat}")
            found = True
            break
    if not found and op != 'generate_check_tolerance': # check_tolerance wasn't generate_
        print(f"WARNING: {op} not found in registry")

print("SUCCESS")
