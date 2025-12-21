#!/usr/bin/env python3
"""
Crystal MCP Server Test Runner
===============================

Unified test runner with multiple test profiles for different scenarios:
- Development (quick feedback)
- Pre-commit (medium coverage)
- Pre-publish (full coverage)
- CI/CD (parallel execution)
- Specific categories

Usage:
    ./run_tests.py --profile dev          # Quick development tests
    ./run_tests.py --profile pre-commit   # Pre-commit checks
    ./run_tests.py --profile full         # Complete test suite
    ./run_tests.py --profile ci           # CI/CD mode
    ./run_tests.py --category protocol    # Specific category
    ./run_tests.py --list-profiles        # Show all profiles
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Dict


class TestRunner:
    """Orchestrates test execution with different profiles"""
    
    def __init__(self):
        # Determine paths relative to this script
        self.script_path = Path(__file__).resolve()
        self.test_suit_dir = self.script_path.parent
        self.root_dir = self.test_suit_dir.parent.parent
        
        # Change working directory to root for consistent relative paths
        os.chdir(self.root_dir)
        
        # Test files relative to root
        self.test_files = {
            "comprehensive": "tests/testsuit/test_mcp_comprehensive.py",
            "operations": "tests/testsuit/test_operation_matrix.py",
            "e2e": "tests/testsuit/test_mcp_e2e.py",
            "all_operations": "tests/testsuit/test_all_operations.py",
            "new_modules": "tests/testsuit/test_new_modules.py",
            "scientific": "tests/testsuit/test_scientific_accuracy.py"
        }
        
        self.profiles = {
            "dev": {
                "name": "Development (Fast Feedback)",
                "description": "Quick tests for development iteration",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py::TestProtocolCompliance",
                    "tests/testsuit/test_mcp_comprehensive.py::TestCrystalGeneration",
                    "-v", "--tb=short", "-x"
                ],
                "estimate": "~30 seconds"
            },
            "pre-commit": {
                "name": "Pre-Commit (Medium Coverage)",
                "description": "Representative tests before committing",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py",
                    "-v", "--tb=short",
                    "-k", "Protocol or Generation or Export or Error"
                ],
                "estimate": "~2 minutes"
            },
            "full": {
                "name": "Full Test Suite (Complete Coverage)",
                "description": "All tests with coverage reporting",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py",
                    "tests/testsuit/test_operation_matrix.py",
                    "tests/testsuit/test_mcp_e2e.py",
                    "tests/testsuit/test_new_modules.py",
                    "-v", "--tb=short",
                    "--cov", "--cov-report=html", "--cov-report=term"
                ],
                "estimate": "~15-20 minutes"
            },
            "publish": {
                "name": "Pre-Publish Validation",
                "description": "Complete validation before publishing - ALL operations",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py",
                    "tests/testsuit/test_operation_matrix.py",
                    "tests/testsuit/test_mcp_e2e.py",
                    "tests/testsuit/test_all_operations.py",
                    "tests/testsuit/test_new_modules.py",
                    "tests/testsuit/test_scientific_accuracy.py",
                    "-v", "--tb=short",
                    "--cov", "--cov-report=html", "--cov-report=term",
                    "--junitxml=test-results.xml",
                    "--durations=20"
                ],
                "estimate": "~25-30 minutes"
            },
            "new-modules": {
                "name": "New Modules Only (December 2024)",
                "description": "Test only the newly added modules",
                "args": [
                    "tests/testsuit/test_new_modules.py",
                    "-v", "--tb=short"
                ],
                "estimate": "~5-10 minutes"
            },
            "scientific": {
                "name": "Scientific Accuracy Tests",
                "description": "Validate structures against reference data",
                "args": [
                    "tests/testsuit/test_scientific_accuracy.py",
                    "tests/testsuit/test_new_modules.py::TestScientificAccuracy",
                    "-v", "--tb=short"
                ],
                "estimate": "~5 minutes"
            },
            "ci": {
                "name": "CI/CD (Parallel)",
                "description": "Optimized for CI/CD with parallel execution",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py",
                    "tests/testsuit/test_operation_matrix.py",
                    "-v", "--tb=short",
                    "-n", "auto",
                    "--cov", "--cov-report=xml",
                    "--junitxml=test-results.xml"
                ],
                "estimate": "~5-8 minutes"
            },
            "protocol": {
                "name": "Protocol Compliance Only",
                "description": "Test MCP and JSON-RPC protocol compliance",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py::TestProtocolCompliance",
                    "tests/testsuit/test_mcp_comprehensive.py::TestToolDiscovery",
                    "tests/testsuit/test_mcp_e2e.py::TestProtocolBasics",
                    "-v"
                ],
                "estimate": "~15 seconds"
            },
            "performance": {
                "name": "Performance Tests",
                "description": "Performance and stress tests",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py::TestPerformance",
                    "tests/testsuit/test_operation_matrix.py::TestStressOperations",
                    "tests/testsuit/test_mcp_e2e.py::TestPerformance",
                    "-v", "--durations=10"
                ],
                "estimate": "~3 minutes"
            },
            "scientific": {
                "name": "Scientific Correctness",
                "description": "Validate physical and chemical accuracy",
                "args": [
                    "tests/testsuit/test_mcp_comprehensive.py::TestScientificCorrectness",
                    "tests/testsuit/test_operation_matrix.py::TestBulkStructures",
                    "-v"
                ],
                "estimate": "~5 minutes"
            }
        }
        
        self.categories = {
            "protocol": "Protocol compliance and tool discovery",
            "generation": "Crystal structure generation",
            "spacegroup": "Space group testing",
            "defects": "Defect generation",
            "surfaces": "Surface and slab generation",
            "2d": "2D materials",
            "transform": "Structure transformations",
            "export": "Export formats",
            "error": "Error handling",
            "performance": "Performance tests",
            "integration": "Integration workflows",
            "scientific": "Scientific correctness"
        }
    
    def verify_environment(self) -> bool:
        """Verify test environment is ready"""
        print(f"🔍 Working Directory: {os.getcwd()}")
        print("🔍 Verifying test environment...")
        
        # Check if server is built
        server_path = Path("dist/index.js")
        if not server_path.exists():
            print(f"❌ Server not built! Expected at {server_path.absolute()}")
            print("Run 'npm run build' first from the project root.")
            return False
        
        print(f"✓ Server found at {server_path}")
        
        # Check test files exist
        missing_files = []
        for name, file in self.test_files.items():
            if not Path(file).exists():
                print(f"⚠️  Test file not found: {file}")
                missing_files.append(file)
        
        if missing_files:
            print(f"❌ Cannot proceed without test files: {', '.join(missing_files)}")
            return False
            
        print("✓ Environment ready\n")
        return True
    
    def list_profiles(self):
        """List all available test profiles"""
        print("\n" + "="*70)
        print("AVAILABLE TEST PROFILES")
        print("="*70 + "\n")
        
        for key, profile in self.profiles.items():
            print(f"📋 {key.upper()}")
            print(f"   {profile['name']}")
            print(f"   {profile['description']}")
            print(f"   Estimated time: {profile['estimate']}")
            print()
    
    def list_categories(self):
        """List all available test categories"""
        print("\n" + "="*70)
        print("AVAILABLE TEST CATEGORIES")
        print("="*70 + "\n")
        
        for key, desc in self.categories.items():
            print(f"🏷️  {key}: {desc}")
        print()
    
    def run_profile(self, profile_name: str) -> int:
        """Run a specific test profile"""
        if profile_name not in self.profiles:
            print(f"❌ Unknown profile: {profile_name}")
            self.list_profiles()
            return 1
        
        profile = self.profiles[profile_name]
        
        print("\n" + "="*70)
        print(f"RUNNING: {profile['name']}")
        print("="*70)
        print(f"{profile['description']}")
        print(f"Estimated time: {profile['estimate']}")
        print("="*70 + "\n")
        
        cmd = ["pytest"] + profile["args"]
        
        print(f"Command: {' '.join(cmd)}\n")
        
        # Ensure tests/testsuit is in PYTHONPATH for cross-imports
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{self.test_suit_dir}:{env.get('PYTHONPATH', '')}"
        
        result = subprocess.run(cmd, env=env)
        return result.returncode
    
    def run_category(self, category: str) -> int:
        """Run tests for a specific category"""
        if category not in self.categories:
            print(f"❌ Unknown category: {category}")
            self.list_categories()
            return 1
        
        print(f"\n🏷️  Running {category} tests...")
        print(f"   {self.categories[category]}\n")
        
        cmd = [
            "pytest",
            "tests/testsuit/test_mcp_comprehensive.py",
            "tests/testsuit/test_operation_matrix.py",
            "tests/testsuit/test_mcp_e2e.py",
            "-v", "--tb=short",
            "-k", category
        ]
        
        # Ensure tests/testsuit is in PYTHONPATH for cross-imports
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{self.test_suit_dir}:{env.get('PYTHONPATH', '')}"
        
        result = subprocess.run(cmd, env=env)
        return result.returncode
    
    def run_custom(self, args: List[str]) -> int:
        """Run pytest with custom arguments"""
        cmd = ["pytest"] + args
        print(f"Running: {' '.join(cmd)}\n")
        result = subprocess.run(cmd)
        return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Crystal MCP Server Test Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --profile dev              # Quick development tests
  %(prog)s --profile full             # Complete test suite
  %(prog)s --profile ci               # CI/CD mode
  %(prog)s --category protocol        # Protocol tests only
  %(prog)s --list-profiles            # Show all profiles
  %(prog)s --custom test_mcp_comprehensive.py::TestPerformance -v
        """
    )
    
    parser.add_argument(
        "--profile", "-p",
        choices=["dev", "pre-commit", "full", "publish", "ci", "protocol", "performance", "scientific"],
        help="Run a predefined test profile"
    )
    
    parser.add_argument(
        "--category", "-c",
        choices=["protocol", "generation", "spacegroup", "defects", "surfaces", "2d", 
                 "transform", "export", "error", "performance", "integration", "scientific"],
        help="Run tests for a specific category"
    )
    
    parser.add_argument(
        "--list-profiles", "-l",
        action="store_true",
        help="List all available test profiles"
    )
    
    parser.add_argument(
        "--list-categories",
        action="store_true",
        help="List all available test categories"
    )
    
    parser.add_argument(
        "--custom",
        nargs=argparse.REMAINDER,
        help="Run pytest with custom arguments"
    )
    
    parser.add_argument(
        "--skip-verify",
        action="store_true",
        help="Skip environment verification"
    )
    
    args = parser.parse_args()
    
    runner = TestRunner()
    
    # Handle list options
    if args.list_profiles:
        runner.list_profiles()
        return 0
    
    if args.list_categories:
        runner.list_categories()
        return 0
    
    # Verify environment unless skipped
    if not args.skip_verify:
        if not runner.verify_environment():
            return 1
    
    # Run tests
    if args.profile:
        return runner.run_profile(args.profile)
    elif args.category:
        return runner.run_category(args.category)
    elif args.custom:
        return runner.run_custom(args.custom)
    else:
        # Default: run pre-commit profile
        print("No profile specified. Running 'pre-commit' profile by default.")
        print("Use --help for more options.\n")
        return runner.run_profile("pre-commit")


if __name__ == "__main__":
    sys.exit(main())
