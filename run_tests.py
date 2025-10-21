#!/usr/bin/env python3
"""
Test runner script for H-wave project.

This script provides a convenient way to run all tests with different configurations
and generate coverage reports.
"""

import sys
import os
import subprocess
import argparse
from pathlib import Path


def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}")
        return False


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="Run H-wave tests")
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    parser.add_argument("--integration", action="store_true", help="Run integration tests only")
    parser.add_argument("--coverage", action="store_true", help="Generate coverage report")
    parser.add_argument("--pytest", action="store_true", help="Use pytest instead of unittest")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--html", action="store_true", help="Generate HTML coverage report")
    
    args = parser.parse_args()
    
    # Change to project root directory
    project_root = Path(__file__).parent
    os.chdir(project_root)
    
    print("🧪 H-wave Test Runner")
    print(f"📁 Working directory: {project_root}")
    
    success = True
    
    if args.pytest:
        # Run with pytest
        cmd = ["python", "-m", "pytest"]
        
        if args.unit:
            cmd.extend(["tests/unit/", "-m", "unit"])
        elif args.integration:
            cmd.extend(["tests/", "-m", "integration"])
        else:
            cmd.append("tests/")
        
        if args.coverage:
            cmd.extend(["--cov=src/hwave", "--cov-report=term-missing"])
            if args.html:
                cmd.append("--cov-report=html:htmlcov")
        
        if args.verbose:
            cmd.append("-v")
        
        success &= run_command(cmd, "pytest tests")
    
    else:
        # Run with unittest
        if args.unit:
            cmd = ["python", "-m", "unittest", "discover", "tests/unit/", "-v"]
            success &= run_command(cmd, "unit tests")
        
        if args.integration or not args.unit:
            cmd = ["python", "-m", "unittest", "discover", "tests/", "-v"]
            success &= run_command(cmd, "integration tests")
    
    # Generate coverage report if requested
    if args.coverage and not args.pytest:
        cmd = ["python", "-m", "coverage", "run", "-m", "unittest", "discover", "tests/"]
        run_command(cmd, "coverage collection")
        
        cmd = ["python", "-m", "coverage", "report", "-m"]
        run_command(cmd, "coverage report")
        
        if args.html:
            cmd = ["python", "-m", "coverage", "html"]
            run_command(cmd, "HTML coverage report")
    
    # Summary
    print(f"\n{'='*60}")
    if success:
        print("🎉 All tests completed successfully!")
    else:
        print("❌ Some tests failed. Check the output above.")
        sys.exit(1)
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
