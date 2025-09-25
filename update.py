#!/usr/bin/env python3
"""
Update script to process dataset.csv and generate all necessary JSON files and images.

This script runs:
1. dataset2electrolyte.py - Processes dataset.csv to create info/electrolytes.json
2. dataset2molecule.py - Processes electrolytes.json to create info/molecules.json and molecule images

Usage: python3 update.py
"""

import subprocess
import sys
from pathlib import Path

def run_script(script_name):
    """Run a Python script and handle errors"""
    print(f"\n{'='*50}")
    print(f"Running {script_name}...")
    print(f"{'='*50}")
    
    try:
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=False, 
                              text=True, 
                              check=True)
        print(f"✅ {script_name} completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running {script_name}: {e}")
        return False
    except FileNotFoundError:
        print(f"❌ Script {script_name} not found!")
        return False

def main():
    """Main update process"""
    print("🚀 Starting dataset update process...")
    
    # Ensure info directory exists
    info_dir = Path("info")
    info_dir.mkdir(exist_ok=True)
    
    molecules_dir = info_dir / "molecules"
    molecules_dir.mkdir(exist_ok=True)
    
    print(f"📁 Ensured directory structure exists: {info_dir}")
    
    # Step 1: Generate electrolytes.json from dataset.csv
    if not run_script("dataset2electrolyte.py"):
        print("❌ Failed to generate electrolytes.json. Aborting.")
        sys.exit(1)
    
    # Step 2: Generate molecules.json and images from electrolytes.json
    if not run_script("dataset2molecule.py"):
        print("❌ Failed to generate molecules.json and images. Aborting.")
        sys.exit(1)
    
    print(f"\n{'='*50}")
    print("🎉 Update process completed successfully!")
    print(f"{'='*50}")
    print("\nGenerated files:")
    print("  📄 info/electrolytes.json")
    print("  📄 info/molecules.json")
    print("  🖼️  info/molecules/ (SVG images)")
    print("  🗺️  info/molecules/smiles_mapping.js")
    print("\n🌐 Your website is ready at http://localhost:8000")

if __name__ == "__main__":
    main()
