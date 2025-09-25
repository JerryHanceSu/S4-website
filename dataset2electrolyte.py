#!/usr/bin/env python3
"""
Convert dataset_og.csv to electrolytes.json for the Electrolyte Landscape visualization.

This script processes the electrolyte dataset and calculates derived values:
- LIC: log10(Ionic Conductivity) 
- LCE: -log10(100-CE)+2
- Salt:Solvent ratio from Salt and Solvent columns
"""

import pandas as pd
import json
import numpy as np
import re
from typing import Dict, Any, Optional

def extract_smiles_from_solvent(solvent_data: str) -> str:
    """Extract SMILES string from solvent data (format: SMILES:ratio)"""
    if not solvent_data or pd.isna(solvent_data):
        return ""
    
    solvent_str = str(solvent_data).strip()
    colon_index = solvent_str.find(':')
    
    if colon_index > 0:
        return solvent_str[:colon_index]
    return solvent_str

def extract_mol_from_salt(salt_data: str) -> Optional[float]:
    """Extract molar concentration from salt data"""
    if not salt_data or pd.isna(salt_data):
        return None
    
    salt_str = str(salt_data).strip()
    colon_index = salt_str.rfind(':')
    
    if colon_index > 0:
        try:
            mol_salt = float(salt_str[colon_index + 1:])
            return mol_salt
        except ValueError:
            return None
    return None

def extract_mol_from_solvent(solvent_data: str) -> Optional[float]:
    """Extract molar ratio from solvent data"""
    if not solvent_data or pd.isna(solvent_data):
        return None
    
    solvent_str = str(solvent_data).strip()
    colon_index = solvent_str.rfind(':')
    
    if colon_index > 0:
        try:
            mol_solvent = float(solvent_str[colon_index + 1:])
            return mol_solvent
        except ValueError:
            return None
    return None

def calculate_salt_solvent_ratio(salt_data: str, solvent_data: str) -> Optional[Dict[str, float]]:
    """Calculate salt:solvent ratio from salt and solvent data"""
    mol_salt = extract_mol_from_salt(salt_data)
    mol_solvent = extract_mol_from_solvent(solvent_data)
    
    if mol_salt and mol_solvent:
        solvent_per_salt = mol_solvent / mol_salt
        return {"salt": 1.0, "solvent": solvent_per_salt}
    
    return None

def format_ratio_for_display(salt_data: str, solvent_data: str) -> str:
    """Format ratio for display (e.g., '1:2.5')"""
    ratio = calculate_salt_solvent_ratio(salt_data, solvent_data)
    if ratio and ratio.get("solvent"):
        return f"1:{ratio['solvent']:.1f}"
    return "(unknown ratio)"

def calculate_lic(ionic_conductivity: float) -> Optional[float]:
    """Calculate LIC = log10(Ionic Conductivity)"""
    if pd.isna(ionic_conductivity) or ionic_conductivity <= 0:
        return None
    return np.log10(ionic_conductivity)

def calculate_lce(ce_percent: float) -> Optional[float]:
    """Calculate LCE = -log10(100-CE)+2"""
    if pd.isna(ce_percent) or ce_percent >= 100 or ce_percent <= 0:
        return None
    return -np.log10(100 - ce_percent) + 2

def process_dataset(input_file: str = "dataset.csv", output_file: str = "info/electrolytes.json"):
    """Process the dataset and generate JSON file for the electrolyte landscape"""
    
    print(f"Loading dataset from {input_file}...")
    df = pd.read_csv(input_file)
    print(f"Loaded {len(df)} rows")
    
    # Process each row
    electrolytes = []
    
    for idx, row in df.iterrows():
        try:
            # Extract basic data
            electrolyte_name = row.get("Electrolyte", "")
            ionic_conductivity = row.get("Ionic Conductivity")
            ce_percent = row.get("CE (%)")
            temperature = row.get("Temperature")
            method = row.get("Method", "")
            current = row.get("Current (mA/cm2)")
            capacity = row.get("Capacity (mAh/cm2)")
            cycle = row.get("Cycle")
            solvent_data = row.get("Solvent", "")
            solvent_class = row.get("Class", "")
            salt_data = row.get("Salt", "")
            reference = row.get("Reference")
            
            # Skip rows with missing critical data
            if pd.isna(ionic_conductivity) or pd.isna(ce_percent):
                print(f"Skipping row {idx + 1}: Missing ionic conductivity or CE data")
                continue
            
            # Calculate derived values
            lic = calculate_lic(ionic_conductivity)
            lce = calculate_lce(ce_percent)
            
            # Extract SMILES and calculate ratio
            smiles = extract_smiles_from_solvent(solvent_data)
            ratio_dict = calculate_salt_solvent_ratio(salt_data, solvent_data)
            ratio_display = format_ratio_for_display(salt_data, solvent_data)
            
            # Create electrolyte entry
            electrolyte_entry = {
                "Electrolyte": str(electrolyte_name) if not pd.isna(electrolyte_name) else "",
                "Ionic Conductivity": float(ionic_conductivity) if not pd.isna(ionic_conductivity) else None,
                "LIC": float(lic) if lic is not None else None,
                "Temperature": float(temperature) if not pd.isna(temperature) else None,
                "CE (%)": float(ce_percent) if not pd.isna(ce_percent) else None,
                "LCE": float(lce) if lce is not None else None,
                "Method": str(method) if not pd.isna(method) else "",
                "Current (mA/cm2)": float(current) if not pd.isna(current) else None,
                "Capacity (mAh/cm2)": float(capacity) if not pd.isna(capacity) else None,
                "Cycle": str(cycle) if not pd.isna(cycle) else "",
                "Solvent": str(solvent_data) if not pd.isna(solvent_data) else "",
                "SMILES": smiles,
                "Class": str(solvent_class) if not pd.isna(solvent_class) else "",
                "Salt": str(salt_data) if not pd.isna(salt_data) else "",
                "SaltSolventRatio": ratio_dict,
                "RatioDisplay": ratio_display,
                "Reference": int(reference) if not pd.isna(reference) else None
            }
            
            electrolytes.append(electrolyte_entry)
            
        except Exception as e:
            print(f"Error processing row {idx + 1}: {e}")
            continue
    
    print(f"Successfully processed {len(electrolytes)} electrolyte entries")
    
    # Save to JSON
    print(f"Saving to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(electrolytes, f, indent=2, ensure_ascii=False)
    
    print(f"Successfully generated {output_file}")
    
    # Print some statistics
    print("\nDataset Statistics:")
    print(f"Total electrolytes: {len(electrolytes)}")
    classes = set(e["Class"] for e in electrolytes if e["Class"])
    print(f"Solvent classes: {sorted(classes)}")
    methods = set(e["Method"] for e in electrolytes if e["Method"])
    print(f"Methods: {sorted(methods)}")
    
    # Check for valid LIC/LCE calculations
    valid_lic = sum(1 for e in electrolytes if e["LIC"] is not None)
    valid_lce = sum(1 for e in electrolytes if e["LCE"] is not None)
    print(f"Valid LIC calculations: {valid_lic}/{len(electrolytes)}")
    print(f"Valid LCE calculations: {valid_lce}/{len(electrolytes)}")

if __name__ == "__main__":
    process_dataset()
