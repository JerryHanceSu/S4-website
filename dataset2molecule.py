#!/usr/bin/env python3
"""
Convert electrolytes.json to molecules.json for the Molecule Landscape visualization.

This script processes the electrolyte dataset and:
- Extracts unique SMILES from solvent data
- Extracts molecule names from electrolyte names (part after "LiFSI")
- Determines fluorinated status by checking for F in SMILES
- Generates 2D embeddings using RDKit Morgan fingerprints + PCA
- Creates high-resolution SVG images of molecules
"""

import json
import re
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
from tqdm import tqdm

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    print("All required packages imported successfully")
except ImportError as e:
    print(f"Missing required package: {e}")
    print("Please install: pip install rdkit-pypi scikit-learn")
    exit(1)

def extract_molecule_name(electrolyte_name: str) -> str:
    """Extract molecule name from electrolyte name (part after LiFSI)"""
    if not electrolyte_name or not isinstance(electrolyte_name, str):
        return "Unknown"
    
    # Pattern to match various formats:
    # "X M LiFSI [NAME]" -> "NAME"
    # "X M LiFSI in [NAME]" -> "NAME"
    # "X LiFSI:Y [NAME]" -> "NAME" 
    # "LiFSI [NAME]" -> "NAME"
    patterns = [
        r'\d+\.?\d*\s*M\s+LiFSI\s+in\s+(.+?)(?:\s*\(|$)',    # "1 M LiFSI in Tetraglyme" -> "Tetraglyme"
        r'\d+\.?\d*\s*M\s+LiFSI\s+(.+?)(?:\s*\(|$)',         # "1 M LiFSI DMB" -> "DMB"
        r'\d+\.?\d*\s+LiFSI:\d+\.?\d*\s+(.+?)(?:\s*\(|$)',   # "1 LiFSI:1.5 DMC" -> "DMC"
        r'LiFSI\s+in\s+(.+?)(?:\s*\(|$)',                    # "LiFSI in NAME" -> "NAME"
        r'LiFSI\s+(.+?)(?:\s*\(|$)',                         # "LiFSI DME" -> "DME"
        r'LiFSI[-:](.+?)(?:\s*\(|$)',                        # "LiFSI-DME" or "LiFSI:DME" -> "DME"
    ]
    
    for pattern in patterns:
        match = re.search(pattern, electrolyte_name, re.IGNORECASE)
        if match:
            name = match.group(1).strip()
            # Clean up common suffixes/prefixes
            name = re.sub(r'\s*\(.*?\)', '', name)  # Remove parentheses content
            name = re.sub(r'\s+', ' ', name)        # Normalize spaces
            return name
    
    # Special case: if format is "X LiFSI:Y NAME", take the last word
    if 'LiFSI:' in electrolyte_name:
        words = electrolyte_name.split()
        if len(words) >= 3:  # At least "X LiFSI:Y NAME"
            return words[-1]  # Return the last word
    
    # Fallback: return the whole name if pattern doesn't match
    return electrolyte_name

def is_fluorinated(smiles: str) -> bool:
    """Check if molecule is fluorinated by examining SMILES string"""
    return 'F' in smiles if smiles else False

def smiles_to_filename(smiles: str) -> str:
    """Convert SMILES to a safe filename using hash"""
    if not smiles:
        return "unknown.svg"
    # Create MD5 hash of SMILES for consistent filename
    hash_obj = hashlib.md5(smiles.encode('utf-8'))
    return f"{hash_obj.hexdigest()}.svg"

def generate_morgan_fingerprint(mol, radius: int = 2, nBits: int = 2048) -> np.ndarray:
    """Generate Morgan fingerprint for a molecule"""
    try:
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        return np.array(fp)
    except Exception as e:
        print(f"Error generating fingerprint: {e}")
        return np.zeros(nBits)

def generate_molecule_image(smiles: str, filename: str, size: Tuple[int, int] = (400, 400)) -> bool:
    """Generate SVG image of molecule using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return False
        
        # Create drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        drawer.SetFontSize(0.8)
        
        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Save SVG
        svg = drawer.GetDrawingText()
        filepath = Path("info/molecules") / filename
        filepath.parent.mkdir(exist_ok=True)
        
        with open(filepath, 'w') as f:
            f.write(svg)
        
        return True
        
    except Exception as e:
        print(f"Error generating image for {smiles}: {e}")
        return False

def create_smiles_mapping(molecules: List[Dict]) -> None:
    """Create JavaScript mapping file for SMILES to filenames"""
    mapping = {}
    for mol in molecules:
        smiles = mol['smiles']
        filename = smiles_to_filename(smiles)
        mapping[smiles] = filename
    
    js_content = f"// Auto-generated SMILES to filename mapping\nwindow.SMILES_TO_FILE = {json.dumps(mapping, indent=2)};\n"
    
    with open("info/molecules/smiles_mapping.js", 'w') as f:
        f.write(js_content)
    
    print(f"Generated SMILES mapping with {len(mapping)} entries")

def process_molecules(input_file: str = "info/electrolytes.json", output_file: str = "info/molecules.json"):
    """Process electrolytes and generate molecule data with embeddings"""
    
    print(f"Loading electrolyte data from {input_file}...")
    with open(input_file, 'r') as f:
        electrolytes = json.load(f)
    
    print(f"Loaded {len(electrolytes)} electrolyte entries")
    
    # Extract unique molecules
    unique_molecules = {}  # smiles -> molecule_data
    
    print("Extracting unique molecules...")
    for electrolyte in tqdm(electrolytes, desc="Processing electrolytes"):
        smiles = electrolyte.get('SMILES', '').strip()
        if not smiles:
            continue
            
        if smiles not in unique_molecules:
            electrolyte_name = electrolyte.get('Electrolyte', '')
            molecule_name = extract_molecule_name(electrolyte_name)
            solvent_class = electrolyte.get('Class', 'Other')
            fluorinated = is_fluorinated(smiles)
            
            unique_molecules[smiles] = {
                'smiles': smiles,
                'name': molecule_name,
                'class': solvent_class,
                'fluorinated': fluorinated
            }
    
    print(f"Found {len(unique_molecules)} unique molecules")
    molecules_list = list(unique_molecules.values())
    
    # Generate molecular fingerprints
    print("Generating molecular fingerprints...")
    fingerprints = []
    valid_molecules = []
    
    for mol_data in tqdm(molecules_list, desc="Computing fingerprints"):
        smiles = mol_data['smiles']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is not None:
            fp = generate_morgan_fingerprint(mol)
            fingerprints.append(fp)
            valid_molecules.append(mol_data)
        else:
            print(f"Skipping invalid SMILES: {smiles}")
    
    if not fingerprints:
        print("No valid molecules found!")
        return
    
    print(f"Generated fingerprints for {len(valid_molecules)} molecules")
    
    # Convert to numpy array and standardize
    fingerprints = np.array(fingerprints)
    scaler = StandardScaler()
    fingerprints_scaled = scaler.fit_transform(fingerprints)
    
    # Apply t-SNE to reduce to 2D
    print("Applying t-SNE for 2D embedding...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(valid_molecules)-1), 
                n_iter=1000, learning_rate=200.0, verbose=1)
    embeddings_2d = tsne.fit_transform(fingerprints_scaled)
    
    print(f"t-SNE embedding completed")
    print(f"Final KL divergence: {tsne.kl_divergence_:.3f}")
    
    # Add coordinates to molecule data
    for i, mol_data in enumerate(valid_molecules):
        mol_data['x'] = float(embeddings_2d[i, 0])
        mol_data['y'] = float(embeddings_2d[i, 1])
    
    # Generate molecule images
    print("Generating molecule images...")
    success_count = 0
    
    for mol_data in tqdm(valid_molecules, desc="Generating images"):
        smiles = mol_data['smiles']
        filename = smiles_to_filename(smiles)
        
        if generate_molecule_image(smiles, filename):
            success_count += 1
    
    print(f"Successfully generated {success_count}/{len(valid_molecules)} molecule images")
    
    # Create SMILES mapping file
    create_smiles_mapping(valid_molecules)
    
    # Save molecules data
    print(f"Saving molecules data to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(valid_molecules, f, indent=2)
    
    print(f"Successfully generated {output_file}")
    
    # Print summary statistics
    print("\nSummary:")
    print(f"Total unique molecules: {len(valid_molecules)}")
    classes = set(mol['class'] for mol in valid_molecules)
    print(f"Solvent classes: {sorted(classes)}")
    fluorinated_count = sum(1 for mol in valid_molecules if mol['fluorinated'])
    print(f"Fluorinated molecules: {fluorinated_count}/{len(valid_molecules)}")
    
    # Print coordinate ranges for verification
    x_coords = [mol['x'] for mol in valid_molecules]
    y_coords = [mol['y'] for mol in valid_molecules]
    print(f"X coordinate range: {min(x_coords):.3f} to {max(x_coords):.3f}")
    print(f"Y coordinate range: {min(y_coords):.3f} to {max(y_coords):.3f}")

if __name__ == "__main__":
    process_molecules()
