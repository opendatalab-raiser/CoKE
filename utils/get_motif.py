import json
from collections import Counter
import os

# Global variables to cache loaded data, avoiding redundant file reads.
_pfam_dict = None
_pfam_descriptions = None

def _load_pfam_data(protein2pfam_path):
    """Loads the protein-to-Pfam mapping data from a JSON file."""
    global _pfam_dict
    if _pfam_dict is None:
        with open(protein2pfam_path, 'r') as file:
            _pfam_dict = json.load(file)

def _load_pfam_descriptions(pfam_descriptions_path):
    """Loads the Pfam descriptions from a JSON file."""
    global _pfam_descriptions
    if _pfam_descriptions is None:
        with open(pfam_descriptions_path, 'r') as file:
            _pfam_descriptions = json.load(file)

def get_motif_pfam(protein_id, protein2pfam_path, pfam_descriptions_path):
    """
    Retrieves Pfam information and its definition for a given protein.
    
    Args:
        protein_id (str): The protein ID.
        protein2pfam_path (str): Path to the interproscan_info.json file.
        pfam_descriptions_path (str): Path to the Pfam descriptions file.
    
    Returns:
        dict: A dictionary mapping Pfam IDs to their definitions, 
              e.g., {"PF04820": "definition content"}.
    """
    _load_pfam_data(protein2pfam_path)
    _load_pfam_descriptions(pfam_descriptions_path)
    
    if protein_id not in _pfam_dict:
        return {}
    
    protein_info = _pfam_dict[protein_id]
    _pfam_dicts = protein_info.get('interproscan_results', {}).get('pfam_id', [])
    pfam_ids = []
    for pfam_dict in _pfam_dicts:
        for key, value in pfam_dict.items():
            pfam_ids.append(key)
    
    result = {}
    for pfam_id in pfam_ids:
        if pfam_id in _pfam_descriptions:
            result[pfam_id] = _pfam_descriptions[pfam_id]['description']
    
    return result

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Get Pfam motif information for a protein.')
    parser.add_argument("--protein_id", type=str, required=True, help='The protein ID, e.g., A8CF74')
    parser.add_argument("--protein2pfam_path", type=str, required=True, help='Path to the protein-Pfam mapping file (InterProScan results JSON file)')
    parser.add_argument("--pfam_descriptions_path", type=str, required=True, help='Path to the Pfam descriptions file (JSON file containing Pfam IDs and descriptions)')
    args = parser.parse_args()
    result = get_motif_pfam(args.protein_id, args.protein2pfam_path, args.pfam_descriptions_path)
    print(json.dumps(result, indent=2))