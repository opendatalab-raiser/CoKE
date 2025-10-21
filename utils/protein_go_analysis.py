import json
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import defaultdict

# Global variable declarations
_go_data = None
_protein_go_dict = None

def _load_go_data(go_info_path):
    """Lazily loads GO data."""
    global _go_data
    if _go_data is None:
        try:
            with open(go_info_path, 'r') as f:
                _go_data = json.load(f)
        except Exception as e:
            print(f"Error loading GO data file: {str(e)}")
            _go_data = None

def _load_protein_go_dict(protein2gopath):
    """Lazily loads the protein-to-GO mapping data."""
    global _protein_go_dict
    if _protein_go_dict is None:
        try:
            _protein_go_dict = {}
            with open(protein2gopath, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    _protein_go_dict[data['protein_id']] = data['GO_id']
        except Exception as e:
            print(f"Error loading protein-to-GO mapping data: {str(e)}")
            _protein_go_dict = None

def get_go_definition(go_id, go_info_path):
    """Gets the definition of a GO term."""
    _load_go_data(go_info_path)
    if _go_data is None:
        return None

    if not go_id.startswith('GO_'):
        go_id = f"GO_{go_id}"
    full_id = f"http://purl.obolibrary.org/obo/{go_id}"
    
    for node in _go_data['graphs'][0]['nodes']:
        if node['id'] == full_id:
            if 'meta' in node and 'definition' in node['meta']:
                return node['meta']['definition']['val']
    return None

def analyze_protein_go(protein_id, protein2gopath, go_info_path):
    """
    Analyzes the GO annotation information for a protein, including GO IDs and definitions.
    
    Args:
        protein_id (str): The protein ID.
        protein2gopath (str): Path to the protein-to-GO mapping file.
    
    Returns:
        dict: A dictionary containing GO information.
    """
    _load_protein_go_dict(protein2gopath)
    if _protein_go_dict is None:
        return {
            "status": "error",
            "message": "Failed to load GO data"
        }

    if protein_id not in _protein_go_dict:
        return {
            "status": "error",
            "message": f"GO annotation for protein {protein_id} not found"
        }
    
    go_ids = _protein_go_dict[protein_id]
    go_info = []
    all_definitions = {}
    
    for go_id in go_ids:
        # Get the GO definition
        definition = get_go_definition(go_id, go_info_path)
        if definition:
            all_definitions[go_id] = definition
        
        go_info.append({
            "go_id": go_id
        })
    
    return {
        "status": "success" if len(go_info) > 0 else "error",
        "protein_id": protein_id,
        "go_annotations": go_info,
        "all_related_definitions": all_definitions
    }

# Example usage
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Analyze GO annotation information for a protein.')
    parser.add_argument('--protein_id', type=str, required=True, help='The protein ID, e.g., A8CF74')
    parser.add_argument('--protein2gopath', type=str, required=True, help='Path to the protein-to-GO mapping file (JSON Lines format)')
    parser.add_argument('--go_info_path', type=str, required=True, help='Path to the GO ontology file (go.json)')
    args = parser.parse_args()
    
    result = analyze_protein_go(args.protein_id, args.protein2gopath, args.go_info_path)
    
    if result["status"] == "success":
        print(f"\nProtein {result['protein_id']} GO annotations:")
        
        for anno in result["go_annotations"]:
            print(f"\nGO ID: {anno['go_id']}")
        
        print("\nAll related GO ID definitions:")
        for go_id, definition in result["all_related_definitions"].items():
            print(f"\nGO:{go_id}")
            print(f"Definition: {definition}")
    else:
        print(result["message"])