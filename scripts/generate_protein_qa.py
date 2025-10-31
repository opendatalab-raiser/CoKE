import os
import json
import lmdb
import argparse
from typing import List, Dict, Optional
from tqdm import tqdm

def get_protein_info(entry_dir: str, protein_id: str) -> Optional[Dict[str, str]]:
    """
    Get protein's function/pathway/subcellular location information from a JSON file.
    
    Args:
        entry_dir: The directory path for the entries.
        protein_id: The protein ID.
        
    Returns:
        Dict: A dictionary containing the three types of information, or None if an error occurs.
    """
    try:
        file_path = os.path.join(entry_dir, f'{protein_id}.json')
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        result = {
            'function': "",
            'pathway': "",
            'subcellular_location': ""
        }
        
        if 'comments' in data and len(data['comments']) > 0:
            for comment in data['comments']:
                # Process FUNCTION
                if comment["commentType"] == "FUNCTION":
                    for text in comment["texts"]:
                        result['function'] += text["value"] + "\n"
                
                # Process PATHWAY
                elif comment["commentType"] == "PATHWAY":
                    for text in comment["texts"]:
                        result['pathway'] += text["value"] + "\n"
                
                # Process SUBCELLULAR LOCATION
                elif comment["commentType"] == "SUBCELLULAR LOCATION":
                    for location in comment.get("subcellularLocations", []):
                        if "location" in location:
                            result['subcellular_location'] += location["location"]["value"] + "\n"
        
        # Remove extra spaces and newlines
        for key in result:
            result[key] = result[key].strip()
            if not result[key]:
                result[key] = None
        
        return result
        
    except Exception as e:
        print(f"Error processing {protein_id}: {str(e)}")
        return None

def generate_qa_pairs(protein_id: str, protein_info: Dict[str, str]) -> List[Dict]:
    """
    Generate QA pairs based on protein information.
    
    Args:
        protein_id: The protein ID.
        protein_info: A dictionary of protein information.
        
    Returns:
        List: A list of QA pairs.
    
    Note: The 'answer' field will be saved as 'ground_truth' in LMDB format
    to match the format expected by integrated_pipeline.py
    """
    qa_pairs = []
    
    # Add function QA pair
    if protein_info['function']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the function of this protein?',
            'answer': protein_info['function'],  # Saved as ground_truth in LMDB
            'question_type': 'function'
        })
    
    # Add pathway QA pair
    if protein_info['pathway']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the pathway of this protein?',
            'answer': protein_info['pathway'],  # Saved as ground_truth in LMDB
            'question_type': 'pathway'
        })
    
    # Add subcellular location QA pair
    if protein_info['subcellular_location']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the subcellular location of this protein?',
            'answer': protein_info['subcellular_location'],  # Saved as ground_truth in LMDB
            'question_type': 'subcellular_location'
        })
    
    return qa_pairs

def save_to_lmdb(qa_pairs: List[Dict], lmdb_path: str):
    """
    Save the QA pairs to an LMDB database.
    
    Args:
        qa_pairs: A list of QA pairs.
        lmdb_path: The path to the LMDB database.
    
    Note: Saves in format [protein_id, question, ground_truth, question_type]
    to match the format expected by integrated_pipeline.py
    """
    # Ensure directory exists
    os.makedirs(os.path.dirname(lmdb_path) if os.path.dirname(lmdb_path) else '.', exist_ok=True)
    
    env = lmdb.open(lmdb_path, map_size=1099511627776)  # 1TB maximum size
    with env.begin(write=True) as txn:
        for idx, qa in enumerate(qa_pairs):
            key = str(idx).encode('utf-8')
            # Save in list format: [protein_id, question, ground_truth, question_type]
            # This format is compatible with get_qa_data function in integrated_pipeline.py
            value = json.dumps([
                qa['protein_id'],
                qa['question'],
                qa['answer'],  # This is the ground_truth answer
                qa['question_type']
            ], ensure_ascii=False).encode('utf-8')
            txn.put(key, value)
    env.close()
    print(f"Saved {len(qa_pairs)} QA pairs to LMDB: {lmdb_path}")

def save_to_json(qa_pairs: List[Dict], json_path: str):
    """
    Save the QA pairs to a JSON file.
    
    Args:
        qa_pairs: A list of QA pairs.
        json_path: The path to the JSON file.
    """
    with open(json_path, 'w') as f:
        json.dump(qa_pairs, f, indent=2)

def read_protein_ids_from_file(file_path: str) -> List[str]:
    """
    Read a list of protein IDs from a text file.
    
    Args:
        file_path: The path to the protein ID file.
        
    Returns:
        List: A list of protein IDs.
    """
    with open(file_path, 'r') as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    return protein_ids

def process_proteins(entry_dir: str, protein_ids: List[str], lmdb_path: str, json_path: str):
    """
    Process the specified list of protein IDs, generate QA pairs, and save them.
    
    Args:
        entry_dir: The directory path for the entries.
        protein_ids: The list of protein IDs to process.
        lmdb_path: The save path for the LMDB database.
        json_path: The save path for the JSON file.
    """
    all_qa_pairs = []
    processed_count = 0
    
    for protein_id in tqdm(protein_ids):
        protein_info = get_protein_info(entry_dir, protein_id)
        if protein_info is None:
            continue
            
        qa_pairs = generate_qa_pairs(protein_id, protein_info)
        all_qa_pairs.extend(qa_pairs)
        processed_count += 1
    
    # Save to LMDB
    save_to_lmdb(all_qa_pairs, lmdb_path)
    
    # Save to JSON
    save_to_json(all_qa_pairs, json_path)
    
    print(f"Successfully processed {len(all_qa_pairs)} QA pairs from {processed_count}/{len(protein_ids)} proteins")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate protein function QA pairs from UniProt Entry files.')
    parser.add_argument('--entry_dir', type=str, required=True, 
                       help='Directory path for UniProt Entries (containing protein JSON files).')
    parser.add_argument('--protein_id_files', type=str, default=None, 
                       help='File path for the list of protein IDs (one ID per line). If not specified, all JSON files in entry_dir will be processed.')
    parser.add_argument('--lmdb_path', type=str, required=True, 
                       help='Output path for the LMDB database.')
    parser.add_argument('--json_path', type=str, required=True, 
                       help='Output path for the JSON file.')
    args = parser.parse_args()
    save_dir = os.path.dirname(args.lmdb_path)
    print(save_dir)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir, exist_ok=True)
    
    if args.protein_id_files:
        # If protein_id_files is specified, process only these IDs
        protein_ids = read_protein_ids_from_file(args.protein_id_files)
    else:
        # Otherwise, process all .json files in the entry directory
        protein_ids = [os.path.splitext(f)[0] for f in os.listdir(args.entry_dir) if f.endswith('.json')]
    
    process_proteins(args.entry_dir, protein_ids, args.lmdb_path, args.json_path)