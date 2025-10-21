from gradio_client import Client
import os
import json
from tqdm import tqdm
import time
import random

def get_protrek_text(sequence, topk=3):
    """
    Retrieves ProTrek text descriptions.
    
    Args:
        sequence (str): The protein sequence.
        topk (int): The number of top results to return.
    
    Returns:
        list: A list containing text descriptions, or an empty list on failure.
    """
    
    try:
        client = Client("http://search-protrek.com/")
        result = client.predict(
            input=sequence,
            nprobe=1000,
            topk=topk,
            input_type="sequence",
            query_type="text",
            subsection_type="Function",
            db="Swiss-Prot",
            api_name="/search"
        )

        full_res = result[-1]["value"]
        data = full_res["data"]

        descriptions = []
        for des, score in data:
            descriptions.append((des, score))
        return descriptions
        
    except Exception as e:
        return []

def get_missing_pids(pids, output_dir):
    """Identifies protein IDs for which output files are missing."""
    missing_pids = []
    for pid in pids:
        if not os.path.exists(f'{output_dir}/{pid}.json'):
            missing_pids.append(pid)
    return missing_pids

def save_protrek_text_to_file():
    """
    Fetches ProTrek descriptions for proteins and saves them to files.
    This script is designed to be stoppable and resumable. If it fails,
    it exits, and can be restarted to continue where it left off.
    """
    import argparse
    from Bio import SeqIO
    
    parser = argparse.ArgumentParser(description='Generate protein function text descriptions using ProTrek.')
    # Supports two modes: 1) Read from a FASTA file 2) Read from PIDs and a JSON file
    parser.add_argument('--input_fasta', type=str, default=None, help='Input FASTA file path (this argument takes precedence).')
    parser.add_argument('--pids_path', type=str, default=None, help='Path to a file with a list of protein IDs (one ID per line).')
    parser.add_argument('--topk', type=int, default=3, help='The number of top results to return (default: 3).')
    parser.add_argument('--protein_seq_info', type=str, default=None, help='Path to the JSON file with protein sequence information (used with pids_path).')
    parser.add_argument('--output_dir', type=str, default='protrek_text', help='Output directory path (default: protrek_text).')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Prioritize using the FASTA file
    if args.input_fasta and os.path.exists(args.input_fasta):
        print(f"Reading sequences from FASTA file: {args.input_fasta}")
        pids = []
        protein_info = {}
        with open(args.input_fasta, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                pids.append(record.id)
                protein_info[record.id] = str(record.seq)
    else:
        # Use the traditional mode: read from PIDs and JSON file
        print(f"Reading sequences from PIDs file and JSON file")
        with open(args.pids_path, 'r') as f:
            pids = f.read().splitlines()
        with open(args.protein_seq_info, 'r') as f:
            protein_info = json.load(f)
    
    # Get the currently missing PIDs
    missing_pids = get_missing_pids(pids, args.output_dir)
    print(f'Currently missing PIDs: {len(missing_pids)}')
    
    if len(missing_pids) == 0:
        print("All PIDs have been processed!")
        return True  # Return True to indicate all tasks are complete
    
    # Process the missing PIDs; exit immediately on failure
    for pid in tqdm(missing_pids, desc="Processing PIDs"):
        try:
            sequence = protein_info[pid]
            descriptions = get_protrek_text(sequence, args.topk)
            
            if descriptions:  # Successfully retrieved descriptions
                with open(f'{args.output_dir}/{pid}.json', 'w') as f:
                    json.dump(descriptions, f)
                print(f"✓ PID {pid} processed successfully")
            else:
                print(f"✗ Failed to get descriptions for PID {pid}, exiting to retry later")
                exit(1)  # Exit immediately with a non-zero status code
                
        except Exception as e:
            print(f"✗ An unexpected error occurred while processing PID {pid}: {str(e)}, exiting to retry later")
            exit(1)  # Exit immediately with a non-zero status code
    
    print("All PIDs in this run were processed successfully!")

if __name__ == "__main__":
    save_protrek_text_to_file()