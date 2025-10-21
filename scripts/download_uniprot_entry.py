import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import requests
from utils.mpr import MultipleProcessRunnerSimplifier
import json

class UniProtDownloader(MultipleProcessRunnerSimplifier):
    def __init__(self, protein_ids, output_dir, **kwargs):
        super().__init__(
            data=protein_ids,
            do=self.download_single,
            save_path=None,  # Each file is saved individually, no aggregation needed
            return_results=False,
            **kwargs
        )
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def download_single(self, process_id, idx, protein_id, _):
        url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            output_path = os.path.join(self.output_dir, f"{protein_id}.json")
            with open(output_path, "w") as f:
                json.dump(response.json(), f, indent=2)
                
        except Exception as e:
            print(f"Error downloading {protein_id}: {str(e)}")

if __name__ == "__main__":
    # Read the list of protein IDs from a file
    import argparse
    parser = argparse.ArgumentParser(description="Download protein entries")
    parser.add_argument("--input_path", type=str, default="data/raw_data/protein_ids_all20250623.txt", help="File path for the list of protein IDs")
    parser.add_argument("--output_dir", type=str, default="data/raw_data/SP_Entry", help="Output directory")
    parser.add_argument("--n_process", type=int, default=64, help="Number of processes")
    parser.add_argument("--max_round", type=int, default=5, help="Maximum number of rounds")
    args = parser.parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    with open(args.input_path) as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    
    # Check for existing protein IDs in the output_dir to avoid re-downloading what is already present from input_path. This enables multi-round downloading.
    for i in range(args.max_round):
        existing_ids = [f.split(".")[0] for f in os.listdir(args.output_dir) if f.endswith(".json")]
        protein_ids = set(protein_ids) - set(existing_ids)
        protein_ids = list(protein_ids)
        print(f"Round {i+1} of downloading, number of protein IDs to download: {len(protein_ids)}")
        if len(protein_ids) == 0:
            print("All protein IDs have been downloaded.")
            break
        downloader = UniProtDownloader(
            protein_ids=protein_ids,
            output_dir=args.output_dir,
            n_process=args.n_process,  # Adjust based on the number of CPU cores
            split_strategy="static",
            log_step=10
        )
        downloader.run()
        
    # print(f"Number of protein IDs to download: {len(protein_ids)}")

    # # Configure download parameters
    # downloader = UniProtDownloader(
    #     protein_ids=protein_ids,
    #     output_dir=args.output_dir,
    #     n_process=args.n_process,  # Adjust based on the number of CPU cores
    #     split_strategy="static",
    #     log_step=10
    # )
    
    # downloader.run()