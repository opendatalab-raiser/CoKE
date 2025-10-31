import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import requests
from utils.mpr import MultipleProcessRunnerSimplifier
from typing import List

class PDBDownloader(MultipleProcessRunnerSimplifier):
    def __init__(self, protein_ids: List[str], output_dir: str, failed_ids_path: str, **kwargs):
        # Create a temporary file path to collect failed IDs
        temp_failed_path = f"{failed_ids_path}.temp"
        
        super().__init__(
            data=protein_ids,
            do=self.download_single,
            save_path=temp_failed_path,  # Set temporary save path
            return_results=True,
            **kwargs
        )
        self.output_dir = output_dir
        self.failed_ids_path = failed_ids_path
        self.temp_failed_path = temp_failed_path
        os.makedirs(self.output_dir, exist_ok=True)
    
    def download_single(self, process_id, idx, protein_id, writer):
        """
        Downloads a single PDB file from AlphaFold.
        """
        url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_id}-F1-model_v6.pdb"
        output_path = os.path.join(self.output_dir, f"AF-{protein_id}-F1-model_v6.pdb")
        
        if os.path.exists(output_path):
            return None
            
        try:
            # Note: Timeout is hardcoded to 30 in the original code
            response = requests.get(url, timeout=30) 
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                return None
            else:
                print(f"Download of protein {protein_id} failed, status code: {response.status_code}")
                if writer:
                    writer.write(f"{protein_id}\n")
                return protein_id
        except Exception as e:
            print(f"Error during download of protein {protein_id}: {str(e)}")
            if writer:
                writer.write(f"{protein_id}\n")
            return protein_id
    
    def save_failed_ids(self, results):
        """
        Reads failed IDs from the temporary file, merges with results, and saves them.
        """
        # Read failed IDs from the temporary file
        failed_ids = []
        if os.path.exists(self.temp_failed_path):
            with open(self.temp_failed_path, 'r') as f:
                failed_ids = [line.strip() for line in f if line.strip()]
            
            # Delete the temporary file
            os.remove(self.temp_failed_path)
        
        # Also collect failed IDs from the returned results
        result_failed_ids = [result for result in results if result is not None]
        
        # Merge failed IDs from both sources and deduplicate
        all_failed_ids = list(set(failed_ids + result_failed_ids))
        
        if all_failed_ids:
            print(f"A total of {len(all_failed_ids)} proteins failed to download")
            # Save the failed IDs to a file
            with open(self.failed_ids_path, 'w') as f:
                for protein_id in all_failed_ids:
                    f.write(f"{protein_id}\n")
            print(f"Saved failed download protein IDs to {self.failed_ids_path}")

def parse_args():
    """Parses command line arguments"""
    import argparse
    parser = argparse.ArgumentParser(description="Download PDB files")
    parser.add_argument("--protein_ids_file", type=str, required=True, help="Path to the protein ID list file")
    parser.add_argument("--output_dir", type=str, default="/zhuangkai/projects/TTS4Protein/data/raw_data/pdb", help="Output directory")
    parser.add_argument("--failed_ids_path", type=str, default="/zhuangkai/projects/TTS4Protein/data/processed_data/protein_wo_pdb.txt", help="Path to save failed IDs")
    parser.add_argument("--n_process", type=int, default=64, help="Number of processes")
    return parser.parse_args()

def main():
    # Read the protein ID list
    args = parse_args()
    protein_ids_file = args.protein_ids_file
    with open(protein_ids_file, 'r') as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    
    print(f"A total of {len(protein_ids)} proteins need to be downloaded")
    
    # Configure download parameters
    downloader = PDBDownloader(
        protein_ids=protein_ids,
        output_dir=args.output_dir,
        failed_ids_path=args.failed_ids_path,
        n_process=args.n_process,
        verbose=True,
        total_only=False,
        log_step=1,
        split_strategy="queue"
    )
    
    results = downloader.run()
    downloader.save_failed_ids(results)

if __name__ == "__main__":
    main()