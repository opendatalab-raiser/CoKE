import os
import subprocess
import tempfile
from typing import Dict, List, Optional
from Bio import SeqIO


class FoldseekSearch:
    def __init__(self, database: str = "foldseek_db/sp", num_threads: int = 64):
        """
        Foldseek search wrapper.
        
        Args:
            database: Path to the Foldseek database.
            num_threads: Number of threads to use for searching.
        """
        self.database = database
        self.num_threads = num_threads
        
        # Check if foldseek is available
        try:
            subprocess.run(['foldseek', 'version'], 
                         capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError("Foldseek is not installed or not in PATH. "
                             "Please install it using: conda install -c bioconda foldseek")
    
    def run_search(self, fasta_file: str, output_file: str, 
                   evalue: float = 0.01, temp_dir: str = "temp_foldseek") -> str:
        """
        Run Foldseek easy-search.
        
        Args:
            fasta_file: Path to the input FASTA file.
            output_file: Path to the output file.
            evalue: E-value threshold.
            temp_dir: Temporary directory for Foldseek.
            
        Returns:
            Path to the output file.
        """
        # Create temporary directory if it doesn't exist
        os.makedirs(temp_dir, exist_ok=True)
        
        # Run foldseek easy-search
        # Format: query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
        cmd = [
            'foldseek', 'easy-search',
            fasta_file,
            self.database,
            output_file,
            temp_dir,
            '--threads', str(self.num_threads),
            '-e', str(evalue),
            '--format-mode', '4',
            '--format-output', 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Foldseek search completed successfully.")
            return output_file
        except subprocess.CalledProcessError as e:
            print(f"Error running Foldseek: {e.stderr}")
            raise
    
    def parse_results(self, result_file: str, fasta_file: str) -> Dict:
        """
        Parse Foldseek search results.
        
        Args:
            result_file: Path to the Foldseek output file.
            fasta_file: Path to the input FASTA file (to get sequences).
            
        Returns:
            A dictionary mapping query IDs to their results.
            Format similar to BLAST results:
            {
                "query_id": {
                    "sequence": "MKTA...",
                    "foldseek_results": [
                        {
                            "ID": "P12345",
                            "Identity%": 45.2,
                            "Coverage%": 89.3,
                            "E-value": "1.0e-60",
                            "Bit Score": 234.5,
                            "Positive%": 45.2  # Same as identity for structure alignment
                        }
                    ]
                }
            }
        """
        # Read sequences from FASTA
        seq_dict = {}
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seq_dict[record.id] = str(record.seq)
        
        # Initialize results dictionary
        results = {query_id: [] for query_id in seq_dict.keys()}
        
        # Parse Foldseek output
        if not os.path.exists(result_file):
            print(f"Warning: Foldseek result file {result_file} does not exist.")
            return {query_id: {"sequence": seq, "foldseek_results": []} 
                    for query_id, seq in seq_dict.items()}
        
        with open(result_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 14:
                    continue
                
                query_id = fields[0]
                target_id = fields[1]
                identity = float(fields[2])
                alnlen = int(fields[3])
                evalue = float(fields[10])
                bits = float(fields[11])
                qlen = int(fields[12])
                tlen = int(fields[13])
                
                # Calculate coverage
                coverage = (alnlen / qlen) * 100 if qlen > 0 else 0
                
                # Extract UniProt ID from target (format may be sp|P12345|NAME or just P12345)
                if '|' in target_id:
                    uniprot_id = target_id.split('|')[1]
                else:
                    uniprot_id = target_id
                
                # Format e-value like BLAST does
                if evalue < 0.001:
                    evalue_str = f"{evalue:.1e}"
                else:
                    evalue_str = str(round(evalue, 4))
                
                result_entry = {
                    "ID": uniprot_id,
                    "Identity%": round(identity, 2),
                    "Coverage%": round(coverage, 2),
                    "E-value": evalue_str,
                    "Bit Score": round(bits, 1),
                    "Positive%": round(identity, 2)  # For structure alignment, use identity
                }
                
                if query_id in results:
                    results[query_id].append(result_entry)
        
        # Format final output
        formatted_results = {}
        for query_id, seq in seq_dict.items():
            formatted_results[query_id] = {
                "sequence": seq,
                "foldseek_results": results.get(query_id, [])
            }
        
        return formatted_results


def run_foldseek_analysis(fasta_file: str, database: str = "foldseek_db/sp",
                         evalue: float = 0.01, num_threads: int = 64,
                         output_file: Optional[str] = None,
                         temp_dir: str = "temp_foldseek") -> Dict:
    """
    High-level function to run Foldseek analysis.
    
    Args:
        fasta_file: Path to the input FASTA file.
        database: Path to the Foldseek database.
        evalue: E-value threshold.
        num_threads: Number of threads.
        output_file: Path to save results (optional).
        temp_dir: Temporary directory.
        
    Returns:
        Dictionary of Foldseek results.
    """
    foldseek = FoldseekSearch(database=database, num_threads=num_threads)
    
    # Create output file if not provided
    if output_file is None:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.m8') as tmp:
            output_file = tmp.name
    
    # Run search
    result_file = foldseek.run_search(fasta_file, output_file, evalue, temp_dir)
    
    # Parse results
    results = foldseek.parse_results(result_file, fasta_file)
    
    return results


if __name__ == "__main__":
    # Test the module
    import sys
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
        results = run_foldseek_analysis(fasta_file)
        print(f"Processed {len(results)} sequences")
        for query_id, data in results.items():
            print(f"{query_id}: {len(data['foldseek_results'])} hits")

