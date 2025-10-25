import os
import subprocess
import tempfile
from typing import Dict, List, Optional


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
    
    def run_search(self, pdb_file: str, output_file: str, 
                   evalue: float = 0.01, temp_dir: str = "temp_foldseek") -> str:
        """
        Run Foldseek easy-search with PDB file.
        
        Args:
            pdb_file: Path to the input PDB file or directory containing PDB files.
            output_file: Path to the output file.
            evalue: E-value threshold.
            temp_dir: Temporary directory for Foldseek.
            
        Returns:
            Path to the output file.
        """
        # Create temporary directory if it doesn't exist
        os.makedirs(temp_dir, exist_ok=True)
        
        # Run foldseek easy-search
        # The default output format is: query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
        cmd = [
            'foldseek', 'easy-search',
            pdb_file,
            self.database,
            output_file,
            temp_dir,
            '--threads', str(self.num_threads),
            '-e', str(evalue)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Foldseek search completed successfully.")
            return output_file
        except subprocess.CalledProcessError as e:
            print(f"Error running Foldseek: {e.stderr}")
            raise
    
    def parse_results(self, result_file: str, pdb_file: str) -> Dict:
        """
        Parse Foldseek search results.
        
        Args:
            result_file: Path to the Foldseek output file.
            pdb_file: Path to the input PDB file or directory (to get sequences).
            
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
        # Extract sequences from PDB files
        seq_dict = self._extract_sequences_from_pdb(pdb_file)
        
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
                if len(fields) < 12:
                    continue
                
                query_id = fields[0]
                target_id = fields[1]
                identity = float(fields[2])
                alnlen = int(fields[3])
                # E-value is in the second-to-last column (index -2)
                evalue = float(fields[-2])
                bits = float(fields[-1])  # Last column is bits
                # For foldseek easy-search, we don't have qlen and tlen in the output
                # We'll calculate coverage based on alignment length
                qlen = alnlen  # Approximation
                tlen = alnlen  # Approximation
                
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

    def _extract_sequences_from_pdb(self, pdb_file: str) -> Dict[str, str]:
        """
        Extract sequences from PDB file(s).
        
        Args:
            pdb_file: Path to PDB file or directory containing PDB files.
            
        Returns:
            Dictionary mapping PDB IDs to sequences.
        """
        seq_dict = {}
        
        if os.path.isfile(pdb_file):
            # Single PDB file
            pdb_files = [pdb_file]
        elif os.path.isdir(pdb_file):
            # Directory containing PDB files
            pdb_files = [os.path.join(pdb_file, f) for f in os.listdir(pdb_file) 
                         if f.endswith('.pdb')]
        else:
            raise ValueError(f"PDB file or directory not found: {pdb_file}")
        
        for pdb_path in pdb_files:
            try:
                # Extract PDB ID from filename
                pdb_id = os.path.basename(pdb_path).replace('.pdb', '')
                
                # Read PDB file and extract sequence
                sequence = self._read_pdb_sequence(pdb_path)
                if sequence:
                    seq_dict[pdb_id] = sequence
                    
            except Exception as e:
                print(f"Warning: Could not extract sequence from {pdb_path}: {str(e)}")
                continue
        
        return seq_dict
    
    def _read_pdb_sequence(self, pdb_path: str) -> str:
        """
        Read amino acid sequence from PDB file.
        
        Args:
            pdb_path: Path to PDB file.
            
        Returns:
            Amino acid sequence string.
        """
        sequence = ""
        atom_lines = []
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    # Extract residue information
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    atom_lines.append((res_num, res_name))
        
        # Sort by residue number and extract sequence
        atom_lines.sort(key=lambda x: x[0])
        
        # Map 3-letter codes to 1-letter codes
        aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        
        for res_num, res_name in atom_lines:
            if res_name in aa_map:
                sequence += aa_map[res_name]
        
        return sequence


def run_foldseek_analysis(pdb_file: str, database: str = "foldseek_db/sp",
                         evalue: float = 0.01, num_threads: int = 64,
                         output_file: Optional[str] = None,
                         temp_dir: str = "temp_foldseek") -> Dict:
    """
    High-level function to run Foldseek analysis with PDB files.
    
    Args:
        pdb_file: Path to the input PDB file or directory containing PDB files.
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
    result_file = foldseek.run_search(pdb_file, output_file, evalue, temp_dir)
    
    # Parse results
    results = foldseek.parse_results(result_file, pdb_file)
    
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

