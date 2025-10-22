import json
import os
import sys
import argparse
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from tqdm import tqdm

# Add paths
root_path = os.path.dirname((os.path.abspath(__file__)))
sys.path.append(root_path)
sys.path.append(os.path.join(root_path, "Models/ProTrek"))

from utils.protein_go_analysis import get_go_definition

class GOIntegrationPipeline:
    def __init__(self,
                 identity_threshold: int = 80,
                 coverage_threshold: int = 80,
                 evalue_threshold: float = 1e-50,
                 topk: int = 2,
                 protrek_threshold: Optional[float] = None,
                 use_protrek: bool = False,
                 use_foldseek: bool = True,
                 foldseek_database: str = "foldseek_db/sp"):
        """
        GO Information Integration Pipeline.

        Args:
            identity_threshold: BLAST identity threshold (0-100).
            coverage_threshold: BLAST coverage threshold (0-100).
            evalue_threshold: BLAST E-value threshold.
            topk: Use the top-k BLAST results (if > 0).
            protrek_threshold: ProTrek score threshold.
            use_protrek: Whether to use second-level ProTrek filtering.
            use_foldseek: Whether to use Foldseek for remote homology search.
            foldseek_database: Path to the Foldseek database.
        """
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.evalue_threshold = evalue_threshold
        self.protrek_threshold = protrek_threshold
        self.use_protrek = use_protrek
        self.use_foldseek = use_foldseek
        self.foldseek_database = foldseek_database
        self.topk = topk
        self.current_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # Two levels above the current file directory
        self.go_info_path = os.path.join(self.current_path, 'data/raw_data/go.json')
        self.protein2go_path = os.path.join(self.current_path, 'data/processed_data/gt_protein2go_sp20250623.json')
        self.pid2seq_path = os.path.join(self.current_path, 'data/processed_data/swissprot_pid2seq.json')

        # Load protein-to-GO mapping data
        self._load_protein_go_dict()
        self._load_pid2seq()

        # If using ProTrek, initialize the model
        if self.use_protrek:
            self._init_protrek_model()

    def _init_protrek_model(self):
        """Initializes the ProTrek model."""
        from model.ProTrek.protrek_trimodal_model import ProTrekTrimodalModel
        import torch

        config = {
            "protein_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/esm2_t33_650M_UR50D",
            "text_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext",
            "structure_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/foldseek_t30_150M",
            "load_protein_pretrained": False,
            "load_text_pretrained": False,
            "from_checkpoint": "Models/ProTrek/weights/ProTrek_650M_UniRef50/ProTrek_650M_UniRef50.pt"
        }

        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.protrek_model = ProTrekTrimodalModel(**config).to(self.device).eval()
        print(f"ProTrek model loaded to device: {self.device}")

    def _load_protein_go_dict(self):
        """Loads protein-to-GO mapping data."""
        self.protein_go_dict = {}
        try:
            with open(self.protein2go_path, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    self.protein_go_dict[data['protein_id']] = data['GO_id']
            print(f"Successfully loaded protein-to-GO mapping data with {len(self.protein_go_dict)} records.")
        except Exception as e:
            print(f"Error loading protein-to-GO mapping data: {str(e)}")
            self.protein_go_dict = {}

    def _load_pid2seq(self):
        """Loads the PID-to-sequence mapping."""
        try:
            with open(self.pid2seq_path, 'r') as f:
                self.pid2seq = json.load(f)
            print(f"Successfully loaded PID-to-sequence mapping data with {len(self.pid2seq)} records.")
        except Exception as e:
            print(f"Error loading PID-to-sequence mapping data: {str(e)}")
            self.pid2seq = {}

    def _get_go_from_uniprot_id(self, uniprot_id: str) -> List[str]:
        """
        Gets GO IDs from a UniProt ID.

        Args:
            uniprot_id: The UniProt ID.

        Returns:
            A list of GO IDs, using the dictionary loaded within the class.
        """
        # Use the dictionary loaded within the class
        return [go_id.split("_")[-1] if "_" in go_id else go_id
                for go_id in self.protein_go_dict.get(uniprot_id, [])]

    def extract_blast_go_ids(self, blast_results: List[Dict], sequence: str) -> List[str]:
        """
        Extracts qualifying GO IDs from BLAST results.

        Args:
            blast_results: A list of BLAST results.
            sequence: The current protein sequence (to avoid self-matching).

        Returns:
            A list of qualifying GO IDs.
        """
        go_ids = []

        if self.topk > 0:
            # Use the top-k strategy
            for result in blast_results[:self.topk]:
                hit_id = result.get('ID', '')
                if self.pid2seq.get(hit_id) == sequence:
                    continue
                go_ids.extend(self._get_go_from_uniprot_id(hit_id))
        else:
            # Use the threshold strategy
            for result in blast_results:
                identity = float(result.get('Identity%', 0))
                coverage = float(result.get('Coverage%', 0))
                evalue = float(result.get('E-value', 1.0))

                # Check if the threshold conditions are met
                if (identity >= self.identity_threshold and
                    coverage >= self.coverage_threshold and
                    evalue <= self.evalue_threshold):

                    # Get the protein_id of this hit
                    hit_id = result.get('ID', '')
                    if self.pid2seq.get(hit_id) == sequence:
                        continue
                    go_ids.extend(self._get_go_from_uniprot_id(hit_id))

        return go_ids

    def extract_foldseek_go_ids(self, foldseek_results: List[Dict], sequence: str) -> List[str]:
        """
        Extracts qualifying GO IDs from Foldseek results.

        Args:
            foldseek_results: A list of Foldseek results.
            sequence: The current protein sequence (to avoid self-matching).

        Returns:
            A list of qualifying GO IDs.
        """
        go_ids = []

        if self.topk > 0:
            # Use the top-k strategy
            for result in foldseek_results[:self.topk]:
                hit_id = result.get('ID', '')
                if self.pid2seq.get(hit_id) == sequence:
                    continue
                go_ids.extend(self._get_go_from_uniprot_id(hit_id))
        else:
            # Use the threshold strategy
            for result in foldseek_results:
                identity = float(result.get('Identity%', 0))
                coverage = float(result.get('Coverage%', 0))
                evalue = float(result.get('E-value', 1.0))

                # Check if the threshold conditions are met
                if (identity >= self.identity_threshold and
                    coverage >= self.coverage_threshold and
                    evalue <= self.evalue_threshold):

                    # Get the protein_id of this hit
                    hit_id = result.get('ID', '')
                    if self.pid2seq.get(hit_id) == sequence:
                        continue
                    go_ids.extend(self._get_go_from_uniprot_id(hit_id))

        return go_ids

    def compare_and_select_by_evalue(self, blast_results: List[Dict], 
                                     foldseek_results: List[Dict],
                                     sequence: str) -> Dict:
        """
        Compare BLAST and Foldseek top-1 results by e-value and select the better one.

        Args:
            blast_results: List of BLAST results.
            foldseek_results: List of Foldseek results.
            sequence: The current protein sequence (to avoid self-matching).

        Returns:
            A dictionary containing:
            - source: "BLAST" or "Foldseek"
            - hit_id: UniProt ID of the selected hit
            - evalue: E-value of the selected hit
            - go_ids: List of GO IDs from the selected hit
        """
        # Find the first non-self hit from BLAST
        blast_hit = None
        blast_evalue = float('inf')
        for result in blast_results:
            hit_id = result.get('ID', '')
            if self.pid2seq.get(hit_id) != sequence:
                blast_hit = result
                try:
                    blast_evalue = float(result.get('E-value', 1.0))
                except (ValueError, TypeError):
                    evalue_str = str(result.get('E-value', '1.0'))
                    blast_evalue = float(evalue_str)
                break

        # Find the first non-self hit from Foldseek
        foldseek_hit = None
        foldseek_evalue = float('inf')
        for result in foldseek_results:
            hit_id = result.get('ID', '')
            if self.pid2seq.get(hit_id) != sequence:
                foldseek_hit = result
                try:
                    foldseek_evalue = float(result.get('E-value', 1.0))
                except (ValueError, TypeError):
                    evalue_str = str(result.get('E-value', '1.0'))
                    foldseek_evalue = float(evalue_str)
                break

        # Compare and select
        if blast_hit is None and foldseek_hit is None:
            return {"source": None, "hit_id": None, "evalue": None, "go_ids": []}
        elif blast_hit is None:
            selected_source = "Foldseek"
            selected_hit = foldseek_hit
            selected_evalue = foldseek_evalue
        elif foldseek_hit is None:
            selected_source = "BLAST"
            selected_hit = blast_hit
            selected_evalue = blast_evalue
        else:
            # Both have hits, compare e-values
            if blast_evalue <= foldseek_evalue:
                selected_source = "BLAST"
                selected_hit = blast_hit
                selected_evalue = blast_evalue
            else:
                selected_source = "Foldseek"
                selected_hit = foldseek_hit
                selected_evalue = foldseek_evalue

        # Get GO IDs from the selected hit
        hit_id = selected_hit.get('ID', '')
        go_ids = self._get_go_from_uniprot_id(hit_id)

        # Format e-value
        if selected_evalue < 0.001:
            evalue_str = f"{selected_evalue:.1e}"
        else:
            evalue_str = str(round(selected_evalue, 4))

        return {
            "source": selected_source,
            "hit_id": hit_id,
            "evalue": evalue_str,
            "go_ids": go_ids
        }

    def first_level_filtering(self, interproscan_info: Dict, blast_info: Dict, 
                             foldseek_info: Dict = None) -> Dict:
        """
        First-level filtering: Combines InterProScan, BLAST, and Foldseek GO information.

        Args:
            interproscan_info: InterProScan results.
            blast_info: BLAST results.
            foldseek_info: Foldseek results (optional).

        Returns:
            A mapping from protein IDs to GO information with sources.
            Format: {
                "protein_id": {
                    "go_ids": ["0006810", "0005524"],
                    "go_sources": {
                        "0006810": {"source": "InterProScan", "evalue": None},
                        "0005524": {"source": "BLAST", "evalue": "1.0e-50"}
                    }
                }
            }
        """
        protein_go_dict = {}

        for protein_id in interproscan_info.keys():
            go_ids = set()
            go_sources = {}

            # Add GO information from InterProScan
            interproscan_gos = interproscan_info[protein_id].get('interproscan_results', {}).get('go_id', [])
            interproscan_gos = [go_id.split(":")[-1] if ":" in go_id else go_id for go_id in interproscan_gos]
            for go_id in interproscan_gos:
                go_ids.add(go_id)
                go_sources[go_id] = {"source": "InterProScan", "evalue": None}

            # Compare BLAST and Foldseek if both are available and use_foldseek is enabled
            if self.use_foldseek and foldseek_info and protein_id in blast_info and protein_id in foldseek_info:
                sequence = blast_info[protein_id]['sequence']
                blast_results = blast_info[protein_id].get('blast_results', [])
                foldseek_results = foldseek_info[protein_id].get('foldseek_results', [])
                
                # Compare and select the best hit by e-value
                selected = self.compare_and_select_by_evalue(blast_results, foldseek_results, sequence)
                
                if selected['source'] is not None:
                    for go_id in selected['go_ids']:
                        clean_go_id = go_id.split(":")[-1] if ":" in go_id else go_id
                        go_ids.add(clean_go_id)
                        go_sources[clean_go_id] = {
                            "source": selected['source'],
                            "evalue": selected['evalue']
                        }
            
            # If not using Foldseek or Foldseek data not available, use BLAST only
            elif protein_id in blast_info:
                sequence = blast_info[protein_id]['sequence']
                blast_results = blast_info[protein_id].get('blast_results', [])
                blast_gos = self.extract_blast_go_ids(blast_results, sequence)
                
                # Get e-value for the first valid hit
                blast_evalue = None
                for result in blast_results[:self.topk if self.topk > 0 else len(blast_results)]:
                    hit_id = result.get('ID', '')
                    if self.pid2seq.get(hit_id) != sequence:
                        blast_evalue = result.get('E-value', None)
                        break
                
                for go_id in blast_gos:
                    clean_go_id = go_id.split(":")[-1] if ":" in go_id else go_id
                    go_ids.add(clean_go_id)
                    if clean_go_id not in go_sources:  # Don't override InterProScan
                        go_sources[clean_go_id] = {
                            "source": "BLAST",
                            "evalue": blast_evalue
                        }

            protein_go_dict[protein_id] = {
                "go_ids": list(go_ids),
                "go_sources": go_sources
            }

        return protein_go_dict

    def calculate_protrek_scores(self, protein_sequences: Dict[str, str],
                                protein_go_dict: Dict[str, List[str]]) -> Dict[str, Dict]:
        """
        Calculates ProTrek scores.

        Args:
            protein_sequences: A dictionary of protein sequences.
            protein_go_dict: Protein-to-GO mapping.

        Returns:
            A dictionary containing GO scores.
        """
        import torch
        results = {}

        for protein_id, go_ids in tqdm(protein_go_dict.items(), desc="Calculating ProTrek scores"):
            if protein_id not in protein_sequences:
                continue

            protein_seq = protein_sequences[protein_id]
            go_scores = {}

            # Get GO definitions
            go_definitions = {}
            for go_id in go_ids:
                definition = get_go_definition(go_id, self.go_info_path)
                if definition:
                    go_definitions[go_id] = definition

            if not go_definitions:
                continue

            try:
                with torch.no_grad():
                    # Calculate protein sequence embeddings
                    seq_emb = self.protrek_model.get_protein_repr([protein_seq])

                    # Calculate text embeddings and similarity scores
                    definitions = list(go_definitions.values())
                    text_embs = self.protrek_model.get_text_repr(definitions)

                    # Calculate similarity scores
                    scores = (seq_emb @ text_embs.T) / self.protrek_model.temperature
                    scores = scores.cpu().numpy().flatten()

                    # Map back to GO IDs
                    for i, go_id in enumerate(go_definitions.keys()):
                        go_scores[go_id] = float(scores[i])

            except Exception as e:
                print(f"Error calculating ProTrek scores for {protein_id}: {str(e)}")
                continue

            results[protein_id] = {
                "protein_id": protein_id,
                "GO_id": go_ids,
                "Clip_score": go_scores
            }

        return results

    def second_level_filtering(self, protrek_results: Dict[str, Dict]) -> Dict[str, List[str]]:
        """
        Second-level filtering: Filters GO terms based on the ProTrek threshold.

        Args:
            protrek_results: ProTrek calculation results.

        Returns:
            The filtered protein-to-GO mapping.
        """
        filtered_results = {}

        for protein_id, data in protrek_results.items():
            clip_scores = data.get('Clip_score', {})
            filtered_gos = []

            for go_id, score in clip_scores.items():
                if score >= self.protrek_threshold:
                    filtered_gos.append(go_id)

            if filtered_gos:
                filtered_results[protein_id] = filtered_gos

        return filtered_results

    def generate_filename(self, base_name: str, is_intermediate: bool = False) -> str:
        """Generates a filename containing parameter information."""
        if self.topk > 0:
            # If using top-k, only include top-k information
            params = f"topk{self.topk}"
        else:
            # Otherwise, use the original parameter combination
            params = f"identity{self.identity_threshold}_coverage{self.coverage_threshold}_evalue{self.evalue_threshold:.0e}"

        if self.use_protrek and self.protrek_threshold is not None:
            params += f"_protrek{self.protrek_threshold}"

        if is_intermediate:
            return f"{base_name}_intermediate_{params}.json"
        else:
            return f"{base_name}_final_{params}.json"

    def run(self, interproscan_info: Dict = None, blast_info: Dict = None,
            foldseek_info: Dict = None,
            interproscan_file: str = None, blast_file: str = None,
            foldseek_file: str = None,
            output_dir: str = "output"):
        """
        Runs the GO integration pipeline.

        Args:
            interproscan_info: A dictionary of InterProScan results.
            blast_info: A dictionary of BLAST results.
            foldseek_info: A dictionary of Foldseek results.
            interproscan_file: Path to the InterProScan results file.
            blast_file: Path to the BLAST results file.
            foldseek_file: Path to the Foldseek results file.
            output_dir: The output directory.
        """
        # Load data
        if interproscan_info is None and interproscan_file:
            with open(interproscan_file, 'r') as f:
                interproscan_info = json.load(f)

        if blast_info is None and blast_file:
            with open(blast_file, 'r') as f:
                blast_info = json.load(f)

        if foldseek_info is None and foldseek_file and self.use_foldseek:
            with open(foldseek_file, 'r') as f:
                foldseek_info = json.load(f)

        if not interproscan_info or not blast_info:
            raise ValueError("Must provide interproscan_info and blast_info data or file paths.")

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        print("Starting first-level filtering...")
        # First-level filtering
        protein_go_dict = self.first_level_filtering(interproscan_info, blast_info, foldseek_info)

        if not self.use_protrek:
            # If not using second-level filtering, save the results directly
            output_file = os.path.join(output_dir, self.generate_filename("go_integration"))
            with open(output_file, 'w') as f:
                for protein_id, go_data in protein_go_dict.items():
                    result = {
                        "protein_id": protein_id, 
                        "GO_id": go_data.get("go_ids", []),
                        "GO_sources": go_data.get("go_sources", {})
                    }
                    f.write(json.dumps(result) + '\n')

            print(f"First-level filtering complete. Results saved to: {output_file}")
            return output_file

        print("Starting second-level filtering...")
        # Extract protein sequences
        protein_sequences = {}
        for protein_id, data in interproscan_info.items():
            protein_sequences[protein_id] = data.get('sequence', '')

        # Convert protein_go_dict to just go_ids for ProTrek
        protein_go_ids_only = {pid: go_data.get("go_ids", []) 
                               for pid, go_data in protein_go_dict.items()}

        # Calculate ProTrek scores
        protrek_results = self.calculate_protrek_scores(protein_sequences, protein_go_ids_only)

        # Save intermediate results
        intermediate_file = os.path.join(output_dir, self.generate_filename("go_integration", is_intermediate=True))
        with open(intermediate_file, 'w') as f:
            for result in protrek_results.values():
                f.write(json.dumps(result) + '\n')

        print(f"ProTrek score calculation complete. Intermediate results saved to: {intermediate_file}")

        # Second-level filtering
        if self.protrek_threshold is not None:
            final_results = self.second_level_filtering(protrek_results)

            # Save final results
            final_file = os.path.join(output_dir, self.generate_filename("go_integration"))
            with open(final_file, 'w') as f:
                for protein_id, go_ids in final_results.items():
                    result = {"protein_id": protein_id, "GO_id": go_ids}
                    f.write(json.dumps(result) + '\n')

            print(f"Second-level filtering complete. Final results saved to: {final_file}")
            return final_file, intermediate_file

        return intermediate_file

def main():
    parser = argparse.ArgumentParser(description="GO Information Integration Pipeline: Integrates InterProScan, BLAST, and Foldseek results, with an optional second-level filtering using ProTrek.")
    parser.add_argument("--interproscan_file", type=str, required=True,
                       help="Path to the InterProScan results file (JSON format).")
    parser.add_argument("--blast_file", type=str, required=True,
                       help="Path to the BLAST results file (JSON format).")
    parser.add_argument("--foldseek_file", type=str, default=None,
                       help="Path to the Foldseek results file (JSON format).")
    parser.add_argument("--identity", type=int, default=80,
                       help="BLAST identity threshold (0-100, default: 80).")
    parser.add_argument("--coverage", type=int, default=80,
                       help="BLAST coverage threshold (0-100, default: 80).")
    parser.add_argument("--evalue", type=float, default=1e-50,
                       help="BLAST E-value threshold (default: 1e-50).")
    parser.add_argument("--topk", type=int, default=2,
                       help="Use the top-k BLAST/Foldseek results (default: 2; set to 0 to use the threshold strategy).")
    parser.add_argument("--protrek_threshold", type=float, default=None,
                       help="ProTrek score threshold (only used if --use_protrek is specified).")
    parser.add_argument("--use_protrek", action="store_true",
                       help="Whether to use second-level ProTrek filtering (requires GPU support).")
    parser.add_argument("--use_foldseek", action="store_true", default=True,
                       help="Whether to use Foldseek for remote homology search (default: True).")
    parser.add_argument("--foldseek_database", type=str, default="foldseek_db/sp",
                       help="Path to the Foldseek database (default: foldseek_db/sp).")
    parser.add_argument("--output_dir", type=str, default="go_integration_results",
                       help="Output directory path (default: go_integration_results).")

    args = parser.parse_args()

    # Create a pipeline instance
    pipeline = GOIntegrationPipeline(
        identity_threshold=args.identity,
        coverage_threshold=args.coverage,
        evalue_threshold=args.evalue,
        topk=args.topk,
        protrek_threshold=args.protrek_threshold,
        use_protrek=args.use_protrek,
        use_foldseek=args.use_foldseek,
        foldseek_database=args.foldseek_database
    )

    # Run the pipeline
    pipeline.run(
        interproscan_file=args.interproscan_file,
        blast_file=args.blast_file,
        foldseek_file=args.foldseek_file,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()