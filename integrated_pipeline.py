import os
import json
import sys
import argparse
from typing import Dict, List, Optional
from pathlib import Path
from tqdm import tqdm

# Add necessary paths
root_path = os.path.dirname(os.path.abspath(__file__))
print(root_path)
sys.path.append(root_path)
sys.path.append(os.path.join(root_path, "Models/ProTrek"))

# Import required modules
from tools.interproscan import InterproScan
from Bio.Blast.Applications import NcbiblastpCommandline
from utils.utils import extract_interproscan_metrics, get_seqnid, extract_blast_metrics, rename_interproscan_keys
from tools.go_integration_pipeline import GOIntegrationPipeline
from utils.generate_protein_prompt import generate_prompt, get_interpro_manager, get_lmdb_connection
from utils.openai_access import call_chatgpt
# Add MPR import
from utils.mpr import MultipleProcessRunnerSimplifier

class IntegratedProteinPipeline:
    def __init__(self,
                 blast_database: str = "uniprot_swissprot",
                 expect_value: float = 0.01,
                 blast_num_threads: int = 256,  # Number of BLAST threads
                 interproscan_path: str = "interproscan/interproscan-5.75-106.0/interproscan.sh",
                 interproscan_libraries: List[str] = None,
                 go_topk: int = 2,
                 selected_info_types: List[str] = None,
                 pfam_descriptions_path: str = None,
                 go_info_path: str = None,
                 interpro_data_path: str = None,
                 lmdb_path: str = None,
                 n_process_prompt: int = 256,  # Number of parallel processes for prompt generation
                 n_process_llm: int = 64,  # Number of parallel processes for LLM answer generation
                 is_enzyme: bool = False,  # Whether the protein is an enzyme
                 skip_protrek_check: bool = False,  # Whether to skip the ProTrek check
                 use_foldseek: bool = True,  # Whether to use Foldseek
                 foldseek_database: str = "foldseek_db/sp",  # Foldseek database path
                 foldseek_num_threads: int = 64,  # Number of Foldseek threads
                 args: argparse.Namespace = None):
        """
        Integrated protein analysis pipeline.

        Args:
            blast_database: BLAST database name.
            expect_value: BLAST E-value threshold.
            blast_num_threads: Number of threads for BLAST.
            interproscan_path: Path to the InterProScan executable.
            interproscan_libraries: List of InterProScan libraries to use.
            go_topk: Top-k parameter for GO integration.
            selected_info_types: Information types to include in the prompt generation.
            pfam_descriptions_path: Path to the Pfam descriptions file.
            go_info_path: Path to the GO information file.
            interpro_data_path: Path to the InterPro data file.
            lmdb_path: Path to the LMDB database.
            n_process_prompt: Number of parallel processes for prompt generation.
            n_process_llm: Number of parallel processes for LLM answer generation.
            is_enzyme: Whether the protein is an enzyme.
            skip_protrek_check: Whether to skip checking for ProTrek results.
            use_foldseek: Whether to use Foldseek for remote homology search.
            foldseek_database: Path to the Foldseek database.
            foldseek_num_threads: Number of threads for Foldseek.
        """
        self.blast_database = blast_database
        self.expect_value = expect_value
        self.blast_num_threads = blast_num_threads
        self.interproscan_path = interproscan_path
        self.interproscan_libraries = interproscan_libraries or [
            "PFAM", "PIRSR", "PROSITE_PROFILES", "SUPERFAMILY", "PRINTS",
            "PANTHER", "CDD", "GENE3D", "NCBIFAM", "SFLM", "MOBIDB_LITE",
            "COILS", "PROSITE_PATTERNS", "FUNFAM", "SMART"
        ]
        self.go_topk = go_topk
        self.selected_info_types = selected_info_types or ['motif', 'go', 'protrek']  # Add protrek
        self.n_process_prompt = n_process_prompt
        self.n_process_llm = n_process_llm
        self.is_enzyme = is_enzyme
        self.skip_protrek_check = skip_protrek_check
        self.use_foldseek = use_foldseek
        self.foldseek_database = foldseek_database
        self.foldseek_num_threads = foldseek_num_threads

        # File path configuration
        self.pfam_descriptions_path = pfam_descriptions_path
        self.go_info_path = go_info_path
        self.interpro_data_path = interpro_data_path
        self.lmdb_path = lmdb_path
        self.interproscan_info_path = args.interproscan_info_path if args else None
        self.blast_info_path = args.blast_info_path if args else None
        self.protrek_dir = args.protrek_dir if args else None

        # Initialize the GO integration pipeline
        self.go_pipeline = GOIntegrationPipeline(
            topk=self.go_topk,
            use_foldseek=self.use_foldseek,
            foldseek_database=self.foldseek_database
        )

        # Initialize the InterPro manager (if needed)
        self.interpro_manager = None
        other_types = [t for t in self.selected_info_types if t not in ['motif', 'go', 'protrek']]
        if other_types and self.interpro_data_path:
            try:
                from utils.generate_protein_prompt import get_interpro_manager
                self.interpro_manager = get_interpro_manager(self.interpro_data_path, None)
            except Exception as e:
                print(f"Failed to initialize InterPro manager: {str(e)}")

    def get_prompt_template(self, selected_info_types=None):
        """
        Gets the prompt template, supporting optional information types and enzyme designation.
        """
        if selected_info_types is None:
            selected_info_types = ['motif', 'go']

        # Select a different base prompt based on the is_enzyme parameter
        if self.is_enzyme:
            from utils.prompts import ENZYME_PROMPT
            base_prompt = ENZYME_PROMPT
        else:
            from utils.prompts import FUNCTION_PROMPT
            base_prompt = FUNCTION_PROMPT

        PROMPT_TEMPLATE = base_prompt + '\n' + """
        input information:

        {%- if 'motif' in selected_info_types and motif_pfam %}

        motif:{% for motif_id, motif_info in motif_pfam.items() %}
        {{motif_id}}: {{motif_info}}
        {% endfor %}
        {%- endif %}

        {%- if 'go' in selected_info_types and go_data.status == 'success' %}

        GO:{% for go_entry in go_data.go_annotations %}
        ▢ GO term{{loop.index}}: {{go_entry.go_id}} (来源: {{go_entry.source}}, E-value: {{go_entry.evalue}})
        • definition: {{ go_data.all_related_definitions.get(go_entry.go_id, 'not found definition') }}
        {% endfor %}
        {%- endif %}

        {%- for info_type in selected_info_types %}
        {%- if info_type == 'protrek' and interpro_descriptions.get('protrek') and not motif_pfam and go_data.status == 'error' %}

        protrek is a tool that can be used to predict the function of a protein.
        To be specific, ProTrek is a tri-modal protein language model that jointly models protein sequence, structure and function (SSF).
        It employs contrastive learning with three core alignment strategies:
        (1) using structure as the supervision signal for AA sequences and vice versa,
        (2) mutual supervision between sequences and functions, and
        (3) mutual supervision between structures and functions.
        This tri-modal alignment training enables ProTrek to tightly associate SSF by bringing genuine sample pairs (sequence-structure, sequence-function, and structure-function) closer together while pushing negative samples farther apart in the latent space.
        These are the protrek output descriptions of the protein, the matching score is the score of the description in the protrek output.
        The higher the matching score, the more relevant the description is to the protein. 15 represents that the description is quite relevant to the protein.
        protrek:{% for description,score in interpro_descriptions.protrek %}
        ▢ Related description {{loop.index}}: {{description}} (matching score: {{score}})
        {% endfor %}
        {%- elif info_type not in ['motif', 'go', 'protrek'] and interpro_descriptions.get(info_type) %}

        {{info_type}}:{% for ipr_id, ipr_info in interpro_descriptions[info_type].items() %}
        ▢ {{ipr_id}}: {{ipr_info.name}}
        • description: {{ipr_info.abstract}}
        {% endfor %}
        {%- endif %}
        {%- endfor %}

        {%- if question %}
        question: \n {{question}}
        {%- endif %}
        """

        return PROMPT_TEMPLATE

    def _check_and_run_protrek(self, input_fasta: str, protrek_dir: str):
        """
        Checks if protrek_dir contains results for all required sequences; if not, runs the ProTrek script.

        Args:
            input_fasta: Path to the input FASTA file.
            protrek_dir: Directory for ProTrek results.
        """
        import subprocess
        from Bio import SeqIO

        # Read all protein IDs from the FASTA file
        protein_ids = []
        with open(input_fasta, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                protein_ids.append(record.id)

        # Check which protein IDs are missing ProTrek results
        missing_ids = []
        os.makedirs(protrek_dir, exist_ok=True)
        for pid in protein_ids:
            protrek_file = os.path.join(protrek_dir, f"{pid}.json")
            if not os.path.exists(protrek_file):
                missing_ids.append(pid)

        if not missing_ids:
            print(f"✓ ProTrek results already exist for all proteins.")
            return

        print(f"⚠ Found {len(missing_ids)} proteins missing ProTrek results. Starting ProTrek script...")

        # Run the run_protrek_text.sh script
        script_path = os.path.join(root_path, "scripts", "run_protrek_text.sh")

        # Call the script, passing necessary parameters
        # Do not use capture_output, so output is displayed directly in the terminal
        try:
            print(f"Running command: bash {script_path} {input_fasta} {protrek_dir}")
            result = subprocess.run(
                ["bash", script_path, input_fasta, protrek_dir],
                check=True
            )
            print("✓ ProTrek script ran successfully.")
        except subprocess.CalledProcessError as e:
            print(f"✗ ProTrek script failed with return code: {e.returncode}")
            raise
        except KeyboardInterrupt:
            print("\n⚠ User interrupted the ProTrek script.")
            print("Hint: You can run the ProTrek script manually later, or it will be automatically detected on the next pipeline run.")
            raise

    def step1_run_blast_and_interproscan(self, input_fasta: str, temp_dir: str = "temp") -> tuple:
        """
        Step 1: Run BLAST, InterProScan, and optionally Foldseek analysis.

        Args:
            input_fasta: Path to the input FASTA file.
            temp_dir: Directory for temporary files.

        Returns:
            A tuple: (interproscan_info, blast_info, foldseek_info).
        """
        print("Step 1: Running BLAST, InterProScan, and Foldseek analysis...")

        # Create temporary directory
        os.makedirs(temp_dir, exist_ok=True)

        # Get sequence dictionary
        seq_dict = get_seqnid(input_fasta)
        print(f"Read {len(seq_dict)} sequences.")

        # Run BLAST
        print("Running BLAST analysis...")
        blast_xml = os.path.join(temp_dir, "blast_results.xml")
        blast_cmd = NcbiblastpCommandline(
            query=input_fasta,
            db=self.blast_database,
            out=blast_xml,
            outfmt=5,  # XML format
            evalue=self.expect_value,
            num_threads=self.blast_num_threads
        )
        blast_cmd()

        # Extract BLAST results
        blast_results = extract_blast_metrics(blast_xml)
        blast_info = {}
        for uid, info in blast_results.items():
            blast_info[uid] = {"sequence": seq_dict[uid], "blast_results": info}

        # Run Foldseek if enabled
        foldseek_info = {}
        if self.use_foldseek:
            print("Running Foldseek analysis...")
            try:
                from tools.foldseek import run_foldseek_analysis
                foldseek_output = os.path.join(temp_dir, "foldseek_results.m8")
                foldseek_info = run_foldseek_analysis(
                    fasta_file=input_fasta,
                    database=self.foldseek_database,
                    evalue=self.expect_value,
                    num_threads=self.foldseek_num_threads,
                    output_file=foldseek_output,
                    temp_dir=os.path.join(temp_dir, "foldseek_tmp")
                )
                print(f"Foldseek analysis complete: Found results for {len(foldseek_info)} proteins.")
            except Exception as e:
                print(f"Warning: Foldseek analysis failed: {str(e)}")
                print("Continuing without Foldseek results...")
                foldseek_info = {}

        # Run InterProScan
        print("Running InterProScan analysis...")
        interproscan_json = os.path.join(temp_dir, "interproscan_results.json")
        interproscan = InterproScan(self.interproscan_path)
        input_args = {
            "fasta_file": input_fasta,
            "goterms": True,
            "pathways": True,
            "save_dir": interproscan_json
        }
        interproscan.run(**input_args)

        # Extract InterProScan results
        interproscan_results = extract_interproscan_metrics(
            interproscan_json,
            librarys=self.interproscan_libraries
        )

        interproscan_info = {}
        for id, seq in seq_dict.items():
            info = interproscan_results[seq]
            info = rename_interproscan_keys(info)
            interproscan_info[id] = {"sequence": seq, "interproscan_results": info}

        # Note: Do not clean up temp files here, as they might be saved to the tool_results directory.
        # Cleanup is handled in the finally block of the run() method.

        print(f"Step 1 complete: Processed {len(interproscan_info)} proteins.")
        return interproscan_info, blast_info, foldseek_info

    def step2_integrate_go_information(self, interproscan_info: Dict, blast_info: Dict, 
                                      foldseek_info: Dict = None) -> Dict:
        """
        Step 2: Integrate GO information.

        Args:
            interproscan_info: InterProScan results.
            blast_info: BLAST results.
            foldseek_info: Foldseek results (optional).

        Returns:
            A dictionary of integrated GO information.
        """
        print("Step 2: Integrating GO information...")

        # Use the GO integration pipeline for first-level filtering
        protein_go_dict = self.go_pipeline.first_level_filtering(interproscan_info, blast_info, foldseek_info)

        print(f"Step 2 complete: Integrated GO information for {len(protein_go_dict)} proteins.")
        return protein_go_dict

    def step3_generate_prompts(self, interproscan_info: Dict, blast_info: Dict,
                              protein_go_dict: Dict) -> Dict:
        """
        Step 3: Generate protein prompts.

        Args:
            interproscan_info: InterProScan results.
            blast_info: BLAST results.
            protein_go_dict: Integrated GO information.

        Returns:
            A dictionary mapping protein IDs to prompts (or QA pairs if lmdb is used).
        """
        print("Step 3: Generating protein prompts...")

        # Create a temporary GO integration file format (for the generate_prompt function)
        temp_go_data = {}
        for protein_id, go_ids in protein_go_dict.items():
            temp_go_data[protein_id] = go_ids

        prompts_data = {}

        if self.lmdb_path:
            # If an lmdb path is provided, process QA data
            from utils.generate_protein_prompt import get_qa_data

            global_index = 0
            for protein_id in tqdm(interproscan_info.keys(), desc="Generating prompts"):
                # Get QA pairs
                qa_pairs = get_qa_data(protein_id, self.lmdb_path)

                for qa_pair in qa_pairs:
                    question = qa_pair['question']
                    ground_truth = qa_pair['ground_truth']
                    question_type = qa_pair['question_type']

                    # Generate prompt (requires modifying generate_prompt to support in-memory data)
                    prompt = self._generate_prompt_from_memory(
                        protein_id, interproscan_info, temp_go_data, question
                    )

                    if prompt:
                        prompts_data[global_index] = {
                            "index": global_index,
                            "protein_id": protein_id,
                            "prompt": prompt,
                            "question": question,
                            "ground_truth": ground_truth,
                            "question_type": question_type
                        }
                        global_index += 1
        else:
            # If no lmdb path, process as before
            for protein_id in tqdm(interproscan_info.keys(), desc="Generating prompts"):
                prompt = self._generate_prompt_from_memory(
                    protein_id, interproscan_info, temp_go_data
                )
                if prompt:
                    prompts_data[protein_id] = prompt

        print(f"Step 3 complete: Generated {len(prompts_data)} prompts.")
        return prompts_data

    def step3_generate_prompts_parallel(self, interproscan_info: Dict, blast_info: Dict,
                                      protein_go_dict: Dict) -> Dict:
        """
        Step 3: Generate protein prompts in parallel.

        Args:
            interproscan_info: InterProScan results.
            blast_info: BLAST results.
            protein_go_dict: Integrated GO information.

        Returns:
            A dictionary mapping protein IDs to prompts (or QA pairs if lmdb is used).
        """
        print("Step 3: Generating protein prompts in parallel...")

        # Create a temporary GO integration file format (for the generate_prompt function)
        temp_go_data = {}
        for protein_id, go_ids in protein_go_dict.items():
            temp_go_data[protein_id] = go_ids

        # Prepare shared data
        self._shared_interproscan_info = interproscan_info
        self._shared_temp_go_data = temp_go_data

        if self.lmdb_path:
            # If an lmdb path is provided, process QA data
            from utils.generate_protein_prompt import get_qa_data

            # Prepare all QA tasks to be processed
            qa_tasks = []
            global_index = 0

            for protein_id in interproscan_info.keys():
                qa_pairs = get_qa_data(protein_id, self.lmdb_path)
                for qa_pair in qa_pairs:
                    qa_tasks.append({
                        'global_index': global_index,
                        'protein_id': protein_id,
                        'question': qa_pair['question'],
                        'ground_truth': qa_pair['ground_truth'],
                        'question_type': qa_pair['question_type']
                    })
                    global_index += 1

            print(f"Preparing to process {len(qa_tasks)} QA tasks in parallel...")

            # Use MPR for parallel processing
            prompts_data = {}

            def process_qa_task(process_id, idx, qa_task, writer):
                try:
                    # Initialize lmdb connection for each process
                    if self.lmdb_path:
                        get_lmdb_connection(self.lmdb_path)

                    protein_id = qa_task['protein_id']
                    question = qa_task['question']
                    ground_truth = qa_task['ground_truth']
                    global_index = qa_task['global_index']
                    question_type = qa_task['question_type']

                    # Generate prompt
                    prompt = self._generate_prompt_from_memory(
                        protein_id, self._shared_interproscan_info,
                        self._shared_temp_go_data, question
                    )

                    if prompt:
                        result = {
                            "index": global_index,
                            "protein_id": protein_id,
                            "prompt": prompt,
                            "question": question,
                            "ground_truth": ground_truth,
                            "question_type": question_type
                        }

                        if writer:
                            writer.write(json.dumps(result) + '\n')

                except Exception as e:
                    print(f"Error processing QA task (protein_id: {qa_task['protein_id']}): {str(e)}")

            # Run parallel processing
            mprs = MultipleProcessRunnerSimplifier(
                data=qa_tasks,
                do=process_qa_task,
                n_process=self.n_process_prompt,
                split_strategy="static",
                return_results=True
            )

            results = mprs.run()

            # Parse results
            for result_line in results:
                try:
                    result = json.loads(result_line)
                    prompts_data[result['index']] = result
                except Exception as e:
                    print(f"Error parsing result: {str(e)}")

        else:
            # If no lmdb path, process by protein ID
            protein_ids = list(interproscan_info.keys())

            print(f"Preparing to process {len(protein_ids)} proteins in parallel...")

            def process_protein(process_id, idx, protein_id, writer):
                try:
                    prompt = self._generate_prompt_from_memory(
                        protein_id, self._shared_interproscan_info,
                        self._shared_temp_go_data
                    )

                    if prompt and writer:
                        result = {
                            "protein_id": protein_id,
                            "prompt": prompt
                        }
                        writer.write(json.dumps(result) + '\n')

                except Exception as e:
                    print(f"Error processing protein (protein_id: {protein_id}): {str(e)}")

            # Run parallel processing
            mprs = MultipleProcessRunnerSimplifier(
                data=protein_ids,
                do=process_protein,
                n_process=self.n_process_prompt,
                split_strategy="static",
                return_results=True
            )

            results = mprs.run()

            # Parse results
            prompts_data = {}
            for result_line in results:
                try:
                    result = json.loads(result_line)
                    prompts_data[result['protein_id']] = result['prompt']
                except Exception as e:
                    print(f"Error parsing result: {str(e)}")

        print(f"Step 3 complete: Generated {len(prompts_data)} prompts in parallel.")
        return prompts_data

    def _generate_prompt_from_memory(self, protein_id: str, interproscan_info: Dict,
                                   protein_go_dict: Dict, question: str = None) -> str:
        """
        Generates a prompt from in-memory data, including complete motif and GO definitions.
        """
        try:
            from utils.protein_go_analysis import get_go_definition
            from jinja2 import Template

            # Get GO analysis results (new format with sources)
            go_data_raw = protein_go_dict.get(protein_id, {})
            go_ids = go_data_raw.get("go_ids", []) if isinstance(go_data_raw, dict) else go_data_raw
            go_sources = go_data_raw.get("go_sources", {}) if isinstance(go_data_raw, dict) else {}
            
            go_annotations = []
            all_related_definitions = {}

            if go_ids:
                for go_id in go_ids:
                    # Ensure correct GO ID format
                    clean_go_id = go_id.split(":")[-1] if ":" in go_id else go_id
                    
                    # Get source and e-value information
                    source_info = go_sources.get(clean_go_id, {})
                    source = source_info.get("source", "Unknown")
                    evalue = source_info.get("evalue", None)
                    
                    go_annotations.append({
                        "go_id": clean_go_id,
                        "source": source,
                        "evalue": evalue if evalue else "N/A"
                    })

                    # Get GO definition
                    if self.go_info_path and os.path.exists(self.go_info_path):
                        definition = get_go_definition(clean_go_id, self.go_info_path)
                        if definition:
                            all_related_definitions[clean_go_id] = definition

            # Get motif information
            motif_pfam = {}
            if self.pfam_descriptions_path and os.path.exists(self.pfam_descriptions_path):
                try:
                    # Extract pfam info from interproscan results
                    interproscan_results = interproscan_info[protein_id].get('interproscan_results', {})
                    pfam_entries = interproscan_results.get('pfam_id', [])

                    # Load pfam descriptions
                    with open(self.pfam_descriptions_path, 'r') as f:
                        pfam_descriptions = json.load(f)

                    # Build motif_pfam dictionary
                    for entry in pfam_entries:
                        for pfam_id, ipr_id in entry.items():
                            if pfam_id and pfam_id in pfam_descriptions:
                                motif_pfam[pfam_id] = pfam_descriptions[pfam_id]['description']

                except Exception as e:
                    print(f"Error getting motif information: {str(e)}")

            # Get InterPro description information
            interpro_descriptions = {}
            other_types = [t for t in self.selected_info_types if t not in ['motif', 'go', 'protrek']]
            if other_types and self.interpro_manager:
                try:
                    interpro_descriptions = self.interpro_manager.get_description(protein_id, other_types)
                except Exception as e:
                    print(f"Error getting InterPro description information: {str(e)}")

            # Get ProTrek description information (only if motif and GO are absent)
            if 'protrek' in self.selected_info_types and not motif_pfam and not go_annotations:
                if self.protrek_dir is not None:
                    protrek_path = os.path.join(self.protrek_dir, f"{protein_id}.json")
                    with open(protrek_path, 'r') as f:
                        protrek_descriptions = json.load(f)
                    if protrek_descriptions:
                        interpro_descriptions['protrek'] = protrek_descriptions
                else:
                    try:
                        from utils.get_protrek_text import get_protrek_text
                        sequence = interproscan_info[protein_id]['sequence']
                        protrek_descriptions = get_protrek_text(sequence, topk=3)
                        if protrek_descriptions:
                            interpro_descriptions['protrek'] = protrek_descriptions
                    except Exception as e:
                        print(f"Error getting ProTrek information: {str(e)}")

            # Prepare template data
            template_data = {
                "protein_id": protein_id,
                "selected_info_types": self.selected_info_types,
                "go_data": {
                    "status": "success" if go_annotations else "error",
                    "go_annotations": go_annotations,
                    "all_related_definitions": all_related_definitions
                },
                "motif_pfam": motif_pfam,
                "interpro_descriptions": interpro_descriptions,
                "question": question
            }

            # Generate prompt using the template
            PROMPT_TEMPLATE = self.get_prompt_template(self.selected_info_types)
            template = Template(PROMPT_TEMPLATE)
            return template.render(**template_data)

        except Exception as e:
            print(f"Error generating prompt (protein_id: {protein_id}): {str(e)}")
            # If an error occurs, return a simplified version of the prompt
            return self._generate_simple_prompt(protein_id, interproscan_info, protein_go_dict, question)

    def _generate_simple_prompt(self, protein_id: str, interproscan_info: Dict,
                               protein_go_dict: Dict, question: str = None) -> str:
        """
        Generates a simplified version of the prompt (as a fallback).
        """
        # Get protein sequence
        sequence = interproscan_info[protein_id].get('sequence', '')

        # Get GO information
        go_ids = protein_go_dict.get(protein_id, [])

        # Get motif information
        interproscan_results = interproscan_info[protein_id].get('interproscan_results', {})
        pfam_entries = interproscan_results.get('pfam_id', [])

        # Simplified prompt generation logic
        prompt_parts = []

        if self.is_enzyme:
            from utils.prompts import ENZYME_PROMPT
            prompt_parts.append(ENZYME_PROMPT)
        else:
            from utils.prompts import FUNCTION_PROMPT
            prompt_parts.append(FUNCTION_PROMPT)

        prompt_parts.append("\ninput information:")

        # Add motif information
        if 'motif' in self.selected_info_types and pfam_entries:
            prompt_parts.append("\nmotif:")
            for entry in pfam_entries:
                for key, value in entry.items():
                    if value:
                        prompt_parts.append(f"{value}: No detailed description")

        # Add GO information
        if 'go' in self.selected_info_types and go_ids:
            prompt_parts.append("\nGO:")
            for i, go_id in enumerate(go_ids[:10], 1):
                prompt_parts.append(f"▢ GO term{i}: {go_id}")
                prompt_parts.append(f"• definition: No detailed definition")

        # Add ProTrek information (only if motif and GO are absent)
        if 'protrek' in self.selected_info_types and not pfam_entries and not go_ids:
            try:
                from utils.get_protrek_text import get_protrek_text
                protrek_descriptions = get_protrek_text(sequence, topk=3)
                if protrek_descriptions:
                    prompt_parts.append("\nprotrek:")
                    for i, (description, score) in enumerate(protrek_descriptions, 1):
                        prompt_parts.append(f"▢ Related description {i}: {description} (matching score: {score})")
            except Exception as e:
                print(f"Error getting ProTrek information: {str(e)}")

        if question:
            prompt_parts.append(f"\nquestion: \n{question}")

        return "\n".join(prompt_parts)

    def step4_generate_llm_answers(self, prompts_data: Dict, save_dir: str) -> None:
        """
        Step 4: Generate LLM answers.

        Args:
            prompts_data: The prompt data.
            save_dir: The directory to save results.
        """
        print("Step 4: Generating LLM answers...")

        # Create save directory
        os.makedirs(save_dir, exist_ok=True)

        if self.lmdb_path:
            # If an lmdb path is provided, process QA data
            for index, qa_item in tqdm(prompts_data.items(), desc="Generating LLM answers"):
                try:
                    protein_id = qa_item['protein_id']
                    prompt = qa_item['prompt']
                    question = qa_item['question']
                    ground_truth = qa_item['ground_truth']
                    question_type = qa_item.get('question_type', None)

                    # Call LLM to generate an answer
                    llm_response = call_chatgpt(prompt)
                    if question_type is None:
                        # Build result data
                        result = {
                            'protein_id': protein_id,
                            'index': index,
                            'question': question,
                            'ground_truth': ground_truth,
                            'llm_answer': llm_response
                        }
                    else:
                        result = {
                            'protein_id': protein_id,
                            'index': index,
                            'question': question,
                            'ground_truth': ground_truth,
                            'llm_answer': llm_response,
                            'question_type': question_type
                        }

                    # Save the file
                    save_path = os.path.join(save_dir, f"{protein_id}_{index}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)

                except Exception as e:
                    print(f"Error processing index {index}: {str(e)}")
        else:
            # If no lmdb path, process as before
            for protein_id, prompt in tqdm(prompts_data.items(), desc="Generating LLM answers"):
                try:
                    # Call LLM to generate an answer
                    llm_response = call_chatgpt(prompt)

                    # Build result data
                    result = {
                        'protein_id': protein_id,
                        'prompt': prompt,
                        'llm_answer': llm_response
                    }

                    # Save the file
                    save_path = os.path.join(save_dir, f"{protein_id}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)

                except Exception as e:
                    print(f"Error processing protein {protein_id}: {str(e)}")

        print(f"Step 4 complete: Results saved to {save_dir}")

    def step4_generate_llm_answers_parallel(self, prompts_data: Dict, save_dir: str) -> None:
        """
        Step 4: Generate LLM answers in parallel.

        Args:
            prompts_data: The prompt data.
            save_dir: The directory to save results.
        """
        print("Step 4: Generating LLM answers in parallel...")

        # Create save directory
        os.makedirs(save_dir, exist_ok=True)

        if self.lmdb_path:
            # If an lmdb path is provided, process QA data
            qa_items = list(prompts_data.items())

            print(f"Preparing to process {len(qa_items)} LLM answer generation tasks in parallel...")

            def process_llm_qa(process_id, idx, qa_item_pair, writer):
                try:
                    index, qa_item = qa_item_pair
                    protein_id = qa_item['protein_id']
                    prompt = qa_item['prompt']
                    question = qa_item['question']
                    ground_truth = qa_item['ground_truth']
                    question_type = qa_item.get('question_type', None)

                    # Call LLM to generate an answer
                    llm_response = call_chatgpt(prompt)
                    if question_type is None:
                        # Build result data
                        result = {
                            'protein_id': protein_id,
                            'index': index,
                            'question': question,
                            'ground_truth': ground_truth,
                            'llm_answer': llm_response
                        }
                    else:
                        result = {
                            'protein_id': protein_id,
                            'index': index,
                            'question': question,
                            'ground_truth': ground_truth,
                            'llm_answer': llm_response,
                            'question_type': question_type
                        }

                    # Save the file
                    save_path = os.path.join(save_dir, f"{protein_id}_{index}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)

                except Exception as e:
                    print(f"Error processing LLM QA task (index: {qa_item_pair[0]}): {str(e)}")

            # Run parallel processing
            mprs = MultipleProcessRunnerSimplifier(
                data=qa_items,
                do=process_llm_qa,
                n_process=self.n_process_llm,
                split_strategy="static"
            )
            mprs.run()

        else:
            # If no lmdb path, process as before
            protein_items = list(prompts_data.items())

            print(f"Preparing to process {len(protein_items)} LLM answer generation tasks in parallel...")

            def process_llm_protein(process_id, idx, protein_item_pair, writer):
                try:
                    protein_id, prompt = protein_item_pair

                    # Call LLM to generate an answer
                    llm_response = call_chatgpt(prompt)

                    # Build result data
                    result = {
                        'protein_id': protein_id,
                        'prompt': prompt,
                        'llm_answer': llm_response
                    }

                    # Save the file
                    save_path = os.path.join(save_dir, f"{protein_id}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)

                except Exception as e:
                    print(f"Error processing LLM protein task (protein_id: {protein_item_pair[0]}): {str(e)}")

            # Run parallel processing
            mprs = MultipleProcessRunnerSimplifier(
                data=protein_items,
                do=process_llm_protein,
                n_process=self.n_process_llm,
                split_strategy="static"
            )
            mprs.run()

        print(f"Step 4 complete: Parallel-generated LLM answers saved to {save_dir}")

    def run(self, input_fasta: str, output_dir: str, temp_dir: str = "temp"):
        """
        Runs the complete workflow.

        Args:
            input_fasta: Path to the input FASTA file.
            output_dir: Output directory.
            temp_dir: Directory for temporary files.
        """
        print(f"Starting the integrated protein analysis pipeline...")
        print(f"Input file: {input_fasta}")
        print(f"Output directory: {output_dir}")
        print(f"Number of parallel processes for prompt generation: {self.n_process_prompt}")
        print(f"Number of parallel processes for LLM answer generation: {self.n_process_llm}")

        # Create output directory structure
        os.makedirs(output_dir, exist_ok=True)
        tool_results_dir = os.path.join(output_dir, "tool_results")
        llm_answers_dir = os.path.join(output_dir, "llm_answers")
        os.makedirs(tool_results_dir, exist_ok=True)
        os.makedirs(llm_answers_dir, exist_ok=True)

        # If protrek_dir is needed, check and run the protrek script
        if self.protrek_dir and 'protrek' in self.selected_info_types and not self.skip_protrek_check:
            print(f"Checking ProTrek directory: {self.protrek_dir}")
            self._check_and_run_protrek(input_fasta, self.protrek_dir)
        elif self.protrek_dir and 'protrek' in self.selected_info_types and self.skip_protrek_check:
            print(f"Skipping ProTrek check (using existing results): {self.protrek_dir}")

        try:
            # Step 1: Run BLAST, InterProScan, and Foldseek
            if self.interproscan_info_path is None or self.blast_info_path is None:
                interproscan_info, blast_info, foldseek_info = self.step1_run_blast_and_interproscan(
                    input_fasta, temp_dir
                )
                # Save intermediate results to the tool_results directory
                interproscan_save_path = os.path.join(tool_results_dir, "interproscan_info.json")
                blast_save_path = os.path.join(tool_results_dir, "blast_info.json")
                foldseek_save_path = os.path.join(tool_results_dir, "foldseek_info.json")
                with open(interproscan_save_path, 'w') as f:
                    json.dump(interproscan_info, f, indent=2, ensure_ascii=False)
                with open(blast_save_path, 'w') as f:
                    json.dump(blast_info, f, indent=2, ensure_ascii=False)
                if foldseek_info:
                    with open(foldseek_save_path, 'w') as f:
                        json.dump(foldseek_info, f, indent=2, ensure_ascii=False)
                print(f"Intermediate results saved to: {tool_results_dir}")
            else:
                interproscan_info = json.load(open(self.interproscan_info_path))
                blast_info = json.load(open(self.blast_info_path))
                foldseek_info = {}

            # Step 2: Integrate GO information
            protein_go_dict = self.step2_integrate_go_information(
                interproscan_info, blast_info, foldseek_info
            )

            # Step 3: Generate prompts in parallel
            prompts_data = self.step3_generate_prompts_parallel(
                interproscan_info, blast_info, protein_go_dict
            )

            # Step 4: Generate LLM answers in parallel
            self.step4_generate_llm_answers_parallel(prompts_data, llm_answers_dir)

            print("Integrated pipeline run complete!")

        except Exception as e:
            print(f"Pipeline run failed: {str(e)}")
            raise
        finally:
            # Clean up temporary directory and shared data
            print(f"Cleaning up temporary directory: {temp_dir}")
            if os.path.exists(temp_dir):
                import shutil
                shutil.rmtree(temp_dir)

            # Clean up shared data
            if hasattr(self, '_shared_interproscan_info'):
                delattr(self, '_shared_interproscan_info')
            if hasattr(self, '_shared_temp_go_data'):
                delattr(self, '_shared_temp_go_data')

def main():
    parser = argparse.ArgumentParser(description="Integrated Protein Analysis Pipeline")
    parser.add_argument("--input_fasta", type=str, default="examples/input.fasta", help="Path to the input FASTA file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
    parser.add_argument("--temp_dir", type=str, default="temp", help="Directory for temporary files.")
    parser.add_argument('--interproscan_info_path', type=str, default=None, help="Path to InterProScan results file (skips BLAST and InterProScan steps if provided).")
    parser.add_argument('--blast_info_path', type=str, default=None, help="Path to BLAST results file (skips BLAST and InterProScan steps if provided).")

    # Parallel processing parameters
    parser.add_argument("--n_process_prompt", type=int, default=256, help="Number of parallel processes for prompt generation.")
    parser.add_argument("--n_process_llm", type=int, default=64, help="Number of parallel processes for LLM answer generation.")

    # BLAST parameters
    parser.add_argument("--blast_database", type=str, default="uniprot_swissprot", help="BLAST database.")
    parser.add_argument("--expect_value", type=float, default=0.01, help="BLAST E-value threshold.")
    parser.add_argument("--blast_num_threads", type=int, default=256, help="Number of threads for BLAST.")

    # InterProScan parameters
    parser.add_argument("--interproscan_path", type=str,
                       default="interproscan/interproscan-5.75-106.0/interproscan.sh",
                       help="Path to the InterProScan executable.")

    # Foldseek parameters
    parser.add_argument("--use_foldseek", action="store_true", default=True,
                       help="Whether to use Foldseek for remote homology search (default: True).")
    parser.add_argument("--foldseek_database", type=str, default="foldseek_db/sp",
                       help="Path to the Foldseek database (default: foldseek_db/sp).")
    parser.add_argument("--foldseek_num_threads", type=int, default=64,
                       help="Number of threads for Foldseek (default: 64).")

    # GO integration parameters
    parser.add_argument("--go_topk", type=int, default=2, help="Top-k parameter for GO integration.")

    # Prompt generation parameters
    parser.add_argument("--selected_info_types", type=str, nargs='+',
                       default=['motif', 'go'], help="Selected information types.")
    parser.add_argument("--pfam_descriptions_path", type=str, default='data/raw_data/all_pfam_descriptions.json', help="Path to the Pfam descriptions file.")
    parser.add_argument("--go_info_path", type=str, default='data/raw_data/go.json', help="Path to the GO information file.")
    parser.add_argument("--interpro_data_path", type=str, default='data/raw_data/interpro_data.json', help="Path to the InterPro data file.")
    parser.add_argument("--lmdb_path", type=str, default=None, help="Path to the LMDB database (for QA dataset generation).")
    parser.add_argument("--protrek_dir", type=str, default=None, help="Path to the ProTrek results directory.")
    parser.add_argument("--is_enzyme", action="store_true", help="Indicates if the protein is an enzyme (determines whether to use ENZYME_PROMPT or FUNCTION_PROMPT).")
    parser.add_argument("--skip_protrek_check", action="store_true", help="Skip the ProTrek check and use existing results directly.")

    args = parser.parse_args()

    # Create a pipeline instance
    pipeline = IntegratedProteinPipeline(
        blast_database=args.blast_database,
        expect_value=args.expect_value,
        blast_num_threads=args.blast_num_threads,
        interproscan_path=args.interproscan_path,
        go_topk=args.go_topk,
        selected_info_types=args.selected_info_types,
        pfam_descriptions_path=args.pfam_descriptions_path,
        go_info_path=args.go_info_path,
        interpro_data_path=args.interpro_data_path,
        lmdb_path=args.lmdb_path,
        n_process_prompt=args.n_process_prompt,
        n_process_llm=args.n_process_llm,
        is_enzyme=args.is_enzyme,
        skip_protrek_check=args.skip_protrek_check,
        use_foldseek=args.use_foldseek,
        foldseek_database=args.foldseek_database,
        foldseek_num_threads=args.foldseek_num_threads,
        args=args
    )

    # Run the pipeline
    pipeline.run(args.input_fasta, args.output_dir, args.temp_dir)

if __name__ == "__main__":
    main()