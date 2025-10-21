import json
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from jinja2 import Template
try:
    from utils.protein_go_analysis import analyze_protein_go
    from utils.prompts import ENZYME_PROMPT, FUNCTION_PROMPT
    from utils.get_motif import get_motif_pfam
except ImportError:
    from protein_go_analysis import analyze_protein_go
    from prompts import ENZYME_PROMPT, FUNCTION_PROMPT
    from get_motif import get_motif_pfam
from tqdm import tqdm

class InterProDescriptionManager:
    """A class to manage InterPro description information to avoid repeatedly reading files."""
    
    def __init__(self, interpro_data_path, interproscan_info_path, protrek_text_dir):
        """
        Initializes and reads all necessary data.
        
        Args:
            interpro_data_path: Path to the interpro_data.json file.
            interproscan_info_path: Path to the interproscan_info.json file.
        """
        self.interpro_data_path = interpro_data_path
        self.interproscan_info_path = interproscan_info_path
        self.protrek_text_dir = protrek_text_dir
        self.interpro_data = None
        self.interproscan_info = None
        self._load_data()
    
    def _load_data(self):
        """Loads data files, executed only once."""
        if self.interpro_data_path and os.path.exists(self.interpro_data_path):
            with open(self.interpro_data_path, 'r') as f:
                self.interpro_data = json.load(f)
        
        if self.interproscan_info_path and os.path.exists(self.interproscan_info_path):
            with open(self.interproscan_info_path, 'r') as f:
                self.interproscan_info = json.load(f)
    
    def get_description(self, protein_id, selected_types=None, protrek_topk=3):
        """
        Retrieves InterPro description information for a protein.
        
        Args:
            protein_id: The protein ID.
            selected_types: A list of information types to retrieve, e.g., ['superfamily', 'panther', 'gene3d', 'protrek'].
            protrek_topk: The number of top results to return from ProTrek.
        
        Returns:
            A dictionary containing description information for each type.
        """
        if selected_types is None:
            selected_types = []
        
        if not self.interpro_data and not self.interproscan_info and 'protrek' not in selected_types:
            return {}
        
        result = {}
        
        # Check if the protein exists (for non-protrek types)
        if self.interproscan_info and protein_id not in self.interproscan_info:
            if 'protrek' not in selected_types:
                return result
            protein_info = None
        else:
            protein_info = self.interproscan_info.get(protein_id) if self.interproscan_info else None
        
        interproscan_results = protein_info.get('interproscan_results', {}) if protein_info else {}
        
        # Iterate over the selected types
        for info_type in selected_types:
            if info_type == 'protrek':
                # Handle protrek type
                if protein_info and 'sequence' in protein_info:
                    try:
                        # from utils.get_protrek_text import get_protrek_text
                        # protrek_descriptions = get_protrek_text(protein_info['sequence'], protrek_topk)
                        
                        with open(f'{self.protrek_text_dir}/{protein_id}.json', 'r') as f:
                            protrek_descriptions = json.load(f)
                        if protrek_descriptions:
                            result['protrek'] = protrek_descriptions[:protrek_topk]
                    except Exception as e:
                        print(f"Error getting ProTrek data for {protein_id}: {e}")
            elif info_type in interproscan_results:
                type_descriptions = {}
                
                # Get all IPR IDs for this type
                for entry in interproscan_results[info_type]:
                    for key, ipr_id in entry.items():
                        if ipr_id and ipr_id in self.interpro_data:
                            type_descriptions[ipr_id] = {
                                'name': self.interpro_data[ipr_id].get('name', ''),
                                'abstract': self.interpro_data[ipr_id].get('abstract', '')
                            }
                
                if type_descriptions:
                    result[info_type] = type_descriptions
        
        return result

# Global variables to cache the InterProDescriptionManager instance and lmdb connection
_interpro_manager = None
_lmdb_db = None
_lmdb_path = None

def get_interpro_manager(interpro_data_path, interproscan_info_path, protrek_text_dir):
    """Gets or creates an InterProDescriptionManager instance."""
    global _interpro_manager
    if _interpro_manager is None:
        _interpro_manager = InterProDescriptionManager(interpro_data_path, interproscan_info_path, protrek_text_dir)
    return _interpro_manager

def get_lmdb_connection(lmdb_path):
    """Gets or creates an lmdb connection."""
    global _lmdb_db, _lmdb_path
    if _lmdb_db is None or _lmdb_path != lmdb_path:
        if _lmdb_db is not None:
            _lmdb_db.close()
        
        if lmdb_path and os.path.exists(lmdb_path):
            import lmdb
            _lmdb_db = lmdb.open(lmdb_path, readonly=True)
            _lmdb_path = lmdb_path
        else:
            _lmdb_db = None
            _lmdb_path = None
    
    return _lmdb_db

def get_prompt_template(selected_info_types=None,lmdb_path=None):
    """
    Gets the prompt template, supporting optional information types.
    
    Args:
        selected_info_types: A list of information types to include, e.g., ['motif', 'go', 'superfamily', 'panther', 'protrek'].
    """
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']  # Default to include motif and go information
    if lmdb_path is None:
        PROMPT_TEMPLATE = ENZYME_PROMPT + "\n"
    else:
        PROMPT_TEMPLATE = FUNCTION_PROMPT + "\n"
    PROMPT_TEMPLATE += """
    input information:

    {%- if 'motif' in selected_info_types and motif_pfam %}

    motif:{% for motif_id, motif_info in motif_pfam.items() %}
    {{motif_id}}: {{motif_info}}
    {% endfor %}
    {%- endif %}

    {%- if 'go' in selected_info_types and go_data.status == 'success' %}

    GO:{% for go_entry in go_data.go_annotations %}
    ▢ GO term{{loop.index}}: {{go_entry.go_id}}
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

    """
    if lmdb_path is not None:
        PROMPT_TEMPLATE += "\n" + "question: \n {{question}}"
    return PROMPT_TEMPLATE

def get_qa_data(protein_id, lmdb_path):
    """
    Retrieves all QA pairs for a specific protein from lmdb.
    
    Args:
        protein_id: The protein ID.
        lmdb_path: The path to the lmdb database.
    
    Returns:
        A list of QA pairs, where each element contains a question and ground_truth.
    """
    if not lmdb_path or not os.path.exists(lmdb_path):
        return []
    
    import json
    
    qa_pairs = []
    
    try:
        db = get_lmdb_connection(lmdb_path)
        if db is None:
            return []
        
        with db.begin() as txn:
            # Iterate through data with numeric indices to find matching protein_id
            cursor = txn.cursor()
            for key, value in cursor:
                try:
                    # Try to decode the key as a number (data with numeric indices)
                    key_str = key.decode('utf-8')
                    if key_str.isdigit():
                        # This is numeric-indexed data, containing protein_id, question, ground_truth
                        data = json.loads(value.decode('utf-8'))
                        if isinstance(data, list) and len(data) == 3:
                            stored_protein_id, question, ground_truth = data[0], data[1], data[2]
                            if stored_protein_id == protein_id:
                                qa_pairs.append({
                                    'question': question,
                                    'ground_truth': ground_truth,
                                })
                        elif isinstance(data, list) and len(data) == 4:
                            stored_protein_id, question, ground_truth, question_type = data[0], data[1], data[2], data[3]
                            if stored_protein_id == protein_id:
                                qa_pairs.append({
                                    'question': question,
                                    'ground_truth': ground_truth,
                                    'question_type': question_type
                                })
                except Exception as e:
                    # Skip this entry if parsing fails
                    continue
    except Exception as e:
        print(f"Error reading lmdb for protein {protein_id}: {e}")
    
    return qa_pairs

def generate_prompt(protein_id, protein2gopath, protrek_text_dir, protein2pfam_path, pfam_descriptions_path, go_info_path, 
                   interpro_data_path=None, interproscan_info_path=None, selected_info_types=None, lmdb_path=None, interpro_manager=None, question=None, protrek_topk=3):
    """
    Generates a protein prompt.
    
    Args:
        selected_info_types: A list of information types to include, e.g., ['motif', 'go', 'superfamily', 'panther', 'protrek'].
        interpro_data_path: Path to the interpro_data.json file.
        interproscan_info_path: Path to the interproscan_info.json file.
        interpro_manager: An InterProDescriptionManager instance, used if provided.
        question: The question text for QA tasks.
        protrek_topk: The number of top results to return from ProTrek.
    """
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']
    
    # Get analysis results
    analysis = analyze_protein_go(protein_id, protein2gopath, go_info_path)
    motif_pfam = get_motif_pfam(protein_id, protein2pfam_path, pfam_descriptions_path)
    
    # Get InterPro and ProTrek description information
    interpro_descriptions = {}
    other_types = [t for t in selected_info_types if t not in ['motif', 'go']]
    if other_types:
        if interpro_manager:
            # Use the provided manager instance
            interpro_descriptions = interpro_manager.get_description(protein_id, other_types, protrek_topk)
        elif interpro_data_path and interproscan_info_path:
            # Use the globally cached manager
            manager = get_interpro_manager(interpro_data_path, interproscan_info_path, protrek_text_dir)
            interpro_descriptions = manager.get_description(protein_id, other_types, protrek_topk)

    # Prepare template data
    template_data = {
        "protein_id": protein_id,
        "selected_info_types": selected_info_types,
        "go_data": {
            "status": analysis["status"],
            "go_annotations": analysis["go_annotations"] if analysis["status"] == "success" else [],
            "all_related_definitions": analysis["all_related_definitions"] if analysis["status"] == "success" else {}
        },
        "motif_pfam": motif_pfam,
        "interpro_descriptions": interpro_descriptions,
        "question": question
    }

    PROMPT_TEMPLATE = get_prompt_template(selected_info_types,lmdb_path)
    template = Template(PROMPT_TEMPLATE)
    return template.render(**template_data)

def save_prompts_parallel(protein_ids, output_path, protein2gopath, protrek_text_dir, protein2pfam_path, pfam_descriptions_path, go_info_path, 
                         interpro_data_path=None, interproscan_info_path=None, selected_info_types=None, lmdb_path=None, n_process=8, protrek_topk=3):
    """Generates and saves protein prompts in parallel."""
    import json
    try:
        from utils.mpr import MultipleProcessRunnerSimplifier
    except ImportError:
        from mpr import MultipleProcessRunnerSimplifier
    
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']
    
    # Create an InterProDescriptionManager instance before starting parallel processing
    interpro_manager = None
    other_types = [t for t in selected_info_types if t not in ['motif', 'go']]
    if other_types and interpro_data_path and interproscan_info_path:
        interpro_manager = InterProDescriptionManager(interpro_data_path, interproscan_info_path, protrek_text_dir)
    
    # Shared variable for tracking the global index
    # if lmdb_path:
    #     import multiprocessing
    #     global_index = multiprocessing.Value('i', 0)  # Shared integer, initial value 0
    #     index_lock = multiprocessing.Lock()  # Lock for synchronized access
    # else:
    #     global_index = None
    #     index_lock = None
    import multiprocessing
    global_index = multiprocessing.Value('i', 0)  # Shared integer, initial value 0
    index_lock = multiprocessing.Lock()  # Lock for synchronized access
    results = {}
    
    def process_protein(process_id, idx, protein_id, writer):
        protein_id = protein_id.strip()
        
        # Initialize lmdb connection for each process
        if lmdb_path:
            get_lmdb_connection(lmdb_path)
        
        if lmdb_path:
            # If lmdb_path is provided, process QA data
            qa_pairs = get_qa_data(protein_id, lmdb_path)
            for qa_pair in qa_pairs:
                question = qa_pair.get('question', None)
                ground_truth = qa_pair.get('ground_truth',None)
                question_type = qa_pair.get('question_type', None)
                prompt = generate_prompt(protein_id, protein2gopath, protrek_text_dir, protein2pfam_path, pfam_descriptions_path, go_info_path,
                                       interpro_data_path, interproscan_info_path, selected_info_types, lmdb_path, interpro_manager, question, protrek_topk)
                if prompt == "":
                    continue
                if writer:
                    # Get and increment the global index
                    with index_lock:
                        current_index = global_index.value
                        global_index.value += 1
                    if question_type is None:
                        result = {
                            "index": current_index,
                            "protein_id": protein_id, 
                            "prompt": prompt, 
                            "question": question, 
                            "ground_truth": ground_truth
                        }
                    else:
                        result = {
                            "index": current_index,
                            "protein_id": protein_id, 
                            "prompt": prompt, 
                            "question": question, 
                            "ground_truth": ground_truth,
                            "question_type": question_type
                        }
                    writer.write(json.dumps(result) + '\n')
        else:
            # If no lmdb_path, process as before
            prompt = generate_prompt(protein_id, protein2gopath, protrek_text_dir, protein2pfam_path, pfam_descriptions_path, go_info_path,
                                   interpro_data_path, interproscan_info_path, selected_info_types, lmdb_path, interpro_manager, protrek_topk=protrek_topk)
            if prompt == "":
                return
            if writer:
                with index_lock:
                    current_index = global_index.value
                    global_index.value += 1
                result = {
                    "index": current_index,
                    "protein_id": protein_id, 
                    "prompt": prompt, 
                    "question_type": "ec_number"
                }
                writer.write(json.dumps(result) + '\n')
    
    # Use MultipleProcessRunnerSimplifier for parallel processing
    runner = MultipleProcessRunnerSimplifier(
        data=protein_ids,
        do=process_protein,
        save_path=output_path + '.tmp',
        n_process=n_process,
        split_strategy="static"
    )
    
    runner.run()
    
    # Clean up the global lmdb connection
    global _lmdb_db
    if _lmdb_db is not None:
        _lmdb_db.close()
        _lmdb_db = None
    
    # if not lmdb_path:
    #     # If no lmdb_path, merge all results into a single dictionary (for backward compatibility)
    #     final_results = {}
    #     with open(output_path + '.tmp', 'r') as f:
    #         for line in f:
    #             if line.strip():  # Ignore empty lines
    #                 final_results.update(json.loads(line))
        
    #     # Save the final result in the correct JSON format
    #     with open(output_path, 'w') as f:
    #         json.dump(final_results, f, indent=2)
    # else:
    #     # If lmdb_path is provided, save directly in jsonl format
    #     import shutil
    #     shutil.move(output_path + '.tmp', output_path)
    import shutil
    shutil.move(output_path + '.tmp', output_path)
    
    # Remove the temporary file (if it still exists)
    if os.path.exists(output_path + '.tmp'):
        os.remove(output_path + '.tmp')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate prompts for protein function analysis')
    parser.add_argument('--protein_path', type=str, required=True, 
                       help='Path to the protein ID list file (one ID per line)')
    parser.add_argument('--protein2pfam_path', type=str, required=True, 
                       help='Path to the protein-Pfam mapping file (InterProScan results JSON file)')
    parser.add_argument('--protrek_text_dir', type=str, default=None, 
                       help='Path to the ProTrek text description directory (if using protrek info type)')
    parser.add_argument('--pfam_descriptions_path', type=str, required=True, 
                       help='Path to the Pfam descriptions file (JSON file with Pfam IDs and descriptions)')
    parser.add_argument('--protein2gopath', type=str, required=True, 
                       help='Path to the protein-GO mapping file (JSON Lines format)')
    parser.add_argument('--go_info_path', type=str, required=True, 
                       help='Path to the GO ontology file (go.json)')
    parser.add_argument('--interpro_data_path', type=str, default=None, 
                       help='Path to the InterPro data file (if using interpro-related info types)')
    parser.add_argument('--interproscan_info_path', type=str, default=None, 
                       help='Path to the InterProScan results file (if using interpro-related info types)')
    parser.add_argument('--lmdb_path', type=str, default=None, 
                       help='Path to the LMDB database (if generating prompts in QA format)')
    parser.add_argument('--output_path', type=str, required=True, 
                       help='Output file path (JSON or JSONL format)')
    parser.add_argument('--selected_info_types', type=str, nargs='+', default=['motif', 'go'], 
                       help='Select the information types to include (default: motif go). Options: motif, go, protrek')
    parser.add_argument('--n_process', type=int, default=8, 
                       help='Number of parallel processes (default: 8)')
    parser.add_argument('--protrek_topk', type=int, default=3, 
                       help='Number of top results to return from ProTrek (default: 3)')
    args = parser.parse_args()
    # Update output_path to include selected_info_types
    args.output_path = args.output_path.replace('.json', '_' + '_'.join(args.selected_info_types) + '.json')
    print(args)

    with open(args.protein_path, 'r') as file:
        protein_ids = file.readlines()
    
    save_prompts_parallel(
        protein_ids=protein_ids,
        output_path=args.output_path,
        n_process=args.n_process,
        protein2gopath=args.protein2gopath,
        protrek_text_dir=args.protrek_text_dir,
        protein2pfam_path=args.protein2pfam_path,
        pfam_descriptions_path=args.pfam_descriptions_path,
        go_info_path=args.go_info_path,
        interpro_data_path=args.interpro_data_path,
        interproscan_info_path=args.interproscan_info_path,
        selected_info_types=args.selected_info_types,
        lmdb_path=args.lmdb_path,
        protrek_topk=args.protrek_topk
    )

    # Test example
    # protein_id = 'A8CF74'
    # prompt = generate_prompt(protein_id, 'data/processed_data/go_integration_final_topk2.json', 
    #                         'data/processed_data/interproscan_info.json', 'data/raw_data/all_pfam_descriptions.json', 
    #                         'data/raw_data/go.json', 'data/raw_data/interpro_data.json', 
    #                         'data/processed_data/interproscan_info.json', 
    #                         ['motif', 'go', 'superfamily', 'panther'])
    # print(prompt)