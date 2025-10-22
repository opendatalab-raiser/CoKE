import os
import json
import sys
import tempfile
import gradio as gr
from typing import Dict, List, Optional
from pathlib import Path
from Bio import SeqIO
from io import StringIO
import queue
import threading
import uuid
from datetime import datetime
import shutil

# Add necessary paths
root_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root_path)
sys.path.append(os.path.join(root_path, "Models/ProTrek"))

# Import required modules
from tools.interproscan import InterproScan
from Bio.Blast.Applications import NcbiblastpCommandline
from utils.utils import extract_interproscan_metrics, get_seqnid, extract_blast_metrics, rename_interproscan_keys
from tools.go_integration_pipeline import GOIntegrationPipeline
from utils.openai_access_demo import call_chatgpt
from utils.prompts import FUNCTION_PROMPT,ENZYME_PROMPT
from utils.get_protrek_text import get_protrek_text
import multiprocessing

# 数据保存配置
USER_DATA_DIR = "user_data"  # 用户数据保存目录
os.makedirs(USER_DATA_DIR, exist_ok=True)

total_cores = multiprocessing.cpu_count()
#obtain available memory
try:
    with open('/proc/meminfo', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('MemAvailable:'):
                available_memory_gb = int(line.split()[1]) / 1024 / 1024
                break
        else:
            # 如果没有MemAvailable，用MemFree + Buffers + Cached估算
            memfree = memcached = buffers = 0
            for line in lines:
                if line.startswith('MemFree:'):
                    memfree = int(line.split()[1])
                elif line.startswith('Cached:'):
                    memcached = int(line.split()[1])
                elif line.startswith('Buffers:'):
                    buffers = int(line.split()[1])
            available_memory_gb = (memfree + memcached + buffers) / 1024 / 1024
except:
    available_memory_gb = 64

# 全局资源限制配置
MEMORY_PER_PROCESS_GB = 8
MAX_CONCURRENT_ANALYSES = min(total_cores, available_memory_gb // MEMORY_PER_PROCESS_GB)
print(f"Total cores: {total_cores}, Available memory: {available_memory_gb}GB, MAX_CONCURRENT_ANALYSES: {MAX_CONCURRENT_ANALYSES}")
analysis_semaphore = threading.Semaphore(MAX_CONCURRENT_ANALYSES)
active_analyses = 0  # 当前活跃的分析数量
analyses_lock = threading.Lock()  # 用于保护计数器

def get_prompt_template(selected_info_types=None, is_enzyme=False):
    """
    Get prompt template, supports optional information types, enzyme types and conversation history
    
    Args:
        selected_info_types: List of information types to include, e.g. ['motif', 'go', 'superfamily', 'panther', 'protrek']
        is_enzyme: Whether it is an enzyme protein
    """
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']  # Default includes motif and go information

    # Select different base prompt based on whether it's an enzyme
    base_prompt = ENZYME_PROMPT if is_enzyme else FUNCTION_PROMPT
    
    PROMPT_TEMPLATE = base_prompt + '\n' + """
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

    {%- if conversation_history %}

    conversation history:
    {%- for item in conversation_history %}
    Q: {{item.question}}
    A: {{item.answer}}
    
    {%- endfor %}
    {%- endif %}

    question: \n {{question}}
    """

    return PROMPT_TEMPLATE

class ProteinAnalysisDemo:
    def __init__(self, cpu_cores=10, session_id=None):  # 限制每个分析任务的CPU核心数
        """
        Protein Analysis Demo Class
        """
        self.blast_database = "uniprot_swissprot"
        self.expect_value = 0.01
        self.interproscan_path = "interproscan/interproscan-5.75-106.0/interproscan.sh"
        self.interproscan_libraries = [
            "PFAM", "PIRSR", "PROSITE_PROFILES", "SUPERFAMILY", "PRINTS", 
            "PANTHER", "CDD", "GENE3D", "NCBIFAM", "SFLM", "MOBIDB_LITE", 
            "COILS", "PROSITE_PATTERNS", "FUNFAM", "SMART"
        ]
        self.go_topk = 2
        self.selected_info_types = ['motif', 'go', 'protrek']
        self.is_enzyme = False
        self.cpu_cores = cpu_cores  # 限制CPU使用
        
        # Session management
        self.session_id = session_id if session_id else str(uuid.uuid4())
        self.session_dir = os.path.join(USER_DATA_DIR, self.session_id)
        os.makedirs(self.session_dir, exist_ok=True)
        
        # File path configuration
        self.pfam_descriptions_path = 'data/raw_data/all_pfam_descriptions.json'
        self.go_info_path = 'data/raw_data/go.json'
        self.interpro_data_path = 'data/raw_data/interpro_data.json'
        
        # Initialize GO integration pipeline
        self.go_pipeline = GOIntegrationPipeline(topk=self.go_topk)
        
        # Conversation state management - 每个实例独立的状态
        self.current_protein_data = None
        self.conversation_history = []
        
        # 加载已有的会话数据（如果存在）
        self._load_session_data()
        
        # Initialize InterPro manager (if needed)
        self.interpro_manager = None
        other_types = [t for t in self.selected_info_types if t not in ['motif', 'go', 'protrek']]
        if other_types and os.path.exists(self.interpro_data_path):
            try:
                from utils.generate_protein_prompt import get_interpro_manager
                self.interpro_manager = get_interpro_manager(self.interpro_data_path, None)
            except Exception as e:
                print(f"Failed to initialize InterPro manager: {str(e)}")
    
    def _load_session_data(self):
        """加载已有的会话数据"""
        session_file = os.path.join(self.session_dir, "session_data.json")
        if os.path.exists(session_file):
            try:
                with open(session_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    self.current_protein_data = data.get('current_protein_data')
                    self.conversation_history = data.get('conversation_history', [])
                print(f"Loaded session data for {self.session_id}")
            except Exception as e:
                print(f"Error loading session data: {str(e)}")
    
    def _save_session_data(self):
        """保存当前会话数据"""
        session_file = os.path.join(self.session_dir, "session_data.json")
        try:
            data = {
                'session_id': self.session_id,
                'timestamp': datetime.now().isoformat(),
                'current_protein_data': self.current_protein_data,
                'conversation_history': self.conversation_history
            }
            with open(session_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"Saved session data for {self.session_id}")
        except Exception as e:
            print(f"Error saving session data: {str(e)}")
    
    def export_session_data(self) -> str:
        """导出会话数据为JSON文件"""
        export_file = os.path.join(self.session_dir, f"export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        try:
            data = {
                'session_id': self.session_id,
                'export_time': datetime.now().isoformat(),
                'sequence': self.current_protein_data.get('sequence') if self.current_protein_data else None,
                'sequence_length': len(self.current_protein_data.get('sequence', '')) if self.current_protein_data else 0,
                'analysis_results': {
                    'protein_id': self.current_protein_data.get('protein_id') if self.current_protein_data else None,
                    'sequence_source': self.current_protein_data.get('sequence_source') if self.current_protein_data else None,
                    'is_enzyme': self.current_protein_data.get('is_enzyme') if self.current_protein_data else False,
                },
                'conversation_history': self.conversation_history,
                'total_questions': len(self.conversation_history)
            }
            with open(export_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            return export_file
        except Exception as e:
            print(f"Error exporting session data: {str(e)}")
            return None
    
    def reset_conversation(self):
        """
        Reset conversation state
        """
        self.current_protein_data = None
        self.conversation_history = []
        self._save_session_data()
    
    def validate_protein_sequence(self, sequence: str) -> bool:
        """
        Validate protein sequence format
        """
        if not sequence:
            return False
        
        # Remove whitespace characters
        sequence = sequence.strip().upper()
        
        # Check if it contains valid amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        sequence_chars = set(sequence.replace('\n', '').replace(' ', ''))
        
        return sequence_chars.issubset(valid_aa) and len(sequence) > 0
    
    def create_temp_fasta(self, sequence: str, seq_id: str = "demo_protein") -> str:
        """
        Create temporary FASTA file
        """
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        temp_file.write(f">{seq_id}\n{sequence}\n")
        temp_file.close()
        return temp_file.name
    
    def run_blast_analysis(self, fasta_file: str, temp_dir: str) -> Dict:
        """
        Run BLAST analysis
        """
        blast_xml = os.path.join(temp_dir, "blast_results.xml")
        
        try:
            blast_cmd = NcbiblastpCommandline(
                query=fasta_file,
                db=self.blast_database,
                out=blast_xml,
                outfmt=5,  # XML format
                evalue=self.expect_value,
                num_threads=self.cpu_cores  # 限制BLAST使用的线程数
            )
            blast_cmd()
            
            # Extract BLAST results
            blast_results = extract_blast_metrics(blast_xml)
            
            # Get sequence dictionary
            seq_dict = get_seqnid(fasta_file)
            
            blast_info = {}
            for uid, info in blast_results.items():
                blast_info[uid] = {"sequence": seq_dict[uid], "blast_results": info}
            
            return blast_info
            
        except Exception as e:
            print(f"BLAST analysis error: {str(e)}")
            return {}
        finally:
            if os.path.exists(blast_xml):
                os.remove(blast_xml)
    
    def run_interproscan_analysis(self, fasta_file: str, temp_dir: str) -> Dict:
        """
        Run InterProScan analysis with limited CPU cores
        """
        interproscan_json = os.path.join(temp_dir, "interproscan_results.json")
        
        try:
            # 使用限制的CPU核心数初始化InterProScan
            interproscan = InterproScan(self.interproscan_path, cpu_cores=self.cpu_cores)
            input_args = {
                "fasta_file": fasta_file,
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
            
            # Get sequence dictionary
            seq_dict = get_seqnid(fasta_file)
            
            interproscan_info = {}
            for id, seq in seq_dict.items():
                info = interproscan_results[seq]
                info = rename_interproscan_keys(info)
                interproscan_info[id] = {"sequence": seq, "interproscan_results": info}
            
            return interproscan_info
            
        except Exception as e:
            print(f"InterProScan analysis error: {str(e)}")
            return {}
        finally:
            if os.path.exists(interproscan_json):
                os.remove(interproscan_json)
    
    def generate_prompt(self, protein_id: str, interproscan_info: Dict, 
                       protein_go_dict: Dict, question: str, is_enzyme: bool = False) -> str:
        """
        Generate prompt from data in memory, including complete motif and GO definitions
        """
        try:
            from utils.protein_go_analysis import get_go_definition
            from jinja2 import Template
            
            # Get GO analysis results
            go_ids = protein_go_dict.get(protein_id, [])
            go_annotations = []
            all_related_definitions = {}
            
            if go_ids:
                for go_id in go_ids:
                    # Ensure GO ID format is correct
                    clean_go_id = go_id.split(":")[-1] if ":" in go_id else go_id
                    go_annotations.append({"go_id": clean_go_id})
                    
                    # Get GO definition
                    if os.path.exists(self.go_info_path):
                        definition = get_go_definition(clean_go_id, self.go_info_path)
                        if definition:
                            all_related_definitions[clean_go_id] = definition
            
            # Get motif information
            motif_pfam = {}
            if os.path.exists(self.pfam_descriptions_path):
                try:
                    # Extract pfam information from interproscan results
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
                interpro_descriptions = self.interpro_manager.get_description(protein_id, other_types)
            
            # Get ProTrek description information (only when motif and go are both absent)
            if 'protrek' in self.selected_info_types and not motif_pfam and not go_annotations:
                try:
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
                "conversation_history": self.conversation_history,
                "question": question
            }
            
            # Use template to generate prompt
            PROMPT_TEMPLATE = get_prompt_template(self.selected_info_types, is_enzyme)
            template = Template(PROMPT_TEMPLATE)
            return template.render(**template_data)
            
        except Exception as e:
            print(f"Error generating prompt (protein_id: {protein_id}): {str(e)}")
            return f"Error generating prompt: {str(e)}"
    
    def process_new_protein(self, sequence_input: str, is_enzyme: bool = False) -> tuple:
        """
        Process new protein sequence analysis with resource control
        """
        global active_analyses
        
        # 显示等待提示
        with analyses_lock:
            current_active = active_analyses
        
        if current_active >= MAX_CONCURRENT_ANALYSES:
            wait_msg = f"⏳ System is at capacity ({current_active}/{MAX_CONCURRENT_ANALYSES} analyses running). Your request is queued and will start automatically when resources are available. Please wait..."
        else:
            wait_msg = f"🔄 Starting analysis... ({current_active + 1}/{MAX_CONCURRENT_ANALYSES} concurrent analyses)"
        
        # 阻塞等待信号量（会自动排队）
        analysis_semaphore.acquire(blocking=True)
        
        # 更新活跃分析计数
        with analyses_lock:
            active_analyses += 1
            current_active = active_analyses
        
        try:
            # Reset conversation state
            self.reset_conversation()
            
            # Determine which sequence input to use
            final_sequence = None
            sequence_source = ""
            
            if sequence_input.strip():
                # Use sequence from textbox input
                if self.validate_protein_sequence(sequence_input):
                    final_sequence = sequence_input.strip().upper().replace('\n', '').replace(' ', '')
                    sequence_source = "From text input"
                else:
                    return None, "Invalid sequence format. Please enter a valid protein sequence"
            else:
                return None, "Please enter a protein sequence"
            
            # Create temporary directory and files
            with tempfile.TemporaryDirectory() as temp_dir:
                try:
                    # Create temporary FASTA file
                    temp_fasta = self.create_temp_fasta(final_sequence, "demo_protein")
                    
                    # Run analysis
                    status_msg = f"Sequence source: {sequence_source}\nSequence length: {len(final_sequence)} amino acids\n\nAnalyzing...\n"
                    
                    # Step 1: BLAST and InterProScan analysis
                    blast_info = self.run_blast_analysis(temp_fasta, temp_dir)
                    interproscan_info = self.run_interproscan_analysis(temp_fasta, temp_dir)
                    
                    if not blast_info or not interproscan_info:
                        return None, status_msg + "Analysis failed: Unable to get BLAST or InterProScan results"
                    
                    # Step 2: Integrate GO information
                    protein_go_dict = self.go_pipeline.first_level_filtering(interproscan_info, blast_info)
                    
                    # Save protein data
                    self.current_protein_data = {
                        "protein_id": "demo_protein",
                        "sequence": final_sequence,
                        "sequence_source": sequence_source,
                        "interproscan_info": interproscan_info,
                        "protein_go_dict": protein_go_dict,
                        "is_enzyme": is_enzyme
                    }
                    
                    # 保存会话数据
                    self._save_session_data()
                    
                    return self.current_protein_data, f"{status_msg}Analysis complete! You can now start asking questions.\n\n"
                    
                except Exception as e:
                    return None, f"Error during analysis: {str(e)}"
                finally:
                    # Clean up temporary files
                    if 'temp_fasta' in locals() and os.path.exists(temp_fasta):
                        os.remove(temp_fasta)
        finally:
            # 更新活跃分析计数并释放信号量
            with analyses_lock:
                active_analyses -= 1
            analysis_semaphore.release()
    
    def ask_question(self, question: str, files=None) -> str:
        """
        Ask question about current protein
        """
        if not self.current_protein_data:
            return "Please analyze a protein sequence first"
        
        if not question.strip() and not files:
            return "Please enter your question or upload an image"
        
        try:
            # Generate prompt
            prompt = self.generate_prompt(
                self.current_protein_data["protein_id"],
                self.current_protein_data["interproscan_info"],
                self.current_protein_data["protein_go_dict"],
                question,
                self.current_protein_data["is_enzyme"]
            )
            
            # Call LLM to generate answer
            llm_response = call_chatgpt(prompt, files=files)
            
            # Save to conversation history
            self.conversation_history.append({
                "question": question,
                "answer": llm_response,
                "timestamp": datetime.now().isoformat(),
                "files": [f.name if hasattr(f, 'name') else str(f) for f in files] if files else None
            })
            
            # 保存会话数据
            self._save_session_data()
            
            return llm_response
            
        except Exception as e:
            return f"Error generating answer: {str(e)}"

def create_demo():
    """
    Create Gradio demo interface with per-session state
    """
    with gr.Blocks(title="Protein Function Analysis Demo") as demo:
        gr.Markdown("# 🧬 Protein Function Analysis Demo")
        gr.Markdown("Enter a protein sequence for analysis, then ask multiple questions")
        # gr.Markdown(f"💡 **System Info**: Supports up to {MAX_CONCURRENT_ANALYSES} concurrent analyses. Each analysis uses 10 CPU cores and takes ~60-90 seconds.")

        # 为每个session创建独立的analyzer实例
        analyzer_state = gr.State(None)

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### 📝 Sequence Input")
                sequence_input = gr.Textbox(
                    label="Protein Sequence",
                    placeholder="Enter protein sequence (single letter amino acid code)...",
                    lines=5,
                    max_lines=10
                )
                
                gr.Markdown("### ⚙️ Analysis Options")
                is_enzyme = gr.Checkbox(
                    label="This is an enzyme protein",
                    value=False,
                    info="Check this to use specialized enzyme analysis template"
                )
                
                analyze_btn = gr.Button("🔍 Analyze Protein", variant="primary", size="lg")
                
                analysis_status = gr.Textbox(
                    label="Analysis Status",
                    lines=5,
                    interactive=False
                )
                
                gr.Markdown("### ❓ Question Area")
                question_input = gr.Textbox(
                    label="Your Question",
                    placeholder="Enter your question about this protein...",
                    lines=3
                )
                
                # Add image upload component
                file_upload = gr.Files(
                    label="Upload Images (Optional)",
                    file_types=["image"],
                    file_count="multiple"
                )
                
                ask_btn = gr.Button("💬 Ask", variant="secondary", size="lg")
                
                clear_btn = gr.Button("🔄 Clear Conversation", variant="secondary")
                
                # export_btn = gr.Button("💾 Export Session Data", variant="secondary")
                # download_file = gr.File(label="Download", visible=False)
            
            with gr.Column(scale=2):
                gr.Markdown("### 💬 Conversation History")
                chatbot = gr.Chatbot(
                    label="Conversation",
                    height=1110,
                    show_copy_button=True
                )
        
        # Examples
        gr.Markdown("### 💡 Examples")
        gr.Examples(
            examples=[
                ["MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV", False],
                ["MGQTKSKIKSKYASYLSFIKILLKRGGVKVSTKNLIKLFQIIEQFCPWFPEQGTLDLKDWKRIGKELKQAGRKGNIIPLTVWNDWAIIKAALEPFQTEEDSVSVSDAPGSCIIDCNENTRKKSQKETESLHCEYVAEPVMAQSTQNVDYNQLQEVIYPETLKLEGKVPELVGPSESKPRGTSRLPAGQVPVTLQPQTQVKENKTQPPVAYQYWPPAELQYRPPLESQYGYPGMPPAPQGRAPYPQPPTRRLNPTAPPSRRGSELHEIIDKSRKEGDTEAWQFPVTLEPMPPGEGAQEGEPPTVEARYKSFSIKMLKDMKEGVKQYGPNSPYMRTLLDSIAHGHRLIPYDWEILAKSSLSPSQFLQFKTWWIDGVQEQVRRNRAANPPVNIDADQLLGIGQNWSTISQQALMQNEAIEQVRAICLRAWEKIQDPGSTCPSFNTVRQGSKEPYPDFVARLQDVAQKSIAIEKARKVIVELMAYENPNPECQSAIKPLKGKVPAGSDVISEYVKACDGMGGAMHKAMLMAQAITGVVLGGQVRTFGGKCYNCGQIGHLKKNCPVLNKQNITIQATTTGREPPDLCPRCKKGKHWASQCRSKFDKNGQPLSGNEQRGQPQAPQQTGAFPIQPFVPHGFQGQQPPLSQVFQGISQLPQYNNCPPPQAAVQQ", False],
                ["MHVPQFISTGALLALLARPAAAHTRMFSVWVNGVDQGDGQNVYIRTPPNTDPIKDLASPALACNVKGGEPVPQFVSASAGDKLTFEWYRVKRGDDIIDPSHSGPITTWIAAFTSPTMDGTGPVWSKIHEEGYDASTKSWAVDKLIANKGMWDFTLPSQLKPGKYMLRQEIVAHHESDATFDKNPKRGAQFYPSCVQVDVKGVGGDAVPDQAFDFNKGYKYSDPGIAFDMYTDFDSYPIPGPPVWDAQDEGCCFIDGVDTTSVKEVVKQIICVLK", True]
            ],
            inputs=[sequence_input, is_enzyme]
        )

        # Contact information
        gr.Markdown("---")
        gr.Markdown("### 📧 Contact Us")
        gr.Markdown("If you have any questions or suggestions, please contact us: **zhuangkai@westlake.edu.cn**")
        
        # Session info display
        session_info = gr.Markdown("", visible=False)
        
        # Event handler functions
        def process_protein(sequence_input, is_enzyme, analyzer):
            # 为每个会话创建新的analyzer实例（如果还没有）
            if analyzer is None:
                analyzer = ProteinAnalysisDemo(cpu_cores=10)  # 每个分析限制10核心
            
            protein_data, status = analyzer.process_new_protein(sequence_input, is_enzyme)
            session_msg = f"\n\n📊 **Session Info**: ID `{analyzer.session_id[:8]}...` | Questions: {len(analyzer.conversation_history)}"
            return status, [], analyzer, session_msg  # 返回analyzer以保持会话状态
        
        def handle_question(question, history, files, analyzer):
            if analyzer is None:
                return history, "", None, analyzer, ""
                
            if not question.strip() and not files:
                return history, "", None, analyzer, ""
            
            response = analyzer.ask_question(question, files)
            
            # Display images in chat history
            if files:
                # Pass file information correctly
                file_messages = []
                for file in files:
                    if hasattr(file, 'name'):  # Uploaded file object
                        file_messages.append(file.name)  # Use filename string directly
                    else:  # Other cases
                        file_messages.append(str(file))  # Convert to string
                
                # Merge question and files into one message
                file_message = "\n".join(file_messages)
                history.append([f"{question}\n{file_message}", response])
            else:
                history.append([question, response])
            
            session_msg = f"\n\n📊 **Session Info**: ID `{analyzer.session_id[:8]}...` | Questions: {len(analyzer.conversation_history)}"
            return history, "", None, analyzer, session_msg
        
        def clear_conversation(analyzer):
            if analyzer:
                analyzer.reset_conversation()
            session_msg = f"\n\n📊 **Session Info**: ID `{analyzer.session_id[:8]}...` | Questions: 0"
            return [], "Conversation cleared", analyzer, session_msg
        
        def export_session(analyzer):
            if analyzer is None:
                return None
            
            export_path = analyzer.export_session_data()
            if export_path:
                return export_path
            return None
        
        # Bind events
        analyze_btn.click(
            fn=process_protein,
            inputs=[sequence_input, is_enzyme, analyzer_state],
            outputs=[analysis_status, chatbot, analyzer_state, session_info]
        )
        
        ask_btn.click(
            fn=handle_question,
            inputs=[question_input, chatbot, file_upload, analyzer_state],
            outputs=[chatbot, question_input, file_upload, analyzer_state, session_info]
        )
        
        question_input.submit(
            fn=handle_question,
            inputs=[question_input, chatbot, file_upload, analyzer_state],
            outputs=[chatbot, question_input, file_upload, analyzer_state, session_info]
        )
        
        clear_btn.click(
            fn=clear_conversation,
            inputs=[analyzer_state],
            outputs=[chatbot, analysis_status, analyzer_state, session_info]
        )
        
        # export_btn.click(
        #     fn=export_session,
        #     inputs=[analyzer_state],
        #     outputs=[download_file]
        # )
    
    return demo

if __name__ == "__main__":
    demo = create_demo()
    demo.queue(  # 启用队列系统以支持并发
        default_concurrency_limit=MAX_CONCURRENT_ANALYSES * 2  # 允许更多请求排队
    )
    demo.launch(
        server_name="0.0.0.0",
        server_port=30003,
        share=False,
        debug=True,
        max_threads=50,  # 限制最大线程数
        root_path="/CoKE"  # 使Gradio正确处理子路径
    )

