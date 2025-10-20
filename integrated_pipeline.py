import os
import json
import sys
import argparse
from typing import Dict, List, Optional
from pathlib import Path
from tqdm import tqdm

# 添加必要的路径
root_path = os.path.dirname(os.path.abspath(__file__))
print(root_path)
sys.path.append(root_path)
sys.path.append(os.path.join(root_path, "Models/ProTrek"))

# 导入所需模块
from tools.interproscan import InterproScan
from Bio.Blast.Applications import NcbiblastpCommandline
from utils.utils import extract_interproscan_metrics, get_seqnid, extract_blast_metrics, rename_interproscan_keys
from tools.go_integration_pipeline import GOIntegrationPipeline
from utils.generate_protein_prompt import generate_prompt, get_interpro_manager, get_lmdb_connection
from utils.openai_access import call_chatgpt
# 添加MPR导入
from utils.mpr import MultipleProcessRunnerSimplifier

class IntegratedProteinPipeline:
    def __init__(self, 
                 blast_database: str = "uniprot_swissprot",
                 expect_value: float = 0.01,
                 blast_num_threads: int = 256,  # BLAST线程数
                 interproscan_path: str = "interproscan/interproscan-5.75-106.0/interproscan.sh",
                 interproscan_libraries: List[str] = None,
                 go_topk: int = 2,
                 selected_info_types: List[str] = None,
                 pfam_descriptions_path: str = None,
                 go_info_path: str = None,
                 interpro_data_path: str = None,
                 lmdb_path: str = None,
                 n_process_prompt: int = 256,  # 生成prompt的并行进程数
                 n_process_llm: int = 64,  # 生成LLM答案的并行进程数
                 is_enzyme: bool = False,  # 是否是酶
                 skip_protrek_check: bool = False,  # 是否跳过ProTrek检查
                 args: argparse.Namespace = None):
        """
        整合蛋白质分析管道
        
        Args:
            blast_database: BLAST数据库名称
            expect_value: BLAST E-value阈值
            blast_num_threads: BLAST线程数
            interproscan_path: InterProScan程序路径
            interproscan_libraries: InterProScan库列表
            go_topk: GO整合的topk参数
            selected_info_types: prompt生成时选择的信息类型
            pfam_descriptions_path: Pfam描述文件路径
            go_info_path: GO信息文件路径
            interpro_data_path: InterPro数据文件路径
            lmdb_path: LMDB数据库路径
            n_process_prompt: 生成prompt的并行进程数
            n_process_llm: 生成LLM答案的并行进程数
            is_enzyme: 是否是酶
            skip_protrek_check: 是否跳过ProTrek检查
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
        self.selected_info_types = selected_info_types or ['motif', 'go', 'protrek']  # 添加protrek
        self.n_process_prompt = n_process_prompt  # 生成prompt的并行进程数
        self.n_process_llm = n_process_llm  # 生成LLM答案的并行进程数
        self.is_enzyme = is_enzyme  # 是否是酶
        self.skip_protrek_check = skip_protrek_check  # 是否跳过ProTrek检查
        
        # 文件路径配置
        self.pfam_descriptions_path = pfam_descriptions_path
        self.go_info_path = go_info_path
        self.interpro_data_path = interpro_data_path
        self.lmdb_path = lmdb_path
        self.interproscan_info_path = args.interproscan_info_path if args else None
        self.blast_info_path = args.blast_info_path if args else None
        self.protrek_dir = args.protrek_dir if args else None
        
        # 初始化GO整合管道
        self.go_pipeline = GOIntegrationPipeline(topk=self.go_topk)
        
        # 初始化InterPro管理器（如果需要）
        self.interpro_manager = None
        other_types = [t for t in self.selected_info_types if t not in ['motif', 'go', 'protrek']]
        if other_types and self.interpro_data_path:
            try:
                from utils.generate_protein_prompt import get_interpro_manager
                self.interpro_manager = get_interpro_manager(self.interpro_data_path, None)
            except Exception as e:
                print(f"初始化InterPro管理器失败: {str(e)}")
    
    def get_prompt_template(self, selected_info_types=None):
        """
        获取prompt模板，支持可选的信息类型和酶类型
        """
        if selected_info_types is None:
            selected_info_types = ['motif', 'go']

        # 根据is_enzyme参数选择不同的基础prompt
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

        {%- if question %}
        question: \n {{question}}
        {%- endif %}
        """

        return PROMPT_TEMPLATE

    def _check_and_run_protrek(self, input_fasta: str, protrek_dir: str):
        """
        检查protrek_dir是否包含所有需要的序列，如果没有则运行protrek脚本
        
        Args:
            input_fasta: 输入FASTA文件路径
            protrek_dir: ProTrek结果目录
        """
        import subprocess
        from Bio import SeqIO
        
        # 读取fasta文件，获取所有的蛋白质ID
        protein_ids = []
        with open(input_fasta, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                protein_ids.append(record.id)
        
        # 检查哪些蛋白质ID缺少protrek结果
        missing_ids = []
        os.makedirs(protrek_dir, exist_ok=True)
        for pid in protein_ids:
            protrek_file = os.path.join(protrek_dir, f"{pid}.json")
            if not os.path.exists(protrek_file):
                missing_ids.append(pid)
        
        if not missing_ids:
            print(f"✓ 所有蛋白质的ProTrek结果都已存在")
            return
        
        print(f"⚠ 发现 {len(missing_ids)} 个蛋白质缺少ProTrek结果，开始运行ProTrek脚本...")
        
        # 运行run_protrek_text.sh脚本
        script_path = os.path.join(root_path, "scripts", "run_protrek_text.sh")
        
        # 调用脚本，传递必要的参数
        # 不使用capture_output，让输出直接显示到终端
        try:
            print(f"运行命令: bash {script_path} {input_fasta} {protrek_dir}")
            result = subprocess.run(
                ["bash", script_path, input_fasta, protrek_dir],
                check=True
            )
            print("✓ ProTrek脚本运行成功")
        except subprocess.CalledProcessError as e:
            print(f"✗ ProTrek脚本运行失败，返回码: {e.returncode}")
            raise
        except KeyboardInterrupt:
            print("\n⚠ 用户中断了ProTrek脚本")
            print("提示: 你可以稍后手动运行ProTrek脚本，或者在下次运行管道时它会自动检测缺失的结果")
            raise
    
    def step1_run_blast_and_interproscan(self, input_fasta: str, temp_dir: str = "temp") -> tuple:
        """
        步骤1: 运行BLAST和InterProScan分析
        
        Args:
            input_fasta: 输入FASTA文件路径
            temp_dir: 临时文件目录
            
        Returns:
            tuple: (interproscan_info, blast_info)
        """
        print("步骤1: 运行BLAST和InterProScan分析...")
        
        # 创建临时目录
        os.makedirs(temp_dir, exist_ok=True)
        
        # 获取序列字典
        seq_dict = get_seqnid(input_fasta)
        print(f"读取到 {len(seq_dict)} 个序列")
        
        # 运行BLAST
        print("运行BLAST分析...")
        blast_xml = os.path.join(temp_dir, "blast_results.xml")
        blast_cmd = NcbiblastpCommandline(
            query=input_fasta,
            db=self.blast_database,
            out=blast_xml,
            outfmt=5,  # XML格式
            evalue=self.expect_value,
            num_threads=self.blast_num_threads
        )
        blast_cmd()
        
        # 提取BLAST结果
        blast_results = extract_blast_metrics(blast_xml)
        blast_info = {}
        for uid, info in blast_results.items():
            blast_info[uid] = {"sequence": seq_dict[uid], "blast_results": info}
        
        # 运行InterProScan
        print("运行InterProScan分析...")
        interproscan_json = os.path.join(temp_dir, "interproscan_results.json")
        interproscan = InterproScan(self.interproscan_path)
        input_args = {
            "fasta_file": input_fasta,
            "goterms": True,
            "pathways": True,
            "save_dir": interproscan_json
        }
        interproscan.run(**input_args)
        
        # 提取InterProScan结果
        interproscan_results = extract_interproscan_metrics(
            interproscan_json, 
            librarys=self.interproscan_libraries
        )
        
        interproscan_info = {}
        for id, seq in seq_dict.items():
            info = interproscan_results[seq]
            info = rename_interproscan_keys(info)
            interproscan_info[id] = {"sequence": seq, "interproscan_results": info}
        
        # 注意：不在这里清理临时文件，因为可能需要保存到tool_results目录
        # 清理工作在run()方法的finally块中完成
        
        print(f"步骤1完成: 处理了 {len(interproscan_info)} 个蛋白质")
        return interproscan_info, blast_info
    
    def step2_integrate_go_information(self, interproscan_info: Dict, blast_info: Dict) -> Dict:
        """
        步骤2: 整合GO信息
        
        Args:
            interproscan_info: InterProScan结果
            blast_info: BLAST结果
            
        Returns:
            Dict: 整合后的GO信息
        """
        print("步骤2: 整合GO信息...")
        
        # 使用GO整合管道进行第一层筛选
        protein_go_dict = self.go_pipeline.first_level_filtering(interproscan_info, blast_info)
        
        print(f"步骤2完成: 为 {len(protein_go_dict)} 个蛋白质整合了GO信息")
        return protein_go_dict
    
    def step3_generate_prompts(self, interproscan_info: Dict, blast_info: Dict, 
                              protein_go_dict: Dict) -> Dict:
        """
        步骤3: 生成蛋白质prompt
        
        Args:
            interproscan_info: InterProScan结果
            blast_info: BLAST结果
            protein_go_dict: 整合的GO信息
            
        Returns:
            Dict: 蛋白质ID到prompt的映射（如果有lmdb则包含QA对）
        """
        print("步骤3: 生成蛋白质prompt...")
        
        # 创建临时的GO整合文件格式（用于generate_prompt函数）
        temp_go_data = {}
        for protein_id, go_ids in protein_go_dict.items():
            temp_go_data[protein_id] = go_ids
        
        prompts_data = {}
        
        if self.lmdb_path:
            # 如果有lmdb路径，处理QA数据
            from utils.generate_protein_prompt import get_qa_data
            
            global_index = 0
            for protein_id in tqdm(interproscan_info.keys(), desc="生成prompts"):
                # 获取QA对
                qa_pairs = get_qa_data(protein_id, self.lmdb_path)
                
                for qa_pair in qa_pairs:
                    question = qa_pair['question']
                    ground_truth = qa_pair['ground_truth']
                    question_type = qa_pair['question_type']
                    
                    # 生成prompt（需要修改generate_prompt函数以支持内存数据）
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
            # 如果没有lmdb路径，按原来的方式处理
            for protein_id in tqdm(interproscan_info.keys(), desc="生成prompts"):
                prompt = self._generate_prompt_from_memory(
                    protein_id, interproscan_info, temp_go_data
                )
                if prompt:
                    prompts_data[protein_id] = prompt
        
        print(f"步骤3完成: 生成了 {len(prompts_data)} 个prompt")
        return prompts_data
    
    def step3_generate_prompts_parallel(self, interproscan_info: Dict, blast_info: Dict, 
                                      protein_go_dict: Dict) -> Dict:
        """
        步骤3: 并行生成蛋白质prompt
        
        Args:
            interproscan_info: InterProScan结果
            blast_info: BLAST结果
            protein_go_dict: 整合的GO信息
            
        Returns:
            Dict: 蛋白质ID到prompt的映射（如果有lmdb则包含QA对）
        """
        print("步骤3: 并行生成蛋白质prompt...")
        
        # 创建临时的GO整合文件格式（用于generate_prompt函数）
        temp_go_data = {}
        for protein_id, go_ids in protein_go_dict.items():
            temp_go_data[protein_id] = go_ids
        
        # 准备共享数据
        self._shared_interproscan_info = interproscan_info
        self._shared_temp_go_data = temp_go_data
        
        if self.lmdb_path:
            # 如果有lmdb路径，处理QA数据
            from utils.generate_protein_prompt import get_qa_data
            
            # 准备所有需要处理的QA任务
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
            
            print(f"准备并行处理 {len(qa_tasks)} 个QA任务...")
            
            # 使用MPR并行处理
            prompts_data = {}
            
            def process_qa_task(process_id, idx, qa_task, writer):
                try:
                    # 为每个进程初始化lmdb连接
                    if self.lmdb_path:
                        get_lmdb_connection(self.lmdb_path)
                    
                    protein_id = qa_task['protein_id']
                    question = qa_task['question']
                    ground_truth = qa_task['ground_truth']
                    global_index = qa_task['global_index']
                    question_type = qa_task['question_type']
                    
                    # 生成prompt
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
                    print(f"处理QA任务时出错 (protein_id: {qa_task['protein_id']}): {str(e)}")
            
            # 运行并行处理
            mprs = MultipleProcessRunnerSimplifier(
                data=qa_tasks,
                do=process_qa_task,
                n_process=self.n_process_prompt,
                split_strategy="static",
                return_results=True
            )
            
            results = mprs.run()
            
            # 解析结果
            for result_line in results:
                try:
                    result = json.loads(result_line)
                    prompts_data[result['index']] = result
                except Exception as e:
                    print(f"解析结果时出错: {str(e)}")
                    
        else:
            # 如果没有lmdb路径，按蛋白质ID处理
            protein_ids = list(interproscan_info.keys())
            
            print(f"准备并行处理 {len(protein_ids)} 个蛋白质...")
            
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
                    print(f"处理蛋白质时出错 (protein_id: {protein_id}): {str(e)}")
            
            # 运行并行处理
            mprs = MultipleProcessRunnerSimplifier(
                data=protein_ids,
                do=process_protein,
                n_process=self.n_process_prompt,
                split_strategy="static",
                return_results=True
            )
            
            results = mprs.run()
            
            # 解析结果
            prompts_data = {}
            for result_line in results:
                try:
                    result = json.loads(result_line)
                    prompts_data[result['protein_id']] = result['prompt']
                except Exception as e:
                    print(f"解析结果时出错: {str(e)}")
        
        print(f"步骤3完成: 并行生成了 {len(prompts_data)} 个prompt")
        return prompts_data

    def _generate_prompt_from_memory(self, protein_id: str, interproscan_info: Dict, 
                                   protein_go_dict: Dict, question: str = None) -> str:
        """
        从内存中的数据生成prompt，包含完整的motif和GO定义
        """
        try:
            from utils.protein_go_analysis import get_go_definition
            from jinja2 import Template
            
            # 获取GO分析结果
            go_ids = protein_go_dict.get(protein_id, [])
            go_annotations = []
            all_related_definitions = {}
            
            if go_ids:
                for go_id in go_ids:
                    # 确保GO ID格式正确
                    clean_go_id = go_id.split(":")[-1] if ":" in go_id else go_id
                    go_annotations.append({"go_id": clean_go_id})
                    
                    # 获取GO定义
                    if self.go_info_path and os.path.exists(self.go_info_path):
                        definition = get_go_definition(clean_go_id, self.go_info_path)
                        if definition:
                            all_related_definitions[clean_go_id] = definition
            
            # 获取motif信息
            motif_pfam = {}
            if self.pfam_descriptions_path and os.path.exists(self.pfam_descriptions_path):
                try:
                    # 从interproscan结果中提取pfam信息
                    interproscan_results = interproscan_info[protein_id].get('interproscan_results', {})
                    pfam_entries = interproscan_results.get('pfam_id', [])
                    
                    # 加载pfam描述
                    with open(self.pfam_descriptions_path, 'r') as f:
                        pfam_descriptions = json.load(f)
                    
                    # 构建motif_pfam字典
                    for entry in pfam_entries:
                        for pfam_id, ipr_id in entry.items():
                            if pfam_id and pfam_id in pfam_descriptions:
                                motif_pfam[pfam_id] = pfam_descriptions[pfam_id]['description']
                                
                except Exception as e:
                    print(f"获取motif信息时出错: {str(e)}")
            
            # 获取InterPro描述信息
            interpro_descriptions = {}
            other_types = [t for t in self.selected_info_types if t not in ['motif', 'go', 'protrek']]
            if other_types and self.interpro_manager:
                try:
                    interpro_descriptions = self.interpro_manager.get_description(protein_id, other_types)
                except Exception as e:
                    print(f"获取InterPro描述信息时出错: {str(e)}")
            
            # 获取ProTrek描述信息（只在motif和go都没有的时候）
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
                        print(f"获取ProTrek信息时出错: {str(e)}")
            
            # 准备模板数据
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
            
            # 使用模板生成prompt
            PROMPT_TEMPLATE = self.get_prompt_template(self.selected_info_types)
            template = Template(PROMPT_TEMPLATE)
            return template.render(**template_data)
            
        except Exception as e:
            print(f"生成prompt时出错 (protein_id: {protein_id}): {str(e)}")
            # 如果出错，返回简化版本的prompt
            return self._generate_simple_prompt(protein_id, interproscan_info, protein_go_dict, question)
    
    def _generate_simple_prompt(self, protein_id: str, interproscan_info: Dict, 
                               protein_go_dict: Dict, question: str = None) -> str:
        """
        生成简化版本的prompt（作为备用）
        """
        # 获取蛋白质序列
        sequence = interproscan_info[protein_id].get('sequence', '')
        
        # 获取GO信息
        go_ids = protein_go_dict.get(protein_id, [])
        
        # 获取motif信息
        interproscan_results = interproscan_info[protein_id].get('interproscan_results', {})
        pfam_entries = interproscan_results.get('pfam_id', [])
        
        # 简化的prompt生成逻辑
        prompt_parts = []
        
        if self.is_enzyme:
            from utils.prompts import ENZYME_PROMPT
            prompt_parts.append(ENZYME_PROMPT)
        else:
            from utils.prompts import FUNCTION_PROMPT
            prompt_parts.append(FUNCTION_PROMPT)
        
        prompt_parts.append("\ninput information:")
        
        # 添加motif信息
        if 'motif' in self.selected_info_types and pfam_entries:
            prompt_parts.append("\nmotif:")
            for entry in pfam_entries:
                for key, value in entry.items():
                    if value:
                        prompt_parts.append(f"{value}: 无详细描述")
        
        # 添加GO信息
        if 'go' in self.selected_info_types and go_ids:
            prompt_parts.append("\nGO:")
            for i, go_id in enumerate(go_ids[:10], 1):
                prompt_parts.append(f"▢ GO term{i}: {go_id}")
                prompt_parts.append(f"• definition: 无详细定义")
        
        # 添加ProTrek信息（只在motif和go都没有的时候）
        if 'protrek' in self.selected_info_types and not pfam_entries and not go_ids:
            try:
                from utils.get_protrek_text import get_protrek_text
                protrek_descriptions = get_protrek_text(sequence, topk=3)
                if protrek_descriptions:
                    prompt_parts.append("\nprotrek:")
                    for i, (description, score) in enumerate(protrek_descriptions, 1):
                        prompt_parts.append(f"▢ Related description {i}: {description} (matching score: {score})")
            except Exception as e:
                print(f"获取ProTrek信息时出错: {str(e)}")
        
        if question:
            prompt_parts.append(f"\nquestion: \n{question}")
        
        return "\n".join(prompt_parts)
    
    def step4_generate_llm_answers(self, prompts_data: Dict, save_dir: str) -> None:
        """
        步骤4: 生成LLM答案
        
        Args:
            prompts_data: prompt数据
            save_dir: 保存目录
        """
        print("步骤4: 生成LLM答案...")
        
        # 创建保存目录
        os.makedirs(save_dir, exist_ok=True)
        
        if self.lmdb_path:
            # 如果有lmdb路径，处理QA数据
            for index, qa_item in tqdm(prompts_data.items(), desc="生成LLM答案"):
                try:
                    protein_id = qa_item['protein_id']
                    prompt = qa_item['prompt']
                    question = qa_item['question']
                    ground_truth = qa_item['ground_truth']
                    question_type = qa_item.get('question_type', None)
        
                    # 调用LLM生成答案
                    llm_response = call_chatgpt(prompt)
                    if question_type is None:
                        # 构建结果数据
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
                    
                    # 保存文件
                    save_path = os.path.join(save_dir, f"{protein_id}_{index}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)
                        
                except Exception as e:
                    print(f"处理索引 {index} 时出错: {str(e)}")
        else:
            # 如果没有lmdb路径，按原来的方式处理
            for protein_id, prompt in tqdm(prompts_data.items(), desc="生成LLM答案"):
                try:
                    # 调用LLM生成答案
                    llm_response = call_chatgpt(prompt)
                    
                    # 构建结果数据
                    result = {
                        'protein_id': protein_id,
                        'prompt': prompt,
                        'llm_answer': llm_response
                    }
                    
                    # 保存文件
                    save_path = os.path.join(save_dir, f"{protein_id}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)
                        
                except Exception as e:
                    print(f"处理蛋白质 {protein_id} 时出错: {str(e)}")
        
        print(f"步骤4完成: 结果已保存到 {save_dir}")

    def step4_generate_llm_answers_parallel(self, prompts_data: Dict, save_dir: str) -> None:
        """
        步骤4: 并行生成LLM答案
        
        Args:
            prompts_data: prompt数据
            save_dir: 保存目录
        """
        print("步骤4: 并行生成LLM答案...")
        
        # 创建保存目录
        os.makedirs(save_dir, exist_ok=True)
        
        if self.lmdb_path:
            # 如果有lmdb路径，处理QA数据
            qa_items = list(prompts_data.items())
            
            print(f"准备并行处理 {len(qa_items)} 个LLM答案生成任务...")
            
            def process_llm_qa(process_id, idx, qa_item_pair, writer):
                try:
                    index, qa_item = qa_item_pair
                    protein_id = qa_item['protein_id']
                    prompt = qa_item['prompt']
                    question = qa_item['question']
                    ground_truth = qa_item['ground_truth']
                    question_type = qa_item.get('question_type', None)
        
                    # 调用LLM生成答案
                    llm_response = call_chatgpt(prompt)
                    if question_type is None:
                        # 构建结果数据
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
                    
                    # 保存文件
                    save_path = os.path.join(save_dir, f"{protein_id}_{index}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)
                        
                except Exception as e:
                    print(f"处理LLM QA任务时出错 (index: {qa_item_pair[0]}): {str(e)}")
            
            # 运行并行处理
            mprs = MultipleProcessRunnerSimplifier(
                data=qa_items,
                do=process_llm_qa,
                n_process=self.n_process_llm,
                split_strategy="static"
            )
            mprs.run()
            
        else:
            # 如果没有lmdb路径，按原来的方式处理
            protein_items = list(prompts_data.items())
            
            print(f"准备并行处理 {len(protein_items)} 个LLM答案生成任务...")
            
            def process_llm_protein(process_id, idx, protein_item_pair, writer):
                try:
                    protein_id, prompt = protein_item_pair
                    
                    # 调用LLM生成答案
                    llm_response = call_chatgpt(prompt)
                    
                    # 构建结果数据
                    result = {
                        'protein_id': protein_id,
                        'prompt': prompt,
                        'llm_answer': llm_response
                    }
                    
                    # 保存文件
                    save_path = os.path.join(save_dir, f"{protein_id}.json")
                    with open(save_path, 'w') as f:
                        json.dump(result, f, indent=2, ensure_ascii=False)
                        
                except Exception as e:
                    print(f"处理LLM蛋白质任务时出错 (protein_id: {protein_item_pair[0]}): {str(e)}")
            
            # 运行并行处理
            mprs = MultipleProcessRunnerSimplifier(
                data=protein_items,
                do=process_llm_protein,
                n_process=self.n_process_llm,
                split_strategy="static"
            )
            mprs.run()
        
        print(f"步骤4完成: 并行生成的LLM答案已保存到 {save_dir}")

    def run(self, input_fasta: str, output_dir: str, temp_dir: str = "temp"):
        """
        运行完整的工作流
        
        Args:
            input_fasta: 输入FASTA文件路径
            output_dir: 输出目录
            temp_dir: 临时文件目录
        """
        print(f"开始运行整合蛋白质分析管道...")
        print(f"输入文件: {input_fasta}")
        print(f"输出目录: {output_dir}")
        print(f"生成prompt并行进程数: {self.n_process_prompt}")
        print(f"生成LLM答案并行进程数: {self.n_process_llm}")
        
        # 创建输出目录结构
        os.makedirs(output_dir, exist_ok=True)
        tool_results_dir = os.path.join(output_dir, "tool_results")
        llm_answers_dir = os.path.join(output_dir, "llm_answers")
        os.makedirs(tool_results_dir, exist_ok=True)
        os.makedirs(llm_answers_dir, exist_ok=True)
        
        # 如果需要protrek_dir，检查并运行protrek脚本
        if self.protrek_dir and 'protrek' in self.selected_info_types and not self.skip_protrek_check:
            print(f"检查ProTrek目录: {self.protrek_dir}")
            self._check_and_run_protrek(input_fasta, self.protrek_dir)
        elif self.protrek_dir and 'protrek' in self.selected_info_types and self.skip_protrek_check:
            print(f"跳过ProTrek检查（使用已有结果）: {self.protrek_dir}")
        
        try:
            # 步骤1: 运行BLAST和InterProScan
            if self.interproscan_info_path is None or self.blast_info_path is None:
                interproscan_info, blast_info = self.step1_run_blast_and_interproscan(
                    input_fasta, temp_dir
                )
                # 保存中间结果到tool_results目录
                interproscan_save_path = os.path.join(tool_results_dir, "interproscan_info.json")
                blast_save_path = os.path.join(tool_results_dir, "blast_info.json")
                with open(interproscan_save_path, 'w') as f:
                    json.dump(interproscan_info, f, indent=2, ensure_ascii=False)
                with open(blast_save_path, 'w') as f:
                    json.dump(blast_info, f, indent=2, ensure_ascii=False)
                print(f"中间结果已保存到: {tool_results_dir}")
            else:
                interproscan_info = json.load(open(self.interproscan_info_path))
                blast_info = json.load(open(self.blast_info_path))
            
            # 步骤2: 整合GO信息
            protein_go_dict = self.step2_integrate_go_information(
                interproscan_info, blast_info
            )
            
            # 步骤3: 并行生成prompt
            prompts_data = self.step3_generate_prompts_parallel(
                interproscan_info, blast_info, protein_go_dict
            )
            
            # 步骤4: 并行生成LLM答案
            self.step4_generate_llm_answers_parallel(prompts_data, llm_answers_dir)
            
            print("整合管道运行完成！")
            
        except Exception as e:
            print(f"管道运行出错: {str(e)}")
            raise
        finally:
            # 清理临时目录和共享数据
            print(f"清理临时目录: {temp_dir}")
            if os.path.exists(temp_dir):
                import shutil
                shutil.rmtree(temp_dir)
            
            # 清理共享数据
            if hasattr(self, '_shared_interproscan_info'):
                delattr(self, '_shared_interproscan_info')
            if hasattr(self, '_shared_temp_go_data'):
                delattr(self, '_shared_temp_go_data')

def main():
    parser = argparse.ArgumentParser(description="整合蛋白质分析管道")
    parser.add_argument("--input_fasta", type=str, default="examples/input.fasta", help="输入FASTA文件路径")
    parser.add_argument("--output_dir", type=str, required=True, help="输出目录")
    parser.add_argument("--temp_dir", type=str, default="temp", help="临时文件目录")
    parser.add_argument('--interproscan_info_path', type=str, default=None, help="InterProScan结果文件路径（如果提供则跳过BLAST和InterProScan步骤）")
    parser.add_argument('--blast_info_path', type=str, default=None, help="BLAST结果文件路径（如果提供则跳过BLAST和InterProScan步骤）")
    
    # 并行处理参数
    parser.add_argument("--n_process_prompt", type=int, default=256, help="生成prompt的并行进程数")
    parser.add_argument("--n_process_llm", type=int, default=64, help="生成LLM答案的并行进程数")
    
    # BLAST参数
    parser.add_argument("--blast_database", type=str, default="uniprot_swissprot", help="BLAST数据库")
    parser.add_argument("--expect_value", type=float, default=0.01, help="BLAST E-value阈值")
    parser.add_argument("--blast_num_threads", type=int, default=256, help="BLAST线程数")
    
    # InterProScan参数
    parser.add_argument("--interproscan_path", type=str, 
                       default="interproscan/interproscan-5.75-106.0/interproscan.sh",
                       help="InterProScan程序路径")
    
    # GO整合参数
    parser.add_argument("--go_topk", type=int, default=2, help="GO整合topk参数")
    
    # Prompt生成参数
    parser.add_argument("--selected_info_types", type=str, nargs='+', 
                       default=['motif', 'go'], help="选择的信息类型")
    parser.add_argument("--pfam_descriptions_path", type=str, default='data/raw_data/all_pfam_descriptions.json', help="Pfam描述文件路径")
    parser.add_argument("--go_info_path", type=str, default='data/raw_data/go.json', help="GO信息文件路径")
    parser.add_argument("--interpro_data_path", type=str, default='data/raw_data/interpro_data.json', help="InterPro数据文件路径")
    parser.add_argument("--lmdb_path", type=str, default=None, help="LMDB数据库路径（用于QA数据集生成）")
    parser.add_argument("--protrek_dir", type=str, default=None, help="ProTrek结果目录路径")
    parser.add_argument("--is_enzyme", action="store_true", help="是否是酶（决定使用ENZYME_PROMPT还是FUNCTION_PROMPT）")
    parser.add_argument("--skip_protrek_check", action="store_true", help="跳过ProTrek检查，直接使用已有结果")
    
    args = parser.parse_args()
    
    # 创建管道实例
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
        args=args
    )
    
    # 运行管道
    pipeline.run(args.input_fasta, args.output_dir, args.temp_dir)

if __name__ == "__main__":
    main() 