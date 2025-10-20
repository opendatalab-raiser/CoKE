import json
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from jinja2 import Template
try:
    from utils.protein_go_analysis import analyze_protein_go
    from utils.prompts import ENZYME_PROMPT_SIMPLE, ENZYME_PROMPT, RELATION_SEMANTIC_PROMPT, FUNCTION_PROMPT
    from utils.get_motif import get_motif_pfam
except ImportError:
    from protein_go_analysis import analyze_protein_go
    from prompts import ENZYME_PROMPT, RELATION_SEMANTIC_PROMPT, FUNCTION_PROMPT
    from get_motif import get_motif_pfam
from tqdm import tqdm

class InterProDescriptionManager:
    """管理InterPro描述信息的类，避免重复读取文件"""
    
    def __init__(self, interpro_data_path, interproscan_info_path, protrek_text_dir):
        """
        初始化时读取所有需要的数据
        
        Args:
            interpro_data_path: interpro_data.json文件路径
            interproscan_info_path: interproscan_info.json文件路径
        """
        self.interpro_data_path = interpro_data_path
        self.interproscan_info_path = interproscan_info_path
        self.protrek_text_dir = protrek_text_dir
        self.interpro_data = None
        self.interproscan_info = None
        self._load_data()
    
    def _load_data(self):
        """加载数据文件，只执行一次"""
        if self.interpro_data_path and os.path.exists(self.interpro_data_path):
            with open(self.interpro_data_path, 'r') as f:
                self.interpro_data = json.load(f)
        
        if self.interproscan_info_path and os.path.exists(self.interproscan_info_path):
            with open(self.interproscan_info_path, 'r') as f:
                self.interproscan_info = json.load(f)
    
    def get_description(self, protein_id, selected_types=None, protrek_topk=3):
        """
        获取蛋白质的InterPro描述信息
        
        Args:
            protein_id: 蛋白质ID
            selected_types: 需要获取的信息类型列表，如['superfamily', 'panther', 'gene3d', 'protrek']
            protrek_topk: ProTrek返回的top结果数量
        
        Returns:
            dict: 包含各类型描述信息的字典
        """
        if selected_types is None:
            selected_types = []
        
        if not self.interpro_data and not self.interproscan_info and 'protrek' not in selected_types:
            return {}
        
        result = {}
        
        # 检查蛋白质是否存在（对于非protrek类型）
        if self.interproscan_info and protein_id not in self.interproscan_info:
            if 'protrek' not in selected_types:
                return result
            protein_info = None
        else:
            protein_info = self.interproscan_info.get(protein_id) if self.interproscan_info else None
        
        interproscan_results = protein_info.get('interproscan_results', {}) if protein_info else {}
        
        # 遍历选定的类型
        for info_type in selected_types:
            if info_type == 'protrek':
                # 处理protrek类型
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
                
                # 获取该类型的所有IPR ID
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

# 全局变量来缓存InterProDescriptionManager实例和lmdb连接
_interpro_manager = None
_lmdb_db = None
_lmdb_path = None

def get_interpro_manager(interpro_data_path, interproscan_info_path, protrek_text_dir):
    """获取或创建InterProDescriptionManager实例"""
    global _interpro_manager
    if _interpro_manager is None:
        _interpro_manager = InterProDescriptionManager(interpro_data_path, interproscan_info_path, protrek_text_dir)
    return _interpro_manager

def get_lmdb_connection(lmdb_path):
    """获取或创建lmdb连接"""
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
    获取prompt模板，支持可选的信息类型
    
    Args:
        selected_info_types: 需要包含的信息类型列表，如['motif', 'go', 'superfamily', 'panther', 'protrek']
    """
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']  # 默认包含motif和go信息
    if lmdb_path is None:
        PROMPT_TEMPLATE = ENZYME_PROMPT_SIMPLE + "\n"
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
    从lmdb中获取指定蛋白质的所有QA对
    
    Args:
        protein_id: 蛋白质ID
        lmdb_path: lmdb数据库路径
    
    Returns:
        list: QA对列表，每个元素包含question和ground_truth
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
            # 遍历数字索引的数据，查找匹配的protein_id
            cursor = txn.cursor()
            for key, value in cursor:
                try:
                    # 尝试将key解码为数字（数字索引的数据）
                    key_str = key.decode('utf-8')
                    if key_str.isdigit():
                        # 这是数字索引的数据，包含protein_id, question, ground_truth
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
                    # 如果解析失败，跳过这个条目
                    continue
    except Exception as e:
        print(f"Error reading lmdb for protein {protein_id}: {e}")
    
    return qa_pairs

def generate_prompt(protein_id, protein2gopath, protrek_text_dir, protein2pfam_path, pfam_descriptions_path, go_info_path, 
                   interpro_data_path=None, interproscan_info_path=None, selected_info_types=None, lmdb_path=None, interpro_manager=None, question=None, protrek_topk=3):
    """
    生成蛋白质prompt
    
    Args:
        selected_info_types: 需要包含的信息类型列表，如['motif', 'go', 'superfamily', 'panther', 'protrek']
        interpro_data_path: interpro_data.json文件路径
        interproscan_info_path: interproscan_info.json文件路径
        interpro_manager: InterProDescriptionManager实例，如果提供则优先使用
        question: 问题文本，用于QA任务
        protrek_topk: ProTrek返回的top结果数量
    """
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']
    
    # 获取分析结果
    analysis = analyze_protein_go(protein_id, protein2gopath, go_info_path)
    motif_pfam = get_motif_pfam(protein_id, protein2pfam_path, pfam_descriptions_path)
    
    # 获取InterPro和ProTrek描述信息
    interpro_descriptions = {}
    other_types = [t for t in selected_info_types if t not in ['motif', 'go']]
    if other_types:
        if interpro_manager:
            # 使用提供的manager实例
            interpro_descriptions = interpro_manager.get_description(protein_id, other_types, protrek_topk)
        elif interpro_data_path and interproscan_info_path:
            # 使用全局缓存的manager
            manager = get_interpro_manager(interpro_data_path, interproscan_info_path, protrek_text_dir)
            interpro_descriptions = manager.get_description(protein_id, other_types, protrek_topk)

    # 准备模板数据
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
    """并行生成和保存protein prompts"""
    import json
    try:
        from utils.mpr import MultipleProcessRunnerSimplifier
    except ImportError:
        from mpr import MultipleProcessRunnerSimplifier
    
    if selected_info_types is None:
        selected_info_types = ['motif', 'go']
    
    # 在并行处理开始前创建InterProDescriptionManager实例
    interpro_manager = None
    other_types = [t for t in selected_info_types if t not in ['motif', 'go']]
    if other_types and interpro_data_path and interproscan_info_path:
        interpro_manager = InterProDescriptionManager(interpro_data_path, interproscan_info_path, protrek_text_dir)
    
    # 用于跟踪全局index的共享变量
    # if lmdb_path:
    #     import multiprocessing
    #     global_index = multiprocessing.Value('i', 0)  # 共享整数，初始值为0
    #     index_lock = multiprocessing.Lock()  # 用于同步访问
    # else:
    #     global_index = None
    #     index_lock = None
    import multiprocessing
    global_index = multiprocessing.Value('i', 0)  # 共享整数，初始值为0
    index_lock = multiprocessing.Lock()  # 用于同步访问
    results = {}
    
    def process_protein(process_id, idx, protein_id, writer):
        protein_id = protein_id.strip()
        
        # 为每个进程初始化lmdb连接
        if lmdb_path:
            get_lmdb_connection(lmdb_path)
        
        if lmdb_path:
            # 如果有lmdb_path，处理QA数据
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
                    # 获取并递增全局index
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
            # 如果没有lmdb_path，按原来的方式处理
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
    
    # 使用MultipleProcessRunnerSimplifier进行并行处理
    runner = MultipleProcessRunnerSimplifier(
        data=protein_ids,
        do=process_protein,
        save_path=output_path + '.tmp',
        n_process=n_process,
        split_strategy="static"
    )
    
    runner.run()
    
    # 清理全局lmdb连接
    global _lmdb_db
    if _lmdb_db is not None:
        _lmdb_db.close()
        _lmdb_db = None
    
    # if not lmdb_path:
    #     # 如果没有lmdb_path，合并所有结果到一个字典（兼容旧格式）
    #     final_results = {}
    #     with open(output_path + '.tmp', 'r') as f:
    #         for line in f:
    #             if line.strip():  # 忽略空行
    #                 final_results.update(json.loads(line))
        
    #     # 保存最终结果为正确的JSON格式
    #     with open(output_path, 'w') as f:
    #         json.dump(final_results, f, indent=2)
    # else:
    #     # 如果有lmdb_path，直接保存为jsonl格式
    #     import shutil
    #     shutil.move(output_path + '.tmp', output_path)
    import shutil
    shutil.move(output_path + '.tmp', output_path)
    
    # 删除临时文件（如果还存在的话）
    if os.path.exists(output_path + '.tmp'):
        os.remove(output_path + '.tmp')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='生成蛋白质功能分析的prompt')
    parser.add_argument('--protein_path', type=str, required=True, 
                       help='蛋白质ID列表文件路径（每行一个ID）')
    parser.add_argument('--protein2pfam_path', type=str, required=True, 
                       help='蛋白质-Pfam映射文件路径（InterProScan结果JSON文件）')
    parser.add_argument('--protrek_text_dir', type=str, default=None, 
                       help='ProTrek文本描述目录路径（如果使用protrek信息类型）')
    parser.add_argument('--pfam_descriptions_path', type=str, required=True, 
                       help='Pfam描述文件路径（包含Pfam ID和描述的JSON文件）')
    parser.add_argument('--protein2gopath', type=str, required=True, 
                       help='蛋白质-GO映射文件路径（JSON Lines格式）')
    parser.add_argument('--go_info_path', type=str, required=True, 
                       help='GO本体文件路径（go.json）')
    parser.add_argument('--interpro_data_path', type=str, default=None, 
                       help='InterPro数据文件路径（如果使用interpro相关信息类型）')
    parser.add_argument('--interproscan_info_path', type=str, default=None, 
                       help='InterProScan结果文件路径（如果使用interpro相关信息类型）')
    parser.add_argument('--lmdb_path', type=str, default=None, 
                       help='LMDB数据库路径（如果生成QA格式的prompt）')
    parser.add_argument('--output_path', type=str, required=True, 
                       help='输出文件路径（JSON或JSONL格式）')
    parser.add_argument('--selected_info_types', type=str, nargs='+', default=['motif', 'go'], 
                       help='选择要包含的信息类型（默认：motif go）。可选：motif, go, protrek')
    parser.add_argument('--n_process', type=int, default=8, 
                       help='并行处理的进程数（默认：8）')
    parser.add_argument('--protrek_topk', type=int, default=3, 
                       help='ProTrek返回的top结果数量（默认：3）')
    args = parser.parse_args()
    #更新output_path，需要包含selected_info_types
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

    # 测试示例
    # protein_id = 'A8CF74'
    # prompt = generate_prompt(protein_id, 'data/processed_data/go_integration_final_topk2.json', 
    #                         'data/processed_data/interproscan_info.json', 'data/raw_data/all_pfam_descriptions.json', 
    #                         'data/raw_data/go.json', 'data/raw_data/interpro_data.json', 
    #                         'data/processed_data/interproscan_info.json', 
    #                         ['motif', 'go', 'superfamily', 'panther'])
    # print(prompt)

        