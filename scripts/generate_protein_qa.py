import os
import json
import lmdb
import argparse
from typing import List, Dict, Optional
from tqdm import tqdm

def get_protein_info(entry_dir: str, protein_id: str) -> Optional[Dict[str, str]]:
    """
    从JSON文件中获取蛋白质的function/pathway/subcellular location信息
    
    Args:
        entry_dir: entry目录路径
        protein_id: 蛋白质ID
        
    Returns:
        Dict: 包含三种类型信息的字典，如果出错返回None
    """
    try:
        file_path = os.path.join(entry_dir, f'{protein_id}.json')
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        result = {
            'function': "",
            'pathway': "",
            'subcellular_location': ""
        }
        
        if 'comments' in data and len(data['comments']) > 0:
            for comment in data['comments']:
                # 处理FUNCTION
                if comment["commentType"] == "FUNCTION":
                    for text in comment["texts"]:
                        result['function'] += text["value"] + "\n"
                
                # 处理PATHWAY
                elif comment["commentType"] == "PATHWAY":
                    for text in comment["texts"]:
                        result['pathway'] += text["value"] + "\n"
                
                # 处理SUBCELLULAR LOCATION
                elif comment["commentType"] == "SUBCELLULAR LOCATION":
                    for location in comment.get("subcellularLocations", []):
                        if "location" in location:
                            result['subcellular_location'] += location["location"]["value"] + "\n"
        
        # 去除多余的空格和换行
        for key in result:
            result[key] = result[key].strip()
            if not result[key]:
                result[key] = None
        
        return result
        
    except Exception as e:
        print(f"Error processing {protein_id}: {str(e)}")
        return None

def generate_qa_pairs(protein_id: str, protein_info: Dict[str, str]) -> List[Dict]:
    """
    根据蛋白质信息生成QA对
    
    Args:
        protein_id: 蛋白质ID
        protein_info: 蛋白质信息字典
        
    Returns:
        List: QA对列表
    """
    qa_pairs = []
    
    # 添加function QA对
    if protein_info['function']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the function of this protein?',
            'answer': protein_info['function'],
            'question_type': 'function'
        })
    
    # 添加pathway QA对
    if protein_info['pathway']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the pathway of this protein?',
            'answer': protein_info['pathway'],
            'question_type': 'pathway'
        })
    
    # 添加subcellular location QA对
    if protein_info['subcellular_location']:
        qa_pairs.append({
            'protein_id': protein_id,
            'question': 'What is the subcellular location of this protein?',
            'answer': protein_info['subcellular_location'],
            'question_type': 'subcellular_location'
        })
    
    return qa_pairs

def save_to_lmdb(qa_pairs: List[Dict], lmdb_path: str):
    """
    将QA对保存到LMDB数据库
    
    Args:
        qa_pairs: QA对列表
        lmdb_path: LMDB数据库路径
    """
    env = lmdb.open(lmdb_path, map_size=1099511627776)  # 1TB最大大小
    with env.begin(write=True) as txn:
        for idx, qa in enumerate(qa_pairs):
            key = str(idx).encode('utf-8')
            value = json.dumps([
                qa['protein_id'],
                qa['question'],
                qa['answer'],
                qa['question_type']
            ]).encode('utf-8')
            txn.put(key, value)

def save_to_json(qa_pairs: List[Dict], json_path: str):
    """
    将QA对保存到JSON文件
    
    Args:
        qa_pairs: QA对列表
        json_path: JSON文件路径
    """
    with open(json_path, 'w') as f:
        json.dump(qa_pairs, f, indent=2)

def read_protein_ids_from_file(file_path: str) -> List[str]:
    """
    从文本文件中读取蛋白质ID列表
    
    Args:
        file_path: 蛋白质ID文件路径
        
    Returns:
        List: 蛋白质ID列表
    """
    with open(file_path, 'r') as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    return protein_ids

def process_proteins(entry_dir: str, protein_ids: List[str], lmdb_path: str, json_path: str):
    """
    处理指定的蛋白质ID列表，生成QA对并保存
    
    Args:
        entry_dir: entry目录路径
        protein_ids: 要处理的蛋白质ID列表
        lmdb_path: LMDB数据库保存路径
        json_path: JSON文件保存路径
    """
    all_qa_pairs = []
    processed_count = 0
    
    for protein_id in tqdm(protein_ids):
        protein_info = get_protein_info(entry_dir, protein_id)
        if protein_info is None:
            continue
            
        qa_pairs = generate_qa_pairs(protein_id, protein_info)
        all_qa_pairs.extend(qa_pairs)
        processed_count += 1
    
    # 保存到LMDB
    save_to_lmdb(all_qa_pairs, lmdb_path)
    
    # 保存到JSON
    save_to_json(all_qa_pairs, json_path)
    
    print(f"Successfully processed {len(all_qa_pairs)} QA pairs from {processed_count}/{len(protein_ids)} proteins")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='从UniProt Entry文件生成蛋白质功能QA对')
    parser.add_argument('--entry_dir', type=str, required=True, 
                       help='UniProt Entry目录路径（包含蛋白质JSON文件）')
    parser.add_argument('--protein_id_files', type=str, default=None, 
                       help='蛋白质ID列表文件路径（每行一个ID）。如果不指定，则处理entry_dir下所有JSON文件')
    parser.add_argument('--lmdb_path', type=str, required=True, 
                       help='LMDB数据库输出路径')
    parser.add_argument('--json_path', type=str, required=True, 
                       help='JSON文件输出路径')
    args = parser.parse_args()
    save_dir = os.path.dirname(args.lmdb_path)
    print(save_dir)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir, exist_ok=True)
    
    if args.protein_id_files:
        # 如果指定了protein_id_files，只处理这些ID
        protein_ids = read_protein_ids_from_file(args.protein_id_files)
    else:
        # 否则处理entry目录下所有.json文件
        protein_ids = [os.path.splitext(f)[0] for f in os.listdir(args.entry_dir) if f.endswith('.json')]
    
    process_proteins(args.entry_dir, protein_ids, args.lmdb_path, args.json_path)