import json
from collections import Counter
import os

_pfam_dict = None
_pfam_descriptions = None

def _load_pfam_data(protein2pfam_path):
    global _pfam_dict
    if _pfam_dict is None:
        with open(protein2pfam_path, 'r') as file:
            _pfam_dict = json.load(file)

def _load_pfam_descriptions(pfam_descriptions_path):
    global _pfam_descriptions
    if _pfam_descriptions is None:
        with open(pfam_descriptions_path, 'r') as file:
            _pfam_descriptions = json.load(file)

def get_motif_pfam(protein_id, protein2pfam_path, pfam_descriptions_path):
    """
    获取指定蛋白质的pfam信息及其定义
    
    参数:
    protein_id: str - 蛋白质ID
    protein2pfam_path: str - interproscan_info.json文件路径
    pfam_descriptions_path: str - pfam描述文件路径
    
    返回:
    dict - pfam_id到定义的映射字典，例如{"PF04820": "definition content"}
    """
    _load_pfam_data(protein2pfam_path)
    _load_pfam_descriptions(pfam_descriptions_path)
    
    if protein_id not in _pfam_dict:
        return {}
    
    protein_info = _pfam_dict[protein_id]
    _pfam_dicts = protein_info.get('interproscan_results', {}).get('pfam_id', [])
    pfam_ids = []
    for pfam_dict in _pfam_dicts:
        for key,value in pfam_dict.items():
            pfam_ids.append(key)
    
    result = {}
    for pfam_id in pfam_ids:
        if pfam_id in _pfam_descriptions:
            result[pfam_id] = _pfam_descriptions[pfam_id]['description']
    
    return result

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='获取蛋白质的Pfam motif信息')
    parser.add_argument("--protein_id", type=str, required=True, help='蛋白质ID，例如：A8CF74')
    parser.add_argument("--protein2pfam_path", type=str, required=True, help='蛋白质-Pfam映射文件路径（InterProScan结果JSON文件）')
    parser.add_argument("--pfam_descriptions_path", type=str, required=True, help='Pfam描述文件路径（包含Pfam ID和描述的JSON文件）')
    args = parser.parse_args()
    result = get_motif_pfam(args.protein_id, args.protein2pfam_path, args.pfam_descriptions_path)
    print(result)