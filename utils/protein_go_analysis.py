import json
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from collections import defaultdict

# 全局变量声明
_go_data = None
_protein_go_dict = None

def _load_go_data(go_info_path):
    """懒加载GO数据"""
    global _go_data
    if _go_data is None:
        try:
            with open(go_info_path, 'r') as f:
                _go_data = json.load(f)
        except Exception as e:
            print(f"加载GO数据文件时发生错误: {str(e)}")
            _go_data = None

def _load_protein_go_dict(protein2gopath):
    """懒加载蛋白质-GO映射数据"""
    global _protein_go_dict
    if _protein_go_dict is None:
        try:
            _protein_go_dict = {}
            with open(protein2gopath, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    _protein_go_dict[data['protein_id']] = data['GO_id']
        except Exception as e:
            print(f"加载蛋白质-GO映射数据时发生错误: {str(e)}")
            _protein_go_dict = None

def get_go_definition(go_id, go_info_path):
    """获取GO term的定义"""
    _load_go_data(go_info_path)
    if _go_data is None:
        return None

    if not go_id.startswith('GO_'):
        go_id = f"GO_{go_id}"
    full_id = f"http://purl.obolibrary.org/obo/{go_id}"
    
    for node in _go_data['graphs'][0]['nodes']:
        if node['id'] == full_id:
            if 'meta' in node and 'definition' in node['meta']:
                return node['meta']['definition']['val']
    return None

def analyze_protein_go(protein_id, protein2gopath, go_info_path):
    """
    分析蛋白质的GO注释信息，包括GO ID和定义
    
    参数：
    protein_id: str - 蛋白质ID
    protein2gopath: str - 蛋白质-GO映射文件路径
    
    返回：
    dict - 包含GO信息的字典
    """
    _load_protein_go_dict(protein2gopath)
    if _protein_go_dict is None:
        return {
            "status": "error",
            "message": "GO数据加载失败"
        }

    if protein_id not in _protein_go_dict:
        return {
            "status": "error",
            "message": f"未找到蛋白质 {protein_id} 的GO注释"
        }
    
    go_ids = _protein_go_dict[protein_id]
    go_info = []
    all_definitions = {}
    
    for go_id in go_ids:
        # 获取GO定义
        definition = get_go_definition(go_id, go_info_path)
        if definition:
            all_definitions[go_id] = definition
        
        go_info.append({
            "go_id": go_id
        })
    
    return {
        "status": "success" if len(go_info) > 0 else "error",
        "protein_id": protein_id,
        "go_annotations": go_info,
        "all_related_definitions": all_definitions
    }

# 使用示例
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='分析蛋白质的GO注释信息')
    parser.add_argument('--protein_id', type=str, required=True, help='蛋白质ID，例如：A8CF74')
    parser.add_argument('--protein2gopath', type=str, required=True, help='蛋白质-GO映射文件路径（JSON Lines格式）')
    parser.add_argument('--go_info_path', type=str, required=True, help='GO本体文件路径（go.json）')
    args = parser.parse_args()
    
    result = analyze_protein_go(args.protein_id, args.protein2gopath, args.go_info_path)
    
    if result["status"] == "success":
        print(f"\nProtein {result['protein_id']} GO annotations:")
        
        for anno in result["go_annotations"]:
            print(f"\nGO ID: {anno['go_id']}")
        
        print("\nAll related GO ID definitions:")
        for go_id, definition in result["all_related_definitions"].items():
            print(f"\nGO:{go_id}")
            print(f"Definition: {definition}")
    else:
        print(result["message"]) 