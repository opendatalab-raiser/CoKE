from gradio_client import Client
import os
import json
from tqdm import tqdm
import time
import random

def get_protrek_text(sequence, topk=3):
    """
    获取ProTrek文本描述
    
    Args:
        sequence: 蛋白质序列
        topk: 返回的top结果数量
    
    Returns:
        list: 包含文本描述的列表，如果失败返回空列表
    """
    
    try:
        client = Client("http://search-protrek.com/")
        result = client.predict(
            input=sequence,
            nprobe=1000,
            topk=topk,
            input_type="sequence",
            query_type="text",
            subsection_type="Function",
            db="Swiss-Prot",
            api_name="/search"
        )

        full_res = result[-1]["value"]
        data = full_res["data"]

        descriptions = []
        for des, score in data:
            descriptions.append((des, score))
        return descriptions
        
    except Exception as e:
        return []

def get_missing_pids(pids, output_dir):
    missing_pids = []
    for pid in pids:
        if not os.path.exists(f'{output_dir}/{pid}.json'):
            missing_pids.append(pid)
    return missing_pids

def save_protrek_text_to_file():
    import argparse
    from Bio import SeqIO
    
    parser = argparse.ArgumentParser(description='使用ProTrek生成蛋白质功能文本描述')
    # 支持两种模式：1) 从FASTA文件读取 2) 从pids和json文件读取
    parser.add_argument('--input_fasta', type=str, default=None, help='输入FASTA文件路径（优先使用此参数）')
    parser.add_argument('--pids_path', type=str, default=None, help='蛋白质ID列表文件路径（每行一个ID）')
    parser.add_argument('--topk', type=int, default=3, help='返回的top结果数量（默认：3）')
    parser.add_argument('--protein_seq_info', type=str, default=None, help='蛋白质序列信息JSON文件路径（与pids_path配合使用）')
    parser.add_argument('--output_dir', type=str, default='protrek_text', help='输出目录路径（默认：protrek_text）')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 优先使用FASTA文件
    if args.input_fasta and os.path.exists(args.input_fasta):
        print(f"从FASTA文件读取序列: {args.input_fasta}")
        pids = []
        protein_info = {}
        with open(args.input_fasta, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                pids.append(record.id)
                print(record.id)
                print(str(record.seq))
                protein_info[record.id] = str(record.seq)
    else:
        # 使用传统模式：从pids和json文件读取
        print(f"从PIDs文件和JSON文件读取序列")
        with open(args.pids_path, 'r') as f:
            pids = f.read().splitlines()
        with open(args.protein_seq_info, 'r') as f:
            protein_info = json.load(f)
    
    # 获取当前缺失的PIDs
    missing_pids = get_missing_pids(pids, args.output_dir)
    print(f'当前缺失的PIDs: {len(missing_pids)}')
    
    if len(missing_pids) == 0:
        print("所有PIDs都已处理完成!")
        return True  # 返回True表示所有任务完成
    
    # 处理缺失的PIDs，一旦失败就立即退出
    for pid in tqdm(missing_pids, desc="处理PIDs"):
        try:
            sequence = protein_info[pid]
            descriptions = get_protrek_text(sequence, args.topk)
            
            if descriptions:  # 成功获取到描述
                with open(f'{args.output_dir}/{pid}.json', 'w') as f:
                    json.dump(descriptions, f)
                print(f"✓ PID {pid} 处理成功")
            else:
                print(f"✗ PID {pid} 获取描述失败，退出程序等待重试")
                exit(1)  # 立即退出，返回非零状态码
                
        except Exception as e:
            print(f"✗ 处理PID {pid} 时发生未预期的错误: {str(e)}，退出程序等待重试")
            exit(1)  # 立即退出，返回非零状态码
    
    print("本次运行所有PIDs都处理成功!")

if __name__ == "__main__":
    success = save_protrek_text_to_file()
