import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import requests
from utils.mpr import MultipleProcessRunnerSimplifier

class PDBDownloader(MultipleProcessRunnerSimplifier):
    def __init__(self, protein_ids, output_dir, failed_ids_path, **kwargs):
        # 创建一个临时文件路径用于收集失败的ID
        temp_failed_path = f"{failed_ids_path}.temp"
        
        super().__init__(
            data=protein_ids,
            do=self.download_single,
            save_path=temp_failed_path,  # 设置临时保存路径
            return_results=True,
            **kwargs
        )
        self.output_dir = output_dir
        self.failed_ids_path = failed_ids_path
        self.temp_failed_path = temp_failed_path
        os.makedirs(self.output_dir, exist_ok=True)
    
    def download_single(self, process_id, idx, protein_id, writer):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_id}-F1-model_v6.pdb"
        output_path = os.path.join(self.output_dir, f"AF-{protein_id}-F1-model_v6.pdb")
        
        if os.path.exists(output_path):
            return None
            
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                return None
            else:
                print(f"下载蛋白质 {protein_id} 失败，状态码: {response.status_code}")
                if writer:
                    writer.write(f"{protein_id}\n")
                return protein_id
        except Exception as e:
            print(f"下载蛋白质 {protein_id} 时出错: {str(e)}")
            if writer:
                writer.write(f"{protein_id}\n")
            return protein_id
    
    def save_failed_ids(self, results):
        # 从临时文件中读取失败的ID
        failed_ids = []
        if os.path.exists(self.temp_failed_path):
            with open(self.temp_failed_path, 'r') as f:
                failed_ids = [line.strip() for line in f if line.strip()]
            
            # 删除临时文件
            os.remove(self.temp_failed_path)
        
        # 也从返回的结果中收集失败的ID
        result_failed_ids = [result for result in results if result is not None]
        
        # 合并两个来源的失败ID并去重
        all_failed_ids = list(set(failed_ids + result_failed_ids))
        
        if all_failed_ids:
            print(f"共有 {len(all_failed_ids)} 个蛋白质下载失败")
            # 保存失败的ID到文件
            with open(self.failed_ids_path, 'w') as f:
                for protein_id in all_failed_ids:
                    f.write(f"{protein_id}\n")
            print(f"已将下载失败的蛋白质ID保存到 {self.failed_ids_path}")

def parse_args():
    """解析命令行参数"""
    import argparse
    parser = argparse.ArgumentParser(description="下载PDB文件")
    parser.add_argument("--protein_ids_file", type=str, required=True, help="蛋白质ID列表文件路径")
    parser.add_argument("--output_dir", type=str, default="/zhuangkai/projects/TTS4Protein/data/raw_data/pdb", help="输出目录")
    parser.add_argument("--failed_ids_path", type=str, default="/zhuangkai/projects/TTS4Protein/data/processed_data/protein_wo_pdb.txt", help="失败ID保存路径")
    parser.add_argument("--n_process", type=int, default=64, help="进程数")
    return parser.parse_args()

def main():
    # 读取蛋白质ID列表
    args = parse_args()
    protein_ids_file = args.protein_ids_file
    with open(protein_ids_file, 'r') as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    
    print(f"共有 {len(protein_ids)} 个蛋白质需要下载")
    
    # 配置下载参数
    downloader = PDBDownloader(
        protein_ids=protein_ids,
        output_dir=args.output_dir,
        failed_ids_path=args.failed_ids_path,
        n_process=args.n_process,
        verbose=True,
        total_only=False,
        log_step=1,
        split_strategy="queue"
    )
    
    results = downloader.run()
    downloader.save_failed_ids(results)

if __name__ == "__main__":
    main()
