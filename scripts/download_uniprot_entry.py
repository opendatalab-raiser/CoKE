import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import requests
from utils.mpr import MultipleProcessRunnerSimplifier
import json

class UniProtDownloader(MultipleProcessRunnerSimplifier):
    def __init__(self, protein_ids, output_dir, **kwargs):
        super().__init__(
            data=protein_ids,
            do=self.download_single,
            save_path=None,  # 每个文件单独保存，不需要聚合
            return_results=False,
            **kwargs
        )
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def download_single(self, process_id, idx, protein_id, _):
        url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            output_path = os.path.join(self.output_dir, f"{protein_id}.json")
            with open(output_path, "w") as f:
                json.dump(response.json(), f, indent=2)
                
        except Exception as e:
            print(f"Error downloading {protein_id}: {str(e)}")

if __name__ == "__main__":
    # 从文件读取蛋白质ID列表
    import argparse
    parser = argparse.ArgumentParser(description="下载蛋白质Entry")
    parser.add_argument("--input_path", type=str, default="data/raw_data/protein_ids_all20250623.txt", help="蛋白质ID列表文件路径")
    parser.add_argument("--output_dir", type=str, default="data/raw_data/SP_Entry", help="输出目录")
    parser.add_argument("--n_process", type=int, default=64, help="进程数")
    parser.add_argument("--max_round", type=int, default=5, help="最大轮数")
    args = parser.parse_args()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    with open(args.input_path) as f:
        protein_ids = [line.strip() for line in f if line.strip()]
    
    #检查output_dir存在的protein_id，和input_path重复的就不要下载了,这里启用多轮下载
    for i in range(args.max_round):
        existing_ids = [f.split(".")[0] for f in os.listdir(args.output_dir) if f.endswith(".json")]
        protein_ids = set(protein_ids) - set(existing_ids)
        protein_ids = list(protein_ids)
        print(f"第{i+1}轮下载，需要下载的蛋白质ID数量: {len(protein_ids)}")
        if len(protein_ids) == 0:
            print("所有蛋白质ID都已下载完成")
            break
        downloader = UniProtDownloader(
            protein_ids=protein_ids,
            output_dir=args.output_dir,
            n_process=args.n_process,  # 根据CPU核心数调整
            split_strategy="static",
            log_step=10
        )
        downloader.run()
        
    # print(f"需要下载的蛋白质ID数量: {len(protein_ids)}")

    # # 配置下载参数
    # downloader = UniProtDownloader(
    #     protein_ids=protein_ids,
    #     output_dir=args.output_dir,
    #     n_process=args.n_process,  # 根据CPU核心数调整
    #     split_strategy="static",
    #     log_step=10
    # )
    
    # downloader.run()