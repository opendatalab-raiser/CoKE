import json
import os
import sys
import argparse
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
from tqdm import tqdm

# 添加路径
root_path = os.path.dirname((os.path.abspath(__file__)))
sys.path.append(root_path)
sys.path.append(os.path.join(root_path, "Models/ProTrek"))

from utils.protein_go_analysis import get_go_definition

class GOIntegrationPipeline:
    def __init__(self, 
                 identity_threshold: int = 80,
                 coverage_threshold: int = 80, 
                 evalue_threshold: float = 1e-50,
                 topk: int = 2,
                 protrek_threshold: Optional[float] = None,
                 use_protrek: bool = False):
        """
        GO信息整合管道
        
        Args:
            identity_threshold: BLAST identity阈值 (0-100)
            coverage_threshold: BLAST coverage阈值 (0-100) 
            evalue_threshold: BLAST E-value阈值
            protrek_threshold: ProTrek分数阈值
            use_protrek: 是否使用第二层ProTrek筛选
        """
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.evalue_threshold = evalue_threshold
        self.protrek_threshold = protrek_threshold
        self.use_protrek = use_protrek
        self.topk = topk
        self.current_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))#当前文件目录的上两层
        self.go_info_path = os.path.join(self.current_path, 'data/raw_data/go.json')
        self.protein2go_path = os.path.join(self.current_path, 'data/processed_data/gt_protein2go_sp20250623.json')
        # self.pid2samepid_path = os.path.join(self.current_path, 'data/processed_data/swissprot_pid2samepid.json')
        self.pid2seq_path = os.path.join(self.current_path, 'data/processed_data/swissprot_pid2seq.json')
        
        # 加载蛋白质-GO映射数据
        self._load_protein_go_dict()
        self._load_pid2seq()
        # self._load_pid2samepid_dict()
        
        # 如果使用protrek，初始化模型
        if self.use_protrek:
            self._init_protrek_model()
    
    def _init_protrek_model(self):
        """初始化ProTrek模型"""
        from model.ProTrek.protrek_trimodal_model import ProTrekTrimodalModel

        config = {
            "protein_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/esm2_t33_650M_UR50D",
            "text_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext",
            "structure_config": "Models/ProTrek/weights/ProTrek_650M_UniRef50/foldseek_t30_150M",
            "load_protein_pretrained": False,
            "load_text_pretrained": False,
            "from_checkpoint": "Models/ProTrek/weights/ProTrek_650M_UniRef50/ProTrek_650M_UniRef50.pt"
        }
        
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.protrek_model = ProTrekTrimodalModel(**config).to(self.device).eval()
        print(f"ProTrek模型已加载到设备: {self.device}")
    
    def _load_protein_go_dict(self):
        """加载蛋白质-GO映射数据"""
        self.protein_go_dict = {}
        try:
            with open(self.protein2go_path, 'r') as f:
                for line in f:
                    data = json.loads(line)
                    self.protein_go_dict[data['protein_id']] = data['GO_id']
            print(f"成功加载蛋白质-GO映射数据，共{len(self.protein_go_dict)}条记录")
        except Exception as e:
            print(f"加载蛋白质-GO映射数据时发生错误: {str(e)}")
            self.protein_go_dict = {}
    
    # def _load_pid2samepid_dict(self):
    #     """加载pid->相同序列的pid的mapping"""
    #     try:
    #         with open(self.pid2samepid_path, 'r') as f:
    #             self.pid2samepid_dict=json.load(f)
    #         print(f"成功加载pid->same pid映射数据，共{len(self.pid2samepid_dict)}条记录")
    #     except Exception as e:
    #         print(f"加载pid->same pid映射数据时发生错误: {str(e)}")
    #         self.pid2samepid_dict = {}
            
    def _load_pid2seq(self):
        """加载pid->序列"""
        try:
            with open(self.pid2seq_path, 'r') as f:
                self.pid2seq=json.load(f)
            print(f"成功加载pid->seq映射数据，共{len(self.pid2seq)}条记录")
        except Exception as e:
            print(f"加载pid->same seq映射数据时发生错误: {str(e)}")
            self.pid2seq = {}
        
    
    def _get_go_from_uniprot_id(self, uniprot_id: str) -> List[str]:
        """
        从Uniprot ID获取GO ID
        
        Args:
            uniprot_id: Uniprot ID
        
        Returns:
            使用类内部加载的字典
        """
        # 使用类内部加载的字典
        return [go_id.split("_")[-1] if "_" in go_id else go_id 
                for go_id in self.protein_go_dict.get(uniprot_id, [])]
    
    def extract_blast_go_ids(self, blast_results: List[Dict],sequence: str) -> List[str]:
        """
        从BLAST结果中提取符合条件的GO ID
        
        Args:
            blast_results: BLAST结果列表
            sequence: 当前蛋白质序列（避免自身匹配）
            
        Returns:
            符合条件的GO ID列表
        """
        go_ids = []
        
        if self.topk > 0:
            # 使用topk策略
            for result in blast_results[:self.topk]:
                hit_id = result.get('ID', '')
                if self.pid2seq[hit_id] ==sequence:
                    continue
                go_ids.extend(self._get_go_from_uniprot_id(hit_id))
        else:
            # 使用阈值策略
            for result in blast_results:
                identity = float(result.get('Identity%', 0))
                coverage = float(result.get('Coverage%', 0))
                evalue = float(result.get('E-value', 1.0))
                
                # 检查是否符合阈值条件
                if (identity >= self.identity_threshold and 
                    coverage >= self.coverage_threshold and 
                    evalue <= self.evalue_threshold):
                    
                    # 获取该hit的protein_id
                    hit_id = result.get('ID', '')
                    if self.pid2seq[hit_id] ==sequence:
                        continue
                    go_ids.extend(self._get_go_from_uniprot_id(hit_id))
        
        return go_ids
    
    def first_level_filtering(self, interproscan_info: Dict, blast_info: Dict) -> Dict[str, List[str]]:
        """
        第一层筛选：合并interproscan和符合条件的blast GO信息
        
        Args:
            interproscan_info: InterProScan结果
            blast_info: BLAST结果
            
        Returns:
            蛋白质ID到GO ID列表的映射
        """
        protein_go_dict = {}
        
        for protein_id in interproscan_info.keys():
            go_ids = set()
            
            # 添加interproscan的GO信息
            interproscan_gos = interproscan_info[protein_id].get('interproscan_results', {}).get('go_id', [])
            interproscan_gos = [go_id.split(":")[-1] if ":" in go_id else go_id for go_id in interproscan_gos]
            if interproscan_gos:
                go_ids.update(interproscan_gos)
            
            #添加符合条件的blast GO信息
            if protein_id in blast_info:
                sequence=blast_info[protein_id]['sequence']
                blast_results = blast_info[protein_id].get('blast_results', [])
                blast_gos = self.extract_blast_go_ids(blast_results,sequence)
                go_ids.update(blast_gos)
            
            protein_go_dict[protein_id] = list(go_ids)
        
        return protein_go_dict
    
    def calculate_protrek_scores(self, protein_sequences: Dict[str, str], 
                                protein_go_dict: Dict[str, List[str]]) -> Dict[str, Dict]:
        """
        计算ProTrek分数
        
        Args:
            protein_sequences: 蛋白质序列字典
            protein_go_dict: 蛋白质GO映射
            
        Returns:
            包含GO分数的字典
        """
        results = {}
        
        for protein_id, go_ids in tqdm(protein_go_dict.items(), desc="计算ProTrek分数"):
            if protein_id not in protein_sequences:
                continue
                
            protein_seq = protein_sequences[protein_id]
            go_scores = {}
            
            # 获取GO定义
            go_definitions = {}
            for go_id in go_ids:
                definition = get_go_definition(go_id,self.go_info_path)
                if definition:
                    go_definitions[go_id] = definition
            
            if not go_definitions:
                continue
            
            try:
                with torch.no_grad():
                    # 计算蛋白质序列嵌入
                    seq_emb = self.protrek_model.get_protein_repr([protein_seq])
                    
                    # 计算文本嵌入和相似度分数
                    definitions = list(go_definitions.values())
                    text_embs = self.protrek_model.get_text_repr(definitions)
                    
                    # 计算相似度分数
                    scores = (seq_emb @ text_embs.T) / self.protrek_model.temperature
                    scores = scores.cpu().numpy().flatten()
                    
                    # 映射回GO ID
                    for i, go_id in enumerate(go_definitions.keys()):
                        go_scores[go_id] = float(scores[i])
                    
            except Exception as e:
                print(f"计算 {protein_id} 的ProTrek分数时出错: {str(e)}")
                continue
            
            results[protein_id] = {
                "protein_id": protein_id,
                "GO_id": go_ids,
                "Clip_score": go_scores
            }
        
        return results
    
    def second_level_filtering(self, protrek_results: Dict[str, Dict]) -> Dict[str, List[str]]:
        """
        第二层筛选：根据ProTrek阈值筛选GO
        
        Args:
            protrek_results: ProTrek计算结果
            
        Returns:
            筛选后的蛋白质GO映射
        """
        filtered_results = {}
        
        for protein_id, data in protrek_results.items():
            clip_scores = data.get('Clip_score', {})
            filtered_gos = []
            
            for go_id, score in clip_scores.items():
                if score >= self.protrek_threshold:
                    filtered_gos.append(go_id)
            
            if filtered_gos:
                filtered_results[protein_id] = filtered_gos
        
        return filtered_results
    
    def generate_filename(self, base_name: str, is_intermediate: bool = False) -> str:
        """生成包含参数信息的文件名"""
        if self.topk > 0:
            # 如果使用topk，则只包含topk信息
            params = f"topk{self.topk}"
        else:
            # 否则使用原有的参数组合
            params = f"identity{self.identity_threshold}_coverage{self.coverage_threshold}_evalue{self.evalue_threshold:.0e}"
        
        if self.use_protrek and self.protrek_threshold is not None:
            params += f"_protrek{self.protrek_threshold}"
        
        if is_intermediate:
            return f"{base_name}_intermediate_{params}.json"
        else:
            return f"{base_name}_final_{params}.json"
    
    def run(self, interproscan_info: Dict = None, blast_info: Dict = None,
            interproscan_file: str = None, blast_file: str = None,
            output_dir: str = "output"):
        """
        运行GO整合管道
        
        Args:
            interproscan_info: InterProScan结果字典
            blast_info: BLAST结果字典  
            interproscan_file: InterProScan结果文件路径
            blast_file: BLAST结果文件路径
            output_dir: 输出目录
        """
        # 加载数据
        if interproscan_info is None and interproscan_file:
            with open(interproscan_file, 'r') as f:
                interproscan_info = json.load(f)
        
        if blast_info is None and blast_file:
            with open(blast_file, 'r') as f:
                blast_info = json.load(f)
        
        if not interproscan_info or not blast_info:
            raise ValueError("必须提供interproscan_info和blast_info数据或文件路径")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        print("开始第一层筛选...")
        # 第一层筛选
        protein_go_dict = self.first_level_filtering(interproscan_info, blast_info)
        
        if not self.use_protrek:
            # 不使用第二层筛选，直接保存结果
            output_file = os.path.join(output_dir, self.generate_filename("go_integration"))
            with open(output_file, 'w') as f:
                for protein_id, go_ids in protein_go_dict.items():
                    result = {"protein_id": protein_id, "GO_id": go_ids}
                    f.write(json.dumps(result) + '\n')
            
            print(f"第一层筛选完成，结果已保存到: {output_file}")
            return output_file
        
        print("开始第二层筛选...")
        # 提取蛋白质序列
        protein_sequences = {}
        for protein_id, data in interproscan_info.items():
            protein_sequences[protein_id] = data.get('sequence', '')
        
        # 计算ProTrek分数
        protrek_results = self.calculate_protrek_scores(protein_sequences, protein_go_dict)
        
        # 保存中间结果
        intermediate_file = os.path.join(output_dir, self.generate_filename("go_integration", is_intermediate=True))
        with open(intermediate_file, 'w') as f:
            for result in protrek_results.values():
                f.write(json.dumps(result) + '\n')
        
        print(f"ProTrek分数计算完成，中间结果已保存到: {intermediate_file}")
        
        # 第二层筛选
        if self.protrek_threshold is not None:
            final_results = self.second_level_filtering(protrek_results)
            
            # 保存最终结果
            final_file = os.path.join(output_dir, self.generate_filename("go_integration"))
            with open(final_file, 'w') as f:
                for protein_id, go_ids in final_results.items():
                    result = {"protein_id": protein_id, "GO_id": go_ids}
                    f.write(json.dumps(result) + '\n')
            
            print(f"第二层筛选完成，最终结果已保存到: {final_file}")
            return final_file, intermediate_file
        
        return intermediate_file

def main():
    parser = argparse.ArgumentParser(description="GO信息整合管道：整合InterProScan和BLAST结果，可选使用ProTrek进行二次筛选")
    parser.add_argument("--interproscan_file", type=str, required=True, 
                       help="InterProScan结果文件路径（JSON格式）")
    parser.add_argument("--blast_file", type=str, required=True, 
                       help="BLAST结果文件路径（JSON格式）")
    parser.add_argument("--identity", type=int, default=80, 
                       help="BLAST identity阈值（0-100，默认：80）")
    parser.add_argument("--coverage", type=int, default=80, 
                       help="BLAST coverage阈值（0-100，默认：80）")
    parser.add_argument("--evalue", type=float, default=1e-50, 
                       help="BLAST E-value阈值（默认：1e-50）")
    parser.add_argument("--topk", type=int, default=2, 
                       help="使用BLAST的topk结果（默认：2，设为0则使用阈值策略）")
    parser.add_argument("--protrek_threshold", type=float, default=None, 
                       help="ProTrek分数阈值（仅当use_protrek为True时使用）")
    parser.add_argument("--use_protrek", action="store_true", 
                       help="是否使用第二层ProTrek筛选（需要GPU支持）")
    parser.add_argument("--output_dir", type=str, default="go_integration_results", 
                       help="输出目录路径（默认：go_integration_results）")
    
    args = parser.parse_args()
    
    # 创建管道实例
    pipeline = GOIntegrationPipeline(
        identity_threshold=args.identity,
        coverage_threshold=args.coverage,
        evalue_threshold=args.evalue,
        topk=args.topk,
        protrek_threshold=args.protrek_threshold,
        use_protrek=args.use_protrek
    )
    
    # 运行管道
    pipeline.run(
        interproscan_file=args.interproscan_file,
        blast_file=args.blast_file,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main() 