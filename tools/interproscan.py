import os
import datetime
import multiprocessing
import math
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

class InterproScan():
    def __init__(self, bash_path, cpu_cores=None, temp_base_dir=None):
        self.bash_path = bash_path
        
        # 获取系统资源信息
        total_cores = multiprocessing.cpu_count()
        
        # 简单获取可用内存（GB），如果没有特殊工具就估算
        try:
            with open('/proc/meminfo', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('MemAvailable:'):
                        available_memory_gb = int(line.split()[1]) / 1024 / 1024  # KB to GB
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
            # 如果读取失败，设置一个保守值
            available_memory_gb = 8.0
        
        # 自动检测CPU核心数，预留2-3个核心给系统
        if cpu_cores is None:
            self.cpu_cores = max(1, total_cores - 3)
        else:
            self.cpu_cores = min(cpu_cores, max(1, total_cores - 2))
        
        # 按每个InterProScan进程12GB内存估算最大并行数
        MEMORY_PER_PROCESS_GB = 12
        max_parallel_by_memory = int(available_memory_gb * 0.8 / MEMORY_PER_PROCESS_GB)
        self.max_parallel_chunks = min(self.cpu_cores, max_parallel_by_memory, 1024)  # 最多1024个并行
        
        # 设置临时目录基路径
        self.temp_base_dir = temp_base_dir if temp_base_dir else "/tmp"
        
        print(f"InterProScan初始化:")
        print(f"  系统: {total_cores}核心, {available_memory_gb:.1f}GB可用内存")
        print(f"  配置: 使用{self.cpu_cores}核心, 每进程约8GB内存")
        print(f"  最大并行块数: {self.max_parallel_chunks}")
    
    def run(self, fasta_file, goterms, pathways, save_dir, parallel_chunks=None) -> dict:
        start_time = datetime.datetime.now()
        
        # 使用更快的临时目录
        temp_dir = f"{self.temp_base_dir}/interproscan_{os.getpid()}_{int(start_time.timestamp())}"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        
        # 检查序列数量决定是否需要并行处理
        seqs = self.read_fasta_to_list(fasta_file)
        seq_count = len(seqs)
        seqtype = self.is_protein_sequence(seqs)
        
        print(f"检测到{seq_count}条序列")
        
        # 简单决定是否需要并行处理
        if seq_count > 20 and parallel_chunks is None:
            # 计算分块数：min(CPU核心数, 序列数/20, 内存限制的并行数)
            optimal_chunks = min(
                self.max_parallel_chunks,  # 考虑CPU和内存限制
                math.ceil(seq_count / 20)  # 每块约20条序列
            )
            if optimal_chunks > 1:
                parallel_chunks = optimal_chunks
                print(f"序列数量较多({seq_count}条)，将分成{parallel_chunks}个块并行处理")
        
        # 确保不超过系统限制
        if parallel_chunks:
            parallel_chunks = min(parallel_chunks, self.max_parallel_chunks)
            print(f"实际使用{parallel_chunks}个并行块 (最大内存需求: {parallel_chunks * 8}GB)")
        
        if parallel_chunks and parallel_chunks > 1 and seq_count > parallel_chunks:
            return self._run_parallel(fasta_file, goterms, pathways, save_dir, parallel_chunks, temp_dir, seqtype)
        else:
            return self._run_single(fasta_file, goterms, pathways, save_dir, temp_dir, seqtype)
    
    def _run_single(self, fasta_file, goterms, pathways, save_dir, temp_dir, seqtype) -> dict:
        """单线程运行InterProScan"""
        start_time = datetime.datetime.now()
        
        if not os.path.exists(os.path.dirname(save_dir)):
            os.makedirs(os.path.dirname(save_dir))
        
        # 使用原始脚本，添加CPU参数
        cmd = f"{self.bash_path} -i {fasta_file} -o {save_dir} -f JSON --cpu {self.cpu_cores}"
        cmd += f" -T {temp_dir}"
        
        # 禁用在线预计算查找，避免网络连接问题
        cmd += " -dp"
        
        # 只运行Pfam分析
        cmd += " -appl Pfam"
        
        if goterms:
            cmd += " -goterms"
        if pathways:
            cmd += " -pa"
        if seqtype:
            cmd += " -t p"
        else:
            cmd += " -t n"
            
        print(f"执行命令: {cmd}")
        
        try:
            result = os.system(cmd)
            end_time = datetime.datetime.now()
            spend_time = (end_time - start_time).total_seconds()
            
            # 清理临时目录
            self._cleanup_temp_dir(temp_dir)
            
            if result == 0 and os.path.exists(save_dir):
                print(f"InterProScan成功完成，耗时 {spend_time:.2f}秒")
                return {"output_dir": save_dir, "duration": spend_time}
            else:
                raise Exception(f"InterProScan运行失败，返回码: {result}")
        
        except Exception as e:
            self._cleanup_temp_dir(temp_dir)
            return {"error": str(e)}
    
    def _run_parallel(self, fasta_file, goterms, pathways, save_dir, chunks, temp_dir, seqtype) -> dict:
        """并行运行InterProScan"""
        start_time = datetime.datetime.now()
        
        print(f"开始并行处理，分成{chunks}个块")
        print(f"预计内存使用: {chunks * 8}GB ({chunks}块 × 8GB/块)")
        
        # 分割FASTA文件
        chunk_files = self._split_fasta(fasta_file, chunks, temp_dir)
        
        if not os.path.exists(os.path.dirname(save_dir)):
            os.makedirs(os.path.dirname(save_dir))
        
        results = []
        errors = []
        
        # 使用线程池并行执行
        with ThreadPoolExecutor(max_workers=chunks) as executor:
            future_to_chunk = {}
            
            for i, chunk_file in enumerate(chunk_files):
                chunk_output = f"{temp_dir}/chunk_{i}_output.json"
                chunk_temp = f"{temp_dir}/chunk_{i}_temp"
                
                future = executor.submit(
                    self._run_chunk, chunk_file, goterms, pathways, 
                    chunk_output, chunk_temp, seqtype, i
                )
                future_to_chunk[future] = (i, chunk_output)
            
            # 收集结果
            for future in as_completed(future_to_chunk):
                chunk_id, chunk_output = future_to_chunk[future]
                try:
                    result = future.result()
                    if 'error' in result:
                        errors.append(f"块{chunk_id}: {result['error']}")
                    else:
                        results.append(chunk_output)
                        print(f"块{chunk_id}完成，耗时{result['duration']:.2f}秒")
                except Exception as exc:
                    errors.append(f"块{chunk_id}异常: {str(exc)}")
        
        end_time = datetime.datetime.now()
        total_time = (end_time - start_time).total_seconds()
        
        if errors:
            self._cleanup_temp_dir(temp_dir)
            return {"error": f"部分块处理失败: {'; '.join(errors)}"}
        
        # 合并结果
        try:
            self._merge_results(results, save_dir)
            self._cleanup_temp_dir(temp_dir)
            print(f"并行处理完成，总耗时{total_time:.2f}秒")
            return {"output_dir": save_dir, "duration": total_time, "chunks_processed": len(results)}
        except Exception as e:
            self._cleanup_temp_dir(temp_dir)
            return {"error": f"合并结果失败: {str(e)}"}
    
    def _run_chunk(self, chunk_file, goterms, pathways, output_file, temp_dir, seqtype, chunk_id):
        """运行单个分块"""
        start_time = datetime.datetime.now()
        
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        
        # 每个块使用单核心避免资源竞争
        cmd = f"{self.bash_path} -i {chunk_file} -o {output_file} -f JSON --cpu 1"
        cmd += f" -T {temp_dir}"
        
        # 禁用在线预计算查找，避免网络连接问题
        cmd += " -dp"
        
        # 只运行Pfam分析
        cmd += " -appl Pfam"
        
        if goterms:
            cmd += " -goterms"
        if pathways:
            cmd += " -pa"
        if seqtype:
            cmd += " -t p"
        else:
            cmd += " -t n"
        
        try:
            result = os.system(cmd)
            end_time = datetime.datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            if result == 0 and os.path.exists(output_file):
                return {"duration": duration}
            else:
                return {"error": f"命令执行失败，返回码: {result}"}
        except Exception as e:
            return {"error": str(e)}
    
    def _split_fasta(self, fasta_file, chunks, temp_dir):
        """将FASTA文件分割成多个块"""
        seqs = self.read_fasta_sequences_with_headers(fasta_file)
        total_seqs = len(seqs)
        chunk_size = math.ceil(total_seqs / chunks)
        
        chunk_files = []
        
        for i in range(chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, total_seqs)
            
            if start_idx >= total_seqs:
                break
            
            chunk_file = f"{temp_dir}/chunk_{i}.fasta"
            with open(chunk_file, 'w') as f:
                for j in range(start_idx, end_idx):
                    header, sequence = seqs[j]
                    f.write(f">{header}\n{sequence}\n")
            
            chunk_files.append(chunk_file)
            print(f"创建块{i}: {end_idx - start_idx}条序列")
        
        return chunk_files
    

    def _merge_results(self, result_files, output_file):
        """
        合并多个JSON结果文件。
        此版本专门为 InterProScan 的标准输出格式优化，
        它会合并所有文件中 "results" 键下的列表内容。
        """
        import json
        
        # 初始化一个最终的字典结构。元数据（如版本号）将从第一个有效文件中获取。
        final_merged_data = {}
        # 创建一个空列表，用于收集所有分块的 "results" 列表内容
        all_sequence_results = []
        
        is_first_file = True

        for result_file in result_files:
            # 确保结果文件存在且不为空，防止因任务失败产生空文件导致错误
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                try:
                    with open(result_file, 'r') as f:
                        data = json.load(f)
                        
                        # 从第一个有效的文件中，复制所有元数据（除了 "results"）
                        if is_first_file:
                            for key, value in data.items():
                                if key != "results":
                                    final_merged_data[key] = value
                            is_first_file = False
                        
                        # 检查 "results" 键是否存在，并且其值是一个列表
                        if "results" in data and isinstance(data["results"], list):
                            # 使用 extend 将当前文件的结果列表追加到总列表中
                            all_sequence_results.extend(data["results"])
                            
                except json.JSONDecodeError:
                    print(f"警告: JSON文件格式错误，跳过文件: {result_file}")
                except Exception as e:
                    print(f"警告: 处理文件 {result_file} 时发生未知错误: {e}")

        # 将合并后的完整结果列表赋值给最终字典的 "results" 键
        final_merged_data["results"] = all_sequence_results
        
        # 将最终合并好的数据写入目标文件
        with open(output_file, 'w') as f:
            json.dump(final_merged_data, f, indent=2)
        
        print(f"成功将 {len(result_files)} 个结果文件合并到 {output_file}")
        print(f"总共合并了 {len(all_sequence_results)} 条序列的分析结果。")
    
    
    def _cleanup_temp_dir(self, temp_dir):
        """清理临时目录"""
        import shutil
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                print(f"已清理临时目录: {temp_dir}")
            except Exception as e:
                print(f"清理临时目录失败: {e}")
    
    def read_fasta_sequences_with_headers(self, file_path):
        """读取FASTA文件，返回(header, sequence)元组列表"""
        sequences = []
        current_header = None
        current_seq = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq)))
                    current_header = line[1:]  # 移除'>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_header is not None:
                sequences.append((current_header, "".join(current_seq)))
        
        return sequences
    
    # TODO
    def is_protein_sequence(self, sequences):
        # sequence = "".join(sequences)
        # ATCG AUCG
        # if len(set(sequence.upper())) > 6:
        #     return True
        # else:
        #     return False
        return True
        
    
    def read_fasta_to_list(self, file_path):
        sequences = []
        current_header = None
        current_seq = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None: 
                        sequences.append("".join(current_seq))
                    current_header = line[1:]  
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_header is not None:
                sequences.append("".join(current_seq))
        
        return sequences


if __name__ == '__main__':
    # 测试代码
    print("=== InterProScan多线程加速测试 ===")
    
    # 初始化InterproScan，简化配置
    interproscan = InterproScan(
        bash_path="interproscan/interproscan-5.75-106.0/interproscan.sh",
        cpu_cores=8,  # 指定使用8个核心，也可以不指定让其自动检测
        temp_base_dir="/tmp"  # 指定快速存储作为临时目录
    )
    
    # 示例：使用已有的测试数据
    try:
        from utils.utils import get_protein_sequence_biopython, tofasta
        import pickle

        uids = []
        seqs = []
        
        # 如果测试文件存在，加载测试数据
        test_file = "your_test_file.fasta"
        if os.path.exists(test_file):
            with open(test_file, "rb") as f:
                datas = pickle.load(f)

            for data in datas[:100]:  # 只取前100条进行测试
                uids.append(data["uniprot_id"])
                seqs.append(data["sequence"])

            fasta_file = "example/protein_go_clean_test.fasta"
            tofasta(fasta_file, uids, seqs)
            
            print(f"准备处理{len(seqs)}条序列")
            
            # 运行InterProScan - 自动判断是否需要并行处理
            result = interproscan.run(
                fasta_file=fasta_file,
                goterms=True,
                pathways=True,
                save_dir="output/interproscan_test.json",
                parallel_chunks=4  # 强制分成4块并行处理，也可以不指定让其自动判断
            )
            
            if 'error' in result:
                print(f"处理出错: {result['error']}")
            else:
                print(f"处理完成！")
                print(f"输出文件: {result['output_dir']}")
                print(f"总耗时: {result['duration']:.2f}秒")
                if 'chunks_processed' in result:
                    print(f"并行处理了{result['chunks_processed']}个块")
        
        else:
            print("测试数据文件不存在，跳过测试")
            
    except ImportError:
        print("缺少依赖模块，请确保utils.utils模块可用")
    except Exception as e:
        print(f"测试过程中发生错误: {e}")


    