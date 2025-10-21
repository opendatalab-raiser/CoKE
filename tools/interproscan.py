import os
import datetime
import multiprocessing
import math
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

class InterproScan():
    def __init__(self, bash_path, cpu_cores=None, temp_base_dir=None):
        self.bash_path = bash_path
        
        # Get system resource information
        total_cores = multiprocessing.cpu_count()
        
        # Simply get available memory (GB), estimate if special tools are not available
        try:
            with open('/proc/meminfo', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('MemAvailable:'):
                        available_memory_gb = int(line.split()[1]) / 1024 / 1024  # KB to GB
                        break
                else:
                    # If MemAvailable is not present, estimate using MemFree + Buffers + Cached
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
            # If reading fails, set a conservative value
            available_memory_gb = 8.0
        
        # Automatically detect the number of CPU cores, reserving 2-3 cores for the system
        if cpu_cores is None:
            self.cpu_cores = max(1, total_cores - 3)
        else:
            self.cpu_cores = min(cpu_cores, max(1, total_cores - 2))
        
        # Estimate the maximum number of parallel processes based on 12GB of memory per InterProScan process
        MEMORY_PER_PROCESS_GB = 12
        max_parallel_by_memory = int(available_memory_gb * 0.8 / MEMORY_PER_PROCESS_GB)
        self.max_parallel_chunks = min(self.cpu_cores, max_parallel_by_memory, 1024)  # Maximum of 1024 parallel processes
        
        # Set the base path for the temporary directory
        self.temp_base_dir = temp_base_dir if temp_base_dir else "/tmp"
        
        print(f"InterProScan initialized:")
        print(f"  System: {total_cores} cores, {available_memory_gb:.1f}GB available memory")
        print(f"  Configuration: Using {self.cpu_cores} cores, approx. 8GB memory per process")
        print(f"  Max parallel chunks: {self.max_parallel_chunks}")
    
    def run(self, fasta_file, goterms, pathways, save_dir, parallel_chunks=None) -> dict:
        start_time = datetime.datetime.now()
        
        # Use a faster temporary directory
        temp_dir = f"{self.temp_base_dir}/interproscan_{os.getpid()}_{int(start_time.timestamp())}"
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        
        # Check the number of sequences to decide if parallel processing is needed
        seqs = self.read_fasta_to_list(fasta_file)
        seq_count = len(seqs)
        seqtype = self.is_protein_sequence(seqs)
        
        print(f"Detected {seq_count} sequences")
        
        # Simply decide whether to use parallel processing
        if seq_count > 20 and parallel_chunks is None:
            # Calculate number of chunks: min(CPU cores, num_sequences/20, memory-limited parallel count)
            optimal_chunks = min(
                self.max_parallel_chunks,  # Consider CPU and memory limits
                math.ceil(seq_count / 20)  # Approximately 20 sequences per chunk
            )
            if optimal_chunks > 1:
                parallel_chunks = optimal_chunks
                print(f"Sequence count is high ({seq_count}), splitting into {parallel_chunks} chunks for parallel processing")
        
        # Ensure it does not exceed system limits
        if parallel_chunks:
            parallel_chunks = min(parallel_chunks, self.max_parallel_chunks)
            print(f"Actually using {parallel_chunks} parallel chunks (Max memory requirement: {parallel_chunks * 8}GB)")
        
        if parallel_chunks and parallel_chunks > 1 and seq_count > parallel_chunks:
            return self._run_parallel(fasta_file, goterms, pathways, save_dir, parallel_chunks, temp_dir, seqtype)
        else:
            return self._run_single(fasta_file, goterms, pathways, save_dir, temp_dir, seqtype)
    
    def _run_single(self, fasta_file, goterms, pathways, save_dir, temp_dir, seqtype) -> dict:
        """Run InterProScan in a single thread."""
        start_time = datetime.datetime.now()
        
        if not os.path.exists(os.path.dirname(save_dir)):
            os.makedirs(os.path.dirname(save_dir))
        
        # Use the original script, add the CPU parameter
        cmd = f"{self.bash_path} -i {fasta_file} -o {save_dir} -f JSON --cpu {self.cpu_cores}"
        cmd += f" -T {temp_dir}"
        
        # Disable online pre-calculated match lookup to avoid network issues
        cmd += " -dp"
        
        # Run only Pfam analysis
        cmd += " -appl Pfam"
        
        if goterms:
            cmd += " -goterms"
        if pathways:
            cmd += " -pa"
        if seqtype:
            cmd += " -t p"
        else:
            cmd += " -t n"
            
        print(f"Executing command: {cmd}")
        
        try:
            result = os.system(cmd)
            end_time = datetime.datetime.now()
            spend_time = (end_time - start_time).total_seconds()
            
            # Clean up the temporary directory
            self._cleanup_temp_dir(temp_dir)
            
            if result == 0 and os.path.exists(save_dir):
                print(f"InterProScan completed successfully in {spend_time:.2f} seconds")
                return {"output_dir": save_dir, "duration": spend_time}
            else:
                raise Exception(f"InterProScan failed with return code: {result}")
        
        except Exception as e:
            self._cleanup_temp_dir(temp_dir)
            return {"error": str(e)}
    
    def _run_parallel(self, fasta_file, goterms, pathways, save_dir, chunks, temp_dir, seqtype) -> dict:
        """Run InterProScan in parallel."""
        start_time = datetime.datetime.now()
        
        print(f"Starting parallel processing, splitting into {chunks} chunks")
        print(f"Estimated memory usage: {chunks * 8}GB ({chunks} chunks × 8GB/chunk)")
        
        # Split the FASTA file
        chunk_files = self._split_fasta(fasta_file, chunks, temp_dir)
        
        if not os.path.exists(os.path.dirname(save_dir)):
            os.makedirs(os.path.dirname(save_dir))
        
        results = []
        errors = []
        
        # Use a thread pool for parallel execution
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
            
            # Collect results
            for future in as_completed(future_to_chunk):
                chunk_id, chunk_output = future_to_chunk[future]
                try:
                    result = future.result()
                    if 'error' in result:
                        errors.append(f"Chunk {chunk_id}: {result['error']}")
                    else:
                        results.append(chunk_output)
                        print(f"Chunk {chunk_id} finished in {result['duration']:.2f} seconds")
                except Exception as exc:
                    errors.append(f"Chunk {chunk_id} exception: {str(exc)}")
        
        end_time = datetime.datetime.now()
        total_time = (end_time - start_time).total_seconds()
        
        if errors:
            self._cleanup_temp_dir(temp_dir)
            return {"error": f"Processing failed for some chunks: {'; '.join(errors)}"}
        
        # Merge results
        try:
            self._merge_results(results, save_dir)
            self._cleanup_temp_dir(temp_dir)
            print(f"Parallel processing completed, total time {total_time:.2f} seconds")
            return {"output_dir": save_dir, "duration": total_time, "chunks_processed": len(results)}
        except Exception as e:
            self._cleanup_temp_dir(temp_dir)
            return {"error": f"Failed to merge results: {str(e)}"}
    
    def _run_chunk(self, chunk_file, goterms, pathways, output_file, temp_dir, seqtype, chunk_id):
        """Run a single chunk."""
        start_time = datetime.datetime.now()
        
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        
        # Use a single core for each chunk to avoid resource contention
        cmd = f"{self.bash_path} -i {chunk_file} -o {output_file} -f JSON --cpu 1"
        cmd += f" -T {temp_dir}"
        
        # Disable online pre-calculated match lookup to avoid network issues
        cmd += " -dp"
        
        # Run only Pfam analysis
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
                return {"error": f"Command execution failed with return code: {result}"}
        except Exception as e:
            return {"error": str(e)}
    
    def _split_fasta(self, fasta_file, chunks, temp_dir):
        """Split a FASTA file into multiple chunks."""
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
            print(f"Created chunk {i}: {end_idx - start_idx} sequences")
        
        return chunk_files
    

    def _merge_results(self, result_files, output_file):
        """
        Merge multiple JSON result files.
        This version is specifically optimized for InterProScan's standard output format,
        merging the list content under the "results" key from all files.
        """
        import json
        
        # Initialize a final dictionary structure. Metadata (like version numbers) will be taken from the first valid file.
        final_merged_data = {}
        # Create an empty list to collect the "results" list content from all chunks
        all_sequence_results = []
        
        is_first_file = True

        for result_file in result_files:
            # Ensure the result file exists and is not empty to prevent errors from empty files due to failed tasks
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                try:
                    with open(result_file, 'r') as f:
                        data = json.load(f)
                        
                        # From the first valid file, copy all metadata (except for "results")
                        if is_first_file:
                            for key, value in data.items():
                                if key != "results":
                                    final_merged_data[key] = value
                            is_first_file = False
                        
                        # Check if the "results" key exists and its value is a list
                        if "results" in data and isinstance(data["results"], list):
                            # Use extend to append the current file's results list to the master list
                            all_sequence_results.extend(data["results"])
                            
                except json.JSONDecodeError:
                    print(f"Warning: Malformed JSON, skipping file: {result_file}")
                except Exception as e:
                    print(f"Warning: An unknown error occurred while processing file {result_file}: {e}")

        # Assign the complete merged results list to the "results" key of the final dictionary
        final_merged_data["results"] = all_sequence_results
        
        # Write the final merged data to the target file
        with open(output_file, 'w') as f:
            json.dump(final_merged_data, f, indent=2)
        
        print(f"Successfully merged {len(result_files)} result files into {output_file}")
        print(f"A total of {len(all_sequence_results)} sequence analysis results were merged.")
    
    
    def _cleanup_temp_dir(self, temp_dir):
        """Clean up the temporary directory."""
        import shutil
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                print(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                print(f"Failed to clean up temporary directory: {e}")
    
    def read_fasta_sequences_with_headers(self, file_path):
        """Read a FASTA file and return a list of (header, sequence) tuples."""
        sequences = []
        current_header = None
        current_seq = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq)))
                    current_header = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_header is not None:
                sequences.append((current_header, "".join(current_seq)))
        
        return sequences
    
    # TODO
    def is_protein_sequence(self, sequences):
        # A simple heuristic to distinguish protein from nucleotide sequences.
        # This is not foolproof and can be improved.
        # sequence = "".join(sequences)
        # Legal protein characters are more than the typical ATCG/AUCG of nucleotides.
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
    # Test code
    print("=== InterProScan Multi-threading Acceleration Test ===")
    
    # Initialize InterproScan, simplified configuration
    interproscan = InterproScan(
        bash_path="interproscan/interproscan-5.75-106.0/interproscan.sh",
        cpu_cores=8,  # Specify using 8 cores, or leave it unspecified for auto-detection
        temp_base_dir="/tmp"  # Specify fast storage as the temporary directory
    )
    
    # Example: Use existing test data
    try:
        from utils.utils import get_protein_sequence_biopython, tofasta
        import pickle

        uids = []
        seqs = []
        
        # If the test file exists, load the test data
        test_file = "your_test_file.fasta"
        if os.path.exists(test_file):
            with open(test_file, "rb") as f:
                datas = pickle.load(f)

            for data in datas[:100]:  # Only take the first 100 entries for testing
                uids.append(data["uniprot_id"])
                seqs.append(data["sequence"])

            fasta_file = "example/protein_go_clean_test.fasta"
            tofasta(fasta_file, uids, seqs)
            
            print(f"Preparing to process {len(seqs)} sequences")
            
            # Run InterProScan - automatically determine if parallel processing is needed
            result = interproscan.run(
                fasta_file=fasta_file,
                goterms=True,
                pathways=True,
                save_dir="output/interproscan_test.json",
                parallel_chunks=4  # Force splitting into 4 chunks for parallel processing, or leave it unspecified for auto-detection
            )
            
            if 'error' in result:
                print(f"Error during processing: {result['error']}")
            else:
                print(f"Processing complete!")
                print(f"Output file: {result['output_dir']}")
                print(f"Total time: {result['duration']:.2f}s")
                if 'chunks_processed' in result:
                    print(f"Processed {result['chunks_processed']} chunks in parallel")
        
        else:
            print("Test data file not found, skipping test")
            
    except ImportError:
        print("Missing dependency, please ensure the utils.utils module is available")
    except Exception as e:
        print(f"An error occurred during testing: {e}")