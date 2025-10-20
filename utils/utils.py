from Bio import ExPASy
from Bio import SeqIO
import json
from Bio.Blast import NCBIXML

def get_protein_sequence_biopython(uniprot_id):
    """
    使用BioPython通过UniProt ID获取蛋白质序列
    
    参数:
        uniprot_id (str): UniProt ID (如P12345)
    
    返回:
        str: 蛋白质序列或错误信息
    """
    try:
        with ExPASy.get_sprot_raw(uniprot_id) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            return str(seq_record.seq)
    except Exception as e:
        return f"Error: {str(e)}"


def extract_interproscan_metrics(file_path, librarys="PFAM"):
    """
    从InterProScan JSON结果中提取蛋白质信息和域信息。    
    参数:
        file_path (str): InterProScan JSON结果文件路径
        librarys (list): 需要提取的域库列表，默认为["PFAM"]
    返回:
        dict: 包含蛋白质序列和对应域信息的字典      
    """
    protein_info = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        data = json.load(file)
    results = data["results"]

    for protein in results:
        sequence = protein["sequence"]
        domain_info = {}
        for library in librarys:
            domain_info[library] = []
        domain_info["GO"] = []

        matches = protein["matches"]
        for match in matches:
            if match["signature"]["signatureLibraryRelease"]["library"] in librarys:
                if match["signature"]["entry"]:
                    domain_info[match["signature"]["signatureLibraryRelease"]["library"]].append({match["signature"]["accession"]: match["signature"]["entry"]["accession"]})
                else:
                    domain_info[match["signature"]["signatureLibraryRelease"]["library"]].append({match["signature"]["accession"]: None})
            
            # 处理GO信息
            if match["signature"]["entry"]:
                if match["signature"]["entry"]["goXRefs"]:
                    for goXRef in match["signature"]["entry"]["goXRefs"]:
                        if goXRef["databaseName"] == "GO":
                            domain_info["GO"].append(goXRef["id"])

        protein_info[sequence] = domain_info

    return protein_info


def get_seqnid(file_path):
    seq_dict = {}
    current_header = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None: 
                    seq_dict[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]  # Take only the first part before whitespace
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header is not None:
            seq_dict[current_header] = "".join(current_seq)
    
    return seq_dict


def tofasta(fasta_path, uids, seqs):
    """
    Write sequences in FASTA format to a file.
    
    Parameters:
    - fasta_path: str, path to the output FASTA file
    - uids: list of str, sequence identifiers (headers)
    - seqs: list of str, corresponding sequences
    """
    if len(uids) != len(seqs):
        raise ValueError("Length of uids and seqs must be equal")
    
    with open(fasta_path, 'w') as f:
        for uid, seq in zip(uids, seqs):
            # Write header line starting with '>' followed by the uid
            f.write(f">{uid}\n")
            # Write sequence (you may want to split long sequences into multiple lines)
            f.write(f"{seq}\n")


def extract_blast_metrics(xml_file):
    """
    从BLAST XML结果中提取以下指标：
    - ID (提取UniProt ID)
    - Identity% (相似度百分比)
    - Coverage (覆盖率)
    - E-value
    - Bit Score
    - Positive% (相似残基百分比)
    """
    with open(xml_file) as f:
        blast_records = NCBIXML.parse(f)
        results = {}
        
        for blast_record in blast_records:
            _results = []
            query_length = blast_record.query_length
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # 提取UniProt ID (格式如 sp|A0A0H2ZM56|ADHE_STRP2)
                    hit_id = alignment.hit_id.split("|")[1] if "|" in alignment.hit_id else alignment.hit_id
                    
                    # 计算关键指标
                    identity_percent = (hsp.identities / hsp.align_length) * 100
                    coverage = (hsp.align_length / query_length) * 100
                    positive_percent = (hsp.positives / hsp.align_length) * 100
                    
                    # 存储结果
                    _results.append({
                        "ID": hit_id,
                        "Identity%": round(identity_percent, 2),
                        "Coverage%": round(coverage, 2),
                        "E-value": f"{hsp.expect:.1e}" if hsp.expect < 0.001 else round(hsp.expect, 4),
                        "Bit Score": round(hsp.bits, 1),
                        "Positive%": round(positive_percent, 2)
                    })
            results[blast_record.query] = _results
        return results


def rename_interproscan_keys(interproscan_results):
    new_results = {}
    for key, value in interproscan_results.items():
        if key == "PFAM":
            new_results["pfam_id"] = value
        elif key == "GO":
            new_results["go_id"] = value
        else:
            new_results[key.lower()] = value
        
    return new_results

