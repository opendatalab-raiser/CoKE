from Bio import ExPASy
from Bio import SeqIO
import json
from Bio.Blast import NCBIXML

def get_protein_sequence_biopython(uniprot_id):
    """
    Gets a protein sequence using its UniProt ID with BioPython.
    
    Args:
        uniprot_id (str): The UniProt ID (e.g., P12345).
    
    Returns:
        str: The protein sequence or an error message.
    """
    try:
        with ExPASy.get_sprot_raw(uniprot_id) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            return str(seq_record.seq)
    except Exception as e:
        return f"Error: {str(e)}"


def extract_interproscan_metrics(file_path, librarys="PFAM"):
    """
    Extracts protein information and domain information from InterProScan JSON results.    
    Args:
        file_path (str): Path to the InterProScan JSON result file.
        librarys (list): A list of domain libraries to extract, defaults to ["PFAM"].
    Returns:
        dict: A dictionary containing protein sequences and their corresponding domain information.      
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
            
            # Process GO information
            if match["signature"]["entry"]:
                if match["signature"]["entry"]["goXRefs"]:
                    for goXRef in match["signature"]["entry"]["goXRefs"]:
                        if goXRef["databaseName"] == "GO":
                            domain_info["GO"].append(goXRef["id"])

        protein_info[sequence] = domain_info

    return protein_info


def get_seqnid(file_path):
    """
    Parses a FASTA file and returns a dictionary of sequence IDs to sequences.
    """
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
    Writes sequences in FASTA format to a file.
    
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
            # Write sequence
            f.write(f"{seq}\n")


def extract_blast_metrics(xml_file):
    """
    Extracts the following metrics from a BLAST XML result file:
    - ID (extracts the UniProt ID)
    - Identity% (percentage of identical matches)
    - Coverage (query coverage percentage)
    - E-value
    - Bit Score
    - Positive% (percentage of similar residues)
    """
    with open(xml_file) as f:
        blast_records = NCBIXML.parse(f)
        results = {}
        
        for blast_record in blast_records:
            _results = []
            query_length = blast_record.query_length
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Extract UniProt ID (e.g., from format sp|A0A0H2ZM56|ADHE_STRP2)
                    hit_id = alignment.hit_id.split("|")[1] if "|" in alignment.hit_id else alignment.hit_id
                    
                    # Calculate key metrics
                    identity_percent = (hsp.identities / hsp.align_length) * 100
                    coverage = (hsp.align_length / query_length) * 100
                    positive_percent = (hsp.positives / hsp.align_length) * 100
                    
                    # Store the results
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
    """
    Renames keys in the InterProScan results dictionary for consistency.
    """
    new_results = {}
    for key, value in interproscan_results.items():
        if key == "PFAM":
            new_results["pfam_id"] = value
        elif key == "GO":
            new_results["go_id"] = value
        else:
            new_results[key.lower()] = value
        
    return new_results