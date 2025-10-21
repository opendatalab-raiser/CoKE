# Context as the Key Engine

This is an integrated pipeline for protein function analysis that combines multiple tools—BLAST, InterProScan, GO annotation, and LLM-based reasoning—to automate the end‑to‑end workflow from sequence to function prediction.

## Features

* **Multi‑tool integration:** Integrates BLAST, InterProScan, ProTrek, and other bioinformatics tools.
* **Parallel processing:** Supports large‑scale parallel execution with configurable process counts for prompt generation and LLM inference.
* **Flexible modes:** Supports both enzyme function prediction and general function prediction.
* **Intermediate artifacts:** Automatically saves BLAST, InterProScan, and other intermediate results for reuse.
* **QA dataset generation:** Supports LMDB‑based question–answer dataset construction.

## Directory Structure

```
Lost_in_tokenization/
├── examples/                    # Examples
│   ├── input.fasta             # Example input FASTA file
│   └── pids.txt                # Example protein ID list
├── tools/                       # Tool modules
│   ├── blast.py                # BLAST helper
│   ├── interproscan.py         # InterProScan helper
│   └── go_integration_pipeline.py  # GO integration pipeline
├── utils/                       # Utilities
│   ├── utils.py                # Common helpers
│   ├── prompts.py              # LLM prompt templates
│   ├── openai_access.py        # OpenAI API access
│   ├── get_protrek_text.py     # Fetch ProTrek text
│   ├── generate_protein_prompt.py  # Protein prompt generation
│   ├── protein_go_analysis.py  # GO analysis
│   └── mpr.py                  # Multiprocess runner
├── scripts/                     # Scripts
│   └── run_protrek_text.sh     # ProTrek batch script
├── data/                        # Data directory
│   └── raw_data/               # Raw data
│       ├── all_pfam_descriptions.json  # Pfam descriptions
│       ├── go.json             # GO information
│       └── interpro_data.json  # InterPro information
├── integrated_pipeline.py       # Main pipeline
├── setup.sh                     # Environment setup
└── README.md                    # This file
```

## Installation

### 1. Environment Setup

Run the setup script to configure the environment:

```bash
bash setup.sh
```

This script will:

* Install Python dependencies
* Configure BLAST databases
* Install InterProScan

### 2. Dependencies

Ensure the following tools are properly installed:

* **BLAST+**: For sequence similarity search
* **InterProScan**: For protein domain and functional site prediction
* **Python 3.8+**: Python 3.8 or later is recommended

### 3. Data Preparation

Place the following files under `data/raw_data/`:

* `all_pfam_descriptions.json`: Pfam domain descriptions
* `go.json`: GO term definitions
* `interpro_data.json`: InterPro database metadata

## Usage

### Basic Usage

#### 1. Protein function prediction (non‑enzyme)

The simplest invocation with default parameters:

```bash
python integrated_pipeline.py \
    --output_dir output/my_analysis
```

This uses `examples/input.fasta` as the input and runs the full analysis pipeline.

#### 2. Enzyme function prediction (EC number)

For enzyme analysis, add `--is_enzyme`:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/enzyme_analysis \
    --is_enzyme
```

#### 3. Use a custom FASTA file

```bash
python integrated_pipeline.py \
    --input_fasta path/to/your/proteins.fasta \
    --output_dir output/custom_analysis
```

### Advanced Usage

#### 1. Custom parallelism

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/parallel_analysis \
    --n_process_prompt 128 \
    --n_process_llm 32 \
    --blast_num_threads 128
```

Parameter notes:

* `--n_process_prompt`: Number of processes for prompt generation (default: 256)
* `--n_process_llm`: Number of processes for LLM answers (default: 64)
* `--blast_num_threads`: Threads used by BLAST (default: 256)

#### 2. Include ProTrek information

If ProTrek is required for function prediction:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/protrek_analysis \
    --protrek_dir output/protrek_results \
    --selected_info_types motif go protrek
```

The pipeline will automatically check whether ProTrek results are complete. If missing, it will run the ProTrek script.

**Note:** ProTrek may take a long time (especially on unstable networks). If you already ran the ProTrek script manually, skip the automatic check with `--skip_protrek_check`:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/protrek_analysis \
    --protrek_dir output/protrek_results \
    --selected_info_types motif go protrek \
    --skip_protrek_check
```

#### 3. Skip BLAST and InterProScan (reuse existing results)

If BLAST and InterProScan have already been run, reuse prior outputs:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/reuse_analysis \
    --interproscan_info_path output/previous/tool_results/interproscan_info.json \
    --blast_info_path output/previous/tool_results/blast_info.json
```

#### 4. Customize information types

Choose which information sources to include:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/custom_info \
    --selected_info_types motif go domain family
```

Available information types:

* `motif`: Pfam domains
* `go`: GO annotations
* `protrek`: ProTrek predictions
* `domain`: InterPro domains
* `family`: Protein family

### Outputs

After the pipeline finishes, the output directory looks like this:

```
output_dir/
├── tool_results/                # Intermediate tool outputs
│   ├── interproscan_info.json  # InterProScan results
│   └── blast_info.json         # BLAST results
└── llm_answers/                 # LLM‑generated answers
    ├── protein_id_1.json       # Result for protein 1
    ├── protein_id_2.json       # Result for protein 2
    └── ...
```

LLM answer format (non‑QA mode):

```json
{
  "protein_id": "P40571",
  "prompt": "Full input prompt...",
  "llm_answer": "LLM‑generated functional analysis..."
}
```

LLM answer format (QA mode):

```json
{
  "protein_id": "P40571",
  "index": 0,
  "question": "What is the function of this protein?",
  "ground_truth": "Reference answer...",
  "llm_answer": "LLM‑generated answer...",
  "question_type": "function"
}
```

## Full Argument List

```bash
python integrated_pipeline.py --help
```

### Required

* `--output_dir`: Path to the output directory

### Optional

**I/O:**

* `--input_fasta`: Input FASTA file (default: `examples/input.fasta`)
* `--temp_dir`: Temporary directory (default: `temp`)

**Skipping steps:**

* `--interproscan_info_path`: InterProScan results file (if set, skip InterProScan)
* `--blast_info_path`: BLAST results file (if set, skip BLAST)

**Parallelism:**

* `--n_process_prompt`: Processes for prompt generation (default: 256)
* `--n_process_llm`: Processes for LLM answers (default: 64)

**BLAST parameters:**

* `--blast_database`: BLAST database name (default: `uniprot_swissprot`)
* `--expect_value`: E‑value threshold (default: `0.01`)
* `--blast_num_threads`: Number of BLAST threads (default: 256)

**InterProScan:**

* `--interproscan_path`: Path to InterProScan executable

**GO integration:**

* `--go_topk`: `topk` for GO integration (default: 2)

**Prompt generation:**

* `--selected_info_types`: List of information types (default: `motif go`)
* `--pfam_descriptions_path`: Path to Pfam description file
* `--go_info_path`: Path to GO info file
* `--interpro_data_path`: Path to InterPro metadata file
* `--lmdb_path`: Path to LMDB database (for QA dataset generation)
* `--protrek_dir`: Directory of ProTrek results
* `--is_enzyme`: Whether the sequence is an enzyme (selects `ENZYME_PROMPT` vs `FUNCTION_PROMPT`)
* `--skip_protrek_check`: Skip ProTrek check and use existing results directly

## Building a QA Dataset

If you already have UniProt entry descriptions, you can create a QA dataset.

### 1. Download UniProt entries

First, download UniProt entry records for your proteins:

```bash
python download_uniprot_entry.py \
    --pids_path examples/pids.txt \
    --output_dir data/uniprot_entries
```

This downloads UniProt entries for all proteins listed in `pids.txt`.

### 2. Generate QA pairs

Then generate QA pairs from the downloaded entries:

```bash
python generate_protein_qa.py \
    --uniprot_entries_dir data/uniprot_entries \
    --output_lmdb data/protein_qa.lmdb
```

### 3. Run the pipeline with the QA dataset

After generating the QA dataset, use it in the pipeline:

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/qa_analysis \
    --lmdb_path data/protein_qa.lmdb
```

This will create multiple QA pairs per protein and have the LLM answer all questions.

## ProTrek Tool

ProTrek is a tri‑modal (sequence–structure–function) protein language model that can predict protein function.

### Run ProTrek separately

```bash
bash scripts/run_protrek_text.sh \
    examples/input.fasta \
    output/protrek_results \
    3 \
    1000
```

Arguments:

* 1st: Input FASTA file
* 2nd: Output directory
* 3rd: `topk` (number of results to return, default 3)
* 4th: Max retry count (default 1000)

The script automatically detects proteins missing ProTrek results and only processes the missing ones. It will retry on network errors.

## FAQ

### Q1: BLAST is slow—what can I do?

A: Increase the number of BLAST threads:

```bash
python integrated_pipeline.py --blast_num_threads 512 ...
```

### Q2: How can I save compute?

A: Reuse prior BLAST/InterProScan results when available:

```bash
python integrated_pipeline.py \
    --interproscan_info_path output/previous/tool_results/interproscan_info.json \
    --blast_info_path output/previous/tool_results/blast_info.json \
    ...
```

### Q3: ProTrek connection failures?

A: ProTrek uses an online API and may fail if the network is unstable. The script retries automatically up to 1000 times (configurable).

### Q4: How do I distinguish enzymes vs non‑enzymes?

A: Use `--is_enzyme`:

* Enzyme: add `--is_enzyme` (uses `ENZYME_PROMPT`, outputs EC numbers)
* Non‑enzyme: omit the flag (uses `FUNCTION_PROMPT`, outputs general function descriptions)

### Q5: How should I set the parallelism?

A: Tune based on your hardware and rate limits:

* `--n_process_prompt`: Can be high (e.g., 256) since it is mostly I/O‑bound
* `--n_process_llm`: Prefer moderate values (e.g., 32–64) due to API rate limits
* `--blast_num_threads`: Set according to CPU core count

## Citation

If you use this pipeline, please cite the relevant tools:

* **BLAST**: Altschul SF, et al. "Gapped BLAST and PSI‑BLAST: a new generation of protein database search programs." *Nucleic Acids Res.* 1997.
* **InterProScan**: Jones P, et al. "InterProScan 5: genome‑scale protein function classification." *Bioinformatics.* 2014.
* **ProTrek**: [ProTrek‑related publication]
* **Gene Ontology**: The Gene Ontology Consortium. "The Gene Ontology resource: enriching a GOld mine." *Nucleic Acids Res.* 2021.

## License

(To be added.)

## Contact

For questions or suggestions, please contact: [zhuangkai@westlake.edu.cn](mailto:zhuangkai@westlake.edu.cn)

## Changelog

### v1.0.0 (2025‑10‑19)

* Initial release
* Support for BLAST, InterProScan, and GO integration
* Support for ProTrek function prediction
* Support for parallel processing
* Support for QA dataset generation
