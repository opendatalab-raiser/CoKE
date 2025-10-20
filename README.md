# 蛋白质功能分析集成管道

这是一个用于蛋白质功能分析的集成管道，整合了BLAST、InterProScan、GO注释和LLM推理等多种工具，可以自动化完成从序列到功能预测的完整流程。

## 功能特点

- **多工具集成**: 整合BLAST、InterProScan、ProTrek等多个生物信息学工具
- **并行处理**: 支持大规模并行处理，可自定义prompt生成和LLM推理的进程数
- **灵活配置**: 支持酶功能预测和通用功能预测两种模式
- **中间结果保存**: 自动保存BLAST、InterProScan等中间结果，方便复用
- **QA数据集生成**: 支持基于LMDB的问答数据集生成

## 目录结构

```
Lost_in_tokenization/
├── examples/                    # 示例文件
│   ├── input.fasta             # 示例输入FASTA文件
│   └── pids.txt                # 示例蛋白质ID列表
├── tools/                       # 工具模块
│   ├── blast.py                # BLAST工具
│   ├── interproscan.py         # InterProScan工具
│   └── go_integration_pipeline.py  # GO整合管道
├── utils/                       # 工具函数
│   ├── utils.py                # 通用工具函数
│   ├── prompts.py              # LLM提示模板
│   ├── openai_access.py        # OpenAI API访问
│   ├── get_protrek_text.py     # ProTrek文本获取
│   ├── generate_protein_prompt.py  # 蛋白质提示生成
│   ├── protein_go_analysis.py  # GO分析
│   └── mpr.py                  # 多进程运行器
├── scripts/                     # 脚本
│   └── run_protrek_text.sh     # ProTrek批量运行脚本
├── data/                        # 数据目录
│   └── raw_data/               # 原始数据
│       ├── all_pfam_descriptions.json  # Pfam描述
│       ├── go.json             # GO信息
│       └── interpro_data.json  # InterPro数据
├── integrated_pipeline.py       # 主管道脚本
├── setup.sh                     # 环境安装脚本
└── README.md                    # 本文件
```

## 安装

### 1. 环境配置

首先运行安装脚本来配置环境：

```bash
bash setup.sh
```

该脚本会：
- 安装Python依赖包
- 配置BLAST数据库
- 安装InterProScan

### 2. 依赖工具

确保以下工具已正确安装：

- **BLAST+**: 用于序列相似性搜索
- **InterProScan**: 用于蛋白质结构域和功能位点预测
- **Python 3.8+**: 推荐使用Python 3.8或更高版本

### 3. 数据准备

将以下数据文件放置在`data/raw_data/`目录下：

- `all_pfam_descriptions.json`: Pfam结构域描述
- `go.json`: GO术语定义
- `interpro_data.json`: InterPro数据库信息

## 使用方法

### 基本用法

#### 1. 蛋白质功能预测（非酶）

最简单的用法，使用默认参数：

```bash
python integrated_pipeline.py \
    --output_dir output/my_analysis
```

这将使用`examples/input.fasta`作为输入，运行完整的分析管道。

#### 2. 酶功能预测（EC号预测）

如果分析的是酶，添加`--is_enzyme`参数：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/enzyme_analysis \
    --is_enzyme
```

#### 3. 使用自定义FASTA文件

```bash
python integrated_pipeline.py \
    --input_fasta path/to/your/proteins.fasta \
    --output_dir output/custom_analysis
```

### 高级用法

#### 1. 自定义并行参数

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/parallel_analysis \
    --n_process_prompt 128 \
    --n_process_llm 32 \
    --blast_num_threads 128
```

参数说明：
- `--n_process_prompt`: 生成prompt的并行进程数（默认256）
- `--n_process_llm`: 生成LLM答案的并行进程数（默认64）
- `--blast_num_threads`: BLAST使用的线程数（默认256）

#### 2. 包含ProTrek信息

如果需要使用ProTrek进行功能预测：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/protrek_analysis \
    --protrek_dir output/protrek_results \
    --selected_info_types motif go protrek
```

管道会自动检查ProTrek结果是否完整，如果缺失会自动运行ProTrek脚本。

**注意**: ProTrek脚本可能需要很长时间运行（特别是在网络不稳定时）。如果你已经手动运行了ProTrek脚本，可以使用`--skip_protrek_check`跳过自动检查：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/protrek_analysis \
    --protrek_dir output/protrek_results \
    --selected_info_types motif go protrek \
    --skip_protrek_check
```

#### 3. 跳过BLAST和InterProScan（使用已有结果）

如果已经运行过BLAST和InterProScan，可以直接使用已有结果：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/reuse_analysis \
    --interproscan_info_path output/previous/tool_results/interproscan_info.json \
    --blast_info_path output/previous/tool_results/blast_info.json
```

#### 4. 自定义信息类型

选择不同的信息来源：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/custom_info \
    --selected_info_types motif go domain family
```

可选的信息类型：
- `motif`: Pfam结构域
- `go`: GO注释
- `protrek`: ProTrek预测
- `domain`: InterPro结构域
- `family`: 蛋白质家族

### 输出结果

管道运行完成后，输出目录结构如下：

```
output_dir/
├── tool_results/                # 工具运行的中间结果
│   ├── interproscan_info.json  # InterProScan结果
│   └── blast_info.json         # BLAST结果
└── llm_answers/                 # LLM生成的答案
    ├── protein_id_1.json       # 蛋白质1的分析结果
    ├── protein_id_2.json       # 蛋白质2的分析结果
    └── ...
```

LLM答案文件格式（非QA模式）：
```json
{
  "protein_id": "P40571",
  "prompt": "完整的输入prompt...",
  "llm_answer": "LLM生成的功能分析..."
}
```

LLM答案文件格式（QA模式）：
```json
{
  "protein_id": "P40571",
  "index": 0,
  "question": "What is the function of this protein?",
  "ground_truth": "标准答案...",
  "llm_answer": "LLM生成的答案...",
  "question_type": "function"
}
```

## 完整参数列表

```bash
python integrated_pipeline.py --help
```

### 必需参数

- `--output_dir`: 输出目录路径

### 可选参数

**输入/输出：**
- `--input_fasta`: 输入FASTA文件路径（默认：examples/input.fasta）
- `--temp_dir`: 临时文件目录（默认：temp）

**跳过步骤：**
- `--interproscan_info_path`: InterProScan结果文件（如提供则跳过InterProScan）
- `--blast_info_path`: BLAST结果文件（如提供则跳过BLAST）

**并行处理：**
- `--n_process_prompt`: 生成prompt的进程数（默认：256）
- `--n_process_llm`: 生成LLM答案的进程数（默认：64）

**BLAST参数：**
- `--blast_database`: BLAST数据库名称（默认：uniprot_swissprot）
- `--expect_value`: E-value阈值（默认：0.01）
- `--blast_num_threads`: BLAST线程数（默认：256）

**InterProScan参数：**
- `--interproscan_path`: InterProScan程序路径

**GO整合：**
- `--go_topk`: GO整合的topk参数（默认：2）

**Prompt生成：**
- `--selected_info_types`: 信息类型列表（默认：motif go）
- `--pfam_descriptions_path`: Pfam描述文件路径
- `--go_info_path`: GO信息文件路径
- `--interpro_data_path`: InterPro数据文件路径
- `--lmdb_path`: LMDB数据库路径（用于QA数据集生成）
- `--protrek_dir`: ProTrek结果目录
- `--is_enzyme`: 是否是酶（决定使用ENZYME_PROMPT还是FUNCTION_PROMPT）
- `--skip_protrek_check`: 跳过ProTrek检查，直接使用已有结果

## 生成QA数据集

如果你已经有UniProt条目的描述信息，可以生成问答数据集：

### 1. 下载UniProt条目信息

首先，下载UniProt数据库中的蛋白质条目信息：

```bash
python download_uniprot_entry.py \
    --pids_path examples/pids.txt \
    --output_dir data/uniprot_entries
```

这将下载`pids.txt`中列出的所有蛋白质的UniProt条目信息。

### 2. 生成QA对

然后使用下载的条目信息生成问答对：

```bash
python generate_protein_qa.py \
    --uniprot_entries_dir data/uniprot_entries \
    --output_lmdb data/protein_qa.lmdb
```

### 3. 使用QA数据集运行管道

生成QA数据集后，可以在管道中使用：

```bash
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/qa_analysis \
    --lmdb_path data/protein_qa.lmdb
```

这将为每个蛋白质生成多个问答对，并用LLM回答所有问题。

## ProTrek工具

ProTrek是一个基于三模态（序列-结构-功能）的蛋白质语言模型，可以预测蛋白质功能。

### 单独运行ProTrek

```bash
bash scripts/run_protrek_text.sh \
    examples/input.fasta \
    output/protrek_results \
    3 \
    1000
```

参数说明：
- 第1个参数：输入FASTA文件
- 第2个参数：输出目录
- 第3个参数：topk（返回前k个结果，默认3）
- 第4个参数：最大重试次数（默认1000）

脚本会自动检查哪些蛋白质缺少ProTrek结果，并只处理缺失的部分。如果网络错误，会自动重试。

## 常见问题

### Q1: BLAST运行很慢怎么办？

A: 可以增加BLAST线程数：
```bash
python integrated_pipeline.py --blast_num_threads 512 ...
```

### Q2: 如何节省计算资源？

A: 如果已经运行过BLAST和InterProScan，可以复用结果：
```bash
python integrated_pipeline.py \
    --interproscan_info_path output/previous/tool_results/interproscan_info.json \
    --blast_info_path output/previous/tool_results/blast_info.json \
    ...
```

### Q3: ProTrek连接失败怎么办？

A: ProTrek使用在线API，如果网络不稳定可能会失败。脚本会自动重试，最多重试1000次（可配置）。

### Q4: 如何区分酶和非酶？

A: 使用`--is_enzyme`参数：
- 酶：`--is_enzyme`（使用ENZYME_PROMPT，输出EC号）
- 非酶：不加该参数（使用FUNCTION_PROMPT，输出通用功能描述）

### Q5: 并行进程数如何设置？

A: 建议根据你的机器配置：
- `--n_process_prompt`: 可以设置较高（如256），因为主要是IO操作
- `--n_process_llm`: 建议设置较低（如32-64），因为涉及API调用有速率限制
- `--blast_num_threads`: 根据CPU核心数设置

## 引用

如果你使用本管道，请引用相关工具：

- **BLAST**: Altschul SF, et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 1997.
- **InterProScan**: Jones P, et al. "InterProScan 5: genome-scale protein function classification." Bioinformatics. 2014.
- **ProTrek**: [ProTrek相关论文]
- **Gene Ontology**: The Gene Ontology Consortium. "The Gene Ontology resource: enriching a GOld mine." Nucleic Acids Res. 2021.

## 许可证



## 联系方式

如有问题或建议，请联系：zhuangkai@westlake.edu.cn

## 更新日志

### v1.0.0 (2025-10-19)
- 初始版本发布
- 支持BLAST、InterProScan、GO整合
- 支持ProTrek功能预测
- 支持并行处理
- 支持QA数据集生成

