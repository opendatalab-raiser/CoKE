#!/bin/bash

# 脚本用法: bash run_protrek_text.sh <input_fasta> <output_dir> [topk] [max_iterations]
# 参数:
#   input_fasta: 输入的FASTA文件路径
#   output_dir: 输出目录路径
#   topk: (可选) ProTrek返回的top结果数量，默认为3
#   max_iterations: (可选) 最大重试次数，默认为1000

# 检查参数
if [ $# -lt 2 ]; then
    echo "用法: bash $0 <input_fasta> <output_dir> [topk] [max_iterations]"
    echo "示例: bash $0 examples/input.fasta output/protrek_results 3 1000"
    exit 1
fi

# 参数设置
INPUT_FASTA="$1"
OUTPUT_DIR="$2"
TOPK="${3:-3}"  # 默认为3
MAX_ITERATIONS="${4:-1000}"  # 默认为1000
ITERATION=1

# 检查输入文件是否存在
if [ ! -f "$INPUT_FASTA" ]; then
    echo "错误: 输入文件 $INPUT_FASTA 不存在"
    exit 1
fi

echo "开始ProTrek文本获取任务"
echo "输入文件: $INPUT_FASTA"
echo "输出目录: $OUTPUT_DIR"
echo "TopK: $TOPK"
echo "最大迭代次数: $MAX_ITERATIONS"

while [ $ITERATION -le $MAX_ITERATIONS ]; do
    echo "===================="
    echo "第 $ITERATION/$MAX_ITERATIONS 次迭代"
    echo "===================="
    
    # 运行Python脚本，直接从FASTA文件读取
    python utils/get_protrek_text.py \
        --input_fasta "$INPUT_FASTA" \
        --topk $TOPK \
        --output_dir "$OUTPUT_DIR"
    
    # 检查退出码
    if [ $? -eq 0 ]; then
        echo "✓ 所有任务完成成功!"
        break
    else
        echo "✗ 本次迭代有失败，准备重试..."
        if [ $ITERATION -lt $MAX_ITERATIONS ]; then
            echo "等待5秒后开始下次迭代..."
            sleep 5
        fi
    fi
    
    ITERATION=$((ITERATION + 1))
done

if [ $ITERATION -gt $MAX_ITERATIONS ]; then
    echo "已达到最大迭代次数，仍有部分任务未完成"
    exit 1
else
    echo "所有任务已成功完成!"
    exit 0
fi