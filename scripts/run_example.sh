#!/bin/bash

# 示例脚本：运行蛋白质功能分析管道
# 使用examples目录中的示例数据

# 设置工作目录为项目根目录
cd "$(dirname "$0")/.."

echo "=========================================="
echo "蛋白质功能分析管道 - 示例运行"
echo "=========================================="

# 示例1: 基本用法 - 非酶功能预测
echo ""
echo "示例1: 基本用法 - 非酶功能预测"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_basic \
    --selected_info_types motif go

# 示例2: 酶功能预测（EC号预测）
echo ""
echo "示例2: 酶功能预测（EC号预测）"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_enzyme \
    --selected_info_types motif go \
    --is_enzyme

# 示例3: 包含ProTrek信息
echo ""
echo "示例3: 包含ProTrek信息"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_protrek \
    --selected_info_types motif go protrek \
    --protrek_dir output/example_protrek/protrek_results

echo ""
echo "=========================================="
echo "所有示例运行完成！"
echo "=========================================="
echo ""
echo "输出结果位于:"
echo "- output/example_basic/llm_answers/"
echo "- output/example_enzyme/llm_answers/"
echo "- output/example_protrek/llm_answers/"
echo ""
echo "中间结果位于:"
echo "- output/example_basic/tool_results/"
echo "- output/example_enzyme/tool_results/"
echo "- output/example_protrek/tool_results/"

