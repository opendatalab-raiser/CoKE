#!/bin/bash

# Example script: Run the protein function analysis pipeline
# Uses example data from the 'examples' directory

# Set the working directory to the project root
cd "$(dirname "$0")/.."

echo "=========================================="
echo "Protein Function Analysis Pipeline - Example Run"
echo "=========================================="

# Example 1: Basic Usage - Non-enzyme Function Prediction
echo ""
echo "Example 1: Basic Usage - Non-enzyme Function Prediction"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_basic \
    --selected_info_types motif go

# Example 2: Enzyme Function Prediction (EC Number Prediction)
echo ""
echo "Example 2: Enzyme Function Prediction (EC Number Prediction)"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_enzyme \
    --selected_info_types motif go \
    --is_enzyme

# Example 3: Including ProTrek Information
echo ""
echo "Example 3: Including ProTrek Information"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_protrek \
    --selected_info_types motif go protrek \
    --protrek_dir output/example_protrek/protrek_results

echo ""
echo "=========================================="
echo "All example runs completed!"
echo "=========================================="
echo ""
echo "Output results are located in:"
echo "- output/example_basic/llm_answers/"
echo "- output/example_enzyme/llm_answers/"
echo "- output/example_protrek/llm_answers/"
echo ""
echo "Intermediate results are located in:"
echo "- output/example_basic/tool_results/"
echo "- output/example_enzyme/tool_results/"
echo "- output/example_protrek/tool_results/"