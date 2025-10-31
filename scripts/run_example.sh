#!/bin/bash

# Example script: Run the protein function analysis pipeline
# Uses example data from the 'examples' directory

# Set the working directory to the project root
cd "$(dirname "$0")/.."

echo "=========================================="
echo "Protein Function Analysis Pipeline - Example Run"
echo "=========================================="

# Activate conda environment
echo ""
echo "Activating conda environment 'bioanalysis'..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bioanalysis
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate bioanalysis environment"
    echo "Please ensure the environment exists. Run 'bash setup.sh' first."
    exit 1
fi
echo "✓ Environment activated successfully"
echo ""

# Set foldseek environment variable (if needed)
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6.0.34" 2>/dev/null || true

# Example 1: Basic Usage - Non-enzyme Function Prediction
echo ""
echo "Example 1: Basic Usage - Non-enzyme Function Prediction"
echo "------------------------------------------"
python integrated_pipeline.py \
    --input_fasta examples/input.fasta \
    --output_dir output/example_basic \
    --selected_info_types motif go

# Example 2: Enzyme Function Prediction
echo ""
echo "Example 2: Enzyme Function Prediction"
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

# Example 4: Using Foldseek for Structure-based Search
echo ""
echo "Example 4: Using Foldseek for Structure-based Search"
echo "------------------------------------------"
# Check if PDB files exist in examples directory
if ls examples/*.pdb 1> /dev/null 2>&1; then
    python integrated_pipeline.py \
        --input_fasta examples/input.fasta \
        --output_dir output/example_foldseek \
        --pdb_dir examples \
        --use_foldseek \
        --foldseek_database foldseek_db/sp \
        --foldseek_num_threads 64 \
        --selected_info_types motif go
    echo "✓ Foldseek example completed"
else
    echo "⚠️  Warning: No PDB files found in examples/ directory"
    echo "   Skipping Foldseek example. To use Foldseek, place PDB files in examples/ directory."
fi

echo ""
echo "=========================================="
echo "All example runs completed!"
echo "=========================================="
echo ""
echo "Output results are located in:"
echo "- output/example_basic/llm_answers/"
echo "- output/example_enzyme/llm_answers/"
echo "- output/example_protrek/llm_answers/"
if [ -d "output/example_foldseek" ]; then
    echo "- output/example_foldseek/llm_answers/"
fi
echo ""
echo "Intermediate results are located in:"
echo "- output/example_basic/tool_results/"
echo "- output/example_enzyme/tool_results/"
echo "- output/example_protrek/tool_results/"
if [ -d "output/example_foldseek" ]; then
    echo "- output/example_foldseek/tool_results/"
fi