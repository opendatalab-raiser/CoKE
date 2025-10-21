#!/bin/bash

# Usage: bash run_protrek_text.sh <input_fasta> <output_dir> [topk] [max_iterations]
# Arguments:
#   input_fasta: Path to the input FASTA file
#   output_dir: Path to the output directory
#   topk: (Optional) Number of top results returned by ProTrek, default is 3
#   max_iterations: (Optional) Maximum number of retry attempts, default is 1000

# Check for arguments
if [ $# -lt 2 ]; then
    echo "Usage: bash $0 <input_fasta> <output_dir> [topk] [max_iterations]"
    echo "Example: bash $0 examples/input.fasta output/protrek_results 3 1000"
    exit 1
fi

# Set parameters
INPUT_FASTA="$1"
OUTPUT_DIR="$2"
TOPK="${3:-3}"  # Default to 3
MAX_ITERATIONS="${4:-1000}"  # Default to 1000
ITERATION=1

# Check if the input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file $INPUT_FASTA not found"
    exit 1
fi

echo "Starting ProTrek text retrieval task"
echo "Input file: $INPUT_FASTA"
echo "Output directory: $OUTPUT_DIR"
echo "TopK: $TOPK"
echo "Max iterations: $MAX_ITERATIONS"

while [ $ITERATION -le $MAX_ITERATIONS ]; do
    echo "===================="
    echo "Iteration $ITERATION/$MAX_ITERATIONS"
    echo "===================="
    
    # Run the Python script, reading directly from the FASTA file
    python utils/get_protrek_text.py \
        --input_fasta "$INPUT_FASTA" \
        --topk $TOPK \
        --output_dir "$OUTPUT_DIR"
    
    # Check the exit code
    if [ $? -eq 0 ]; then
        echo "✓ All tasks completed successfully!"
        break
    else
        echo "✗ Failures occurred in this iteration, preparing to retry..."
        if [ $ITERATION -lt $MAX_ITERATIONS ]; then
            echo "Waiting 5 seconds before starting the next iteration..."
            sleep 5
        fi
    fi
    
    ITERATION=$((ITERATION + 1))
done

if [ $ITERATION -gt $MAX_ITERATIONS ]; then
    echo "Maximum number of iterations reached, but some tasks remain unfinished."
    exit 1
else
    echo "All tasks have been completed successfully!"
    exit 0
fi