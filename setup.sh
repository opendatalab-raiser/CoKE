#!/bin/bash

# Define InterProScan version
IPS_VERSION="5.75-106.0"
CONDA_ENV_NAME="bioanalysis"

IPS_DIR="interproscan-${IPS_VERSION}"
IPS_TAR="interproscan-${IPS_VERSION}-64-bit.tar.gz"
IPS_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPS_VERSION}/${IPS_TAR}"


# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first"
    exit 1
fi

# Check if conda environment already exists
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "Found existing conda environment '${CONDA_ENV_NAME}', using it..."
else
    # Create conda environment with Java 11
    echo "Creating conda environment '${CONDA_ENV_NAME}' with OpenJDK 11..."
    conda create -y -n ${CONDA_ENV_NAME} openjdk=11 python
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create conda environment"
        exit 1
    fi
fi

# Activate conda environment
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${CONDA_ENV_NAME}

# # Install Java 11 if not present in the existing environment
# if ! java -version 2>&1 | grep -q "openjdk.*11"; then
#     echo "Installing OpenJDK 11 in existing environment..."
#     conda install -y openjdk=11
# fi

# # Create installation directory
# echo "Setting up InterProScan ${IPS_VERSION}..."
# mkdir -p interproscan
# cd interproscan || exit 1

# # Download InterProScan and checksum
# echo "Downloading InterProScan..."
# wget -nc "${IPS_URL}"
# wget -nc "${IPS_URL}.md5"

# # Verify MD5 checksum
# echo "Verifying download integrity..."
# if ! md5sum -c "${IPS_TAR}.md5"; then
#     echo "ERROR: MD5 checksum verification failed!"
#     echo "The downloaded file may be corrupted. Please try downloading again."
#     exit 1
# fi

# # Extract package
# echo "Extracting InterProScan..."
# tar -xzf "${IPS_TAR}"

# # Verify Java installation in conda env
# echo "Checking Java environment in conda env..."
# JAVA_VER=$(java -version 2>&1 | head -n 1 | awk -F '"' '{print $2}')

# if [[ "$JAVA_VER" =~ ^11\. ]]; then
#     echo "Found compatible Java version in conda env: $JAVA_VER"
# else
#     echo "Error: Java version in conda env is not 11.x (found: $JAVA_VER)"
#     exit 1
# fi

# # Run setup
# echo "Running InterProScan setup..."
# cd "${IPS_DIR}" || exit 1
# python setup.py -f interproscan.properties

# echo ""
# echo "InterProScan installation completed in conda environment '${CONDA_ENV_NAME}'!"
# echo "To use InterProScan, first activate the conda environment:"
# echo "conda activate ${CONDA_ENV_NAME}"
# echo "Then add InterProScan to your PATH:"
# echo "export PATH=\$PATH:$(pwd)"
# echo "You may also need to set INTERPROSCAN_HOME=$(pwd)"

# cd ../..

# # install biopython for blast
# echo "Installing Biopython for BLAST support..."
# pip install biopython

# echo "Biopython installation completed."

# # Install BLAST from bioconda
# echo "Installing BLAST from bioconda..."
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda install -c bioconda blast=2.16.0 -y

# mkdir -p blast_db
# cd blast_db || exit 1

# echo "Downloading UniProt SwissProt database..."
# wget --quiet --show-progress -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# if [ -f "uniprot_sprot.fasta.gz" ]; then
#     echo "Decompressing database..."
#     # gunzip -k uniprot_sprot.fasta.gz
#     gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta
    
#     if [ -f "uniprot_sprot.fasta" ]; then
#         echo "Creating BLAST database..."
#         makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_swissprot -parse_seqids -title "UniProt SwissProt"
        
#         # Verify database creation
#         if [ -f "uniprot_swissprot.phr" ]; then
#             echo "BLAST database created successfully."
#             echo "You can now use it with: blastp -db uniprot_swissprot -query your_file.fasta"
#         else
#             echo "Error: BLAST database files not created!" >&2
#             exit 1
#         fi
#     else
#         echo "Error: Failed to decompress database!" >&2
#         exit 1
#     fi
# else
#     echo "Error: Failed to download database!" >&2
#     exit 1
# fi

# # Get the absolute path for the BLAST database
# BLASTDB_PATH=$(pwd)
# export BLASTDB=$BLASTDB_PATH
# echo "BLASTDB environment variable set for current session: $BLASTDB"

# # Function to add BLASTDB to shell config file
# setup_blastdb_env() {
#     local shell_config_file=$1
#     if [ -f "$shell_config_file" ]; then
#         if grep -q "export BLASTDB=" "$shell_config_file"; then
#             echo "BLASTDB variable is already set in $shell_config_file."
#         else
#             echo "Adding BLASTDB to $shell_config_file..."
#             echo "" >> "$shell_config_file"
#             echo "# Set BLASTDB for local BLAST searches" >> "$shell_config_file"
#             echo "export BLASTDB=\"$BLASTDB_PATH\"" >> "$shell_config_file"
#             echo "Successfully added BLASTDB to $shell_config_file."
#             echo "Please run 'source $shell_config_file' or restart your terminal to apply changes."
#         fi
#     else
#         echo "Could not find $shell_config_file. Please add the following line to your shell's startup file:"
#         echo "export BLASTDB=\"$BLASTDB_PATH\""
#     fi
# }

# # Detect shell and update config file
# if [ -n "$BASH_VERSION" ]; then
#     setup_blastdb_env "$HOME/.bashrc"
# elif [ -n "$ZSH_VERSION" ]; then
#     setup_blastdb_env "$HOME/.zshrc"
# else
#     echo "Could not automatically detect your shell's config file."
#     echo "Please add the following line to your shell's startup file (e.g., .bash_profile, .profile):"
#     echo "export BLASTDB=\"$BLASTDB_PATH\""
# fi

# # install python packages
# echo "Installing required Python packages..."
# pip install -i https://pypi.tuna.tsinghua.edu.cn/simple \
#     fastmcp \
#     biopython \
#     numpy \
#     pandas \
#     jinja2 \
#     openai \
#     requests \
#     gradio-client \
#     jsonlines \
#     tqdm \
#     gradio \
#     lmdb

# Install and configure Foldseek
echo "Installing Foldseek from bioconda..."
# First install updated libstdcxx-ng to avoid library version issues
conda install -c conda-forge libstdcxx-ng -y
conda install -c conda-forge -c bioconda foldseek -y

if [ $? -ne 0 ]; then
    echo "Error: Failed to install Foldseek"
    exit 1
fi

echo "Foldseek installation completed."

# Create Foldseek database directory
mkdir -p foldseek_db
cd foldseek_db || exit 1

echo "Downloading AlphaFold Swiss-Prot database for Foldseek..."
echo "This may take a while depending on your network connection..."

# Check if clash proxy is available
if curl -s --connect-timeout 2 http://127.0.0.1:7890 > /dev/null 2>&1; then
    echo "Detected clash proxy, using it for download..."
    export http_proxy=http://127.0.0.1:7890
    export https_proxy=http://127.0.0.1:7890
    export HTTP_PROXY=http://127.0.0.1:7890
    export HTTPS_PROXY=http://127.0.0.1:7890
fi

# Use LD_PRELOAD to ensure foldseek uses the updated libstdc++ from conda
# This fixes the GLIBCXX_3.4.29 and CXXABI_1.3.13 version issues
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6.0.34"
echo "Using preloaded library: $LD_PRELOAD"

foldseek databases Alphafold/Swiss-Prot sp tmp
unset LD_PRELOAD

if [ $? -ne 0 ]; then
    echo "Error: Failed to download Foldseek database!"
    echo "You can try downloading manually later by running:"
    echo "foldseek databases Alphafold/Swiss-Prot foldseek_db/sp tmp"
    cd ..
else
    echo "Foldseek database downloaded successfully."
    cd ..
    
    # Get the absolute path for the Foldseek database
    FOLDSEEK_DB_PATH=$(pwd)/foldseek_db/sp
    export FOLDSEEK_DB=$FOLDSEEK_DB_PATH
    echo "FOLDSEEK_DB environment variable set for current session: $FOLDSEEK_DB"
    
    # Function to add FOLDSEEK_DB to shell config file
    setup_foldseek_env() {
        local shell_config_file=$1
        if [ -f "$shell_config_file" ]; then
            if grep -q "export FOLDSEEK_DB=" "$shell_config_file"; then
                echo "FOLDSEEK_DB variable is already set in $shell_config_file."
            else
                echo "Adding FOLDSEEK_DB to $shell_config_file..."
                echo "" >> "$shell_config_file"
                echo "# Set FOLDSEEK_DB for Foldseek searches" >> "$shell_config_file"
                echo "export FOLDSEEK_DB=\"$FOLDSEEK_DB_PATH\"" >> "$shell_config_file"
                echo "Successfully added FOLDSEEK_DB to $shell_config_file."
                echo "Please run 'source $shell_config_file' or restart your terminal to apply changes."
            fi
        else
            echo "Could not find $shell_config_file. Please add the following line to your shell's startup file:"
            echo "export FOLDSEEK_DB=\"$FOLDSEEK_DB_PATH\""
        fi
    }
    
    # Detect shell and update config file
    if [ -n "$BASH_VERSION" ]; then
        setup_foldseek_env "$HOME/.bashrc"
    elif [ -n "$ZSH_VERSION" ]; then
        setup_foldseek_env "$HOME/.zshrc"
    else
        echo "Could not automatically detect your shell's config file."
        echo "Please add the following line to your shell's startup file (e.g., .bash_profile, .profile):"
        echo "export FOLDSEEK_DB=\"$FOLDSEEK_DB_PATH\""
    fi
fi

echo ""
echo "===================================="
echo "Setup Complete!"
echo "===================================="
echo "Summary:"
echo "- InterProScan installed in conda environment '${CONDA_ENV_NAME}'"
echo "- BLAST database: $(pwd)/blast_db"
echo "- Foldseek database: $(pwd)/foldseek_db/sp"
echo ""
echo "To activate the environment:"
echo "conda activate ${CONDA_ENV_NAME}"
echo ""
echo "Environment variables set:"
echo "export BLASTDB=\"$(pwd)/blast_db\""
echo "export FOLDSEEK_DB=\"$(pwd)/foldseek_db/sp\""