#!/bin/bash

# ==============================================================================
# SLURM Directives: Confirmed
# ==============================================================================
#SBATCH --job-name=emu_disney_water_analysis 
#SBATCH --account=tgoldwire                 
#SBATCH --qos=duttonc                       
#SBATCH --partition=hpg-turin               
#SBATCH --cpus-per-task=8                   
#SBATCH --mem=32G                           
#SBATCH --time=24:00:00                     
#SBATCH --output=/blue/duttonc/tgoldwire/emu/disney_water/logs/emu_water_%j.log 
#SBATCH --mail-type=END,FAIL                
#SBATCH --mail-user=tgoldwire@ufl.edu       

# ==============================================================================
# Setup Conda Environment (Keeping as-is)
# ==============================================================================
module load conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate/create conda environment and dependencies
if ! conda info --envs | grep -q "^emu_env[[:space:]]"; then
    echo "Creating emu_env..." >&2
    conda create -y -n emu_env python=3.9
fi
conda activate emu_env

# Check and install Emu if needed
if ! command -v emu &> /dev/null; then
    echo "Installing emu and dependencies..." >&2
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install -y emu
    pip install osfclient
fi

# ==============================================================================
# Setup Emu Database (Keeping as-is)
# ==============================================================================
EMU_DB_DIR="/blue/duttonc/duttonc/databases/emu_db"
if [ ! -d "$EMU_DB_DIR" ] || [ ! -f "$EMU_DB_DIR/species_taxid.fasta" ]; then
    echo "Downloading Emu database from OSF..." >&2
    mkdir -p "$EMU_DB_DIR"
    export EMU_DATABASE_DIR="$EMU_DB_DIR"
    cd "$EMU_DB_DIR"
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
    tar -xvf emu.tar
    rm emu.tar
    cd -
fi

# ==============================================================================
# Define Input and Output Paths: Confirmed
# ==============================================================================
INPUT_ROOT="/blue/duttonc/DisneyWater/pass"
OUTPUT_ROOT="/blue/duttonc/tgoldwire/emu/disney_water"
# 🌟 CRITICAL FIX: Define the fixed thread count based on SBATCH header
EMU_THREADS=8 

# Setup directories
mkdir -p "$OUTPUT_ROOT/logs" 
mkdir -p "$OUTPUT_ROOT"

# Define the standard output file pattern for Emu counts
EXPECTED_OUTFILE_SUFFIX="_counts.tsv"

# ==============================================================================
# Main Emu Analysis Loop: Barcode 1-24 and 'unclassified'
# ==============================================================================
DIRS_TO_PROCESS=()
for i in $(seq 1 24); do 
    # Use printf to zero-pad single-digit barcodes (e.g., barcode01)
    DIRS_TO_PROCESS+=("barcode$(printf "%02d" "$i")")
done
DIRS_TO_PROCESS+=("unclassified")


for SAMPLE_DIR_NAME in "${DIRS_TO_PROCESS[@]}"; do
    SAMPLE_DIR="$INPUT_ROOT/$SAMPLE_DIR_NAME"
    OUTDIR="$OUTPUT_ROOT/$SAMPLE_DIR_NAME"
    
    # Define Merged File Location/Name (Critical for emu combine-outputs)
    MERGED_FASTQ_GZ="$OUTDIR/${SAMPLE_DIR_NAME}.fastq.gz" 
    
    # The file emu combine-outputs is looking for:
    OUTFILE="$OUTDIR/${SAMPLE_DIR_NAME}.fastq${EXPECTED_OUTFILE_SUFFIX}" 

    echo "Checking $SAMPLE_DIR_NAME, OUTFILE=$OUTFILE..."

    # Check and skip if already processed
    if [ -f "$OUTFILE" ]; then
        echo "Sample $SAMPLE_DIR_NAME already processed at $OUTFILE. Skipping."
        continue
    fi
    
    # Check if sample folder exists
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "Sample folder $SAMPLE_DIR not found, skipping."
        continue
    fi

    mkdir -p "$OUTDIR"

    # --- File Merging Logic ---
    
    # Check for gzipped files
    if ls "$SAMPLE_DIR"/*.fastq.gz 1> /dev/null 2>&1; then
        echo "Merging gzipped FASTQ files for $SAMPLE_DIR_NAME..."
        # Use zcat for gzipped files, piping to a new gzip stream (faster than cat)
        zcat "$SAMPLE_DIR"/*.fastq.gz | gzip > "$MERGED_FASTQ_GZ"
    # Check for uncompressed files
    elif ls "$SAMPLE_DIR"/*.fastq 1> /dev/null 2>&1; then
        echo "Merging uncompressed FASTQ files for $SAMPLE_DIR_NAME and compressing..."
        # Use cat for uncompressed files, then pipe to gzip
        cat "$SAMPLE_DIR"/*.fastq | gzip > "$MERGED_FASTQ_GZ"
    else
        echo "No FASTQ or FASTQ.GZ files found for $SAMPLE_DIR_NAME, skipping."
        continue
    fi
    
    # --- Run Emu abundance ---
    if [ -f "$MERGED_FASTQ_GZ" ]; then
        echo "Running emu abundance for $SAMPLE_DIR_NAME..."
        # 🌟 CRITICAL FIX: Using the fixed EMU_THREADS variable instead of SLURM_CPUS_PER_TASK
        emu abundance "$MERGED_FASTQ_GZ" \
            --db "$EMU_DB_DIR" \
            --threads "$EMU_THREADS" \
            --output-dir "$OUTDIR" \
            --keep-counts
            
        # Check if emu abundance ran successfully
        if [ $? -ne 0 ]; then
            echo "ERROR: emu abundance failed for $SAMPLE_DIR_NAME."
            # Don't delete merged file if Emu failed, to allow for inspection
            continue
        fi
    
        # CLEANUP: Delete the large merged FASTQ file
        rm -f "$MERGED_FASTQ_GZ"
    else
        echo "Error: Merged FASTQ file not created for $SAMPLE_DIR_NAME."
    fi
done

# ==============================================================================
# Combine Outputs 
# ==============================================================================
echo "Combining all Emu outputs..."
emu combine-outputs "$OUTPUT_ROOT" species --counts

# Also generate CSV for spreadsheet use
TSV_COMBINED="$OUTPUT_ROOT/emu-combined-species.tsv"
CSV_COMBINED="$OUTPUT_ROOT/emu-combined-species.csv"
if [ -f "$TSV_COMBINED" ]; then
    awk 'BEGIN{FS=OFS="\t"}{print}' "$TSV_COMBINED" | sed 's/\t/,/g' > "$CSV_COMBINED"
    echo "Combined CSV file: $CSV_COMBINED"
else
    echo "ERROR: Combined TSV file was not generated by emu combine-outputs."
fi

echo "SUCCESS: All samples processed."
