#!/bin/bash

# Exit on error
set -e

# ------------ Step 1: Create directories if not exist ------------
echo "[INFO] Creating folders if needed..."
mkdir -p reads
mkdir -p ref_genome

# ------------ Step 2: Download ERR024407 if not present ------------
cd reads
if [ ! -f ERR024407_1.fastq.gz ]; then
    echo "[INFO] Downloading ERR024407_1.fastq.gz..."
    wget -O ERR024407_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR024/ERR024407/ERR024407_1.fastq.gz
else
    echo "[SKIP] ERR024407_1.fastq.gz already exists."
fi

# ------------ Step 3: Convert to full FASTA ------------
if [ ! -f ecoli_subset.fasta ]; then
    echo "[INFO] Converting FASTQ to full FASTA..."
    zcat ERR024407_1.fastq.gz | awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' > ecoli_subset.fasta
else
    echo "[SKIP] ecoli_subset.fasta already exists."
fi

# ------------ Step 4: Generate smaller subsets ------------
if [ ! -f ecoli_1000.fasta ]; then
    echo "[INFO] Creating 1000-read subset..."
    head -n 2000 ecoli_subset.fasta > ecoli_1000.fasta
else
    echo "[SKIP] ecoli_1000.fasta already exists."
fi

if [ ! -f ecoli_5000.fasta ]; then
    echo "[INFO] Creating 5000-read subset..."
    head -n 10000 ecoli_subset.fasta > ecoli_5000.fasta
else
    echo "[SKIP] ecoli_5000.fasta already exists."
fi

cd ..

# ------------ Step 5: Download reference genome ------------
cd ref_genome
if [ ! -f ecoli_k12_reference.fasta ]; then
    if [ ! -f GCF_000005845.2_ASM584v2_genomic.fna.gz ]; then
        echo "[INFO] Downloading reference genome..."
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
    else
        echo "[SKIP] GCF_000005845.2_ASM584v2_genomic.fna.gz already exists."
    fi

    echo "[INFO] Extracting and renaming reference genome..."
    gunzip -f GCF_000005845.2_ASM584v2_genomic.fna.gz
    mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_k12_reference.fasta
else
    echo "[SKIP] ecoli_k12_reference.fasta already exists."
fi
cd ..

echo "[DONE] All required files are in 'reads/' and 'ref_genome/'."
