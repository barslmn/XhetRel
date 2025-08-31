#!/bin/bash

# This script is designed to run the pipeline in colab environment.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Script Input ---
# The script expects the VCF directory path as the first argument.
vcf_dir=$1

if [ -z "$vcf_dir" ];
    then
        echo "Error: No input directory provided. Please provide the path to the VCF files."
        exit 1
fi

echo "✅ VCF Directory: $vcf_dir"

# --- Tool Installation Functions ---
install_bcftools() {
    if [ ! -f "/usr/local/bin/bcftools" ]; then
        apt-get update -qq 2>/dev/null
        apt-get install -y libgsl27 -qq 2>/dev/null
        wget -q http://barissalman.com/public/xhetrel/bcftools-ubuntu22-precompiled.tar.gz -O /tmp/bcftools.tar.gz
        tar -xzf /tmp/bcftools.tar.gz -C /usr/local/
        ln -sf /usr/local/bcftools/bin/bcftools /usr/local/bin/bcftools
        rm /tmp/bcftools.tar.gz
    fi
}

install_tools() {
    echo "Installing tools..."
    install_bcftools
    # The repo is already cloned by the Python cell, so we just need the other tools.
    curl -L -o nextflow-24.10.8-dist https://github.com/nextflow-io/nextflow/releases/download/v24.10.8/nextflow-24.10.8-dist 2>/dev/null
    chmod +x nextflow-24.10.8-dist
    pip install -q git+https://github.com/MultiQC/MultiQC
    apt-get update -y -qq >/dev/null
    apt-get install -y vcftools -qq >/dev/null
}

check_tools() {
    missing_tools="0"
    if ! command -v vcftools &> /dev/null; then missing_tools="1"; fi
    if ! command -v bcftools &> /dev/null; then missing_tools="1"; fi
    if ! command -v multiqc &> /dev/null; then missing_tools="1"; fi
    if [ ! -f "./nextflow-24.10.8-dist" ]; then missing_tools="1"; fi
    if [ $missing_tools = "1" ]; then
        install_tools
    fi
}

# --- Analysis Function ---
run_analysis() {
    echo "🚀 Starting analysis..."
    NXF_ANSI_LOG=false ./nextflow-24.10.8-dist XhetRel/main.nf --input_dir "$vcf_dir" --output_dir "/content"
    echo "✅ Analysis complete! Report generated."
}

# --- Main Execution ---
start=$(date +%s)

echo "🛠️ Checking for required tools..."
check_tools
echo "✅ Tools are ready."

run_analysis

end=$(date +%s)
runtime=$((end-start))
echo "⏱️ Total runtime: $runtime seconds."
