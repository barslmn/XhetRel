#!/usr/bin/env nextflow

params.input_dir = './'
params.output_dir = './'

 // Works for all three assemblies;start position is from GRCh38, end position is from T2T.
 // We are just taking the smallest matching region all three assemblies in order
 // to escape the PAR regions in any of them.
 // https://en.wikipedia.org/wiki/Pseudoautosomal_region
params.nonPARregion="chrX:2781479-153925834"

// Parametric filtering parameters (for Xhet analysis)
params.vaf_threshold      = params.vaf_threshold ?: 0.25
params.dp_threshold       = params.dp_threshold ?: 20
params.gq_threshold       = params.gq_threshold ?: 0
params.apply_pass_filter  = params.apply_pass_filter ?: true

nextflow.enable.dsl=2

process prepareVCFs {
    container 'docker.io/staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${name}.preprocessed.vcf.gz"), path("${name}.preprocessed.vcf.gz.csi")

    shell:
    '''
    #!/bin/bash
    set -e
    set -o pipefail

    # Define the final output name for clarity
    FINAL_VCF="!{name}.preprocessed.vcf.gz"
    INPUT_VCF="!{vcf}"

    bcftools view "$INPUT_VCF" -Oz -o "$FINAL_VCF"
    bcftools index "$FINAL_VCF"
    # --- 2. Unconditionally Add Custom VCF Headers ---
    echo "  - Ensuring custom headers are present in $FINAL_VCF..."
    OLD_HEADER="old_header"
    NEW_HEADER="new_header"
    TEMP_VCF="temp_vcf"
    # Define the header lines to be added
    HEADER_LINES='##FILTER=<ID=likely_sequence_context_artefAct,Description="likey sequence context artefact">
##INFO=<ID=SBS,Number=.,Type=Float,Description="SBS">
##INFO=<ID=SBP,Number=.,Type=Float,Description="SBP">
##INFO=<ID=SBF,Number=.,Type=Float,Description="SBF">
##INFO=<ID=SBT,Number=.,Type=Float,Description="SBT">
##INFO=<ID=SBR,Number=.,Type=Float,Description="SBR">
##INFO=<ID=DPHQ,Number=.,Type=Float,Description="DPHQ">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=DETECTOR,Number=.,Type=Float,Description="DETECTOR">'

    CHR_MAP="chr_map"
    cat > "$CHR_MAP" <<EOL
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
EOL

    # Extract the header from the input VCF
    bcftools view -h "$FINAL_VCF" > "$OLD_HEADER"

    # Construct the new header by inserting the custom lines before the final header row
    echo "$(sed \\$d "$OLD_HEADER"; echo "$HEADER_LINES"; sed -n \\$p "$OLD_HEADER")" > "$NEW_HEADER"

    # Reheader the VCF file in-place
    bcftools reheader -h "$NEW_HEADER" "$FINAL_VCF" | bcftools annotate --rename-chrs "$CHR_MAP" | bcftools view -Oz -o "$TEMP_VCF" && mv "$TEMP_VCF" "$FINAL_VCF"

    # Clean up temporary header files
    rm -f "$OLD_HEADER" "$NEW_HEADER" "$CHR_MAP"
    echo "  - Custom header check complete."

    # --- 3. Check if AD annotation transfer is needed (Robust Check) ---
    TRANSFER_NEEDED=false
    # First, ensure INFO/AD exists in the header, otherwise there's nothing to transfer.
    if bcftools view -h "$FINAL_VCF" | grep -q '##INFO=<ID=AD,'; then
        # Case 1: Header correctly says FORMAT/AD is missing. Transfer is needed.
        if ! bcftools view -h "$FINAL_VCF" | grep -q '##FORMAT=<ID=AD,'; then
            echo "Condition met: INFO/AD found and FORMAT/AD not declared in header."
            TRANSFER_NEEDED=true
        else
            # Case 2: Header claims FORMAT/AD exists. We must verify the data records.
            echo "Header claims FORMAT/AD exists. Verifying first data record..."
            FIRST_FORMAT_FIELD=$(set +o pipefail; bcftools view -H "$FINAL_VCF" | head -n 1 | cut -f 9)

            # Transfer if the first record exists but its FORMAT field lacks AD.
            if [[ -n "$FIRST_FORMAT_FIELD" ]] && ! (echo "$FIRST_FORMAT_FIELD" | tr ':' '\n' | grep -qw "AD"); then
                echo "Condition met: Header declares FORMAT/AD, but it is missing from the data records (e.g., '$FIRST_FORMAT_FIELD')."
                TRANSFER_NEEDED=true
            fi
        fi
    fi

    if [[ "$TRANSFER_NEEDED" == "true" ]]; then
        echo "Proceeding with AD annotation transfer..."

        # If the header has a misleading FORMAT/AD line, remove it first to prevent conflicts.
        if bcftools view -h "$FINAL_VCF" | grep -q '##FORMAT=<ID=AD,'; then
            echo "  - Misleading FORMAT/AD line found in header. Removing it before annotation..."
            TEMP_VCF_NO_AD_HEADER="temp_vcf_no_ad"
            MODIFIED_HEADER="modified_header"

            # Create a new header file by excluding the old FORMAT/AD line
            bcftools view -h "$FINAL_VCF" | grep -v '##FORMAT=<ID=AD,' | sed "s/##INFO=<ID=AF,Number=1/##INFO=<ID=AF,Number=A/" > "$MODIFIED_HEADER"

            # Reheader the VCF to remove the line
            bcftools reheader -h "$MODIFIED_HEADER" "$FINAL_VCF" | bcftools view -Oz -o "$TEMP_VCF_NO_AD_HEADER"
            mv "$TEMP_VCF_NO_AD_HEADER" "$FINAL_VCF"

            rm -f "$MODIFIED_HEADER"
            echo "  - Misleading FORMAT/AD line removed."
        fi

        ANNOT_FILE="annot_file"
        HEADER_FILE="header_file"
        OUTPUT_VCF="temp_vcf"

        # Perform the Annotation Transfer
        echo "  - Extracting INFO/AD to a temporary annotation file..."
        bcftools query -f '%CHROM\t%POS\t%INFO/AD\n' "$FINAL_VCF" | bgzip -c > "${ANNOT_FILE}.gz"

        echo "  - Indexing the temporary annotation file..."
        tabix -s1 -b2 -e2 "${ANNOT_FILE}.gz"

        echo "  - Creating a temporary header for the new FORMAT/AD field..."
        echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">' > "$HEADER_FILE"

        echo "  - Annotating the VCF to add FORMAT/AD..."
        bcftools annotate -a "${ANNOT_FILE}.gz" -h "$HEADER_FILE" -c CHROM,POS,FORMAT/AD -o "$OUTPUT_VCF" -Oz "$FINAL_VCF"

        echo "  - Overwriting VCF with the AD-transferred version..."
        mv "$OUTPUT_VCF" "$FINAL_VCF"

        echo "  - Cleaning up temporary files..."
        rm -f "${ANNOT_FILE}.gz" "${ANNOT_FILE}.gz.tbi" "$HEADER_FILE"

        echo "--- Transfer complete. ---"
    else
        echo "Skipping AD Transfer: VCF file does not require AD annotation transfer."
        echo "--- No changes made for AD transfer. ---"
    fi

    # --- 5. Final Indexing ---
    bcftools index -f "$FINAL_VCF"
    '''
}

process initialstat {
    container 'docker.io/staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    path "${name}.initial.stats"

    shell:
    '''
    bcftools stats !{vcf} > "!{name}.initial.stats"
    '''
}

// Filter for X heterozygosity analysis - uses parametric filtering
process XhetFilter {
    container 'docker.io/staphb/bcftools:1.20'

    input:
    tuple val(name), path(vcf), path(index)

    output:
    tuple val(name), path("${name}.xhet.filtered.vcf.gz"), path("${name}.xhet.filtered.vcf.gz.csi")

    shell:
    def filter_expr = "FORMAT/VAF >= ${params.vaf_threshold} && FORMAT/DP >= ${params.dp_threshold} && FORMAT/GQ >= ${params.gq_threshold}"
    def pass_filter = params.apply_pass_filter ? "-f .,PASS" : ""

    """
    bcftools +fill-tags !{vcf} -- -t FORMAT/VAF | bcftools annotate -x FORMAT/DP | bcftools +fill-tags -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))' |
    bcftools view -i '${filter_expr}' ${pass_filter} -Oz -o !{name}.xhet.filtered.vcf.gz
    bcftools index !{name}.xhet.filtered.vcf.gz
    """
}

process xhetStat {
    container 'docker.io/staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    path "${name}.xhet.filtered.stats"

    shell:
    '''
    bcftools stats !{vcf} > "!{name}.xhet.filtered.stats"
    '''
}

// Filter for relatedness analysis - uses default filters (VAF=0.25, DP=20, PASS)
process relatednessFilter {
    container 'docker.io/staphb/bcftools:1.20'

    input:
    tuple val(name), path(vcf), path(index)

    output:
    tuple val(name), path("${name}.relatedness.filtered.vcf.gz"), path("${name}.relatedness.filtered.vcf.gz.csi")

    shell:
    """
    bcftools +fill-tags !{vcf} -- -t FORMAT/VAF | bcftools annotate -x FORMAT/DP | bcftools +fill-tags -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))' |
    bcftools view -i 'FORMAT/VAF >= 0.25 && FORMAT/DP >= 20' -f .,PASS -Oz -o !{name}.relatedness.filtered.vcf.gz
    bcftools index !{name}.relatedness.filtered.vcf.gz
    """
}

process relatednessStat {
    container 'docker.io/staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    path "${name}.relatedness.filtered.stats"

    shell:
    '''
    bcftools stats !{vcf} > "!{name}.relatedness.filtered.stats"
    '''
}

process mergeVCFs {
    container 'docker.io/staphb/bcftools:1.20'
    input:
    path(vcfs)

    output:
    tuple path('merged/merged.vcf.gz'), path('merged/merged.vcf.gz.tbi')

    shell:
    '''
    mkdir merged/
    bcftools merge *.vcf.gz -Ou | \
    bcftools +fill-tags -- -t INFO/F_MISSING | \
    bcftools view -i 'F_MISSING<0.2' -Oz -o merged/merged.vcf.gz
    bcftools index -t merged/merged.vcf.gz
    '''
}

process relatedness {
    container 'quay.io/biocontainers/vcftools:0.1.17--pl5321h077b44d_0'
    input:
    tuple path(vcf), path(index)

    output:
    path "merged.relatedness2"

    shell:
    '''
    vcftools --gzvcf !{vcf} --relatedness2 --out merged
    '''
}

process xhetyaml{
    output:
    path "Xhet_mqc.yaml"

    shell:
    '''
    echo 'id: "xhet_section"
    section_name: "X heterozygosity"
    description: "This plot shows the fraction of heterozygous variants to total variants in X chromosome."
    plot_type: "scatter"
    pconfig:
      id: "Xhet_plot"
      title: "X heterozygosity plot"
      xlab: "samples"
      ylab: "Xhet"
    data:' > Xhet_mqc.yaml
    '''
}

process xhet{
    container 'docker.io/staphb/bcftools:1.20'
    input:
    path(xhetyaml)
    tuple val(name), path(vcf), path(index)

    output:
    path "Xhet_mqc.yaml"

    shell:
    '''
    bcftools query "!{vcf}" -r "!{params.nonPARregion}" -f "[%SAMPLE\t%GT]\n" |
    sort |
    uniq -c |
    grep '0/1\\|1/1' |
    awk -v name=!{name} 'NR==1{numhet=$1}
        {sum+=$1;}
        END{printf "  %s: {x: %s, y: %s}\\n", name, name, numhet/sum;}' >> "!{xhetyaml}"
    '''
}

process multiqc {
    container 'docker.io/multiqc/multiqc:v1.23'
    publishDir "${params.output_dir}", mode: 'move'
    input:
    path(multiqc_files)

    output:
    path("*.html")

    shell:
    """
    cat <<EOF > XhetRel_mqc.yaml
XhetRel:
  parametric_filter_vaf_threshold: "${params.vaf_threshold}"
  parametric_filter_dp_threshold: "${params.dp_threshold}"
  parametric_filter_gq_threshold: "${params.gq_threshold}"
  parametric_filter_apply_pass_filter: "${params.apply_pass_filter}"
  relatedness_filter_vaf_threshold: "0.25"
  relatedness_filter_dp_threshold: "20"
  relatedness_filter_apply_pass_filter: "true"
  nonPAR_region: "${params.nonPARregion}"
EOF
    multiqc --flat --force --fullnames --dirs --outdir . . --filename XhetRel_flat_multiqc_report.html
    multiqc --interactive --force --fullnames --dirs --outdir . . --filename XhetRel_interactive_multiqc_report.html
    """
}

workflow {
    vcf_files = Channel.fromPath("${params.input_dir}/*{.vcf,.vcf.gz,.bcf}")
    | map { vcf ->
        [
            vcf.baseName.split('.vcf')[0],
            vcf
        ]
    }

    // Prepare VCFs
    vcf_files = prepareVCFs(vcf_files)
    
    // Get initial stats
    initial_stats = initialstat(vcf_files)
    
    // Split into two branches for parallel filtering
    // Branch 1: Parametric filtering for Xhet analysis
    xhet_filtered = XhetFilter(vcf_files)
    xhet_stats = xhetStat(xhet_filtered)
    
    // Branch 2: Default filtering for relatedness analysis
    relatedness_filtered = relatednessFilter(vcf_files)
    relatedness_stats = relatednessStat(relatedness_filtered)
    
    // Merge filtered VCFs for relatedness analysis
    relatedness_filtered
        .map { name, vcf, index -> [vcf, index] }
        .collect()
        .set { relatedness_vcfs_wo_names }
    
    merged_vcf = mergeVCFs(relatedness_vcfs_wo_names)
    relatedness(merged_vcf)
    
    // Xhet analysis on parametric filtered VCFs
    xhetyaml()
    xhet(xhetyaml.out, xhet_filtered)
    
    // Collect all outputs for MultiQC
    Channel.empty()
        .mix(relatedness.out)
        .mix(xhetyaml.out)
        .mix(initial_stats)
        .mix(xhet_stats)
        .mix(relatedness_stats)
        .collect()
        .set { multiqc_files }

    multiqc(multiqc_files)
}
