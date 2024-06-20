#!/usr/bin/env nextflow

params.input_dir = './'
params.output_dir = './'

 // Works for all three assemblies;start position is from GRCh38, end position is from T2T.
 // We are just taking the smallest matching region all three assemblies in order
 // to escape the PAR regions in any of them.
 // https://en.wikipedia.org/wiki/Pseudoautosomal_region
params.nonPARregion="chrX:2781479-153925834"

nextflow.enable.dsl=2

process prepareVCFs {
    container 'staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${name}.vcf.gz"), path("${name}.vcf.gz.csi")

    shell:
    '''
    if [[ $(od -t x1 "!{vcf}" | sed 1q | cut -d " " -f2,3) == "1f 8b" ]]; then
      :
    else
      echo "File !{vcf} is plain text."
      bgzip -c !{vcf}> !{name}.vcf.gz;
    fi

    bcftools index !{name}.vcf.gz;
    '''
}

process initialstat {
    container 'staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    path "${name}.stats"

    shell:
    '''
    bcftools stats !{vcf} > "!{name}.stats"
    '''
}
process filteredstat {
    container 'staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    path "${name}.filtered.stats"

    shell:
    '''
    bcftools stats !{vcf} > "!{name}.filtered.stats"
    '''
}

process filterVCFs {
    container 'staphb/bcftools:1.20'
    input:
    tuple val(name), path(vcf), path(index)

    output:
    tuple val(name), path("${name}.filtered.vcf.gz"), path("${name}.filtered.vcf.gz.csi")

    shell:
    '''
    bcftools +fill-tags !{vcf} -- -t FORMAT/VAF |
    bcftools view \
        -i 'FORMAT/VAF >= 0.25 && FORMAT/DP >= 20' \
        -f .,PASS \
        -Oz -o !{name}.filtered.vcf.gz
    bcftools index !{name}.filtered.vcf.gz
    '''
}


process mergeVCFs {
    container 'staphb/bcftools:1.20'
    input:
    path(vcfs)

    output:
    tuple path('merged/merged.vcf.gz'), path('merged/merged.vcf.gz.tbi')

    shell:
    '''
    mkdir merged/
    bcftools merge *.vcf.gz -Ou | \
    bcftools +fill-tags -- -t INFO/F_MISSING | \
    bcftools view -i 'F_MISSING<.8' -Oz -o merged/merged.vcf.gz
    bcftools index -t merged/merged.vcf.gz
    '''
}


process relatedness {
    container 'biocontainers/vcftools:v0.1.16-1-deb_cv1'
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
    container 'staphb/bcftools:1.20'
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
    container 'multiqc/multiqc:dev'
    publishDir "${params.output_dir}", mode: 'move'
    input:
    path(multiqc_files)

    output:
    path("*.html")

    shell:
    '''
    multiqc --flat --force --fullnames --dirs --outdir . . --filename XhetRel_flat_multiqc_report.html
    multiqc --interactive --force --fullnames --dirs --outdir . . --filename XhetRel_interactive_multiqc_report.html
    '''
}

workflow {
    vcf_files = Channel.fromPath("${params.input_dir}/*{.vcf,.vcf.gz}")
    | map { vcf ->
        [
            vcf.baseName.split('.vcf')[0],
            vcf
        ]
    }

    // vcf_files.view()
    vcf_files = prepareVCFs(vcf_files)
    filtered_vcfs = filterVCFs(vcf_files)

    // filtered_vcfs.view()
    filtered_vcfs .map {
        name, vcf, index -> [vcf, index]
        }
        .collect()
        .set {filtered_vcfs_wo_names}

    merged_vcf = mergeVCFs(filtered_vcfs_wo_names)
    relatedness(merged_vcf)
    xhetyaml()
    xhet(xhetyaml.out, filtered_vcfs)
    initial_stats = initialstat(vcf_files)
    filtered_stats = filteredstat(filtered_vcfs)

    Channel.empty()
        .mix(relatedness.out)
        .mix(xhetyaml.out)
        .mix(initial_stats)
        .mix(filtered_stats)
        .collect()
        .set{multiqc_files}

    // multiqc_files.view()
    multiqc(multiqc_files)
}
