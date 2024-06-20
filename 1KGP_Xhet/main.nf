#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.source = "RefSeq"
params.assembly = "GRCh38"
params.version = "latest"
params.feature = "exon"
params.annot_url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
params.exons = "${params.source}.${params.assembly}.${params.version}.bed.gz"
params.samples_info = "igsr_samples.tsv"
params.samples_url = 'https://www.internationalgenome.org/api/beta/sample/_search/igsr_samples.tsv'
params.samples_data = 'json=%7B%22fields%22%3A%5B%22name%22%2C%22sex%22%2C%22populations.code%22%5D%7D'
params.base_url = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/raw_calls_old'
params.genotypes = "genotypes.tsv.gz"
params.bedgraph = "genotype_counts.bedGraph"

process createbedfile {
    output:
    path("${params.exons}")

    shell:
    '''
    wget -q -O- "!{params.annot_url}" |
    zcat |
    grep  'NC_000023' |
    awk -F"\t" -v feature="!{params.feature}" '
        /^#!/ {print}
        /^##/ {next}
        $3 ~ feature {
        sub(/^NC_[0]+/, "chr");
        sub(/^chr23/, "chrX");
        split($1, chrom, ".");
        split($9, info, "ID=");
        split(info[2], id, ";");
        printf "%s\t%s\t%s\t%s\n", chrom[1], $4, $5, id[1]}' |
    grep "NM" |
    gzip -c > "!{params.exons}" &&
    echo "INFO: BED file created at !{params.exons}" ||
    echo "ERROR: An error occurred while creating BED file at !{params.exons}"
    '''
}

process getsamplesheet {
    output:
    path("${params.samples_info}")

    shell:
    '''
    curl '!{params.samples_url}' --data-raw '!{params.samples_data}' -o "!{params.samples_info}" &&
    echo "INFO: Samples info file downloaded at !{params.samples_info}" ||
    echo "ERROR: An error occurred while downloading samples info file at !{params.samples_info}"
    # TODO: this is for testing remove afterwards
    sed -i 10q '!{params.samples_info}'
    '''
}

process getfiltervcf {
    memory '1 GB'
    cpus 1
    errorStrategy 'ignore'
    input:
    tuple val(sample), val(sex), val(population_code)
    path exons

    output:
    path "${sample}.bcf{,.csi}"

    shell:
    '''
    url="!{params.base_url}/!{population_code}/Sample_!{sample}/analysis/!{sample}.haplotypeCalls.er.raw.vcf.gz"

    bcftools view "$url" \
        -r chrX:2781479-153925834 \
        -T "!{exons}" \\
        -v snps \
        -f .,PASS \
        -i 'GT=="0/1" && FORMAT/DP>=20' \
        -Ob -o !{sample}.bcf
    bcftools index !{sample}.bcf
    '''
}

process mergebcf {
    input:
    path(bcfs)

    output:
    path 'merged/merged.bcf'

    shell:
    '''
    mkdir merged/
    bcftools merge *.bcf -Ob -o merged/merged.bcf &&
    echo "INFO: Merged BCF file created at merged.bcf" ||
    echo "ERROR: An error occurred while merging BCF files"
    '''
}

process getgenotypes {
    input:
    path 'merged.bcf'

    output:
    path(params.genotypes)

    shell:
    '''
    bcftools query \
      merged.bcf \
      -f '%CHROM\t%POS[\t%GT]\n' | \
      gzip -c > !{params.genotypes} &&
      echo "INFO: Genotypes file created at !{params.genotypes}" ||
      echo "ERROR: An error occurred while extracting genotypes"
    '''
}

process createbedgraph {
    input:
    path(params.genotypes)

    output:
    path(params.bedgraph)

    shell:
    '''
    zcat !{params.genotypes} | awk -F"\t" 'BEGIN {
        print "track type=bedGraph visibility=full color=255,153,0 altColor=0,0,0"
    }
    {
      count = 0
      for (i = 3; i <= NF; i++) {
        if ($i == "0/1") count++
      }
      printf "%s\t%s\t%s\t%s\n" $1, $2, $2, count
    }' > !{params.bedgraph} &&
    echo "INFO: BedGraph file created at !{params.bedgraph}" ||
    echo "ERROR: An error occurred while creating BedGraph file"
    '''
}

process creategenebed {
    output:
    path "genes.bed"

    shell:
    '''
    wget -q -O- "!{params.annot_url}" |
    zcat |
    grep  'NC_000023' |
    awk -F"\t" -v feature="gene" '
        /^#!/ {print}
        /^##/ {next}
        $3 ~ feature {
        sub(/^NC_[0]+/, "chr");
        sub(/^chr23/, "chrX");
        split($1, chrom, ".");
        split($9, info, "ID=");
        split(info[2], id, ";");
        printf "%s\t%s\t%s\t%s\n", chrom[1], $4, $5, id[1]}' |
    gzip -c > "genes.bed.gz" &&
    echo "INFO: BED file created at genes.bed" ||
    echo "ERROR: An error occurred while creating BED file at genes.bed"
    '''
}

process getgenecounts {

    input:
    path(genotype_counts), path(genes_bed)

    output:
    path "genotype_counts_by_gene.txt"

    shell:
    '''
    bedtools intersect \
        -a !{genotype_counts} \
        -b !{genes_bed} -wb |
        awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $8}' |
        datamash -s -g 5 sum 4 |
        sort -k 2,2n > genotype_counts_per_gene_by_sum.txt

    bedtools intersect \
        -a !{genotype_counts} \
        -b !{genes_bed} -wb |
        awk '{printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $8}' |
        datamash -s -g 5 count 4 |
        sort -k 2,2n > genotype_counts_per_gene_by_count.txt

    tail -n 20 genotype_counts_per_gene_by_sum.txt | sort > top_20_genes_by_sum_sorted_by_name.txt
    zcat genes.bed.gz | awk '{printf "%s\t%s\t%s\t%s\n", $4,$1,$2,$3}' | sort -k 1,1 > genes_sorted_by_name.bed
    join top_20_genes_by_sum_sorted_by_name.txt genes_sorted_by_name.bed | awk '{printf "%s\t%s\t%s\t%s\n", $3, $4, $5, $1}'> top_20_genes_by_sum_with_positions.bed

    tail -n 20 genotype_counts_per_gene_by_count.txt | sort > top_20_genes_by_count_sorted_by_name.txt
    zcat genes.bed.gz | awk '{printf "%s\t%s\t%s\t%s\n", $4,$1,$2,$3}' | sort -k 1,1 > genes_sorted_by_name.bed
    join top_20_genes_by_count_sorted_by_name.txt genes_sorted_by_name.bed | awk '{printf "%s\t%s\t%s\t%s\n", $3, $4, $5, $1}'> top_20_genes_by_count_with_positions.bed
    '''
}

workflow {
    getsamplesheet()
    samplesheet = getsamplesheet.out
    samplesheet
        .splitCsv(header:["sample", "sex", "population_code"],
            sep:"\t",
            skip:1
        )
        .filter { row -> row.sex == "male" && row.population_code }
        | set {samples}
    createbedfile()
    bcfs = getfiltervcf(samples, createbedfile.out)
    mergebcf(bcfs.collect())
    getgenotypes(mergebcf.out)
    createbedgraph(getgenotypes.out)
    creategenebed()
    getgenecounts(
        createbedgraph.out
        creategenebed.out
    )
}
