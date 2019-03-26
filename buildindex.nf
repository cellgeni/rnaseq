#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                               build genome indexes
========================================================================================
 #### Authors
 Vladimir Kiselev @wikiselev <vk6@sanger.ac.uk>
 Stijn van Dongen <svd@sanger.ac.uk>
 Original development by SciLifeLabs
----------------------------------------------------------------------------------------
*/



params.star_overhang = '74'
params.dna = params.genome ? params.genomes[ params.genome ].dna ?: false : false
params.cdna = params.genome ? params.genomes[ params.genome ].cdna ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.splicesites = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.outdir = './results'
params.name = false


params.aligner = 'star'
if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'salmon'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'salmon'"
}

if (params.aligner == 'hisat2') {
    // exit 1, "Not yet supported: hisat2"
}


[ ["dna", params.dna], ["gtf", params.gtf], ["cdna", params.cdna] ].each {
  if (!it[1]) {
    exit 1, "Please set a ${it[0]} reference file in the config file"
  }
}

Channel.fromPath(params.dna).ifEmpty { exit 1, "dna fasta file not found" } .into { ch_dna_star; ch_dna_hisat2 }
Channel.fromPath(params.gtf).ifEmpty { exit 1, "gtf annot file not found" } .into { ch_gtf_star; ch_gtf_bed; ch_gtf_hisat2_splice; ch_gtf_hisat2_index }
Channel.fromPath(params.cdna).ifEmpty { exit 1, "cdna annot file not found" } .into { ch_cdna_salmon; ch_cdna_transgene }


//  Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Data Type']    = 'Paired-End'
summary['Genome']       = params.genome
summary['DNA file']     = params.dna
summary['CDNA file']    = params.cdna
summary['GTF file']     = params.gtf
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


process makeSTARindex {
    tag "$fasta"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    when:
    params.aligner == 'star'

    input:
    file fasta from ch_dna_star
    file gtf from ch_gtf_star

    output:
    file "star"

    script:
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang ${params.star_overhang} \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta
    """
}

// process makeBED12 {
//     tag "$gtf"
//     publishDir "${params.outdir}/reference_genome", mode: 'copy'

//     input:
//     file gtf from ch_gtf_bed

//     output:
//     file "${gtf.baseName}.bed"

//     script:
//     """
//     gtf2bed $gtf > ${gtf.baseName}.bed
//     """
// }


process makeSalmonIndex {
    tag "$fasta"
    publishDir "${params.outdir}/reference_genome", mode: 'link'

    when:
    params.aligner == 'salmon'

    input:
    file fasta from ch_cdna_salmon

    output:
    file "salmon"

    script:
    """
    mkdir salmon
    salmon index        \\
        -t $fasta       \\
        -p ${task.cpus} \\
        -i salmon
    """
}

process makeTransGeneMatrix {
    tag "$fasta"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    when:
    params.aligner == 'salmon'

    input:
    file fasta from ch_cdna_transgene

    output:
    file "trans_gene*.txt"

    shell:
    '''
    perl -ne 'if (/^>(\\w+)(?:\\.\\d+)\\s+.*?gene:(\\w+)/){print "$1\\t$2\\n"}elsif(/^>(ERCC\\S+)/){print"$1\\t$1-gene\\n"}' \\
      !{fasta} > trans_gene.txt
    if (( $(grep -c ENST trans_gene.txt) < 1000 )); then
       echo "Not enough Ensembl transcripts. This test is present because this script section is ugly.'
       echo 'Currently it makes an effort to recognise ERCC information.'
       echo 'If you want to run a gencode genome, update this section, make this file aware of gencode/Ensembl distinction'
       false
    fi
    '''
}


process hisat2_splicesites {

    label 'hisat2_build'
    tag "$gtf"

    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    when:
    params.aligner == 'hisat2'

    input:
    file gtf from ch_gtf_hisat2_splice

    output:
    file "${gtf.baseName}.hisat2_splice_sites.txt" into ch_hisat2_index_splice

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
    """
}

process hisat2_index {

    label 'hisat2_build'
    tag "$fasta"

    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    when:
    params.aligner == 'hisat2'

    input:
    file fasta from ch_dna_hisat2
    file indexing_splicesites from ch_hisat2_index_splice
    file gtf from ch_gtf_hisat2_index

    output:
    file "${fasta.baseName}.*.ht2"

    script:
    if( task.memory == null ){
        log.info "[HISAT2 index build] Available memory not known - defaulting to 0. Specify process memory requirements to change this."
        avail_mem = 0
    } else {
        log.info "[HISAT2 index build] Available memory: ${task.memory}"
        avail_mem = task.memory.toGiga()
    }
    if( avail_mem > params.hisatBuildMemory ){
        log.info "[HISAT2 index build] Over ${params.hisatBuildMemory} GB available, so using splice sites and exons in HISAT2 index"
        extract_exons = "hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt"
        ss = "--ss $indexing_splicesites"
        exon = "--exon ${gtf.baseName}.hisat2_exons.txt"
    } else {
        log.info "[HISAT2 index build] Less than ${params.hisatBuildMemory} GB available, so NOT using splice sites and exons in HISAT2 index."
        log.info "[HISAT2 index build] Use --hisatBuildMemory [small number] to skip this check."
        extract_exons = ''
        ss = ''
        exon = ''
    }
    """
    $extract_exons
    hisat2-build -p ${task.cpus} $ss $exon $fasta ${fasta.baseName}.hisat2_index
    """
}


