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

if (params.aligner == 'hisat2' || params.aligner == 'salmon') {
    exit 1, "Not yet supported: hisat2 and salmon"
}


if (!params.dna) { exit 1, "Please set a fasta reference file in the config file" }
if (!params.gtf) { exit 1, "Please set a gtf reference file in the config file"   }

Channel.fromPath(params.dna).ifEmpty { exit 1, "dna fasta file not found" } .set { dna_fa_channel }
Channel.fromPath(params.gtf).ifEmpty { exit 1, "gtf annot file not found" } .into { gtf_channel; gtf_makeBED12 }


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
summary['GTF file']     = params.gtf
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


process makeSTARindex {
    tag "$fasta"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file fasta from dna_fa_channel
    file gtf from gtf_channel

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


if(params.aligner == 'salmon') {

    process makeSalmonIndex {
        tag "$fasta"
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file fasta from Channel.fromPath(params.cdna)

        output:
        file "salmon"

        script:
        """
        mkdir salmon
        salmon index \\
            -t $fasta \\
            -i salmon
        """
    }

    process makeTransGeneMatrix {
        tag "$fasta"
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file fasta from Channel.fromPath(params.cdna)

        output:
        file "trans_gene.txt"

        script:
        """
        grep '>' $fasta \\
        | awk '{ print \$1, \$4 }' \\
        | sed 's/>//' \\
        | sed 's/gene://' \\
            > trans_gene.txt
        """
    }
}



if (params.aligner == 'hisat2') {

    process makeHisatSplicesites {
        tag "$gtf"
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file gtf from gtf_makeHisatSplicesites

        output:
        file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

        script:
        """
        hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
        """
    }

    process makeHISATindex {
        tag "$fasta"
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file fasta from fasta
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf_makeHISATindex

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
}



process makeBED12 {
    tag "$gtf"
    publishDir "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file gtf from gtf_makeBED12

    output:
    file "${gtf.baseName}.bed"

    script:
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}


