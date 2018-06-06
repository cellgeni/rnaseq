#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         B U L K - R N A S E Q    P I P E L I N E
========================================================================================
 Cellular Genetics bulk-RNA-Seq analysis pipeline, Wellcome Sanger Institute
 #### Homepage / Documentation
 https://github.com/cellgeni/RNAseq
 #### Authors
 Vladimir Kiselev @wikiselev <vk6@sanger.ac.uk>
 Original development by SciLifeLabs
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     Bulk-RNA-Seq pipeline v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run cellgeni/RNAseq -profile farm3

    Mandatory arguments:
      -profile                      Hardware config to use. farm3 / farm4 / openstack / docker / aws

    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behaviour

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --star_index                  Path to STAR index
      --star_overhang               sjdbOverhang parameter for building a STAR index (has to be (read_length - 1))
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF file
      --bed12                       Path to bed12 file
      --downloadFasta               If no STAR / Fasta reference is supplied, a URL can be supplied to download a Fasta file at the start of the pipeline.
      --downloadGTF                 If no GTF reference is supplied, a URL can be supplied to download a Fasta file at the start of the pipeline.
      --saveReference               Save the generated reference files the the Results directory.
      --saveAlignedIntermediates    Save the BAM files from the Aligment step  - not done by default

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


/*
 * utils is shared between projects. Include it in the PATH so scripts are found.
 */
// env.PATH = "$baseDir/utils:$PATH"


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '1.5'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.scratch = false
params.name = false
params.project = false
params.genome = 'GRCh38'
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false
params.star_overhang = '74'
params.dna = params.genome ? params.genomes[ params.genome ].dna ?: false : false
params.cdna = params.genome ? params.genomes[ params.genome ].cdna ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.splicesites = false
params.download_hisat2index = false
params.download_fasta = false
params.download_gtf = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.outdir = './results'
params.email = false
params.plaintext_email = false


// Choose aligner
params.aligner = 'star'
if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'salmon'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'salmon'"
}

// Validate inputs
if( params.star_index && params.aligner == 'star' ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}
if ( params.hisat2_index && params.aligner == 'hisat2' ){
    hs2_indices = Channel
        .fromPath("${params.hisat2_index}*")
        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat2_index}" }
}
if ( params.salmon_index && params.aligner == 'salmon' ){
    salmon_index = Channel
        .fromPath(params.salmon_index)
        .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
}
if ( params.salmon_trans_gene && params.aligner == 'salmon' ){
    salmon_trans_gene = Channel
        .fromPath(params.salmon_trans_gene)
        .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_trans_gene}" }
}
if ( params.dna ){
    f = file(params.dna)
    if( !f.exists() ) exit 1, "Fasta file not found: ${params.dna}"
}
if ( params.cdna ){
    f = file(params.cdna)
    if( !f.exists() ) exit 1, "Fasta file not found: ${params.cdna}"
}
if ( ( params.aligner == 'hisat2' && !params.download_hisat2index ) && !params.download_fasta ){
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeBED12;
              gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
}
else if ( !params.download_gtf ){
    exit 1, "No GTF annotation specified!"
}
if( params.bed12 ){
    bed12 = Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "BED12 annotation file not found: ${params.bed12}" }
        .into {bed_rseqc; bed_genebody_coverage}
}
if( params.aligner == 'hisat2' && params.splicesites ){
    Channel
        .fromPath(params.bed12)
        .ifEmpty { exit 1, "HISAT2 splice sites file not found: $alignment_splicesites" }
        .into { indexing_splicesites; alignment_splicesites }
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Header log info
log.info "========================================="
log.info "         RNASeq pipeline v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Data Type']    = 'Paired-End'
summary['Genome']       = params.genome
if(params.aligner == 'star'){
    summary['Aligner'] = "STAR"
    if(params.star_index)          summary['STAR Index']   = params.star_index
    else if(params.dna)          summary['Fasta Ref']    = params.dna
    else if(params.download_fasta) summary['Fasta URL']    = params.download_fasta
}
if(params.aligner == 'salmon'){
    summary['Aligner'] = "Salmon"
    summary['Salmon Index']   = params.salmon_index
    if(params.salmon_index)          summary['Salmon Index']   = params.salmon_index
    else if(params.dna)          summary['Fasta Ref']    = params.dna
    else if(params.download_fasta) summary['Fasta URL']    = params.download_fasta
} 
if(params.aligner == 'hisat2') {
    summary['Aligner'] = "HISAT2"
    if(params.hisat2_index)        summary['HISAT2 Index'] = params.hisat2_index
    else if(params.download_hisat2index) summary['HISAT2 Index'] = params.download_hisat2index
    else if(params.dna)          summary['Fasta Ref']    = params.dna
    else if(params.download_fasta) summary['Fasta URL']    = params.download_fasta
    if(params.splicesites)         summary['Splice Sites'] = params.splicesites
}
if(params.gtf)                 summary['GTF Annotation']  = params.gtf
else if(params.download_gtf)   summary['GTF URL']         = params.download_gtf
if(params.bed12)               summary['BED Annotation']  = params.bed12
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * PREPROCESSING - Build STAR index
 */
if(params.aligner == 'star' && !params.star_index){
    process makeSTARindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from Channel.fromPath(params.dna)
        file gtf from Channel.fromPath(params.gtf)

        output:
        file "star" into star_index

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
}

/*
 * PREPROCESSING - Build Salmon index
 */
if(params.aligner == 'salmon' && !params.salmon_index){
    process makeSalmonIndex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from Channel.fromPath(params.cdna)

        output:
        file "salmon" into salmon_index

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
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from Channel.fromPath(params.cdna)

        output:
        file "trans_gene.txt" into salmon_trans_gene

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

/*
 * PREPROCESSING - Build HISAT2 splice sites file
 */
if(params.aligner == 'hisat2' && !params.splicesites){
    process makeHisatSplicesites {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeHisatSplicesites

        output:
        file "${gtf.baseName}.hisat2_splice_sites.txt" into indexing_splicesites, alignment_splicesites

        script:
        """
        hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
        """
    }
}
/*
 * PREPROCESSING - Build HISAT2 index
 */
if(params.aligner == 'hisat2' && !params.hisat2_index && !params.download_hisat2index && fasta){
    process makeHISATindex {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta
        file indexing_splicesites from indexing_splicesites
        file gtf from gtf_makeHISATindex

        output:
        file "${fasta.baseName}.*.ht2" into hs2_indices

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
/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: // This script is bundled with the pipeline, in RNAseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}


