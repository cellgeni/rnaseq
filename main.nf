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

    nextflow run cellgeni/RNAseq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile farm3

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
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

    Trimming options
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Presets:
      --pico                        Sets trimming and standedness settings for the SMARTer Stranded Total RNA-Seq Kit - Pico Input kit. Equivalent to: --forward_stranded --clip_r1 3 --three_prime_clip_r2 3

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --sampleLevel                 Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}


/*
 * utils is shared between projects. Include it in the PATH so scripts are found.
 */
env.PATH = "$baseDir/utils:$PATH"


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
params.reads = 'cram/*.cram'
params.outdir = './results'
params.email = false
params.plaintext_email = false
params.mdsplot_header = "$baseDir/assets/mdsplot_header.txt"
params.heatmap_header = "$baseDir/assets/heatmap_header.txt"
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"

mdsplot_header = file(params.mdsplot_header)
heatmap_header = file(params.heatmap_header)
biotypes_header = file(params.biotypes_header)
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
params.sampleLevel = false

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
params.pico = false
if (params.pico){
  clip_r1 = 3
  clip_r2 = 0
  three_prime_clip_r1 = 0
  three_prime_clip_r2 = 3
  forward_stranded = true
  reverse_stranded = false
  unstranded = false
}

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
summary['Reads']        = params.reads
summary['Data Type']    = 'Paired-End'
summary['Genome']       = params.genome
if( params.pico ) summary['Library Prep'] = "SMARTer Stranded Total RNA-Seq Kit - Pico Input"
summary['Strandedness'] = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
summary['Trim R1'] = clip_r1
summary['Trim R2'] = clip_r2
summary["Trim 3' R1"] = three_prime_clip_r1
summary["Trim 3' R2"] = three_prime_clip_r2
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

/*
 * Create a channel for input sample ids
 */
sample_list = Channel.fromPath('samples.txt')

process irods {
    tag "${sample}"

    maxForks 29 

    input: 
        val sample from sample_list.flatMap{ it.readLines() }
    output: 
        set val(sample), file('*.cram') optional true into cram_files
    script:
    """
    irods.sh ${sample}
    """
}

process crams_to_fastq {
    tag "${sample}"

    scratch true

    input: 
        set val(sample), file(crams) from cram_files
    output: 
        file "${sample}.fastq.gz" into fastqs
    script:

    def avail_mem = task.memory == null ? '' : "${ ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    """
    samtools merge -f ${sample}.cram ${crams}
    # check that the size of the cram file is >0.5Mb
    minimumsize=500000
    actualsize=\$(wc -c <"${sample}.cram")
    if [ \$actualsize -ge \$minimumsize ]; then
        samtools sort \\
            -n \\
            -@ ${task.cpus} \\
            -m ${avail_mem} \\
            ${cram} | \\
            samtools fastq \\
                -N \\
                -@ ${task.cpus} \\
                -1 ${sample}_1.fastq.gz \\
                -2 ${sample}_2.fastq.gz \\
                -
    fi
    """
}

/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}
if(params.aligner == 'star'){
    hisat_stdout = Channel.from(false)
    salmon_stdout = Channel.from(false)
    process star {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }

        input:
        file reads from fastqs
        file index from star_index.collect()
        file gtf from gtf_star.collect()

        output:
        set file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix
        """
    }
    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { logs, bams -> check_log(logs) }
        .flatMap {  logs, bams -> bams }
    .into { bam_count; bam_rseqc; bam_preseq; bam_markduplicates; bam_featurecounts; bam_stringtieFPKM; bam_geneBodyCoverage }
}

if(params.aligner == 'salmon'){
    hisat_stdout = Channel.from(false)
    star_log = Channel.from(false)
    process salmon {
        tag "$prefix"
        publishDir "${params.outdir}/Salmon", mode: 'copy'

        input:
        file reads from fastqs
        file index from salmon_index.collect()
        file trans_gene from salmon_trans_gene.collect()

        output:
        file "${prefix}.quant.sf" into salmon_trans
        file "${prefix}.quant.genes.sf" into salmon_genes

        script:
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        salmon quant \\
            -i $index \\
            -l ISR \\
            -p ${task.cpus} \\
            --seqBias \\
            --gcBias \\
            --posBias \\
            -q \\
            -o . \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -g ${trans_gene} \\
            --useVBOpt \\
            --numBootstraps 100
        mv quant.sf ${prefix}.quant.sf
        mv quant.genes.sf ${prefix}.quant.genes.sf
        """
    }
}

/*
 * STEP 3 - align with HISAT2
 */
if(params.aligner == 'hisat2'){
    star_log = Channel.from(false)
    salmon_stdout = Channel.from(false)
    process hisat2Align {
        tag "$prefix"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                else params.saveAlignedIntermediates ? filename : null
            }

        input:
        file reads from fastqs
        file hs2_indices from hs2_indices.collect()
        file alignment_splicesites from alignment_splicesites.collect()

        output:
        file "${prefix}.bam" into hisat2_bam
        file "${prefix}.hisat2_summary.txt" into alignment_logs
        file '.command.log' into hisat_stdout

        script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/
        prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
        def rnastrandness = ''
        if (forward_stranded && !unstranded){
            rnastrandness = '--rna-strandness FR'
        } else if (reverse_stranded && !unstranded){
            rnastrandness = '--rna-strandness RF'
        }
        """
        hisat2 -x $index_base \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                $rnastrandness \\
                --known-splicesite-infile $alignment_splicesites \\
                --no-mixed \\
                --no-discordant \\
                -p ${task.cpus} \\
                --met-stderr \\
                --new-summary \\
                --summary-file ${prefix}.hisat2_summary.txt \\
                | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
        hisat2 --version
        """
    }

    process hisat2_sortOutput {
        tag "${hisat2_bam.baseName}"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename -> params.saveAlignedIntermediates ? "aligned_sorted/$filename" : null }

        input:
        file hisat2_bam

        output:
        file "${hisat2_bam.baseName}.sorted.bam" into bam_count, bam_rseqc, bam_preseq, bam_markduplicates, bam_featurecounts, bam_stringtieFPKM, bam_geneBodyCoverage

        script:
        def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"
        """
        samtools sort \\
            $hisat2_bam \\
            -@ ${task.cpus} $avail_mem \\
            -o ${hisat2_bam.baseName}.sorted.bam
        """
    }
}

/*
 * STEP 8 Feature counts
 */

if(params.aligner != 'salmon'){
    process featureCounts {
        tag "${bam_featurecounts.baseName - '.sorted'}"
        publishDir "${params.outdir}/featureCounts", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_biotype_counts_mqc.txt") > 0) "biotype_counts/$filename"
                else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
                else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
                else "$filename"
            }

        input:
        file bam_featurecounts
        file gtf from gtf_featureCounts.collect()
        file biotypes_header

        output:
        file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
        file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
        file "${bam_featurecounts.baseName}_biotype_counts_mqc.txt" into featureCounts_biotype
        file '.command.log' into featurecounts_stdout

        script:
        def featureCounts_direction = 0
        if (forward_stranded && !unstranded) {
            featureCounts_direction = 1
        } else if (reverse_stranded && !unstranded){
            featureCounts_direction = 2
        }
        """
        featureCounts -a $gtf -g gene_id -o ${bam_featurecounts.baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
        featureCounts -a $gtf -g gene_biotype -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
        cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n 7 > tmp_file
        cat $biotypes_header tmp_file >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt
        """
    }

/*
 * STEP 9 - Merge featurecounts
 */
    process merge_featureCounts {
        tag "${input_files[0].baseName - '.sorted'}"
        publishDir "${params.outdir}/featureCounts", mode: 'copy'

        input:
        file input_files from featureCounts_to_merge.collect()

        output:
        file 'merged_gene_counts.txt'

        script:
        """
        merge_featurecounts.py -o merged_gene_counts.txt -i $input_files
        """
    }
}

if(params.aligner == 'salmon'){
    
    process mergeSalmonCounts {
        tag "${input_trans[0].baseName - '.quant.sf'}"
        publishDir "${params.outdir}/mergedCounts", mode: 'copy'

        input:
        file input_trans from salmon_trans.collect()
        file input_genes from salmon_genes.collect()

        output:
        file 'merged_*'

        script:
        """
        merge_salmon.R $input_trans trans
        merge_salmon.R $input_genes genes
        """
    }

}
