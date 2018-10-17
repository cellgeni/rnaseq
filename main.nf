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
 Stijn van Dongen <svd@sanger.ac.uk>
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

    Other options:
      --outdir                      The output directory where the results will be saved
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

params.samplefile = false
params.studyid = -1
params.fastqdir = false
params.outdir = './results'
params.fcextra = ""                          // feature counts extra parameters; currently for testing
params.singleend = false

// Configurable variables
params.scratch = false
params.runtag  = "NF"                        // use runtag as primary tag identifying the run; e.g. studyid
params.name = false
params.project = false
params.genome = 'GRCh38'
params.forward_stranded = false
params.reverse_stranded = false
params.skip_qc = false
params.mito_name = 'MT'
params.skip_multiqc = false
params.skip_fastqc = false
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
params.splicesites = false
params.download_hisat2index = false
params.download_fasta = false
params.download_gtf = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"

biotypes_header = file(params.biotypes_header)
output_docs = file("$baseDir/docs/output.md")

gene_biotype = params.gtf.matches(".*gencode.*") ? "gene_type" : "gene_biotype"

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded


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
    if( !f.exists() ) exit 1, "cdna file not found: ${params.cdna}"
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
summary['Sample file']  = params.samplefile
summary['Data Type']    = 'Paired-End'
summary['Genome']       = params.genome
summary['Biotype tag']  = gene_biotype
summary['Strandedness'] = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
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
 * Create a channel for input sample ids
 */
sample_list = Channel.fromPath(params.samplefile)

if (params.studyid > 0) {
    ch_fastqs_dir = Channel.empty()
    process irods {
        tag "${samplename}"

        input: 
            val samplename from sample_list.flatMap{ it.readLines() }
        output: 
            set val(samplename), file('*.cram') optional true into ch_cram_files
        script:
        """
        bash -euo pipefail irods.sh -t ${params.studyid} -s ${samplename}
        """
    }
} else if (params.fastqdir) {
    ch_cram_files = Channel.empty()
    if (params.singleend) {
      process get_fastq_files_single {
          tag "$samplename"

          input:
              val samplename from sample_list.flatMap{ it.readLines() }
          output:
              set val(samplename), file("${samplename}.fastq.gz") optional true into ch_fastqs_dir
          script:
          """
          name=${params.fastqdir}/${samplename}.fastq.gz
          if [[ ! -e \$name ]]; then
            echo "Count file \$name not found"
            false
          else
            ln -s \$name .
          fi
          """
      }
    }
    else {
      process get_fastq_files {
          tag "${samplename}"

          input:
              val samplename from sample_list.flatMap{ it.readLines() }
          output:
              set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_dir
          script:
          """
          list=( \$(ls ${params.fastqdir}/${samplename}_{1,2}.fastq.gz) )
          if [[ 2 == \${#list[@]} ]]; then
            ln -s \${list[0]} .
            ln -s \${list[1]} .
          else
            echo "Count mismatch sample ${samplename} found (\${list[@]})"
            false
          fi
          """
      }
    }
} else {
  exit 1, "Need --fastqdir <dirname> or --studyid <ID> option"
}


process crams_to_fastq {
    tag "${samplename}"

    if (params.scratch) {
       scratch true
    }

    input: 
        set val(samplename), file(crams) from ch_cram_files
    output: 
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_cram
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    """
    samtools merge -@ ${task.cpus} -f ${samplename}.cram ${crams}

    # check that the size of the cram file is >0.5Mb
    minimumsize=500000
    actualsize=\$(wc -c <"${samplename}.cram")

    f1=${samplename}_1.fastq.gz
    f2=${samplename}_2.fastq.gz

    if [ \$actualsize -ge \$minimumsize ]; then
                              # -O {stdout} -u {no compression}
                              # -N {always append /1 and /2 to the read name}
                              # -F 0x900 (bit 1, 8, filter secondary and supplementary reads)
      samtools collate    \\
          -O -u           \\
          -@ ${task.cpus} \\
          ${samplename}.cram pfx-${samplename} | \\
      samtools fastq      \\
          -N              \\
          -F 0x900        \\
          -@ ${task.cpus} \\
          -1 \$f1 -2 \$f2 \\
          -
    fi
    """
}


ch_fastqs_cram
  .mix(ch_fastqs_dir)
  .into{ ch_reads; ch_fastqc }


process fastqc {
    tag "$samplename"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(samplename), file(reads) from ch_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -q $reads
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
        tag "$samplename"
        publishDir "${params.outdir}", mode: 'copy',
            saveAs: { filename ->
                if (filename ==~ /.*\.ReadsPerGene\.out\.tab/) "STARcounts/$filename"
                else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
                else null
            }

        input:
        set val(samplename), file(reads) from ch_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()

        output:
        set val(samplename), file("*Log.final.out"), file ('*.bam') into star_aligned
        file "*.SJ.out.tab"
        file "*.out" into ch_alignment_logs
        file "*.ReadsPerGene.out.tab"

        script:
                  // TODO featurecounts resorts the BAM file; SortedByName is not a STAR option though.
                  // --outSAMunmapped Within: In case someone wants the BAM files.
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads --readFilesCommand zcat \\
            --runThreadN ${task.cpus} \\
            --twopassMode Basic \\
            --outWigType bedGraph \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --runDirPerm All_RWX \\
            --quantMode GeneCounts \\
            --outFileNamePrefix ${samplename}.
        """
    }
    // Filter removes all 'aligned' channels that fail the check
    star_aligned
        .filter { name, logs, bams -> check_log(logs) }
        .map    { name, logs, bams -> [name, bams] }
    .into { ch_featurecounts; ch_indexbam }
}

if(params.aligner == 'salmon'){
    hisat_stdout = Channel.from(false)
    ch_alignment_logs = Channel.from(false)
    process salmon {
        tag "$samplename"
        publishDir "${params.outdir}/Salmon", mode: 'copy'

        input:
        set val(samplename), file(reads) from ch_reads
        file index from salmon_index.collect()
        file trans_gene from salmon_trans_gene.collect()

        output:
        file "${prefix}.quant.sf" into salmon_trans
        file "${prefix}.quant.genes.sf" into salmon_genes

        script:
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
        mv quant.sf ${samplename}.quant.sf
        mv quant.genes.sf ${samplename}.quant.genes.sf
        """

        // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
        // Include the row names so merger can check identity.
        // The merge step will concatenate the rows and re-transpose to obtain final result.
    }
}

/*
 * STEP 3 - align with HISAT2
 */
if(params.aligner == 'hisat2'){
    salmon_stdout = Channel.from(false)
    process hisat2Align {
        tag "$samplename"
        publishDir "${params.outdir}/HISAT2", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/$filename"
                else null
            }

        input:
        set val(samplename), file(reads) from ch_reads
        file hs2_indices from hs2_indices.collect()
        file alignment_splicesites from alignment_splicesites.collect()

        output:
        file "${samplename}.bam" into hisat2_bam
        file "${samplename}.hisat2_summary.txt" into ch_alignment_logs
        file '.command.log' into hisat_stdout

        script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/
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
                --summary-file ${samplename}.hisat2_summary.txt \\
                | samtools view -bS -F 4 -F 8 -F 256 - > ${samplename}.bam
        hisat2 --version
        """
    }

    process hisat2_sortOutput {
        tag "${hisat2_bam.baseName}"

        input:
        file hisat2_bam

        output:
        set val($samplename), file("${hisat2_bam.baseName}.sorted.bam") into ch_featurecounts, ch_indexbam

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


if(params.aligner != 'salmon') {
    process featureCounts {
        tag "${samplename}"
        publishDir "${params.outdir}/featureCounts", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".biotype_counts_mqc.txt") > 0) "biotype_counts/$filename"
                else if (filename.indexOf(".gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
                else if (filename.indexOf(".gene.featureCounts.txt") > 0) "gene_counts/$filename"
                else "$filename"
            }

        input:
        set val(samplename), file(thebam) from ch_featurecounts
        file gtf from gtf_featureCounts.collect()
        file biotypes_header

        output:
        file "*.gene.featureCounts.txt" into featureCounts_to_merge
        file "*.gene.featureCounts.txt.summary"
        file "*.biotype_counts*mqc.{txt,tsv}" into ch_fc_biotype      // captures _counts_gs_mqc.tsv as well

        script:
        def extraparams = params.fcextra.toString() - ~/^dummy/
        def featureCounts_direction = 0
        def pairedend = params.singleend ? "" : "-p"
        if (forward_stranded && !unstranded) {
            featureCounts_direction = 1
        } else if (reverse_stranded && !unstranded){
            featureCounts_direction = 2
        }
        """
        featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
          -o ${samplename}.gene.featureCounts.txt $pairedend      \\
          -s $featureCounts_direction ${extraparams} $thebam

        featureCounts -T ${task.cpus} -a $gtf -g ${gene_biotype}  \\
          -o ${samplename}.biotype.featureCounts.txt $pairedend   \\
          -s $featureCounts_direction ${extraparams} $thebam

        cut -f 1,7 ${samplename}.biotype.featureCounts.txt | tail -n 7 > tmp_file
        cat $biotypes_header tmp_file >> ${samplename}.biotype_counts_mqc.txt
        # mqc_features_stat.py ${samplename}.biotype_counts_mqc.txt -s $samplename -f rRNA -o ${samplename}.biotype_counts_gs_mqc.tsv
        # TODO this is pbb a multiqc script, activate when ready.
        """
    }

    process indexbam {
        tag "${samplename}"
        publishDir "${params.outdir}/mapsummary", mode: 'copy'

        input:
        set val(samplename), file(thebam) from ch_indexbam

        output:
        set val(samplename), file("*.idxstats") into ch_mapsummary

        """
        samtools index $thebam
        samtools idxstats $thebam > ${samplename}.idxstats
        """
    }
    
    process mapsummary {
        tag "${samplename}"
        publishDir "${params.outdir}/mapsummary", mode: 'copy'

        input:
        set val(samplename), file(thestats) from ch_mapsummary

        output:
        file "*_mqc.txt" into ch_mapsummary_qc

        script:
        def mito_name = params.mito_name
        """
        python $baseDir/bin/mito.py -m ${mito_name} -t $thestats > ${samplename}_mqc.txt
        """
    }

    process merge_featureCounts {
          // TODO: ideally we pass the samplename in the channel. Not sure how to do this given below channel.collect().
        tag "$outputname"
        publishDir "${params.outdir}/featureCounts", mode: 'copy'

        input:
        file input_files from featureCounts_to_merge.collect()

        output:
        file '*-fc-genecounts.txt'

        script:
        def outputname = "${params.runtag}-fc-genecounts.txt"
        """
        python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
          --rm-suffix .gene.featureCounts.txt                             \\
          -c -1 --skip-comments --header                                  \\
          -o $outputname -i $input_files
        """
    }
}

if(params.aligner == 'salmon'){
    
    process mergeSalmonCounts {
        tag "${input_trans[0].baseName - '_1.quant.sf'}"
        publishDir "${params.outdir}/mergedCounts", mode: 'copy'

        input:
        file input_trans from salmon_trans.collect()
        file input_genes from salmon_genes.collect()

        output:
        file '*counts.txt'

        script:
        """
        python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
          --rm-suffix _1.quant.genes.sf                                   \\
          -c -1 --skip-comments --header                                  \\
          -o ${params.runtag}-salmon-genecounts.txt -i $input_genes
        python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
          --rm-suffix _1.quant.genes.sf                                   \\
          -c -1 --skip-comments --header                                  \\
          -o ${params.runtag}-salmon-transcounts.txt -i $input_trans
        """
    }

}


process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('mapsummary/*') from ch_mapsummary_qc.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from ch_fc_biotype.collect()
    file ('alignment/*') from ch_alignment_logs.collect()
/*
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from genebody_coverage_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect()
    file ('stringtie/stringtie_log*') from stringtie_log.collect()
    file ('sample_correlation_results/*') from sample_correlation_results.collect().ifEmpty([]) // If the Edge-R is not run create an Empty
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()
*/
    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    # multiqc . -f $rtitle $rfilename \\
    #     -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m fastqc
    multiqc . -f $rtitle $rfilename \\
        -m custom_content -m featureCounts -m star -m fastqc
    """
}


