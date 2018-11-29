#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         B U L K - R N A S E Q    P I P E L I N E
========================================================================================
 Cellular Genetics bulk-RNA-Seq analysis pipeline, Wellcome Sanger Institute
 Documentation:   https://github.com/cellgeni/rnaseq
 Authors:
    Stijn van Dongen <svd@sanger.ac.uk>
    Vladimir Kiselev @wikiselev <vk6@sanger.ac.uk>
    Original development by SciLifeLabs
----------------------------------------------------------------------------------------
*/

params.run_star     = true
params.run_qc       = true
params.run_multiqc  = true
params.run_fastqc   = true
params.run_rnaseq   = true     // feature counts; featureCounts with one of STAR, hisat2, salmon
params.run_mixcr    = false
params.run_hisat2   = true
params.run_salmon   = true
params.save_bam     = false

params.outdir = 'results'
params.runtag = "cgirnaseq"    // use runtag as primary tag identifying the run; e.g. studyid


def helpMessage() {
    log.info"""
    =========================================
     Bulk-RNA-Seq pipeline v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run cellgeni/RNAseq --studyid 01234 --samplefile sample-study01234.txt --genome GRCh38 -profile farm3

    Mandatory arguments:
      -profile                      Hardware config to use. farm3 / farm4 / docker
      --genome                      Name of iGenomes reference
      --samplefile                  File with sample IDs. If input is from --fastqdir,
                                    files names will be expected as <sampleid>.fastq.gz

    Important options:
      --studyid                     Unless input is from --fastqdir, specify a studyID.
                                    Input will be queried in IRODS using ids from samplefile.
      --runtag                      Will be included e.g. in count matrix names. Suggested: studyid
      --outdir        [$params.outdir]

    Modes:
      --outdir                      The output directory where the results will be saved
      --runtag                      Tag for outputs
      --run_star      [$params.run_star]
      --run_hisat2    [$params.run_hisat2]
      --run_salmon    [$params.run_salmon]
      --run_mixcr     [$params.run_mixcr]
      --run_fastqc    [$params.run_fastqc]
      --run_multiqc   [$params.run_multiqc]
      --save_bam      [$params.save_bam]
      --run_qc        [$params.run_qc]      Set to false to avoid fastqc, multiqc
      --run_rnaseq    [$params.run_rnaseq]      Set to false avoid star, hisat2, salmon, featureCounts

    Strandedness:
      --forward_stranded            The library is forward stranded
      --reverse_stranded            The library is reverse stranded
      --unstranded                  The default behaviour
    """.stripIndent()
}


version = '1.7'

params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.samplefile = false
params.studyid = -1
params.fastqdir = false
params.fcextra = ""                          // feature counts extra parameters; currently for testing
params.singleend = false


params.name = false
params.genome = 'GRCh38'
params.forward_stranded = false
params.reverse_stranded = false

params.mito_name = 'MT'
params.unstranded = false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false
params.star_overhang = '74'
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.download_hisat2index = false
params.hisatBuildMemory = 200 // Required amount of memory in GB to build HISAT2 index with splice sites
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"

biotypes_header = file(params.biotypes_header)
output_docs = file("$baseDir/docs/output.md")

gene_biotype = params.gtf.matches(".*gencode.*") ? "gene_type" : "gene_biotype"

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded



if (params.studyid < 0 && !params.fastqdir) {
  exit 1, "Need --fastqdir <dirname> or --studyid <ID> option"
}

params.aligner = 'star'
if (params.aligner != 'star' && params.aligner != 'hisat2' && params.aligner != 'salmon'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'salmon'"
}

if (params.run_hisat2 && !params.hisat2_index) {
    exit 1, "No hisat2 index"
}

ch_star_index = params.run_star
    ? Channel.fromPath(params.star_index)
      .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    : Channel.empty()

ch_salmon_index = params.run_salmon
    ? Channel.fromPath(params.salmon_index)
       .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
    : Channel.empty()

ch_salmon_trans_gene = params.run_salmon
    ?  Channel.fromPath(params.salmon_trans_gene)
       .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_trans_gene}" }
    : Channel.empty()

ch_hs2_indices = params.run_hisat2
    ? Channel.fromPath("${params.hisat2_index}/*.?.ht2")
      .ifEmpty { exit 1, "HISAT2 index files not found from $params.hisat2_index" }
    : Channel.empty()

ch_hisat2_splicesites = params.run_hisat2
  ? Channel.fromPath("${params.hisat2_index}/*hisat2_splice_sites.txt")
    .ifEmpty { exit 1, "HISAT2 splice sites file not found in $params.hisat2_index" }
  : Channel.empty()


Channel.fromPath(params.gtf)
  .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  .into { ch_gtf_star; ch_gtf_featureCounts; }


// Header log info
log.info "========================================="
log.info "         RNASeq pipeline v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']           = workflow.runName
summary['Sample file']        = params.samplefile
summary['Data Type']          = 'Paired-End'
summary['Genome']             = params.genome
summary['Biotype tag']        = gene_biotype
summary['Strandedness']       = ( unstranded ? 'None' : forward_stranded ? 'Forward' : reverse_stranded ? 'Reverse' : 'None' )
summary['STAR Index']         = params.star_index
summary['Salmon Index']       = params.salmon_index
summary['HISAT2 Index']       = params.hisat2_index
summary['GTF Annotation']     = params.gtf
summary['BED Annotation']     = params.bed12
summary['Output dir']         = params.outdir
summary['Working dir']        = workflow.workDir
summary['Container']          = workflow.container
summary['Current home']       = "$HOME"
summary['Current path']       = "$PWD"
summary['Script dir']         = workflow.projectDir
summary['Config Profile']     = workflow.profile
if(workflow.revision)
  summary['Pipeline Release'] = workflow.revision

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


nf_required_version = '0.30.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error """
====================================================
  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.
  Pipeline execution will continue, but things may break.
  Please run `nextflow self-update` to update Nextflow.
============================================================
"""
}



/*
 * The channel for input sample ids. This can either refer to file names in a directory with
 * fastq files, or to Sanger sample IDs in IRODS.
*/

Channel.fromPath(params.samplefile)
.into { sample_list_irods; sample_list_dirpe; sample_list_dirse }


process irods {
    tag "${samplename}"

    when:
      params.studyid > 0

    input: 
        val samplename from sample_list_irods.flatMap{ it.readLines() }

    output: 
        set val(samplename), file('*.cram') optional true into ch_cram_files
        file('*.lostcause.txt') optional true into ch_lostcause_irods

    script:
    """
    if bash -euo pipefail irods.sh -t ${params.studyid} -s ${samplename}; then
      true
    else
      stat=\$?
      tag='UNKNOWN'
      if [[ \$stat == 64 ]]; then tag='nofiles'; fi
      echo -e "${samplename}\\tirods\\t\$tag" > ${samplename}.lostcause.txt
    fi
    """
}


process get_fastq_files_single {
    tag "$samplename"

    when:
    params.fastqdir && params.singleend

    input:
        val samplename from sample_list_dirse.flatMap{ it.readLines() }
    output:
        set val(samplename), file("${samplename}.fastq.gz") optional true into ch_fastqs_dirse
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


process get_fastq_files {
    tag "${samplename}"

    when:
    params.fastqdir && !params.singleend

    input:
        val samplename from sample_list_dirpe.flatMap{ it.readLines() }
    output:
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_dirpe
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


process crams_to_fastq {
    tag "${samplename}"

//  if (params.scratch) {    // This is tricky; need to get job requirements correct to ensure space exists.
//     scratch true          // At the moment we don't use this. Perhaps with a retry regime ... but a lot of fuss to solve.
//  }                        // I've left it as a reminder it's an option (svd).

    input: 
        set val(samplename), file(crams) from ch_cram_files
    output: 
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_irods
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


ch_fastqs_irods
  .mix(ch_fastqs_dirpe, ch_fastqs_dirse)
  .into{ ch_rnaseq; ch_fastqc; ch_mixcr }

ch_rnaseq
  .until{ ! params.run_rnaseq  }
  .into { ch_hisat2; ch_star; ch_salmon }


process fastqc {
    tag "$samplename"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.run_qc && params.run_fastqc

    input:
    set val(samplename), file(reads) from ch_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_multiqc_fastqc

    script:
    """
    fastqc -t ${task.cpus} -q $reads
    """
}

process mixcr {
    tag {samplename}
    publishDir "${params.outdir}/mixcr/", mode: 'copy',
      saveAs: { filename -> "${samplename.md5()[-1..-2]}/$samplename/$filename" }

    when:
    params.run_mixcr

    input:
    set val(samplename), file(reads) from ch_mixcr

    output:
    file("*full_clones.txt")
    file("*.clones.clna")
    file("*.vdjca")

    script:
    """
    mixcr align --species hsa -t ${task.cpus} $reads ${samplename}.alignments.vdjca
    mixcr assemble -t ${task.cpus} ${samplename}.alignments.vdjca ${samplename}.clones.clna
    mixcr exportClones ${samplename}.clones.clna ${samplename}.full_clones.txt
    """
}



// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false

n_star_lowmapping = 0

def star_filter(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        n_star_lowmapping++
        return false
    } else {
        return true
    }
}


// Currently this prefers star if both star and hisat2 are run, otherwise takes hisat2
// This is for processes that we only want to run for one aligner, not both;
// For example when publishing bams, or pushing bams to multiqc.

def pick_aligner(String aligner) {
    return  aligner == 'star' || (!params.run_star && aligner == 'hisat2')
      ? true
      : false
}


process star {
    tag "$samplename"
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename ==~ /.*\.ReadsPerGene\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }

    when:
    params.run_star

    input:
    set val(samplename), file(reads) from ch_star
    file index from ch_star_index.collect()
    file gtf from ch_gtf_star.collect()

    output:
    set val(samplename), file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.ReadsPerGene.out.tab" into ch_merge_starcounts
    file "*.out" into ch_alignment_logs_star
    file "*.SJ.out.tab"

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
  ch_star_accept = Channel.create()
  ch_star_reject = Channel.create()

  star_aligned
      .choice(ch_star_accept, ch_star_reject)
          { namelogsbams -> star_filter(namelogsbams[1]) ? 0 : 1 }

  ch_star_accept
  .map    { name, logs, bams -> ["star", name, bams] }
  .into   { ch_fc_star; ch_bam_star }

  ch_star_reject
  .map    { it -> "${it[0]}\tSTAR\tlowmapping\n" }
  .mix(ch_lostcause_irods)
  .set    { ch_lostcause }


process salmon {
    tag "$samplename"
    publishDir "${params.outdir}/Salmon", mode: 'copy'

    when:
    params.run_salmon

    input:
    set val(samplename), file(reads) from ch_salmon
    file index from ch_salmon_index.collect()
    file trans_gene from ch_salmon_trans_gene.collect()

    output:
    file "${samplename}.quant.sf" into ch_salmon_trans
    file "${samplename}.quant.genes.sf" into ch_salmon_genes

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


process hisat2Align {

    tag "$samplename"

    publishDir "${params.outdir}/HISAT2", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "logs/hisat2/$filename"
            else null
        }

    when:
    params.run_hisat2

    input:
    set val(samplename), file(reads) from ch_hisat2
    file indices from ch_hs2_indices.collect()
    file alignment_splicesites from ch_hisat2_splicesites.collect()

    output:
    set val(samplename), file("${samplename}.bam") into ch_hisat2_bam
    file "${samplename}.hisat2_summary.txt" into ch_alignment_logs_hisat2
    file '.command.log' into hisat_stdout

    script:
    def index_base = indices[0].toString() - ~/.\d.ht2/
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


// TODO; any reason not to merge this with the process above?
process hisat2_sortOutput {

    tag "$samplename"

    input:
    set val(samplename), file(hisat2_bam) from ch_hisat2_bam

    output:
    set val("hisat2"), val(samplename), file("${samplename}.hs2.sorted.bam") into ch_fc_hisat2, ch_bam_hisat2

    script:
    def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"

    """
    samtools sort \\
        $hisat2_bam \\
        -@ ${task.cpus} $avail_mem \\
        -o ${samplename}.hs2.sorted.bam
    """
}


ch_fc_hisat2
  .mix(ch_fc_star)
  .set{ ch_featurecounts }

ch_bam_hisat2
  .mix(ch_bam_star)
  .into{ ch_indexbam; ch_publishbam }


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
    set val(aligner), val(samplename), file(thebam) from ch_featurecounts
    file gtf from ch_gtf_featureCounts.collect()
    file biotypes_header

    output:
    set val(aligner), file("*.gene.fc.txt") into ch_merge_fc
    set val(aligner), file("*.gene.fc.txt.summary") into ch_multiqc_fc
    set val(aligner), file("*.biotype_counts*mqc.txt") into ch_multiqc_fcbiotype

    script:
    def extraparams = params.fcextra.toString() - ~/^dummy/
    def fc_direction = 0
    def tag = "${samplename}.${aligner}"

    def pairedend = params.singleend ? "" : "-p"
    if (forward_stranded && !unstranded) {
        fc_direction = 1
    } else if (reverse_stranded && !unstranded){
        fc_direction = 2
    }
    """
    featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
      -o ${tag}.gene.fc.txt $pairedend                        \\
      -s $fc_direction ${extraparams} $thebam

    featureCounts -T ${task.cpus} -a $gtf -g ${gene_biotype}  \\
      -o ${tag}.biotype.fc.txt $pairedend                     \\
      -s $fc_direction ${extraparams} $thebam

    cut -f 1,7 ${tag}.biotype.fc.txt |                        \\
        tail -n +3 | cat $biotypes_header - >> ${tag}.biotype_counts_mqc.txt

    # Below works, but will require a new NF process for python.
    # rRNA: ribosomal RNA.
    # python3 $baseDir/bin/mqc_features_stat.py ${samplename}.biotype_counts_mqc.txt -s $samplename -f rRNA -o ${samplename}.biotype_counts_gs_mqc.tsv
    """
}


// Note: we lose information about the used aligner currently.
process indexbam {
    tag "${samplename}"

    when:
    pick_aligner(aligner)

    input:
    set val(aligner), val(samplename), file(thebam) from ch_indexbam

    output:
    set val(samplename), file("*.idxstats") into ch_mapsummary

    script:
    """
    samtools index $thebam
    samtools idxstats $thebam > ${samplename}.idxstats
    """
}


ch_publishbam
  .until{ !params.save_bam }
  .subscribe {
      aligner     = it[0]
      samplename  = it[1]
      thebam      = it[2]
      bamname     = thebam.toString()
      if (pick_aligner(aligner)) {
        dir = "${params.outdir}/${aligner}-bams"
        thebam.copyTo("$dir/${samplename.md5()[0..1]}/${samplename}-${aligner}.bam")
      }
  }


process mapsummary {
    tag "${samplename}"
    publishDir "${params.outdir}/mapsummary", mode: 'copy'

    input:
    set val(samplename), file(thestats) from ch_mapsummary

    output:
    file "*_mqc.txt" into ch_multiqc_mapsum

    script:
    def mito_name = params.mito_name
    """
    python $baseDir/bin/mito.py -m ${mito_name} -t $thestats > ${samplename}_mqc.txt
    """
}


// https://www.biostars.org/p/218995/
// column 1: gene ID
// column 2: counts for unstranded RNA-seq
// column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
// column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

process merge_starcounts {

    tag "${input_files[0]}"
    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file input_files from ch_merge_starcounts.collect()

    output:
    file '*-star-genecounts.txt'

    script:
    def outputname = "${params.runtag}-star-genecounts.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      -c 1 --rm-suffix .ReadsPerGene.out.tab                          \\
      -o $outputname -i $input_files
    """
}


ch_merge_fc
  .transpose()
  .groupTuple()
  .collectFile { id, files -> [ id, files.collect{ it.toString() }.join('\n') + '\n' ] }
  .set{ ch_merge_fc_byaligner }


process merge_featureCounts {
    tag "$aligner"
    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file metafile from ch_merge_fc_byaligner

    output:
    file '*-fc-genecounts.txt'

    script:
    aligner = metafile.baseName   // not strictly necessary
    outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .gene.featureCounts.txt                             \\
      -c -1 --skip-comments --header                                  \\
      -o $outputname -i \$(cat $metafile)
    """
}


process merge_salmoncounts {
    tag "${input_trans[0]}"
    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file input_trans from ch_salmon_trans.collect()
    file input_genes from ch_salmon_genes.collect()

    output:
    file '*counts.txt'

    script:
    def outtransname = "${params.runtag}-salmon-transcounts.txt"
    def outgenesname = "${params.runtag}-salmon-genecounts.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix _1.quant.genes.sf                                   \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenesname -i $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix _1.quant.sf                                         \\
      -c -1 --skip-comments --header                                  \\
      -o $outtransname -i $input_trans
    """
}

process lostcause {

    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file (inputs) from ch_lostcause.collect().ifEmpty([])

    output:
    file('*.lostcause_mqc.txt') into ch_multiqc_lostcause

    script:
    def outputname = "${params.runtag}.${workflow.runName}.lostcause_mqc.txt"
    """
    echo -e "# plot_type: 'table'\n# section_name: 'Lost samples'" > $outputname
    echo -e "Sample\tProcess\tMessage" >> $outputname
    cat $inputs | sort >> $outputname
    """
}


ch_multiqc_fc
  .filter{ pick_aligner(it[0]) }
  .map { it[1] }
  .set{ ch_multiqc_fc_hisat2 }

ch_multiqc_fcbiotype
  .filter{ pick_aligner(it[0]) }
  .map{ it[1] }
  .set{ ch_multiqc_fcbiotype_hisat2 }

process multiqc {

    publishDir "${params.outdir}", mode: 'link',
      saveAs: {filename ->
          if (filename.indexOf("multiqc.html") > 0) "combined/$filename"
          else if (filename.indexOf("_data") > 0) "$filename"
          else null
      }

    when:
    params.run_qc && params.run_multiqc

    input:
    file ('lostcause/*') from ch_multiqc_lostcause.collect().ifEmpty([])
    file (fastqc:'fastqc/*') from ch_multiqc_fastqc.collect().ifEmpty([])
    file ('mapsummary/*') from ch_multiqc_mapsum.collect().ifEmpty([])
    file ('featureCounts/*') from ch_multiqc_fc_hisat2.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from ch_multiqc_fcbiotype_hisat2.collect().ifEmpty([])
    file ('star/*') from ch_alignment_logs_star.collect().ifEmpty([])
    file ('hisat2/*') from ch_alignment_logs_hisat2.collect().ifEmpty([])

    output:
    file "*_multiqc.html"
    file "*_data"

    script:
    def filename = "${params.runtag}_multiqc.html"
    def reporttitle = "${params.runtag} (cellgeni/rnaseq)"
    """
    multiqc . -f --title "$reporttitle" --filename "$filename" -m custom_content -m featureCounts -m star -m fastqc
    """
}

process workflow_manifest {

    publishDir "${params.outdir}/combined", mode: 'copy'

    output:
    file('*.manifest.txt')

    script:
                      // Manifest's pipeline version: $workflow.manifest.version
                      // ^ option in nextflow config section, perhaps for later.
    def versionfile = "${params.runtag}.manifest.txt"
    """
cat <<EOF > $versionfile
Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
EOF
    """
}


workflow.onComplete {

    summary = [:]
    summary['star low mapping'] = n_star_lowmapping

    log.info "========================================="
    log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
    log.info "========================================="
}


