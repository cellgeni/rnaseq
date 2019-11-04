#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         B U L K - R N A S E Q    P I P E L I N E
========================================================================================
 Cellular Genetics bulk-RNA-Seq analysis pipeline, Wellcome Sanger Institute
 Documentation: https://github.com/cellgeni/rnaseq
 This pipeline was forked from https://github.com/nf-core/rnaseq
 Post-fork development at Wellcome Sanger Institute:
    Stijn van Dongen <svd@sanger.ac.uk>
    Vladimir Kiselev @wikiselev <vk6@sanger.ac.uk>
 Original development by SciLifeLabs
    Phil Ewels @ewels <phil.ewels@scilifelab.se>
    Rickard Hammarén @Hammarn  <rickard.hammaren@scilifelab.se>
----------------------------------------------------------------------------------------
*/

params.run_star     = true
params.run_qc       = true
params.run_multiqc  = true
params.run_fastqc   = true
params.run_rnaseq   = true     // feature counts; featureCounts with one of STAR, hisat2, salmon
params.run_mixcr    = false
params.run_bracer   = false
params.run_tracer   = false
params.run_hisat2   = true
params.run_salmon   = true
params.save_bam     = false
params.min_reads    = 500
params.min_pct_aln  = 5
params.pe_suffix_pattern  = '_{1,2}.fastq.gz'
params.se_suffix    = '.fastq.gz'

params.dropirodsqc  = false
dropqc              = ""

if (params.dropirodsqc) {
  dropqc = '-Q'
}


params.outdir = 'results'
params.runtag = "cgirnaseq"    // use runtag as primary tag identifying the run; e.g. studyid

params.bracer_genometag = 'Hsap'

params.tracer_genometag = 'Hsap'
params.tracer_outputtag = 'summary'
params.tracer_assemble_publish = false

params.fastqc_publish = false

n_numreads = 0

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
                                    files names will be expected as <sampleid>.fastq.gz (single end)
                                    or <sampleid>_1.fastq.gz and <sampleid>_2.fastq.gz

    Important options:
      --studyid                     Unless input is from --fastqdir, specify a studyID.
                                    Input will be queried in IRODS using ids from samplefile.
      --fastqdir                  ! If specified, use the absolute path of the directory.
                                    Input will be fastq files, with the base name specified by
                                    sample IDs from the sample file (see under --samplefile).
      --bamdir                    ! If specified, use the absolute path of the directory.
      --runtag                      Will be included e.g. in count matrix names. Suggested: studyid
      --outdir        [$params.outdir]

    Modes:
      --outdir                      The output directory where the results will be saved
      --runtag                      Tag for outputs
      --run_star      [$params.run_star]
      --run_hisat2    [$params.run_hisat2]
      --run_salmon    [$params.run_salmon]
      --run_mixcr     [$params.run_mixcr]
      --run_bracer    [$params.run_bracer]
      --run_tracer    [$params.run_tracer]
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
params.bamdir = false
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
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"

biotypes_header = file(params.biotypes_header)
output_docs = file("$baseDir/docs/output.md")

gene_biotype = params.gtf.matches(".*gencode.*") ? "gene_type" : "gene_biotype"

forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded



if (params.studyid < 0 && !params.fastqdir && !params.bamdir) {
  exit 1, "Need --fastqdir <dirname> or --studyid <ID> or --bamdir option"
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
summary['Salmon t/g Map']     = params.salmon_trans_gene
summary['HISAT2 Index']       = params.hisat2_index
summary['GTF Annotation']     = params.gtf
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
.into { sample_list_irods; sample_list_fastqpe; sample_list_fastqse; sample_list_bampe }


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
    if bash -euo pipefail irods.sh -N ${task.cpus} -t ${params.studyid} -s ${samplename} ${dropqc}; then
      true
    else
      stat=\$?
      if [[ \$stat == 64 ]];
        then tag='nofiles';
        echo -e "${samplename}\\tirods\\t\$tag" > ${samplename}.lostcause.txt
      else          
        tag='UNKNOWN'
        echo -e "${samplename}\\tirods\\t\$tag" > ${samplename}.lostcause.txt
        exit \$stat
      fi
    fi
    """
}


process get_fastq_files_single {
    tag "$samplename"
    errorStrategy 'terminate'

    when:
    params.fastqdir && params.singleend

    input:
        val samplename from sample_list_fastqse.flatMap{ it.readLines() }
    output:
        set val(samplename), file("${samplename}.fastq.gz") optional true into ch_fastqs_dirse
        file('numreads.txt') optional true into ch_numreads_fastq_se
    script:
    """
    name=${params.fastqdir}/${samplename}{$params.se_suffix}
    if [[ ! -e \$name ]]; then
      echo "Fastq file \$name not found"
      false
    else
      ln -s \$name ${samplename}.fastq.gz
      echo \$(( \$(zcat \$fname | wc -l) / 4)) > numreads.txt
    fi
    """
}


process get_fastq_files_from_bam {

    tag "${samplename}"
    errorStrategy 'terminate'

    when:
    params.bamdir

    input:
        val samplename from sample_list_bampe.flatMap{ it.readLines() }
    output:
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_bams_dirpe
        file('numreads.txt') optional true into ch_numreads_bam

    shell:
    '''
    bam="!{params.bamdir}/!{samplename}.bam"
    f1="!{samplename}_1.fastq.gz"
    f2="!{samplename}_2.fastq.gz"
    if [[ -e $bam  ]]; then
      samtools fastq -N -F 0x900 -@ !{task.cpus} -1 $f1 -2 $f2 $bam
      echo  $(( $(zcat $f1 | wc -l) / 2)) > numreads.txt
    else
      echo "File $bam not found"
      false
    fi
    '''
}


process get_fastq_files {
    tag "${samplename}"
    errorStrategy 'terminate'

    when:
    params.fastqdir && !params.singleend

    input:
        val samplename from sample_list_fastqpe.flatMap{ it.readLines() }
    output:
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_dirpe
        file('numreads.txt') optional true into ch_numreads_fastq
    script:
    """
    list=( \$(ls ${params.fastqdir}/${samplename}${params.pe_suffix_pattern}) )
    if [[ 2 == \${#list[@]} ]]; then
      f1=\${list[0]}
      f2=\${list[1]}
      ln -s \$f1 ${samplename}_1.fastq.gz
      ln -s \$f2 ${samplename}_2.fastq.gz
      echo  \$(( \$(zcat \$f1 | wc -l) / 2)) > numreads.txt
      # TODO: we could do the same for f2 and introduce check. #shouldWe?
    else
      echo "File count error sample ${samplename} found (\${list[@]})"
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
        file('*.lostcause.txt') optional true into ch_lostcause_cram
        file('numreads.txt') optional true into ch_numreads_crams
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    def cramfile = "${samplename}.cram"
    """
    samtools merge -@ ${task.cpus} -f $cramfile ${crams}

    f1=${samplename}_1.fastq.gz
    f2=${samplename}_2.fastq.gz

    numreads=\$(samtools view -c -F 0x900 $cramfile)
    if (( numreads >= ${params.min_reads} )); then
                              # -O {stdout} -u {no compression}
                              # -N {always append /1 and /2 to the read name}
                              # -F 0x900 (bit 1, 8, filter secondary and supplementary reads)
      echo -n \$numreads > numreads.txt
      samtools collate    \\
          -O -u           \\
          -@ ${task.cpus} \\
          $cramfile pfx-${samplename} | \\
      samtools fastq      \\
          -N              \\
          -F 0x900        \\
          -@ ${task.cpus} \\
          -1 \$f1 -2 \$f2 \\
          -
      sync \$f1 \$f2          # this line and next to tackle k8s weirdness (see k8s)
      sleep 1
    else
      echo -e "${samplename}\\tcram\\tlowreads" > ${samplename}.lostcause.txt
    fi
    """
}


ch_fastqs_irods
  .mix(ch_fastqs_dirpe, ch_bams_dirpe, ch_fastqs_dirse)
  .into{ ch_rnaseq; ch_fastqc; ch_mixcr; ch_bracer; ch_tracer }

ch_rnaseq
  .until{ ! params.run_rnaseq  }
  .into { ch_hisat2; ch_star; ch_salmon }


process fastqc {
    tag "$samplename"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> !params.fastqc_publish ? null : filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

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
    tag "$samplename"
    publishDir "${params.outdir}/mixcr/", mode: 'copy',
      saveAs: { filename -> "${samplename.md5()[0..1]}/$samplename/$filename" }

    when:
    params.run_mixcr

    input:
    set val(samplename), file(reads) from ch_mixcr

    output:
    file("*full_clones.txt")
    file("*.clones.clns")
    file("*.vdjca")

    script:
    """
    mixcr align --species hsa -t ${task.cpus} $reads ${samplename}.alignments.vdjca
    mixcr assemble -t ${task.cpus} ${samplename}.alignments.vdjca ${samplename}.clones.clns
    mixcr exportClones ${samplename}.clones.clns ${samplename}.full_clones.txt
    """
}


    // TODO: we unzip fastqs here, so we spend unnecessary time zipping/unzipping.
    // Options are (1) make zipping optional (2) extract a pure bracer pipeline
    // (3) ...
process bracer_assemble {
    tag "$samplename"

    when:
    params.run_bracer

    input:
    set val(samplename), file(reads) from ch_bracer

    output:
    file('out_asm/out-*') into ch_bracer_summarise

    script:
    spec = params.bracer_genometag
    f1gz = reads[0]
    f2gz = reads[1]
    """
          # on k8s weird errors happen: gzip: Immunodeficiency7112625_1.fastq.gz: unexpected end of file
    sleep 60
    f1=\$(readlink -f $f1gz)
    f2=\$(readlink -f $f2gz)
    touch \$f1
    touch \$f2
    zcat  \$f1 > f1
    zcat  \$f2 > f2
          # output created in out_asm/out-${samplename} 
    bracer assemble -p ${task.cpus} -s $spec out-${samplename} out_asm f1 f2
    """
}


process bracer_summarise {
    tag "bracer summarise"
    publishDir "${params.outdir}/combined", mode: 'copy'      // TODO: meaningful tag, e.g. studyID.

    input:
    file('in_asm/*') from ch_bracer_summarise.collect()

    output:
    file('in_asm/filtered_BCR_summary')

    script:
    spec = params.bracer_genometag
    """
          # all the output directories of the form out-{samplename} are subdirectories of in_asm.
    bracer summarise -s $spec in_asm
    """
}


process tracer_assemble {
    tag "$samplename"

    publishDir "${params.outdir}/tracer_assemble", mode: 'copy',
      saveAs: { filename -> params.tracer_assemble_publish ? filename : null }

    when:
    params.run_tracer

    input:
    set val(samplename), file(reads) from ch_tracer

    output:
    file('out_asm/out-*') into ch_tracer_summarise

    shell:
    spec = params.tracer_genometag
    f1gz = reads[0]
    f2gz = reads[1]
    '''
          # ? output created in out_asm/out-${samplename} 
    zcat  !{f1gz} > f1
    zcat  !{f2gz} > f2
    tracer assemble --loci A B D G -p !{task.cpus} -s !{spec} -c /tracer/docker_helper_files/docker_tracer.conf f1 f2 out-!{samplename} out_asm
    '''
}



process tracer_summarise {
    tag "tracer summarise"
    publishDir "${params.outdir}/combined", mode: 'copy',
      saveAs: { filename -> "tracer_summary_${params.tracer_outputtag}" }

    input:
    file('in_asm/*') from ch_tracer_summarise.collect()

    output:
    file('in_asm/filtered_TCRAB_summary')

    shell:
    spec = params.tracer_genometag
    '''
          # all the output directories of the form out-{samplename} are subdirectories of in_asm.
    tracer summarise -p !{task.cpus} -s !{spec} -c /tracer/docker_helper_files/docker_tracer.conf in_asm

    '''
}





// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false

n_star_lowmapping = 0

def star_filter(logs) {
    def percent_aligned = 100
    logs.eachLine { line ->
        matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        if (matcher.matches()) {
            percent_aligned = matcher[0][1]
        }
    }
    if (percent_aligned.toFloat() < params.min_pct_aln.toFloat()) {
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
          { namelogbam -> star_filter(namelogbam[1]) ? 0 : 1 }

  ch_star_accept
  .map    { name, log, bam -> ["star", name, bam] }
  .into   { ch_fc_star; ch_bam_star }

              // { it -> [text: "${it[0]}\tSTAR\tlowmapping\n"] }
              // ^ This channel output will be merged with the it.text from a file
              // in other channels. This is slightly hacky. Dangersign.
              // This in pursuit of keeping track of where we lose samples.
              // If this is to be rejigged, then it is probably easiest to
              // implement the alignment check in shell code, and use file-based
              // logic similarly as in process irods and process crams_to_fastq.
  ch_star_reject
  .map    { it -> [text: "${it[0]}\tSTAR\tlowmapping\n"] }
  .set    { ch_lostcause_star }


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
    file "my_outs/${samplename}" into ch_alignment_logs_salmon

    script:
    """
    salmon quant \\
        -i $index \\
        -l ISR \\
        -p ${task.cpus} \\
        --seqBias \\
        --gcBias \\
        --posBias \\
        --no-version-check \\
        -q \\
        -o . \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -g ${trans_gene} \\
        --useVBOpt \\
        --numBootstraps 100
    mv quant.sf ${samplename}.quant.sf
    mv quant.genes.sf ${samplename}.quant.genes.sf
    mkdir -p my_outs/${samplename}/libParams
    mkdir -p my_outs/${samplename}/aux_info
    ln -f aux_info/meta_info.json my_outs/${samplename}/aux_info/meta_info.json
    ln -f libParams/flenDist.txt  my_outs/${samplename}/libParams/flenDist.txt
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}


n_hisat2_lowmapping = 0

def hisat2_filter(logs) {
    def percent_aligned = 100
    logs.eachLine { line ->
        // if ((matcher = line =~ /Overall alignment rate: \s*([\d\.]+)%/)) {
        //     percent_aligned = matcher[0][1]          
        // }
        matcher = line =~ /Aligned concordantly 1 time: \s*(\d+)\s+\((.*)%\)/
        if (matcher.matches()) {
          percent_aligned = matcher[0][2]
        }
    }
    if (percent_aligned.toFloat() < params.min_pct_aln.toFloat()) {
        n_hisat2_lowmapping++
        return false
    } else {
        return true
    }
}

process hisat2_align {

    tag "$samplename"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".hisat2_summary.txt") > 0) "HISAT2logs/$filename"
            else null
        }

    when:
    params.run_hisat2

    input:
    set val(samplename), file(reads) from ch_hisat2
    file indices from ch_hs2_indices.collect()
    file alignment_splicesites from ch_hisat2_splicesites.collect()

    output:
    set val(samplename), file('*.hisat2_summary.txt'), file("${samplename}.h2.bam") into ch_hisat2_aligned
    file "*.hisat2_summary.txt" into ch_alignment_logs_hisat2

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
            -S ${samplename}.h2.bam
    """
}

  ch_hisat2_accept = Channel.create()
  ch_hisat2_reject = Channel.create()

  ch_hisat2_aligned
      .choice(ch_hisat2_accept, ch_hisat2_reject)
          { namelogbams -> hisat2_filter(namelogbams[1]) ? 0 : 1 }

  ch_hisat2_accept
  .map    { name, log, bam -> [name, bam] }
  .set    { ch_hisat2_bam }

                  // This is slightly hacky. Dangersign. See previous dangersign.
  ch_hisat2_reject
  .map    { it -> [text: "${it[0]}\thisat2\tlowmapping\n"] }
  .mix(ch_lostcause_irods, ch_lostcause_cram, ch_lostcause_star)
  .set    { ch_lostcause }


// TODO; any reason not to merge this with the process above?
process hisat2_sort {

    tag "$samplename"

    input:
    set val(samplename), file(hisat2_bam) from ch_hisat2_bam

    output:
    set val("hisat2"), val(samplename), file("${samplename}.h2sort.bam") into ch_fc_hisat2, ch_bam_hisat2

    script:
    // def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"
    def avail_mem = task.memory == null ? '' : "-m ${ sprintf "%.0f", 0.9 * ( task.memory.toBytes() - 1000000000 ) / task.cpus}"

    """
                    # TODO: document flags
                    # 4: read unmapped
                    # 256: not primary alignment
                    # -u uncompressed (quicker/less CPU if streaming), implies bam format.
    samtools view -u -bS -F 4 -F 256 $hisat2_bam |
    samtools sort \\
        -@ ${task.cpus} $avail_mem \\
        -o ${samplename}.h2sort.bam \\
        -
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
    outfile = "${tag}.gene.fc.txt"
    """
    featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
      -o ${outfile} $pairedend                                \\
      -s $fc_direction ${extraparams} $thebam

    cut -f 1,7 ${outfile} > reduced.${outfile}   #  This
    mv reduced.${outfile} ${outfile}             #  reduces the file size from ~ 30M to ~1M

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

    tag "$metafile"
    publishDir "${params.outdir}/combined", mode: 'link'
    label 'merge_feature'

    input:
    file metafile from ch_merge_starcounts.map { it.toString() }.collectFile(name: 'star.meta', newLine: true)

    output:
    file '*-star-genecounts.txt'

    script:
    def outputname = "${params.runtag}-star-genecounts.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      -c 1 --rm-suffix .ReadsPerGene.out.tab                          \\
      -o $outputname -I $metafile
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
    label 'merge_feature'

    input:
    file metafile from ch_merge_fc_byaligner

    output:
    file '*-fc-genecounts.txt'

    shell:
    suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    aligner = metafile.baseName   // not strictly necessary
    outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    thesuffix  = suffix[aligner] ?: '.txt'
    '''
    python3 !{workflow.projectDir}/bin/merge_featurecounts.py        \\
      --rm-suffix !{thesuffix}                                       \\
      -c 1 --skip-comments --header                                  \\
      -o !{outputname} -I !{metafile}
    '''
}

process merge_salmoncounts {
    tag "${input_trans}/${input_genes}"
    publishDir "${params.outdir}/combined", mode: 'link'
    label 'merge_feature'

    input:
    file input_trans from ch_salmon_trans.map { it.toString() }.collectFile(name: 'trans.meta', newLine: true)
    file input_genes from ch_salmon_genes.map { it.toString() }.collectFile(name: 'genes.meta', newLine: true)

    output:
    set file('*counts.txt'), file('*tpm.txt')

    script:
    def outtranscount = "${params.runtag}-salmon-transcounts.txt"
    def outgenescount = "${params.runtag}-salmon-genecounts.txt"
    def outtranstpm   = "${params.runtag}-salmon-transtpm.txt"
    def outgenestpm   = "${params.runtag}-salmon-genetpm.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenescount -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -1 --skip-comments --header                                  \\
      -o $outtranscount -I $input_trans

    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -2 --skip-comments --header                                  \\
      -o $outgenestpm -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -2 --skip-comments --header                                  \\
      -o $outtranstpm -I $input_trans
    """
}

process lostcause {

    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file (inputs) from ch_lostcause.collectFile{ ['lostcause.txt', it.text] }

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
  .set{ ch_multiqc_fc_aligner }

ch_multiqc_fcbiotype
  .filter{ pick_aligner(it[0]) }
  .map{ it[1] }
  .set{ ch_multiqc_fcbiotype_aligner }

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
    file ('featureCounts/*') from ch_multiqc_fc_aligner.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from ch_multiqc_fcbiotype_aligner.collect().ifEmpty([])
    file ('star/*') from ch_alignment_logs_star.collect().ifEmpty([])
    file ('hisat2/*') from ch_alignment_logs_hisat2.collect().ifEmpty([])
    file ('salmon/*') from ch_alignment_logs_salmon.collect().ifEmpty([])

    output:
    file "*_multiqc.html"
    file "*_data"

    script:
    def filename = "${params.runtag}_multiqc.html"
    def reporttitle = "${params.runtag} (cellgeni/rnaseq)"
    """
    multiqc . -f --title "$reporttitle" --filename "$filename" -m custom_content -m featureCounts -m star -m fastqc -m salmon
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

                      // The try below is to work around zero-length-string
                      // errors found sometimes with -resume. A diagnosis is needed for
                      // that. TODO.
ch_numreads_crams
  .mix(ch_numreads_fastq, ch_numreads_fastq_se, ch_numreads_bam)
  .map { try { a = it.text.trim().toBigInteger() } catch(e) { a = 0 }; a }
  .sum()
  .subscribe{ n_numreads = it }

workflow.onComplete {

    summary = [:]
    if (params.run_star)   summary['star low mapping'] = n_star_lowmapping
    if (params.run_hisat2) summary['hisat2 low mapping'] = n_hisat2_lowmapping
    summary['total read count'] = n_numreads

    log.info "========================================="
    log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
    log.info "========================================="
}


