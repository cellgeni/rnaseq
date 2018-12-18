# RNAseq pipeline
[![Build Status](https://travis-ci.org/cellgeni/RNAseq.svg?branch=devel)](https://travis-ci.org/cellgeni/RNAseq)

### Introduction

RNAseq is a bioinformatics analysis pipeline used for RNA sequencing data at
the [Cellular Genetics
program](http://www.sanger.ac.uk/science/programmes/cellular-genetics) at [the
Wellcome Sanger Institute](http://www.sanger.ac.uk/).

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics
workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads
and performs extensive quality-control on the results. It can run three
aligners: STAR, hisat2, and salmon. It will create count matrices using
featureCounts for both STAR and hisat2, it will create a merged count matrix
for salmon, and additionally it will create a count matrix from gene counts
produced by STAR itself. It is possible to run any subset of aligners.

This pipeline is primarily used with an LSF cluster and an OpenStack private
cloud. However, the pipeline should be able to run on any system that Nextflow
supports. See the [installation docs](docs/installation.md) for more
information.

### Documentation
The RNAseq pipeline comes with documentation about the pipeline, found in the
`docs/` directory.

### Credits
The original pipeline was developed at the [National Genomics
Infrastructure](https://portal.scilifelab.se/genomics/) at
[SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden by Phil Ewels
([@ewels](https://github.com/ewels)) and Rickard Hammarén
([@Hammarn](https://github.com/Hammarn)).

Our pipeline has diverged in a few ways:

- Integration of IRODS input
- 

### Diagram
This is the flow chart of the pipeline. mixcr is an additional mode that can
be enabled.

