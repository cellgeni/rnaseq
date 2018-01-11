# CGI bulk-RNAseq pipeline
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)

### Introduction

CGI bulk-RNAseq is a bioinformatics analysis pipeline used for RNA sequencing data at the [Cellular Genetics program](https://ngisweden.scilifelab.se/)
at [the Wellcome Sanger Institute](http://www.sanger.ac.uk/), Sweden.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline is primarily used with an LSF cluster and an OpenStack private cloud. However, the pipeline should be able to run on any system that Nextflow supports. See the [installation docs](docs/installation.md) for more information.

### Documentation
The CGI bulk-RNAseq pipeline comes with documentation about the pipeline, found in the `docs/` directory.

### Credits
The original pipeline was developed at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).
