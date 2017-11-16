# ExomeSeq - Best Practice

Exome sequencing best practive analysis pipeline. Started on 
November 2017. This pipeline was developed based on 
[GATK's best practices](https://software.broadinstitute.org/gatk/best-practices/), 
which takes a set of FASTQ files and performs:

- alignment (BWA MEM)
- recalibration (GATK)
- variant calling (GATK)
- variant recalibration (GATK)
- variant evaluation (SnpEff)

## Requirements

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [trimmomatic](trimmomatic)
- [bwa](http://bio-bwa.sourceforge.net)
- [sambamba](http://lomereiter.github.io/sambamba/)
- [samtools](https://samtools.github.io)
- [gatk](https://software.broadinstitute.org/gatk/)
- [snpeff](http://snpeff.sourceforge.net)
- [multiqc](http://multiqc.info)

## Homepage / Documentation

[gitlab.com/andremrsantos/](https://gitlab.com/andremrsantos/exomeseq)

## Authors

- André M. Ribeiro dos Santos <andremrsantos@gmail.com>

### Pipeline overview

1.  **FastQC** for raw sequencing quality control.
2.  **Trimmomatic** to remove low quality reads.
3.  Alignment
    1. **BWA MEM** alignment and sorting.
    2. **Sambamba** Mark Duplicates
    3. **Samtools** Flagstats and Stats
4.  **GATK** Variant Calling
    1. **Base Recalibration** 
    2. **Haplotype Caller**
    3. **Call Pool Genotype**
5.  **GATK** Variant Recalibration
6.  MultiQC