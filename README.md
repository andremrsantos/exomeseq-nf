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
    4. **GATK Base Recalibration**
    5. **Picard** Alignment, InsertSize, and Hs Metrics
4.  **GATK** Variant Calling
    1. **Haplotype Caller**
    2. **Call Genotype**
5.  **GATK** Variant Recalibration
    1. Separate SNPs and Indels
    2. Recali