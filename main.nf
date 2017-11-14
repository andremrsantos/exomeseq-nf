#! /usr/bin/env nextflow


/*
## TITLE

TODO: short description

### Homepage / Documentation

TODO: URL

### Authors

- Phil Ewels <phil.ewels@scilifelab.se>
- Chuan Wang <chuan.wang@scilifelab.se>
- Rickard Hammarén <rickard.hammaren@scilifelab.se>

### Pipeline overview
TODO: pipeline overview
1.   FastQC for raw sequencing reads quality control
2.   Trim Galore! for adapter trimming
3.1. Bowtie 1 alignment against miRBase mature miRNA
3.2. Post-alignment processing of miRBase mature miRNA counts
3.3. edgeR analysis on miRBase mature miRNA counts
  - TMM normalization and a table of top expression mature miRNA
    - MDS plot clustering samples
    - Heatmap of sample similarities
4.1. Bowtie 1 alignment against miRBase hairpin for the unaligned reads in step 3
4.2. Post-alignment processing of miRBase hairpin counts
4.3. edgeR analysis on miRBase hairpin counts
    - TMM normalization and a table of top expression hairpin
    - MDS plot clustering samples
    - Heatmap of sample similarities
5.1. Bowtie 2 alignment against host reference genome
5.2. Post-alignment processing of Bowtie 2
6.   NGI-Visualization of Bowtie 2 alignment statistics
7.   MultiQC
*/

def helpMessage() {
  log.info"""
  =========================================
  Exome-Seq: Best Practice v${version}
  =========================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run lghm/ExomeSeq --reads '*.fastq.gz' --genome /data/rsc/b37/human_g1k_v37.fasta
  Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes).
    --genome                      Genome reference fasta path
  Trimming options
    --length             [36/int] Discard reads that became shorter than length [int] because of either quality or adapter trimming. Default: 18
    --leading            [10/int] Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
    --trailing           [10/int] Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
    --slidingSize         [5/int] Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
    --slidingCutoff      [15/int] Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
  Other options:
    --help                        Print this help text
    --outdir                      The output directory where the results will be saved
    --cpus                [4/int] The number of cpus to reserve for multithread jobs
    --memory            ["10 GB"] The memory size to researve
    --time                        The maximum execution time
  """.stripIndent()
}

// Pipeline version
version = "0.1.0"

// Show help message
params.help = false
if (params.help) {
  helpMessage()
  exit 0
}

// Configuration
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)

// Required params
params.reads = ""
params.genome = ""

// Custom trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15

// Header log info
log.info "====================================="
log.info " Exome-Seq: Best Practice v${version}"
log.info "====================================="
def summary = [:]
summary['Reads']           = params.reads
summary['Genome']          = params.genome
summary['Trim Min Lenght'] = params.length
summary['Trim Leading']    = params.leading
summary['Trim Trailing']   = params.trailing
summary["Trim Sliding Window Size"] = params.slidingSize
summary["Trim Sliding Window Cutoff"] = params.slidingCutoff
summary['Output dir']     = params.outdir
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "====================================="

// Generate reads pairs
genome = Channel
  .fromPath(params.genome)
  .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
reads_trimming = Channel
  .fromFilePairs( params.reads, size: 2)
  .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}.\n" }

// Step 1. TrimGalore
process trimomatic {
  publishDir "${params.outdir}/trimmomatic", mode: "copy",
    overwrite: false,
    saveAs: { filename -> 
      if (filename.indexOf("trimmomatic.log") > 0) "logs/$filename"
      else if (filename.indexOf("trim.fq.gz") > 0) filename }

  input:
  set val(name), file(reads) from reads_trimming

  output:
  set val(name), file(reads), file("*_R{1,2}.trim.fq.gz") into reads_fastqc
  set val(name), file("*trim.fq.gz") into trimmed_reads
  file "*trimmomatic.log" into trimgalore_results, trimgalore_logs

  script:
  lead   = params.leading > 0 ? "LEADING:${params.leading}" : ""
  trail  = params.trailing > 0 ? "TRAILING:${params.trailing}" : ""
  slide  = (params.slidingCutoff > 0 && params.slidingSize > 0) ? "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" : ""
  minlen = params.length > 0 ? "MINLEN:${params.length}" : ""
  """
  trimmomatic PE -threads ${params.cpus} \
  $reads \
  ${name}_R1.trim.fq.gz ${name}_R1.unpair.trim.fq.gz \
  ${name}_R2.trim.fq.gz ${name}_R2.unpair.trim.fq.gz \
  $lead $trail $slide $minlen 2> ${name}.trimmomatic.log
  """
}

// Step 2. FastQC
process fastqc {
  publishDir "${params.outdir}/fastqc", mode: "copy"
    saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

  input:
  set val(name), file(reads), file(trimmed) from reads_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc -t ${params.cpus} -q $reads $trimmed
  """
}

// Step 3.1 BWA Align
process bwamem {
  publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  file "*.bam" into aligned_reads

  script:
  """
  bwa mem -M -t ${params.cpus} \
    -R \"@RG\tID:${name}\tSM:${name}\tPL:illumina\" \
    ${genome} ${reads} | \
    sambamba view -S -f bam /dev/stdin | \
    sambamba sort -o /dev/stdout /dev/stdin | \
    sambamba markdup /dev/stdin ${name}.bam
  """
}

// Step 3.2 Samtools Flagstat and Stat
process samtools_flagstat {
    publishDir "${params.outdir}/alignment/stats", mode: 'copy'

    input:
    file bam from aligned_reads

    output:
    file "${bam.baseName}_flagstat.txt" into flagstat_results
    file "${bam.baseName}_stats.txt" into samtools_stats_results

    script:
    """
    samtools flagstat $bam > ${bam.baseName}_flagstat.txt
    samtools stats $bam > ${bam.baseName}_stats.txt
    """
}

// Step X. MultiQC
process multiqc {
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  input:
  file multiqc_config
  file (fastqc:'fastqc/*') from fastqc_results.collect()
  file ('trimgalore/*') from trimgalore_results.collect()
  file ('flagstat') from flagstat_results
  file ('stats') from samtools_stats_results

  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data" into multiqc_data

  script:
  """
  multiqc -f --config $multiqc_config .
  """
}