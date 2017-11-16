#! /usr/bin/env nextflow run -resume

// Pipeline Parameters
// Configuration
params.multiqc_config = "$baseDir/resource/multiqc_config.yaml"
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

def helpMessage() {
  log.info """
  =========================================
  Exome-Seq: Best Practice v${version}
  =========================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run exomeseq/main.nf \
    --reads '*.fastq.gz' \
    --genome human_g1k_v37.fasta
  
  Mandatory arguments:
    --reads         Path to input data (must be surrounded with quotes).
    --genome        Genome reference fasta path
  
  Trimming options
    --length        Minimal read lenght. Default: ${params.length}.
    --leading       Cut bases off the start of a read whose quality is below. 
                    Default: ${params.leading}.
    --trailing      Cut bases off the end of a read whose quality is below. 
                    Default: ${params.trailing}.
    --slidingSize   In a slidding window cutoff, sets window size.
                    Default: ${params.slidingSize}.
    --slidingCutoff In a slidding window cutoff, sets window quality threshold.
                    Default: ${params.slidindCutoff}.
  
  Other options:
    --help         Print this help text
    --outdir       The output directory where the results will be saved
    --cpus         The number of cpus to reserve for multithread jobs
    --memory       The memory size to researve
    --time         The maximum execution time
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
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "====================================="

// Check genome reference
genome = 
  if (params.genome == "") exit(1, "Required reference genome")
  else file(params.genome)
// Generate reads pairs
Channel
  .fromFilePairs( params.reads, size: 2)
  .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}.\n" }
  .into { reads_trimming; reads_fastqc }

// Step 1. FastQC
process fastqc {
    publishDir "${params.outdir}/logs/fastqc", mode: "copy",
    saveAs: {
	filename ->
	filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
    }

  input:
  set val(name), file(reads) from reads_fastqc

  output:
  file "*_fastqc.{zip,html}" into fastqc_results

  script:
  """
  fastqc -t ${params.cpus} -q $reads
  """
}

// Step 2. TrimGalore
process trimomatic {
    publishDir "${params.outdir}", mode: "copy",
    overwrite: false,
    saveAs: {
      filename -> 
      if (filename.indexOf("trimmomatic.log") > 0) "logs/$filename"
      else if (filename.indexOf("fq.gz") > 0) "reads/$filename"
      else ""
    }

    input:
    set val(name), file(reads) from reads_trimming

    output:
    set val(name), file("*R{1,2}.trim.fq.gz") into trimmed_reads
    file "*trimmomatic.log" into trimmomatic_results

    script:
    lead   = params.leading > 0  ? "LEADING:${params.leading}" : ""
    trail  = params.trailing > 0 ? "TRAILING:${params.trailing}" : ""
    slide  = 
      if (params.slidingCutoff > 0 && params.slidingSize > 0) 
        "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" 
      else ""
    minlen = params.length > 0 ? "MINLEN:${params.length}" : ""
    """
    trimmomatic PE -threads ${params.cpus} \
      $reads \
      ${name}_R1.trim.fq.gz ${name}_R1.unpaired.fq.gz \
      ${name}_R2.trim.fq.gz ${name}_R2.unpaired.fq.gz \
      $lead $trail $slide $minlen \
      2> ${name}.trimmomatic.log
  """
}

// Step 3.1 BWA Align
process bwamem {
    publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

    input:
    set val(name), file(reads) from trimmed_reads

    output:
    set val(name), file("*.bam") into aligned_reads

    script:
    """
    bwa mem -M -t ${params.cpus} \
      -R \"@RG\tID:${name}\tSM:${name}\tPL:illumina\" \
      ${genome} ${reads} | \
      sambamba view -S -f bam /dev/stdin | \
      sambamba sort -m ${params.memory} -o ${name}.bam /dev/stdin
    """
}

// Step 3.2 Sambamba mark duplicates
process markdup {
    publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

    input:
    set val(name), file(bam) from aligned_reads

    output:
    set val(name), file("*.mkd.bam") into aligned_mkd_reads

    script:
    """
    sambamba markdup -t ${params.cpus} ${bam} ${name}.mkd.bam
    """
}

// Step 3.3 Samtools Flagstat and Stat
process samtools_flagstat {
    publishDir "${params.outdir}/logs", mode: 'copy'

    input:
    set val(name), file(bam) from aligned_mkd_reads

    output:
    file "${name}_flagstat.txt" into flagstat_results
    file "${name}_stats.txt" into stats_results

    script:
    """
    samtools flagstat $bam > ${name}_flagstat.txt
    samtools stats $bam > ${name}_stats.txt
    """
}

// Step 6. MultiQC
process multiqc {
  publishDir "${params.outdir}/logs/MultiQC", mode: 'copy'

  input:
  file multiqc_config
  file ('logs/fastqc/*') from fastqc_results.collect()
  file ('logs/trimmomatic*') from trimmomatic_results.collect()
  file ('flagstat') from flagstat_results
  file ('stats') from stats_results

  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data" into multiqc_data

  script:
  """
  multiqc -f --config $multiqc_config .
  """
}
