#! /usr/bin/env nextflow run -resume

// Pipeline Parameters
// Configuration
params.multiqc_config = "$baseDir/resource/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)

// Required params
params.reads  = false
params.genome = false
params.target = false
params.bait   = false

// Custom trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15

// GATK parameters
params.dbsnp = false
params.mills = false
params.kgp3  = false

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
    --target        Targeted regions interval
    --bait          Bait regions interval
  
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

// Check required arguments
if (!params.reads || !params.genome || !params.target || !params.bait) {
  msg = "Parameters '--reads', '--genome', '--bait' and '--target' are required to run."
  exit(1, msg)
}
genome = file(params.genome)
target = file(params.target)
bait   = file(params.bait)

// Exit if not find target
if (!target.exists() || !bait.exists())
  exit(1, "Could not find target intervals at `${target}` or bait at `${bait}`.")

// Parse GATK params
ref_dir = genome.getParent()
ref_ver = ref_dir.getBaseName()
if (!params.dbsnp) 
  dbsnp = file("${ref_dir}/dbsnp_138.${ref_ver}.vcf.gz")
else 
  dbsnp = file(params.dbsnp)
if (!params.mills) 
  mills = file("${ref_dir}/Mills_and_1000G_gold_standard.indels.${ref_ver}.vcf.gz")
else 
  mills = file(params.mills)

if (!params.kgp3) 
  kgp3 = file("${ref_dir}/1000G_phase3_v4_20130502.sites.vcf.gz")
else
  kgp3 = file(params.kgp3)


// Header log info
log.info "====================================="
log.info " Exome-Seq: Best Practice v${version}"
log.info "====================================="
def summary = [:]
summary['Reads']           = params.reads
summary['Genome']          = genome
summary['Target Interval'] = target
summary['dbSNP']           = dbsnp
summary['Mills Indels']    = mills
summary['1000G phase 3']   = kgp3
summary['Trim Min Lenght'] = params.length
summary['Trim Leading']    = params.leading
summary['Trim Trailing']   = params.trailing
summary["Trim Sliding Window Size"]   = params.slidingSize
summary["Trim Sliding Window Cutoff"] = params.slidingCutoff
summary['Output dir']      = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

// Generate reads pairs
Channel
  .fromFilePairs( params.reads, size: 2)
  .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}.\n" }
  .into { reads_trimming; reads_fastqc }

// Step 1. FastQC
process fastqc {
  publishDir "${params.outdir}/reports/fastqc", mode: "copy",
    saveAs: { filename ->
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

// Step 2. Trimmomatic
process trimomatic {
  publishDir "${params.outdir}", mode: "copy", overwrite: false,
    saveAs: { filename -> 
      if (filename.indexOf("trimmomatic.log") > 0) "reports/$filename"
      else if (filename.indexOf("fq.gz") > 0) "seq/$filename"
      else ""
    }

  input:
  set val(name), file(reads) from reads_trimming

  output:
  set val(name), file("*R{1,2}.trim.fq.gz") into reads_align
  file "*trimmomatic.log" into trimmomatic_results

  script:
  lead   = params.leading > 0  ? "LEADING:${params.leading}" : ""
  trail  = params.trailing > 0 ? "TRAILING:${params.trailing}" : ""
  if (params.slidingCutoff > 0 && params.slidingSize > 0) 
    slide  = "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" 
  else 
    slide  = ""
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
  set val(name), file("*.bam"), file("*.bam.bai") into reads_align

  script:
  """
  bwa mem -M -t ${params.cpus} \
    -R \"@RG\tID:${name}\tSM:${name}\tPL:illumina\" \
    ${genome} ${reads} | \
    sambamba view -S -f bam /dev/stdin | \
    sambamba sort -m ${params.memory} -o ${name}.bam /dev/stdin
  """
}

// Step 3.2 Picard Mark Duplicates
process markdup {
  publishDir "${params.outdir}", mode: "copy", overwrite: false,
    saveAs: { fn -> 
      fn.indexOf("metrics") > 0 ? "reports/$fn" : "alignment/$fn"
    }

  input:
  set val(name), file(bam), file(bam_idx) from markdup_alignment

  output:
  set val(name), file("*.mkd.bam"), file("*mkd.bam.bai") into recal_alignment
  
  script:
  """
  picard MarkDuplicates \
    I=${bam} \
    O=${name}.mkd.bam \
    M=${name}.markdup_metrics \
    ASSUME_SORTED=true
  samtools index ${name}.mkd.bam
  """
}

// Step 3.3 Base Recalibration
process base_recalibration {
  publishDir "${params.outdir}/alignment", mode: "copy", overwrite: false

  input:
  set val(name), file(bam), file(bam_idx) from recal_alignment

  output:
  set val(name), file("*recal.bam"), file("*recal.bam.bai") into 
    align_metrics, align_varcall

  script:
  """
  gatk -T BaseRecalibrator \
    -R ${genome} \
    -L ${target} \
    -I ${bam} \
    -o ${name}.recal.table \
    -knownSites ${dbsnp} -knownSites ${mills} -knownSites ${kgp3} \
    -nct ${params.cpus}
  gatk -T PrintReads \
    -R ${genome} \
    -I ${bam} \
    -o ${name}.recal.bam \
    -BQSR ${name}.recal.table \
    -nct ${params.cpus}
  samtools index ${name}.recal.bam
  """
}


// Step 3.4 HsMetrics
process hs_metrics {
  publishDir "${params.outdir}/reports", mode: 'copy'

  input:
  set val(name), file(bam), file(bam_idx) from align_metrics

  output:
  file("*{metrics,pdf}") into hsmetric_results

  script:
  """
  picard CollectAlignmentSummaryMetrics I=${bam} O=${name}.align_metrics \
    REFERENCE_SEQUENCE=${genome} \
    ASSUME_SORTED=true \
    METRIC_ACCUMULATION_LEVEL="ALL_READS"
    VALIDATION_STRINGENCY=SILENT \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
  
  picarc CollectInsertSizeMetrics I=${bam} O=${name}.insert_metrics \
    HISTOGRAM_FILE=${name} \
    DEVIATIONS=10.0 \
    MINIMUM_PCT=0.05 \
    ASSUME_SORTED=true \
    METRIC_ACCUMULATION_LEVEL="ALL_READS"
    VALIDATION_STRINGENCY=SILENT \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
  
  picard CollectHsMetrics I=${bam} O=${name}.hs_metrics \
    BAIT_INTERVALS=${bait} \
    TARGET_INTERVALS=${target} \
    METRIC_ACCUMULATION_LEVEL="ALL_READS" \
    VERBOSITY=INFO \
    VALIDATION_STRINGENCY=SILENT \
    QUIET=false \
    COMPRESSION_LEVEL=5 \
    MAX_RECORDS_IN_RAM=500000 \ 
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false
  """
}


// Step 4.1 Haplotype Caller
process haplotype_call {
  publishDir "${params.outdir}/var/gvcf", mode: "copy", overwrite: false

  input:
  set val(name), file(bam), file(bam_idx) from align_varcall
  file index from align_varcall_idx

  output:
  file("*.gvcf") into varcall

  script:
  """
  gatk -T HaplotypeCaller \
    -R ${genome} \
    -I ${bam} \
    -L ${target} \
    -o ${name}.gvcf \
    --dbsnp ${dbsnp} \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    --emitRefConfidence GVCF \
    --genotyping_mode DISCOVERY \
    -stand_call_conf 30 \
    -nct ${params.cpus}
  """
}

// Step 4.2 genotype call
process genotype_call {
  publishDir "${params.outdir}/var/", mode: "copy", overwrite: false

  input:
  file (gvcfs) from varcall.collect()

  output:
  file("exoseq.vcf") into sample_variant

  script:
  vars = gvcfs.collect({ var -> "--variant $var"}).join(" ")
  """
  echo $gvcfs
  gatk -T GenotypeGVCFs \
    -R ${genome} \
    -L ${target} \
    -o exoseq.vcf \
    -D ${dbsnp} \
    ${vars}
  """ 
}

// Step 6. MultiQC
process multiqc {
  publishDir "${params.outdir}/reports/MultiQC", mode: 'copy'

  input:
  file multiqc_config
  file (fastqc) from fastqc_results.collect()
  file (trim) from trimmomatic_results.collect()
  file (flagstat) from flagstat_results.collect()
  file (stat) from stats_results.collect()
  file (metrics) from hsmetric_results.collect()

  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data" into multiqc_data

  script:
  """
  multiqc -f --config $multiqc_config .
  """
}