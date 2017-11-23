#! /usr/bin/env nextflow run -resume

def helpMessage() {
  log.info """
  =========================================
  Exome-Seq: Best Practice v${version}
  =========================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run exomeseq/main.nf --reads '*.fastq.gz' --genome human_g1k_v37.fasta

  Mandatory options:
    --reads         Path to input data (must be surrounded with quotes).
    --genome        Genome reference fasta path
    --target        Targeted regions interval
    --bait          Bait regions interval
  Reference Databases:
    --snpEff        snpEff genome version [default: GRCh37.75]
    --dbsnp         dbSNP reference database [default: b37/dbSNP_150.b37.vcf.gz]
    --mills         Mills indels gold standard [default: b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
    --kgp3          1000 Genomes High Confidence SNPS [default: b37/1000G_phase1.snps.high_confidence.b37.vcf.gz]
    --omni          1000 Genomes Omni Reference set [default: b37/1000G_omni2.5.b37.vcf.gz]
    --hapmap        HapMap Reference set [default: b37/hapmap_3.3.b37.vcf.gz]
    --axiom         Axiom Exome Reference [default: b37/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz]
    --dbnsfp        dbNSFP Annotation Database [default: dbnsfp/dbNSFP2.9.3.txt.gz]
    --clinvar       Clinvar Annotation Database [default: clinvar/hg19/clinvar_20171029.vcf.gz]
  Trimming options:
    --length        Minimal read lenght. [default: ${params.length}]
    --leading       Cut bases off the start of a read whose quality is below. [default: ${params.leading}]
    --trailing      Cut bases off the end of a read whose quality is below. [default: ${params.trailing}]
    --slidingSize   In a slidding window cutoff, sets window size. [default: ${params.slidingSize}]
    --slidingCutoff In a slidding window cutoff, sets window quality threshold. [default: ${params.slidindCutoff}]
  Other options:
    --help          Print this help text
    --project       Name of the running project
    --cpus          The number of cpus to reserve for multithread jobs
    --memory        The memory size to researve
    --outdir        The output directory where the results will be saved
    --time          The maximum execution time
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

// Pipeline Parameters
// Configuration
params.multiqc_config = "$baseDir/resource/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)

// Required params
params.reads  = false
params.genome = false
params.target = false
params.bait   = false

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

// Custom trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15

// GATK parameters
params.dbsnp = false
params.mills = false
params.kgp3 = false
params.omni = false
params.hapmap = false
params.axiom = false
params.dbnsfp = false
params.clinvar = false

// Parse GATK params
ref_dir = genome.getParent()
ref_ver = ref_dir.getBaseName()
if (!params.dbsnp) 
  dbsnp = file("${ref_dir}/dbsnp_150.${ref_ver}.vcf.gz")
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

if (!params.hapmap)
  hapmap = file("${ref_dir}/hapmap_3.3.${ref_ver}.vcf.gz")
else
  hapmap = file(params.hapmap)

if (!params.omni)
  omni = file("${ref_dir}/1000G_omni2.5.${ref_ver}.vcf.gz")
else
  omni = file(params.omni)

if (!params.axiom)
  axiom = file("${ref_dir}/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz")
else
  axiom = file(params.axiom)

if (!params.dbnsfp)
  dbnsfp = file("/data/rsc/dbnsfp/dbNSFP2.9.3.txt.gz")
else
  dbnsfp = file(params.dbnsfp)

if (!params.clinvar)
  clinvar = file("/data/rsc/clinvar/hg19/clinvar_20171029.vcf.gz")
else
  clinvar = file(params.clinvar)

// Header log info
log.info "====================================="
log.info " Exome-Seq: Best Practice v${version}"
log.info "====================================="
def summary = [:]
summary['Reads']           = params.reads
summary['Genome']          = genome
summary['Target Interval'] = target
summary['Target Baits']    = bait
summary['Mutation References'] = ""
summary['dbSNP']           = dbsnp
summary['Mills Indels']    = mills
summary['Hapmap']          = hapmap
summary['1000G phase 3']   = kgp3
summary['1000G Omni']      = omni
summary['Axiom Exome Plus'] = axiom
summary['dbNSFP']          = dbnsfp
summary['clinvar']         = clinvar
summary['Trimming Options'] = ""
summary['Trim Min Lenght'] = params.length
summary['Trim Leading']    = params.leading
summary['Trim Trailing']   = params.trailing
summary["Trim Sliding Window Size"]   = params.slidingSize
summary["Trim Sliding Window Cutoff"] = params.slidingCutoff
summary["Global Options"]  = ""
summary['Output dir']      = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

// Generate reads pairs
Channel
  .fromFilePairs( params.reads, size: 2)
  .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}.\n" }
  .into { reads_trimming; reads_fastqc; reads_count }

sample_count = 0
reads_count.subscribe { sample_count += 1 }

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
  set val(name), file(reads) from reads_align

  output:
  set val(name), file("*.bam"), file("*.bam.bai") into markdup_alignment

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
  file("*metrics") into markdup_results
  
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
  set val(name), file("*recal.bam"), file("*recal.bam.bai") into align_metrics, align_varcall

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
  
  picard CollectInsertSizeMetrics I=${bam} O=${name}.insert_metrics \
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
  file("*.vcf") into sample_variant

  script:
  vars = gvcfs.collect({ var -> "--variant $var"}).join(" ")
  """
  gatk -T GenotypeGVCFs \
    -R ${genome} \
    -L ${target} \
    -o ${params.project}.raw.vcf \
    -D ${dbsnp} \
    ${vars}
  """ 
}

// Step 5.1 Separate SNPs and INDELS
process selectVariant {
  input:
  file(raw_vcf) from sample_variant

  output:
  set file("*snps.vcf"), file("*snps.vcf.idx") into raw_snps
  set file("*indels.vcf"), file("*indels.vcf.idx") into raw_indels

  script:
  """
  gatk -T SelectVariants \
    -R ${genome} \
    --variant ${raw_vcf} \
    --out ${params.project}_snps.vcf \
    --selectTypeToInclude SNP
  gatk -T SelectVariants \
    -R ${genome} \
    --variant ${raw_vcf} \
    --out ${params.project}_indels.vcf \
    --selectTypeToInclude INDEL \
    --selectTypeToInclude MIXED \
    --selectTypeToInclude MNP
  """
}

// Step 5.2 Recalibrate SNPs
process recalibrateSNPs {
  errorStrategy 'retry'
  maxRetries 3
  
  publishDir "${params.outdir}/var/", mode: "copy", overwrite: false

  input:
  set file(raw_snp), file(raw_snp_idx) from raw_snps

  output:
  set file("*_snps.recal.vcf"), file("*_snps.recal.vcf.idx") into recalibrated_snps

  script:
  inbreed = (sample_count > 10) ? "-an InbreedingCoeff": ""
  """
  gatk -T VariantRecalibrator \
    -R ${genome} \
    -input ${raw_snp} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${kgp3} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=5.0 ${dbsnp} \
    -mode SNP \
    -an QD -an FS -an SOR -an MQ -an MQRankSum \
    -an ReadPosRankSum ${inbreed} \
    -allPoly \
    --maxGaussians 6 \
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 \
    -recalFile ${params.project}_snps.recal \
    -tranchesFile ${params.project}_snps.tranches \
    --num_threads ${params.cpus}
  
  gatk -T ApplyRecalibration \
    -R ${genome} \
    -L ${target} \
    -input ${raw_snp} \
    -recalFile ${params.project}_snps.recal \
    -tranchesFile ${params.project}_snps.tranches \
    -o ${params.project}_snps.recal.vcf \
    -ts_filter_level 99.5 \
    -mode SNP
  """
}

// Step 5.3 Recalibrate INDELS
process recalibrateIndels {
  errorStrategy 'retry'
  maxRetries 3

  publishDir "${params.outdir}/var/", mode: "copy", overwrite: false

  input:
  set file(raw_indel), file(raw_indel_idx) from raw_indels

  output:
  set file("*_indels.flt.vcf"), file("*_indels.flt.vcf.idx") into recalibrated_indels

  script:
  inbreed = (sample_count > 10) ? "-an InbreedingCoeff": ""
  """
  COUNT=\$(cat ${raw_indel} | grep -v '#' | wc -l | awk '{ print \$1 }')
  if [ \$COUNT -gt 10000 ]; then
    gatk -T VariantRecalibrator \
      -R ${genome} \
      -input ${raw_indel} \
      -resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
      -resource:axiomPoly,known=false,training=true,truth=true,prior=10.0 ${axiom} \
      -resource:dbsnp150,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
      -mode INDEL \
      -an QD -an FS -an SOR -an MQRankSum \
      -an ReadPosRankSum ${inbreed} \
      -allPoly \
      --maxGaussians 4 \
      -tranche 100.0 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
      -recalFile ${params.project}_indels.recal \
      -tranchesFile ${params.project}_indels.tranches \
      --num_threads ${params.cpus}
    gatk -T ApplyRecalibration \
      -R ${genome} \
      -L ${target} \
      -input ${raw_indel} \
      -recalFile ${params.project}_indels.recal \
      -tranchesFile ${params.project}_indels.tranches \
      -o ${params.project}_indels.recal.vcf \
      -ts_filter_level 95.0 \
      -mode SNP
  else
    gatk -T VariantFiltration \
      -R ${genome} \
      --variant ${raw_indel} \
      --out ${params.project}_indels.flt.vcf \
      --filterName GATKStandardQD \
      --filterExpression "QD < 2.0" \
      --filterName GATKStandardReadPosRankSum \
      --filterExpression "ReadPosRankSum < -20.0" \
      --filterName GATKStandardFS \
      --filterExpression "FS > 200.0"
  fi
  """
}

// Step 5.4 Merge Recalibrated
process mergeVariant {
  publishDir "${params.outdir}/var/", mode: "copy", overwrite: false

  input:
  set file(snp), file(snp_idx) from recalibrated_snps
  set file(indel), file(indel_idx) from recalibrated_indels

  output:
  set file("*.vcf"), file("*.vcf.idx") into variants

  script:
  """
  gatk -T CombineVariants \
    -R ${genome} \
    -o ${params.project}.final.vcf \
    --variant:snps ${snp} \
    --variant:indels ${indel} \
    -genotypeMergeOptions PRIORITIZE \
    -priority snps,indels
  """
}

// Step 5.5 SnpEff
process snpeff {
  publishDir "${params.outdir}", mode: "copy", overwrite: false,
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "report/$fn" : "var/$fn" }

  input:
  set file(var), file(var_idx) from variants

  output:
  file("*ann.vcf") into annotated_variants
  file("*snpEff.{csv,html}") into snpeff_stats

  script:
  """
  snpEff ann GRCh37.75 \
      -stats ${params.project}.snpEff.html \
      -csvStats ${params.project}.snpEff.csv \
      -q -canon -lof ${var} | \
    snpSift varType - | \
    snpSift annotate ${clinvar} - | \
    snpSift dbnsfp -v -db ${dbnsfp} - \
    > ${params.project}.ann.vcf
  """
}

// Step 6. MultiQC
process multiqc {
  publishDir "${params.outdir}/reports/MultiQC", mode: 'copy'
  errorStrategy "ignore"  

  input:
  file multiqc_config
  file (fastqc)  from fastqc_results.collect()
  file (trim)    from trimmomatic_results.collect()
  file (markdup) from markdup_results.collect()
  file (metrics) from hsmetric_results.collect()
  file (snpeff)  from snpeff_stats.collect()

  output:
  file "*multiqc_report.html" into multiqc_report
  file "*_data" into multiqc_data

  script:
  """
  multiqc -f --config $multiqc_config .
  """
}