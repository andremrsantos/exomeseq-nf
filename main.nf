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
    --genomeVersion Reference Genome versions. Supports either b37, or hg19. [default: b37]     
    --snpEff        snpEff genome version [default: GRCh37.75]
    --dbsnp         dbSNP reference database [default: b37/dbSNP_150.b37.vcf.gz]
    --mills         Mills indels gold standard [default: b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
    --kgsnp         1000 Genomes High Confidence SNPS [default: b37/1000G_phase1.snps.high_confidence.b37.vcf.gz]
    --kgindel       1000 Genomes High Confidence Indels [default: b37/1000G_phase1.indels.b37.vcf.gz]
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

def required(String... args) {
  args.each { arg ->
    if (!params[arg]) exit(1, "The required parameter --${arg} is missing.")
  }
}

def fetchReference(String arg) {
  ref = params.references[params.genomeVersion]

  if (params[arg]) 
    return(file(params[arg]))
  else if (ref && ref[arg])
    return(ref[arg])
  else
    return false
}

// Pipeline version
version = "1.0.1"

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
params.reads = false
params.genome = false
params.target = false
params.bait = false

required("reads", "genome", "target")
genome = file(params.genome)
target = file(params.target)
bait   = ""
if (params.bait) {
  bait = file(params.bait)
}
// Exit if not find target
if (!target.exists())
  exit(1, "Could not find target intervals at `${target}`.")

// Custom trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15

// GATK parameters
params.genomeVersion = "b37"
params.snpEff = false
params.dbsnp = false
params.dbsnp_all = false
params.mills = false
params.kgsnp = false
params.kgindel = false
params.omni = false
params.hapmap = false
params.axiom = false
params.dbnsfp = false
params.clinvar = false

// Parse GATK params
snpeff = fetchReference("snpEff")
dbsnp = file(fetchReference("dbsnp"))
dbsnp_all = file(fetchReference("dbsnp_all"))
mills = file(fetchReference("mills"))
kgsnp = file(fetchReference("kgsnp"))
kgindel = file(fetchReference("kgindel"))
omni = file(fetchReference("omni"))
hapmap = file(fetchReference("hapmap"))
axiom = file(fetchReference("axiom"))
dbnsfp = file(fetchReference("dbnsfp"))
clinvar = file(fetchReference("clinvar"))

// Header log info
log.info "====================================="
log.info " Exome-Seq: Best Practice v${version}"
log.info "====================================="
def summary = [:]
summary['Reads']           = params.reads
summary['Genome']          = genome
summary['Target Interval'] = target
summary['Target Baits']    = bait
summary['References'] = ""
summary['dbSNP Common']    = dbsnp
summary['dbSNP All']       = dbsnp_all
summary['Mills Indels']    = mills
summary['Hapmap']          = hapmap
summary['1000G phase 3 Snps'] = kgsnp
summary['1000G phase 3 INDEL'] = kgindel
summary['1000G Omni']      = omni
summary['Axiom Exome Plus'] = axiom
summary['dbNSFP']          = dbnsfp
summary['clinvar']         = clinvar
summary['Trimming'] = ""
summary['Trim Min Lenght'] = params.length
summary['Trim Leading']    = params.leading
summary['Trim Trailing']   = params.trailing
summary["Trim Sliding Window Size"]   = params.slidingSize
summary["Trim Sliding Window Cutoff"] = params.slidingCutoff
summary["Global"]  = ""
summary['Output dir']      = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

// Generate reads pairs
Channel
  .fromFilePairs( params.reads, size: 2)
  .ifEmpty { exit(1, "Cannot find any reads matching: ${params.reads}.\n") }
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
  publishDir "${params.outdir}", mode: "copy",
    saveAs: { filename -> (filename.indexOf("log") > 0)? "reports/$filename" : "seq/$filename"}

  input:
  set val(name), file(reads) from reads_trimming

  output:
  set val(name), file("*R{1,2}.trim.fq.gz") into reads_align
  file "*trimmomatic.log" into trimmomatic_results

  script:
  lead   = params.leading > 0  ? "LEADING:${params.leading}" : ""
  trail  = params.trailing > 0 ? "TRAILING:${params.trailing}" : ""
  slide  = (params.slidingCutoff > 0 && params.slidingSize > 0) ? "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" : ""
  minlen = params.length > 0 ? "MINLEN:${params.length}" : ""
  """
  trimmomatic PE -threads ${params.cpus} \
    ${reads} \
    ${name}_R1.trim.fq.gz ${name}_R1.unpaired.fq.gz \
    ${name}_R2.trim.fq.gz ${name}_R2.unpaired.fq.gz \
    $lead $trail $slide $minlen \
    2> ${name}.trimmomatic.log
  """
}

// Step 3.1 BWA Align
process bwamem {
  input:
  set val(name), file(reads) from reads_align

  output:
  set val(name), file("*.bam"), file("*.bam.bai") into markdup_alignment

  script:
  """
  bwa mem -t ${params.cpus} \
    -R \"@RG\tID:${name}\tSM:${name}\tPL:illumina\" \
    ${genome} ${reads} | \
    sambamba view -S -f bam /dev/stdin | \
    sambamba sort -o ${name}.raw.bam /dev/stdin
  """
}

// Step 3.2 Picard Mark Duplicates
process markdup {
  publishDir "${params.outdir}/reports", mode: "copy",
    saveAs: { fn -> fn.indexOf("metrics") > 0 ? fn : "" }

  input:
  set val(name), file(bam), file(bam_idx) from markdup_alignment

  output:
  set val(name), file("*.mkd.bam"), file("*mkd.bam.bai") into recal_alignment
  file("*metrics") into markdup_results
  
  script:
  """
  gatk MarkDuplicates \
    --INPUT ${bam} \
    --OUTPUT ${name}.mkd.bam \
    --METRICS_FILE ${name}.markdup_metrics \
    --ASSUME_SORTED
  samtools index ${name}.mkd.bam
  """
}

// Step 3.3 Base Recalibration
process base_recalibration {
  publishDir "${params.outdir}/alignment", mode: "copy"

  input:
  set val(name), file(bam), file(bam_idx) from recal_alignment

  output:
  set val(name), file("*.bam"), file("*.bam.bai") into align_metrics, align_varcall

  script:
  """
  gatk BaseRecalibrator \
    --reference ${genome} \
    --input ${bam} \
    --output ${name}.recal.table \
    --known-sites ${dbsnp} \
    --known-sites ${mills} \
    --known-sites ${kgsnp} \
    --known-sites ${kgindel}
  gatk ApplyBQSR \
    --reference ${genome} \
    --input ${bam} \
    --output ${name}.bam \
    --bqsr-recal-file ${name}.recal.table
  samtools index ${name}.bam
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
  bait
  """
  gatk CollectAlignmentSummaryMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.align_metrics \
    --REFERENCE_SEQUENCE ${genome} \
    --ASSUME_SORTED
  
  gatk CollectInsertSizeMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.insert_metrics \
    --Histogram_FILE ${name} \
    --DEVIATIONS 10.0 \
    --MINIMUM_PCT 0.05 \
    --ASSUME_SORTED
  
  if [ ${bait} -ne "" ]; then
    gatk CollectHsMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.hs_metrics \
    --BAIT_INTERVALS ${bait} \
    --TARGET_INTERVALS ${target}
  fi
  """
}

// Step 4.1 Haplotype Caller
process haplotype_call {
  publishDir "${params.outdir}/var/gvcf", mode: "copy"

  input:
  set val(name), file(bam), file(bam_idx) from align_varcall

  output:
  file("*.gvcf") into varcall
  file("*.gvcf.idx") into varcall_index

  script:
  """
  gatk HaplotypeCaller \
    --reference ${genome} \
    --intervals ${target} \
    --input ${bam} \
    --output ${name}.gvcf \
    --emit-ref-confidence GVCF \
    --dbsnp ${dbsnp} 
  """
}

// Step 4.2 genotype call
process genotype_call {
  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  file (gvcfs) from varcall.collect()
  file (index) from varcall_index.collect()

  output:
  file("*.vcf") into sample_variant

  script:
  vars = gvcfs.collect({ var -> "--variant $var"}).join(" ")
  """
  gatk GenotypeGVCFs \
    --reference ${genome} \
    --intervals ${target} \
    --dbsnp ${dbsnp} \
    --output ${params.project}.raw.vcf \
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
  gatk SelectVariants \
    --reference ${genome} \
    --variant ${raw_vcf} \
    --output ${params.project}_snps.vcf \
    --select-type-to-include SNP
  gatk SelectVariants \
    --reference ${genome} \
    --variant ${raw_vcf} \
    --output ${params.project}_indels.vcf \
    --select-type-to-include INDEL \
    --select-type-to-include MIXED \
    --select-type-to-include MNP
  """
}

// Step 5.2 Recalibrate SNPs
process recalibrateSNPs {
  errorStrategy 'retry'
  maxRetries 3
  
  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set file(raw_snp), file(raw_snp_idx) from raw_snps

  output:
  set file("*_snps.recal.vcf"), file("*_snps.recal.vcf.idx") into recalibrated_snps

  script:
  inbreed = (sample_count > 10) ? "-an InbreedingCoeff": ""
  """
  gatk VariantRecalibrator \
    --reference ${genome} \
    --variant ${raw_snp} \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:${hapmap} \
    --resource omni,known=false,training=true,truth=false,prior=12.0:${omni} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:${kgsnp} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:${dbsnp} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum ${inbreed} \
    -an FS -an SOR \
    -mode SNP \
    --output ${params.project}_snps.recal \
    --tranches-file ${params.project}_snps.tranches \
    --rscript-file ${params.project}_snps.plots.R \
    --max-gaussians 6
  gatk ApplyVQSR \
    --reference ${genome} \
    --intervals ${target} \
    --variant ${raw_snp} \
    --output ${params.project}_snps.recal.vcf \
    --recal-file ${params.project}_snps.recal \
    --tranches-file ${params.project}_snps.tranches \
    --truth-sensitivity-filter-level  99.0 -mode SNP
  """
}

// Step 5.3 Recalibrate INDELS
process recalibrateIndels {
  errorStrategy 'retry'
  maxRetries 3

  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set file(raw_indel), file(raw_indel_idx) from raw_indels

  output:
  set file("*_indels.{recal,flt}.vcf"), file("*_indels.{recal,flt}.vcf.idx") into recalibrated_indels

  script:
  inbreed  = (sample_count > 10) ? "-an InbreedingCoeff": ""
  axiomrsc = (axiom == "") ? "--resource axiomPoly,known=false,training=true,truth=false,prior=10.0:${axiom}" : ""
  """
  COUNT=\$(cat ${raw_indel} | grep -v '#' | wc -l | awk '{ print \$1 }')
  if [ \$COUNT -gt 10000 ]; then
    gatk VariantRecalibrator \
      --reference ${genome} \
      --variant ${raw_indel} \
      --resource mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
      ${axiomrsc} \
      --resource mills,known=false,training=true,truth=false,prior=8.0 ${kgindel} \
      --resource dbsnp150,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
      -an QD -an MQRankSum -an ReadPosRankSum ${inbreed} \
      -an FS -an SOR \
      -mode INDEL \
      -tranche 100.0 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0
      --output ${params.project}_indels.recal \
      --tranches-file ${params.project}_indels.tranches \
      --rscript-file ${params.project}_indels.plots.R \
      --max-gaussians 4
    gatk ApplyVQSR \
      --reference ${genome} \
      --intervals ${target} \
      --variant ${raw_indel} \
      --output ${params.project}_indels.recal.vcf \
      --recal-file ${params.project}_indels.recal \
      --tranches-file ${params.project}_indels.tranches \
      --truth-sensitivity-filter-level  95.0 -mode INDEL
  else
    gatk VariantFiltration \
      --reference ${genome} \
      --variant ${raw_indel} \
      --output ${params.project}_indels.flt.vcf \
      --filter-name "GATKStandardQD" \
      --filter-expression "QD < 2.0" \
      --filter-name "GATKStandardReadPosRankSum" \
      --filter-expression "ReadPosRankSum < -20.0" \
      --filter-name "GATKStandardFS" \
      --filter-expression "FS > 200.0"
  fi
  """
}

// Step 5.4 Merge Recalibrated
process mergeVariant {
  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set file(snp), file(snp_idx) from recalibrated_snps
  set file(indel), file(indel_idx) from recalibrated_indels

  output:
  set file("*.vcf"), file("*.vcf.idx") into variants

  script:
  """
  gatk MergeVcfs \
    --INPUT ${snp} --INPUT ${indel} \
    --OUTPUT ${params.project}.vcf
  """
}

// Step 5.5 SnpEff
process snpeff {
  publishDir "${params.outdir}", mode: "copy",
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "reports/$fn" : "var/$fn" }

  input:
  set file(var), file(var_idx) from variants

  output:
  file("*ann.vcf") into annotated_variants
  file("*snpEff.{csv,html}") into snpeff_stats

  script:
  """
  snpEff ann ${snpeff} \
    -stats ${params.project}.snpEff.html \
    -csvStats ${params.project}.snpEff.csv \
    -lof -canon -strict -no-intergenic -noInteraction \
    ${var} | \
  snpSift annotate -noId -info CLNDN,CLNHGVS,CLNSIG ${clinvar} - | \
  snpSift dbnsfp -v -db ${dbnsfp} - > ${params.project}.ann.vcf
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