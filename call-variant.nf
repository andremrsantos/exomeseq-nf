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

// Pipeline Parameters
version = 0.1
// Required params
params.align = false
params.genome = false
params.target = false

required("align", "genome", "target") 
genome = file(params.genome)
target = file(params.target)

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
params.clinical = false
params.lrg_genome = false
params.lrg_transcript = false

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
if (params.clinical) {
    snpeff = fetchReference("lrg_genome")
    lrg_transcript = file(fetchReference("lrg_transcript"))
}

// Header log info
log.info "====================================="
log.info " Exome-Seq: Best Practice v${version}"
log.info "====================================="
def summary = [:]
summary['Alignments']      = params.align
summary['Genome']          = genome
summary['Target Interval'] = target
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
summary["Global"]  = ""
summary['Output dir']      = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

align_varcall = Channel
    .fromFilePairs(params.align + "{,.bai}")
    .map { it -> tuple(it[0], it[1][0], it[1][1]) }
    .ifEmpty { exit(1, "Cannot find any alignment matching: ${params.align}.\n") }

// Step 4.1 Haplotype Caller
process haplotype_call {
  publishDir "${params.outdir}/var/gvcf", mode: "copy"

  input:
  set val(name), file(bam), file(bai) from align_varcall

  output:
  set val(name), file("*.vcf"), file("*.vcf.idx") into varcall

  script:
  """
  gatk -T HaplotypeCaller \
    -R ${genome} \
    -I ${bam} \
    -L ${target} \
    -D ${dbsnp_all} \
    -o ${name}.vcf \
    -stand_call_conf 30 \
    -nct ${params.cpus}
  """
}

// Step 5.1 Separate SNPs and INDELS
process selectVariant {
  input:
  set val(name), file(raw_vcf), file(raw_idx) from varcall

  output:
  set val(name), file("*snps.vcf"), file("*snps.vcf.idx") into raw_snps
  set val(name), file("*indels.vcf"), file("*indels.vcf.idx") into raw_indels

  script:
  """
  gatk -T SelectVariants \
    -R ${genome} \
    --variant ${raw_vcf} \
    --out ${name}_snps.vcf \
    --selectTypeToInclude SNP
  gatk -T SelectVariants \
    -R ${genome} \
    --variant ${raw_vcf} \
    --out ${name}_indels.vcf \
    --selectTypeToInclude INDEL \
    --selectTypeToInclude MIXED \
    --selectTypeToInclude MNP
  """
}

// Step 5.2 Recalibrate SNPs
process recalibrateSNPs {
  errorStrategy 'retry'
  maxRetries 3
  
  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set val(name), file(raw_snp), file(raw_snp_idx) from raw_snps

  output:
  set val(name), file("*_snps.recal.vcf"), file("*_snps.recal.vcf.idx") into recalibrated_snps

  script:
  """
  gatk -T VariantRecalibrator \
    -R ${genome} \
    -input ${raw_snp} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${kgsnp} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=5.0 ${dbsnp} \
    -mode SNP \
    -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -allPoly \
    --maxGaussians 6 \
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 \
    -recalFile ${name}_snps.recal \
    -tranchesFile ${name}_snps.tranches \
    --num_threads ${params.cpus}
  
  gatk -T ApplyRecalibration \
    -R ${genome} \
    -L ${target} \
    -input ${raw_snp} \
    -recalFile ${name}_snps.recal \
    -tranchesFile ${name}_snps.tranches \
    -o ${name}_snps.recal.vcf \
    -ts_filter_level 99.5 \
    -mode SNP
  """
}

// Step 5.3 Recalibrate INDELS
process recalibrateIndels {
  errorStrategy 'retry'
  maxRetries 3

  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set val(name), file(raw_indel), file(raw_indel_idx) from raw_indels

  output:
  set val(name), file("*_indels.{recal,flt}.vcf"), file("*_indels.{recal,flt}.vcf.idx") into recalibrated_indels

  script:
  axiomrsc = (axiom == "") ? "-resource:axiomPoly,known=false,training=true,truth=true,prior=10.0 ${axiom}" : ""
  """
  COUNT=\$(cat ${raw_indel} | grep -v '#' | wc -l | awk '{ print \$1 }')
  if [ \$COUNT -gt 10000 ]; then
    gatk -T VariantRecalibrator \
      -R ${genome} \
      -input ${raw_indel} \
      -resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
      ${axiomrsc} \
      -resource:mills,known=false,training=true,truth=true,prior=8.0 ${kgindel} \
      -resource:dbsnp150,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
      -mode INDEL \
      -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
      -allPoly \
      --maxGaussians 4 \
      -tranche 100.0 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
      -recalFile ${name}_indels.recal \
      -tranchesFile ${name}_indels.tranches \
      --num_threads ${params.cpus}
    gatk -T ApplyRecalibration \
      -R ${genome} \
      -L ${target} \
      -input ${raw_indel} \
      -recalFile ${name}_indels.recal \
      -tranchesFile ${name}_indels.tranches \
      -o ${name}_indels.recal.vcf \
      -ts_filter_level 95.0 \
      -mode INDEL
  else
    gatk -T VariantFiltration \
      -R ${genome} \
      --variant ${raw_indel} \
      --out ${name}_indels.flt.vcf \
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
  publishDir "${params.outdir}/var/", mode: "copy"

  input:
  set val(name), file(snp), file(snp_idx) from recalibrated_snps
  set val(name), file(indel), file(indel_idx) from recalibrated_indels

  output:
  set val(name), file("*.vcf"), file("*.vcf.idx") into variants

  script:
  """
  gatk -T CombineVariants \
    -R ${genome} \
    -o ${name}.vcf \
    --variant:snps ${snp} \
    --variant:indels ${indel} \
    -genotypeMergeOptions PRIORITIZE \
    -priority snps,indels
  """
}

// Step 5.5 SnpEff
process snpeff {
  publishDir "${params.outdir}", mode: "copy",
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "reports/$fn" : "var/$fn" }

  input:
  set val(name), file(var), file(var_idx) from variants

  output:
  file("*ann.vcf") into annotated_variants
  file("*snpEff.{csv,html}") into snpeff_stats

  script:
  lrg = (params.clinical) ? "" : "-onlyTr ${lrg_transcript}" 
  """
  snpEff ann ${snpeff} \
    -stats ${name}.snpEff.html \
    -csvStats ${name}.snpEff.csv \
    -q -canon -lof ${var} ${lrg} | \
  snpSift varType - | \
  snpSift annotate ${clinvar} - | \
  snpSift dbnsfp -v -db ${dbnsfp} - > ${name}.ann.vcf
  """
}