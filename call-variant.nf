#! /usr/bin/env nextflow run -resume

def helpMessage() {
    log.info """
    =================================================
    Exome-Seq: Best Practice v${version}/Call Variant
    =================================================
    Required options:
    --aligns        Path to alignment files (must be surrounded with quotes).
    --genome        Genome reference fasta path
    --target        Targeted regions interval
  
    Reference Databases:
    --genomeVersion Reference Genome versions. Supports either b37, or hg19. [default: b37]     
    --dbsnp         dbSNP reference database [default: b37/dbSNP_150.b37.vcf.gz]
    --mills         Mills indels gold standard [default: b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
    --kgsnp         1000 Genomes High Confidence SNPS [default: b37/1000G_phase1.snps.high_confidence.b37.vcf.gz]
    --kgindel       1000 Genomes High Confidence Indels [default: b37/1000G_phase1.indels.b37.vcf.gz]
    --omni          1000 Genomes Omni Reference set [default: b37/1000G_omni2.5.b37.vcf.gz]
    --hapmap        HapMap Reference set [default: b37/hapmap_3.3.b37.vcf.gz]
    --axiom         Axiom Exome Reference [default: b37/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz]

    Other options:
    --help          Print this help text
    --project       Name of the running project
    --cpus          The number of cpus to reserve for multithread jobs
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

  if (params[arg]) return(file(params[arg]))
  else if (ref && ref[arg]) return(ref[arg])
  else return false
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
// Required params
params.aligns = false
params.genome = false
params.target = false

params.bysample = false
params.byproject = true

required("aligns", "genome", "target") 
genome = file(params.genome)
target = file(params.target)

// GATK parameters
params.genomeVersion = "b37"
params.dbsnp = false
dbsnp = file(fetchReference("dbsnp"))
params.mills = false
mills = file(fetchReference("mills"))
params.kgsnp = false
kgsnp = file(fetchReference("kgsnp"))
params.kgindel = false
kgindel = file(fetchReference("kgindel"))
params.omni = false
omni = file(fetchReference("omni"))
params.hapmap = false
hapmap = file(fetchReference("hapmap"))
params.axiom = false
axiom = file(fetchReference("axiom"))

// Header log info
log.info "================================================="
log.info "Exome-Seq: Best Practice v${version}/Variant Call"
log.info "================================================="
def summary = [:]
summary['Alignments'] = params.aligns
summary['Genome'] = genome
summary['Target Interval'] = target
summary['References'] = ""
summary['dbSNP Common'] = dbsnp
summary['Mills Indels'] = mills
summary['Hapmap'] = hapmap
summary['1000G phase 3 Snps'] = kgsnp
summary['1000G phase 3 INDEL'] = kgindel
summary['1000G Omni'] = omni
summary['Axiom Exome'] = axiom
summary["Global"] = ""
summary["Project"] = params.project
summary['Output dir'] = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

align_varcall = Channel
    .fromFilePairs(params.aligns + "{,.bai}")
    .map { it -> [it[0], file(it[1][0]), file(it[1][1])] }
    .ifEmpty { exit(1, "Cannot find any alignment matching: ${params.align}.\n") }

// Step 4.1 Haplotype Caller
process haplotype_call {
    publishDir "${params.outdir}/var/gvcf", mode: "copy"

    input:
    set val(name), file(bam), file(bam_idx) from align_varcall

    output:
    set val(name), file("*.gvcf"), file("*.gvcf.idx") into varcall_by_sample, varcall_by_project

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

process genotype_by_sample {
    publishDir "${params.outdir}/var/", mode: "copy"
    
    input:
    set val(name), file(gvcf), file(idx) from varcall_by_sample
    
    output:
    set val(name), val(1), file("*.vcf"), file("*vcf.idx") into genotyped_by_sample

    when:
    params.bysample

    script:
    """
    gatk GenotypeGVCFs \
    --reference ${genome} \
    --intervals ${target} \
    --variant ${gvcf} \
    --dbsnp ${dbsnp} \
    --output ${name}.raw.vcf
    """ 
}

process genotype_by_project {
    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    set val(names), file(gvcfs), file(idxs) from varcall_by_project.collect()

    output:
    set val(params.project), val(count), file("*.vcf"), file("*vcf.idx") into genotyped_by_project

    when:
    params.byproject

    script:
    count = (names instanceof List)? names.size() : 1 
    variants = gvcfs.collect({ v -> "--variant ${v}"}).join(" ")
    """
    echo ${count}
    gatk GenotypeGVCFs \
    --reference ${genome} \
    --intervals ${target} \
    ${variants} \
    --dbsnp ${dbsnp} \
    --output ${params.project}.raw.vcf
    """ 
}

// Step 5.1 Separate SNPs and INDELS
process selectVariant {
    input:
    set val(name), val(count), file(raw_vcf), file(idx) from genotyped_by_sample.mix(genotyped_by_project)
    
    output:
    set val(name), val(count), file("*snps.vcf"), file("*snps.vcf.idx") into raw_snps
    set val(name), val(count), file("*indels.vcf"), file("*indels.vcf.idx") into raw_indels

    script:
    """
    gatk SelectVariants \
    --reference ${genome} \
    --variant ${raw_vcf} \
    --output ${name}_snps.vcf \
    --select-type-to-include SNP

    gatk SelectVariants \
    --reference ${genome} \
    --variant ${raw_vcf} \
    --output ${name}_indels.vcf \
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
    set val(name), val(count), file(raw_snp), file(idx) from raw_snps

    output:
    set val(name), file("*_snps.recal.vcf"), file("*_snps.recal.vcf.idx") into recalibrated_snps

    script:
    inbreed = (count > 10)? "-an InbreedingCoeff" : ""
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
    --output ${name}_snps.recal \
    --tranches-file ${name}_snps.tranches \
    --rscript-file ${name}_snps.plots.R \
    --max-gaussians 6

    gatk ApplyVQSR \
    --reference ${genome} \
    --intervals ${target} \
    --variant ${raw_snp} \
    --output ${name}_snps.recal.vcf \
    --recal-file ${name}_snps.recal \
    --tranches-file ${name}_snps.tranches \
    --truth-sensitivity-filter-level  99.0 -mode SNP
    """
    }

// Step 5.3 Recalibrate INDELS
process recalibrateIndels {
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    set val(name), val(count), file(raw_indel), file(idx) from raw_indels

    output:
    set val(name), file("*_indels.{recal,flt}.vcf"), file("*_indels.{recal,flt}.vcf.idx") into recalibrated_indels

    script:
    inbreed = (count > 10)? "-an InbreedingCoeff" : ""
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
        --output ${name}_indels.recal \
        --tranches-file ${name}_indels.tranches \
        --rscript-file ${name}_indels.plots.R \
        --max-gaussians 4

        gatk ApplyVQSR \
        --reference ${genome} \
        --intervals ${target} \
        --variant ${raw_indel} \
        --output ${name}_indels.recal.vcf \
        --recal-file ${name}_indels.recal \
        --tranches-file ${name}_indels.tranches \
        --truth-sensitivity-filter-level  95.0 -mode INDEL
    else
        gatk VariantFiltration \
        --reference ${genome} \
        --variant ${raw_indel} \
        --output ${name}_indels.flt.vcf \
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
    set val(name), file(snp), file(snp_idx) from recalibrated_snps
    set val(_name), file(indel), file(indel_idx) from recalibrated_indels

    output:
    set val(name), file("*.vcf"), file("*.vcf.idx") into variants

    script:
    """
    gatk MergeVcfs \
    --INPUT ${snp} \
    --INPUT ${indel} \
    --OUTPUT ${name}.vcf
    """
}