#! /usr/bin/env nextflow run -resume

// Pipeline version
version = "1.1.0"

// Define script parameters
params.help = false
// Required params
params.aligns = false
params.gvcfs  = false
params.genome = false
params.target = false
// Define call strategies
params.bysample = false
params.byproject = true
// Reference Database parameters
params.genomeVersion = "b37"
params.dbsnp   = false
params.mills   = false
params.kgsnp   = false
params.kgindel = false
params.omni    = false
params.hapmap  = false
params.axiom   = false

// Check requested help message
if (params.help) {
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
    exit(0)
}

def fetchReference(String arg) {
  ref = params.references[params.genomeVersion]

  if (params[arg]) return(file(params[arg]))
  else if (ref && ref[arg]) return(ref[arg])
  else return false
}

// Check parameters
if (!params.genome) {
    exit(1, "Missing the required parameters --genome.")
}
if (!params.aligns && !params.gvcfs) {
    exit(1, "You must provide --aligns or --gvcfs to proceed.")
}
if (params.aligns && !params.target) {
    exit(1, "--target must be provided on --aligns mode.")
}

// Fetch required parameters
genome  = file(params.genome)
// Fetch database parameters
dbsnp   = file(fetchReference("dbsnp"))
mills   = file(fetchReference("mills"))
kgsnp   = file(fetchReference("kgsnp"))
kgindel = file(fetchReference("kgindel"))
omni    = file(fetchReference("omni"))
hapmap  = file(fetchReference("hapmap"))
axiom   = file(fetchReference("axiom"))

// Log process
log.info """
================================================
ExomeSeq: Best Practice v${version}/Variant Call
================================================
Alignment:        ${params.aligns}
GVCFs:            ${params.gvcfs}
Genome:           ${genome}
Target:           ${params.target}
== Databases
dbSNP:            ${dbsnp}
Mills:            ${mills}
1000G.p3 SNPs:    ${kgsnp}
1000G.p3 Indels:  ${kgindel}
HapMap:           ${hapmap}
1000G Omni:       ${omni}
Axiom Exome:      ${axiom}
== Global
Project:          ${params.project}
Output:           ${params.outdir}
===============================================
"""

// Check input mode --aligns or --gvcfs
if (params.aligns) {
	target = file(params.target)
    alignment = Channel
		.fromPath(params.aligns)
		.ifEmpty({ exit(1, "Could not find any file matching: ${params.align}.\n") })
    alignment_index = Channel.fromPath(params.aligns + ".bai")

    process haplotype_call {
		publishDir "${params.outdir}/var/gvcf", mode: "copy"

		input:
		file(align) from alignment
		file(index) from alignment_index

		output:
		val name into sample_names
		file "*.gvcf"  into sample_varcall, project_varcall
		file "*.gvcf.idx" into sample_varcall_index, project_varcall_index

    	script:
		name = align.baseName
    	"""
        gatk HaplotypeCaller \
        --reference ${genome} \
        --intervals ${target} \
        --input ${align} \
        --output ${name}.gvcf \
        --emit-ref-confidence GVCF \
        --dbsnp ${dbsnp}
        """
    }
} else if (params.gvcfs) {
    sample_names = Channel
	.fromPath(params.gvcfs)
	.map { it -> it.baseName }
    Channel
	.fromPath(params.gvcfs)
	.into { sample_varcall; project_varcall }
    Channel
	.fromPath(params.gvcfs + ".idx")
	.into { sample_varcall_index; project_varcall_index }
}

process genotype_by_sample {
    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    val(name) from sample_names
    file(gvcf) from sample_varcall
    file(idx) from sample_varcall_index

    output:
    val(name) into sample_genotype_name
    val(1) into sample_count
    file("*.vcf") into sample_genotype
    file("*.idx") into sample_genotype_index

    when:
    params.bysample

    script:
    """
    gatk GenotypeGVCFs \
    --reference ${genome} \
    --variant ${gvcf} \
    --dbsnp ${dbsnp} \
    --output ${name}.raw.vcf
    """
}

process genotype_by_project {
    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    file(gvcfs) from project_varcall.collect()
    file(idxs) from project_varcall_index.collect()

    output:
    val(name) into project_genotype_name
    val(count) into project_count
    file("*raw.vcf") into project_genotype
    file("*raw.vcf.idx") into project_genotype_index

    when:
    params.byproject

    script:
    name  = params.project
    count = (gvcfs instanceof List)? names.size() : 1
    vars  = gvcfs.collect({ v -> "--variant ${v}"}).join(" ")
    """
    gatk CombineGVCFs --reference ${genome} --output cohort.gvcf ${vars}
    gatk GenotypeGVCFs \
    --reference ${genome} \
    --variant cohort.gvcf \
    --dbsnp ${dbsnp} \
    --output ${name}.raw.vcf
    """
}

// Step 5.1 Separate SNPs and INDELS
process selectVariant {
    input:
    val(name) from sample_genotype_name.concat(project_genotype_name)
    val(count) from sample_count.concat(project_count)
    file(vcf) from sample_genotype.concat(project_genotype)
    file(idx) from sample_genotype_index.concat(project_genotype_index)

    output:
    val name  into snp_names, indel_names
    val count into snp_count, indel_count
    file "*snps.vcf" into snp_raw
    file "*snps.vcf.idx" into snp_raw_idx
    file "*indels.vcf" into indel_raw
    file "*indels.vcf.idx" into indel_raw_idx

    script:
    """
    gatk SelectVariants \
    --reference ${genome} \
    --variant ${vcf} \
    --output ${name}_snps.vcf \
    --select-type-to-include SNP

    gatk SelectVariants \
    --reference ${genome} \
    --variant ${vcf} \
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
    val(name) from snp_names
    val(count) from snp_count
    file(raw_snp) from snp_raw
    file(raw_idx) from snp_raw_idx

    output:
    val(name) into vcf_names
    file("*snps.recal.vcf") into snp_recal
    file("*snps.recal.vcf.idx") into snp_recal_idx

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
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
    -mode SNP \
    --output ${name}_snps.recal \
    --tranches-file ${name}_snps.tranches \
    --rscript-file ${name}_snps.plots.R \
    --max-gaussians 6

    gatk ApplyVQSR \
    --reference ${genome} \
    --variant ${raw_snp} \
    --output ${name}_snps.recal.vcf \
    --recal-file ${name}_snps.recal \
    --tranches-file ${name}_snps.tranches \
    --truth-sensitivity-filter-level 99.0 -mode SNP
    """
}

// Step 5.3 Recalibrate INDELS
process recalibrateIndels {
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    val(name) from indel_names
    val(count) from indel_count
    file(raw_indel) from indel_raw
    file(raw_idx) from indel_raw_idx

    output:
    file("*indels.recal.vcf") into indel_recal
    file("*indels.recal.vcf.idx") into indel_recal_idx

    script:
    inbreed = (count > 10)? "-an InbreedingCoeff" : ""
    axiomrsc = (axiom == "") ? "--resource axiomPoly,known=false,training=true,truth=false,prior=10.0:${axiom}" : ""
    """
    gatk VariantRecalibrator \
    --reference ${genome} \
    --variant ${raw_indel} \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    ${axiomrsc} \
    --resource kgp,known=false,training=true,truth=false,prior=8.0:${kgindel} \
    --resource dbsnp150,known=true,training=false,truth=false,prior=2.0:${dbsnp} \
    -an QD -an MQRankSum -an ReadPosRankSum ${inbreed} \
    -an FS -an SOR \
    -mode INDEL \
    -tranche 100.0 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    --output ${name}_indels.recal \
    --tranches-file ${name}_indels.tranches \
    --rscript-file ${name}_indels.plots.R \
    --max-gaussians 4

    gatk ApplyVQSR \
    --reference ${genome} \
    --variant ${raw_indel} \
    --output ${name}_indels.recal.vcf \
    --recal-file ${name}_indels.recal \
    --tranches-file ${name}_indels.tranches \
    --truth-sensitivity-filter-level  95.0 -mode INDEL
    """
}

// Step 5.4 Merge Recalibrated
process mergeVariant {
    publishDir "${params.outdir}/var/", mode: "copy"

    input:
    val(name) from vcf_names
    file(snp) from snp_recal
    file(snp_idx) from snp_recal_idx
    file(indel) from indel_recal
    file(indel_idx) from indel_recal_idx

    output:
    file("*.vcf") into vcf
    file("*.vcf.idx") into vcf_idx

    script:
    """
    gatk MergeVcfs \
    --INPUT ${snp} \
    --INPUT ${indel} \
    --OUTPUT ${name}.vcf
    """
}
