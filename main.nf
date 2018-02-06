#! /usr/bin/env nextflow run -resume

// Pipeline version
version = "1.0.1"

// Show help message
params.help = false
// Required params
params.reads = false
params.genome = false
params.target = false
params.bait = false
// Trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15
// GATK recalibration parameters
params.genomeVersion = "b37"
params.dbsnp = false
params.mills = false
params.kgsnp = false
params.kgindel = false
params.omni = false
params.hapmap = false
params.axiom = false
params.dbnsfp = false
params.clinvar = false
params.kg = false
params.exac = false
params.annotations = false
// Clinical
params.clinical = false
params.lrg = false
//Configuration
params.bysample = false
params.byproject = true
params.multiqc = "$baseDir/resource/multiqc_config.yaml"

// Check and process parameters
if (params.help) {
    log.info """
    ==========================================
    Exome-Seq: Best Practice v${version}/Align
    ==========================================
    Required options:
    --reads         Path to input data (must be surrounded with quotes).
    --genome        Genome reference fasta path
    --target        Targeted regions interval
    --bait          Bait regions interval
    Reference Databases:
    --genomeVersion Reference Genome versions.
    --dbsnp         dbSNP reference database.
    --mills         Mills indels gold standard.
    --kgsnp         1000 Genomes High Confidence SNPS.
    --kgindel       1000 Genomes High Confidence Indels.
    --omni          1000 Genomes Omni Reference set.
    --hapmap        HapMap Reference set.
    --axiom         Axiom Exome Reference.
    --snpEff        snpEff genome version.
    --dbnsfp        dbNSFP Annotation Database.
    --clinvar       Clinvar Annotation Database.
    --kg            1000 Genomes Annotation Database.
    --exac          ExAC Annotation Database.
    --clinical      Clinical Annotation Flag.
    --lrg           LRG Transcript list.
    Trimming options:
    --length        Minimal read lenght.
    --leading       Remove leading bases whose quality is bellow threshold.
    --trailing      Remove trailing bases whose quality is bellow threshold.
    --slidingSize   Slidding window size.
    --slidingCutoff Slidding window quality threshold.
    Other options:
    --help          Print this help text
    --cpus          The number of cpus to reserve for multithread jobs
    --project       Name of the running project
    --outdir        The output directory where the results will be saved
    """.stripIndent()
}

def fetchReference(String arg) {
    ref = params.references[params.genomeVersion]

    if (params[arg]) return(file(params[arg]))
    else if (ref && ref[arg]) return(ref[arg])
    else return false
}

if (!params.reads || !params.genome || !params.target) {
    exit(1, "Missing required parameters: --reads, --genome, or --target")
}
genome = file(params.genome)
target = file(params.target)
bait = (params.bait) ? file(params.bait) : target
// GATK recalibration parameters
dbsnp = file(fetchReference("dbsnp"))
mills = file(fetchReference("mills"))
kgsnp = file(fetchReference("kgsnp"))
kgindel = file(fetchReference("kgindel"))
omni    = file(fetchReference("omni"))
hapmap  = file(fetchReference("hapmap"))
axiom   = file(fetchReference("axiom"))

snpeff      = fetchReference("snpEff")
dbnsfp      = file(fetchReference("dbnsfp"))
clinvar     = file(fetchReference("clinvar"))
kg          = file(fetchReference("kg"))
exac        = file(fetchReference("exac"))
annotations = fetchReference("annotations")
lrg         = file(fetchReference("lrg"))

// Header log info
log.info """
==========================================
Exome-Seq: Best Practice v${version}/Align
==========================================
Reads:                 ${params.reads}
Genome:                ${genome}
Target:                ${target}
Bait:                  ${bait}
-- Reference
dbSNP:                 ${dbsnp}
Mills Indels:          ${mills}
1000G.p3 Snps:         ${kgsnp}
1000G.p3 Indels:       ${kgindel}
HapMap:                ${hapmap}
1000G Omni:            ${omni}
Axiom Exome:           ${axiom}
snpEff:                ${snpeff}
dbNSFP:                ${dbnsfp}
annotations:           ${annotations}
clinvar:               ${clinvar}
1000G.p3:              ${kg}
EXAc:                  ${exac}
Clinical:              ${params.clinical}
LRG:                   ${lrg}
-- Trimming
Min. Length:           ${params.length}
Leading:               ${params.leading}
Trailing:              ${params.trailing}
Sld. Window size:      ${params.slidingSize}
Sld. Window threshold: ${params.slidingCutoff}
-- Global
Project:               ${params.project}
Output:                ${params.outdir}
MultiQC:               ${params.multiqc}
==========================================
"""

// Generate input channel
Channel
    .fromFilePairs(params.reads, size: 2)
    .ifEmpty { exit(1, "Cannot find any reads matching: ${params.reads}.\n") }
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
    publishDir "${params.outdir}", mode: "copy",
    saveAs: { filename ->
	(filename.indexOf("log") > 0)? "reports/$filename" : "seq/$filename"
    }

    input:
    set val(name), file(reads) from reads_trimming

    output:
    set val(name), file("*R{1,2}.trim.fq.gz") into reads_align
    file "*trimmomatic.log" into trimmomatic_results

    script:
    lead  = (params.leading > 0) ? "LEADING:${params.leading}" : ""
    trail = (params.trailing > 0) ? "TRAILING:${params.trailing}" : ""
    slide = (params.slidingCutoff > 0 && params.slidingSize > 0) ? "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" : ""
    len   = (params.length > 0) ? "MINLEN:${params.length}" : ""
    """
    trimmomatic PE -threads ${params.cpus} ${reads} \
        ${name}_R1.trim.fq.gz ${name}_R1.unpaired.fq.gz \
        ${name}_R2.trim.fq.gz ${name}_R2.unpaired.fq.gz \
        $lead $trail $slide $len \
    2> ${name}.trimmomatic.log
    """
}

// Step 3.1 BWA Align
process bwamem {
    publishDir "${params.outdir}/reports", mode: "copy",
    saveAs: { filename -> (filename.indexOf("log") > 0)? filename : "" }

    input:
    set val(name), file(reads) from reads_align

    output:
    set val(name), file("*.bam"), file("*.bam.bai") into recal_alignment
    file "*.samblaster.log" into markdup_results

    script:
    rgid = "@RG\tID:${name}\tSM:${name}\tPL:illumina"
    """
    bwa mem -t ${params.cpus} -R "${rgid}" ${genome} ${reads} | \
    samblaster 2> ${name}.samblaster.log | \
    sambamba view -S -f bam /dev/stdin | \
    sambamba sort -o ${name}.raw.bam /dev/stdin
    """
}

// Step 3.3 Base Recalibration
process base_recalibration {
    publishDir "${params.outdir}/alignment", mode: "copy"

    input:
    set val(name), file(bam), file(bam_idx) from recal_alignment

    output:
    set val(name), file("*.bam"), file("*.bam.bai") into align_hs, align_ins, align_st, align_varcall
    file "*.bam" into alignment
    file "*.bam.bai" into alignment_index

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
process align_metrics {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_idx) from align_st

    output:
    file("*{metrics,pdf}") into align_metric_results

    script:
    """
    gatk CollectAlignmentSummaryMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.align_metrics \
    --REFERENCE_SEQUENCE ${genome} \
    --ASSUME_SORTED
    """
}

// Step 3.4 HsMetrics
process insert_metrics {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_idx) from align_ins

    output:
    file("*{metrics,pdf}") into insert_metric_results

    script:
    """
    gatk CollectInsertSizeMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.insert_metrics \
    --Histogram_FILE ${name} \
    --DEVIATIONS 10.0 \
    --MINIMUM_PCT 0.05 \
    --ASSUME_SORTED
    """
}

// Step 3.4 HsMetrics
process hs_metrics {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_idx) from align_hs

    output:
    file("*{metrics,pdf}") into hs_metric_results

    script:
    bait
    """
    gatk CollectHsMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.hs_metrics \
    --BAIT_INTERVALS ${bait} \
    --TARGET_INTERVALS ${target}
    """
}

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
    --intervals ${target} \
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
    --intervals ${target} \
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
    -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum ${inbreed} \
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
    --truth-sensitivity-filter-level 99.5 -mode SNP
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
    axiomrsc = (axiom == "") ? "" : "--resource axiom,known=false,training=true,truth=false,prior=10.0:${axiom}"
    """
    gatk VariantRecalibrator \
    --reference ${genome} \
    --variant ${raw_indel} \
    --resource mills,known=false,training=true,truth=true,prior=12.0:${mills} \
    ${axiomrsc} \
    --resource kgp,known=false,training=true,truth=false,prior=8.0:${kgindel} \
    --resource dbsnp150,known=true,training=false,truth=false,prior=2.0:${dbsnp} \
    -an DP -an QD -an MQRankSum -an ReadPosRankSum ${inbreed} \
    -an FS -an SOR \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
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
    --truth-sensitivity-filter-level  99.0 -mode INDEL
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
    val name into variant_names
    file("*.vcf") into variants

    script:
    """
    gatk MergeVcfs \
    --INPUT ${snp} \
    --INPUT ${indel} \
    --OUTPUT ${name}.vcf
    """
}

// Step 5.5 SnpEff
process snpeff {
    publishDir "${params.outdir}", mode: "copy",
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "reports/$fn" : "var/$fn" }

    input:
    val name from variant_names
    file var from variants

    output:
    file("*ann.vcf") into annotated_variants
    file("*snpEff.{csv,html}") into snpeff_stats

    script:
    clinical = (params.clinical)? "-onlyTr ${lrg}" : ""
    annotation_list = annotations.join(",")
    """
    snpEff ann ${snpeff} \
        -stats ${name}.snpEff.html \
        -csvStats ${name}.snpEff.csv \
        -canon -no-intergenic -noInteraction \
        ${clinical} ${var} | \
    snpSift annotate -noId -name KG -info AF ${kg} - | \
    snpSift annotate -noId -name EXAC -info AF ${exac} - | \
    snpSift annotate -noId -info CLNDN,CLNHGVS,CLNSIG ${clinvar} - | \
    snpSift dbnsfp \
        -f ${annotation_list} \
        -db ${dbnsfp} - > ${name}.ann.vcf
    """
}

// Step 6. MultiQC
process multiqc {
    publishDir "${params.outdir}/reports/MultiQC", mode: 'copy'
    errorStrategy "ignore"

    input:
    file (fastqc)  from fastqc_results.collect()
    file (trim)    from trimmomatic_results.collect()
    file (markdup) from markdup_results.collect()
    file (align_metrics) from align_metric_results.collect()
    file (insert_metrics) from insert_metric_results.collect()
    file (hs_metric_results) from hs_metric_results.collect()
    file (snpeff)  from snpeff_stats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_data

    script:
    """
    multiqc -f --config ${params.multiqc} .
    """
}
