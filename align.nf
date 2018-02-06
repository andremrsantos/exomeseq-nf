#! /usr/bin/env nextflow run -resume

// Pipeline version
version = "1.1.0"

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
    Trimming options:
    --length        Minimal read lenght.
    --leading       Remove leading bases whose quality is bellow threshold.
    --trailing      Remove trailing bases whose quality is bellow threshold.
    --slidingSize   Slidding window size.
    --slidingCutoff Slidding window quality threshold.
    Other options:
    --help          Print this help text
    --project       Name of the running project
    --cpus          The number of cpus to reserve for multithread jobs
    --outdir        The output directory where the results will be saved
    --time          The maximum execution time
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
-- Trimming
Min. Length:           ${params.length}
Leading:               ${params.leading}
Trailing:              ${params.trailing}
Sld. Window size:      ${params.slidingSize}
Sld. Window threshold: ${params.slidingCutoff}
-- Global
Project:               ${params.project}
Output:                ${params.outdir}
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
