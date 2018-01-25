#! /usr/bin/env nextflow run -resume

def helpMessage() {
    log.info """
    ==========================================
    Exome-Seq: Best Practice v${version}/Align
    ==========================================
    Required options:
    --reads             Path to input data (must be surrounded with quotes).
    --genome            Genome reference fasta path
    --target            Targeted regions interval
    --bait              Bait regions interval

    Reference Databases:
    --genomeVersion     Reference Genome versions. Supports either b37, or hg19. [default: b37] 
    --dbsnp             dbSNP reference database [default: b37/dbSNP_150.b37.vcf.gz]
    --mills             Mills indels gold standard [default: b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
    --kgsnp             1000 Genomes High Confidence SNPS [default: b37/1000G_phase1.snps.high_confidence.b37.vcf.gz]
    --kgindel           1000 Genomes High Confidence Indels [default: b37/1000G_phase1.indels.b37.vcf.gz]
    
    Trimming options:
    --length            Minimal read lenght. [default: ${params.length}]
    --leading           Cut bases off the start of a read whose quality is below. [default: ${params.leading}]
    --trailing          Cut bases off the end of a read whose quality is below. [default: ${params.trailing}]
    --slidingSize       In a slidding window cutoff, sets window size. [default: ${params.slidingSize}]
    --slidingCutoff     In a slidding window cutoff, sets window quality threshold. [default: ${params.slidindCutoff}]
    
    Other options:
    --help              Print this help text
    --project           Name of the running project
    --cpus              The number of cpus to reserve for multithread jobs
    --outdir            The output directory where the results will be saved
    --time              The maximum execution time
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
params.reads = false
params.genome = false
params.target = false
params.bait = false

required("reads", "genome", "target")
genome = file(params.genome)
target = file(params.target)
bait = ""
if (params.bait) {
    bait = file(params.bait)
}

// Custom trimming options
params.length = 36
params.leading = 10
params.trailing = 10
params.slidingSize = 5
params.slidingCutoff = 15

// GATK recalibration parameters
params.genomeVersion = "b37"
params.dbsnp = false
dbsnp = file(fetchReference("dbsnp"))
params.mills = false
mills = file(fetchReference("mills"))
params.kgsnp = false
kgsnp = file(fetchReference("kgsnp"))
params.kgindel = false
kgindel = file(fetchReference("kgindel"))

// Header log info
log.info "=========================================="
log.info "Exome-Seq: Best Practice v${version}/Align"
log.info "=========================================="
def summary = [:]
summary['Reads'] = params.reads
summary['Genome'] = genome
summary['Target Interval'] = target
summary['Target Baits'] = bait
summary['Reference'] = ""
summary['dbSNP Common'] = dbsnp
summary['Mills Indels'] = mills
summary['1000G phase 3 Snps'] = kgsnp
summary['1000G phase 3 INDEL'] = kgindel
summary['Trimming'] = ""
summary['Trim Min Lenght'] = params.length
summary['Trim Leading'] = params.leading
summary['Trim Trailing'] = params.trailing
summary["Trim Sliding Window Size"] = params.slidingSize
summary["Trim Sliding Window Cutoff"] = params.slidingCutoff
summary["Global"] = ""
summary["Project"] = params.project
summary['Output dir'] = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

// Generate reads pairs
Channel
    .fromFilePairs(params.reads, size: 2)
    .ifEmpty { exit(1, "Cannot find any reads matching: ${params.reads}.\n") }
    .into { reads_trimming; reads_fastqc; sample_counts }

sample_count = 0
sample_counts.subscribe { sample_count += 1 }

// Step 1. FastQC
process fastqc {
    publishDir "${params.outdir}/reports/fastqc", 
    mode: "copy",
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
    lead = params.leading > 0 ? "LEADING:${params.leading}" : ""
    trail = params.trailing > 0 ? "TRAILING:${params.trailing}" : ""
    slide = (params.slidingCutoff > 0 && params.slidingSize > 0) ? "SLIDINGWINDOW:${params.slidingSize}:${params.slidingCutoff}" : ""
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
    samblaster 2&> ${name}.samblaster.log | \
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
    gatk CollectAlignmentSummaryMetrics 
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

    gatk CollectHsMetrics \
    --INPUT ${bam} \
    --OUTPUT ${name}.hs_metrics \
    --BAIT_INTERVALS ${bait} \
    --TARGET_INTERVALS ${target}
    """
}