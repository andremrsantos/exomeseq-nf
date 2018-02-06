#! /usr/bin/env nextflow run -resume

// Pipeline version
version = "1.1.0"

// Show help message
params.help = false
// Required params
params.variants = false
// GATK parameters
params.genomeVersion = "b37"
params.dbnsfp = false
params.clinvar = false
params.kg = false
params.exac = false
params.annotations = false
// Clinical
params.clinical = false
params.lrg = false

// Check and Process parameters
if (params.help) {
    log.info """
    =================================================
    Exome-Seq: Best Practice v${version}/Annotate
    =================================================
    Required options:
    --variants      Path to variant VCF files (must be surrounded with quotes).
    Reference Databases:
    --genomeVersion Reference Genome versions. Supports either b37, or hg19.
    --snpEff        snpEff genome version.
    --dbnsfp        dbNSFP Annotation Database.
    --clinvar       Clinvar Annotation Database.
    --kg            1000 Genomes Annotation Database.
    --exac          ExAC Annotation Database.
    --clinical      Clinical Annotation Flag.
    --lrg           LRG Transcript list.
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

// Required params
params.variants = false
if (!params.variants) {
    exit(1, "The required parameter `--variants` is missing.")
}
// GATK parameters
snpeff      = fetchReference("snpEff")
dbnsfp      = file(fetchReference("dbnsfp"))
clinvar     = file(fetchReference("clinvar"))
kg          = file(fetchReference("kg"))
exac        = file(fetchReference("exac"))
annotations = fetchReference("annotations")
lrg         = file(fetchReference("lrg"))

// Log process
log.info """
================================================
ExomeSeq: Best Practice v${version}/Annotate
================================================
Variants:         ${params.variants}
== Databases
snpEff:           ${snpeff}
dbNSFP:           ${dbnsfp}
annotations:      ${annotations}
clinvar:          ${clinvar}
1000G.p3:         ${kg}
EXAc:             ${exac}
Clinical:         ${params.clinical}
LRG:              ${lrg}
== Global
Project:          ${params.project}
Output:           ${params.outdir}
===============================================
"""

// Set input channels
variant_name = Channel
    .fromPath(params.variants)
    .map { it -> it.baseName }
variants = Channel
    .fromPath(params.variants)
    .ifEmpty { exit(1, "Cannot find any file matching: ${params.variants}.") }

// Step 5.5 SnpEff
process snpeff {
    publishDir "${params.outdir}", mode: "copy",
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "reports/$fn" : "var/$fn" }

    input:
    val name from variant_name
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
