#! /usr/bin/env nextflow run -resume

def helpMessage() {
    log.info """
    =================================================
    Exome-Seq: Best Practice v${version}/Call Variant
    =================================================
    Required options:
    --variants      Path to variants files (must be surrounded with quotes).
  
    Reference Databases:
    --genomeVersion Reference Genome versions. Supports either b37, or hg19. [default: b37]     
    --snpEff        snpEff genome version [default: GRCh37.75]
    --dbnsfp        dbNSFP Annotation Database [default: dbnsfp/dbNSFP2.9.3.txt.gz]
    --clinvar       Clinvar Annotation Database [default: clinvar/hg19/clinvar_20171029.vcf.gz]
    --kg            1000 Genomes Annotation Database [default: 1kg/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz]
    --exac          ExAC Annotation Database [default: exac/ExAC.r1.sites.vep.vcf.gz]
    --clinical      Clinical Annotation Flag
    --lrg           LRG Transcript list [default: LRG/LRG_ensembl.txt]

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

// Pipeline version
version = "1.0.1"

// Show help message
params.help = false
if (params.help) {
  helpMessage()
  exit 0
}

// Required params
params.variants = false
if (!params.variants) {
    exit(1, "The required parameter --variants is missing.")
}

// GATK parameters
params.genomeVersion = "b37"
snpeff = fetchReference("snpEff")
params.dbnsfp = false
dbnsfp = file(fetchReference("dbnsfp"))
params.clinvar = false
clinvar = file(fetchReference("clinvar"))
params.kg = false
kg = file(fetchReference("kg"))
params.exac = false
exac = file(fetchReference("exac"))
params.annotations = false
annotations = fetchReference("annotations")
// Clinical
params.clinical = false
params.lrg = false
lrg = file(fetchReference("lrg"))

// Header log info
log.info "================================================="
log.info "Exome-Seq: Best Practice v${version}/Variant Call"
log.info "================================================="
def summary = [:]
summary['Variants'] = params.variants
summary['References'] = ""
summary['dbNSFP'] = dbnsfp
summary['clinvar'] = clinvar
summary['1000G'] = kg
summary['ExAC'] = exac
summary['LRG'] = lrg
summary["Global"] = ""
summary["Project"] = params.project
summary['Output dir'] = params.outdir
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "====================================="

variants = Channel
    .fromFilePairs(params.variants + "{,.idx}")
    .map { it -> [it[0], file(it[1][0]), file(it[1][1])] }    
    .ifEmpty { exit(1, "Cannot find any variant matching: ${params.align}.\n") }

// Step 5.5 SnpEff
process snpeff {
    publishDir "${params.outdir}", mode: "copy",
    saveAs: { fn -> (fn.indexOf("snpEff") > 0) ? "reports/$fn" : "var/$fn" }

    input:
    set val(name), file(var), file(idx) from variants

    output:
    file("*ann.vcf") into annotated_variants
    file("*snpEff.{csv,html}") into snpeff_stats

    script:
    clin = (params.clinical)? "--onlyTr ${lrg}" : ""
    anns = annotations.join(",")
    """
    snpEff ann ${snpeff} \
        -stats ${name}.snpEff.html \
        -csvStats ${name}.snpEff.csv \
        -lof -canon -strict \
        -no-intergenic -noInteraction \
        ${clin} \
        ${var} | \
    snpSift annotate -noId -name KG -info AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF ${kg} - | \
    snpSift annotate -noId -name EXAC -info AF ${exac} - | \
    snpSift annotate -noId -info CLNDN,CLNHGVS,CLNSIG ${clinvar} - | \
    snpSift dbnsfp -f ${anns} -db ${dbnsfp} - > ${name}.ann.vcf
    """
}