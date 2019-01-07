#!/usr/bin/env nextflow
/*
=========================================================================
                            METHFLOW PIPELINE
=========================================================================
  DNA Methylation (BS-Seq) Analysis Pipeline. Started December 2018.
  # Homepage / Documentation
    https://github.com/rlebron-bioinfo/methflow
  # Authors
    Ricardo Lebrón <rlebron@go.ugr.es>
-------------------------------------------------------------------------
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */

params.indir = false


// Validate inputs

if( params.indir ){
  indir = Channel
    .fromPath(params.indir, checkIfExists: true)
    .ifEmpty { exit 1, "Input directory not found: ${params.indir}" }
} else {
  exit 1, "No input directory specified!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

log.info """========================================================================
MethFlow - diff_meth: DNA Methylation (BS-Seq) Analysis Pipeline v${params.version}
========================================================================"""

def summary = [:]
summary['Pipeline Name']  = 'MethFlow - tools'
summary['Pipeline Version'] = params.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Input dir'] = params.indir
summary['Flatten'] = params.flatten ? 'Yes' : 'No'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = System.getProperty("user.name")
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================================================"

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented

try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/*
 * STEP 0 - Get Input Directories
 */

process getInputDirectories {
    input:
    file indir from indir.collect()

    output:
    file "flat_input/*" into methylation_profiles

    script:
    if (params.comprehensive) {
      """
      mkdir flat_input
      cp `find -L $indir -name \"*CG.output\"` flat_input/
      cp `find -L $indir -name \"*CHG.output\"` flat_input/
      cp `find -L $indir -name \"*CHH.output\"` flat_input/
      """
    } else {
      """
      mkdir flat_input
      cp `find -L $indir -name \"*CG.output\"` flat_input/
      """
    }
}

/*
 * STEP 1 - Convert to MethylExtract
 */

process convertToMethylExtract {
    publishDir "${params.outdir}/MethylExtract", mode: 'copy'

    input:
    file infile from methylation_profiles.flatten()
    output:
    file "*.output" into methylextract

    script:
    """
    bedMethylToMe $infile ${infile%.bed}.output
    """
}
