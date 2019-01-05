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
params.comparisons = false
params.groups = false


// Validate inputs
if( params.indir ){
    indir = Channel
        .fromPath(params.indir, checkIfExists: true)
        .ifEmpty { exit 1, "Input directory not found: ${params.indir}" }
        .into { indir_1; indir_2 }
} else {
  exit 1, "No input directory specified!"
}
if ( params.comparisons ){
    comparisons = Channel
        .fromPath(params.comparisons, checkIfExists: true)
        .ifEmpty { exit 1, "Comparisons file not found: ${params.comparisons}" }
        .into { comparisons_1; comparisons_2 }
}
else {
  exit 1, "No comparisons file specified!"
}
if ( params.groups ){
    groups = Channel
        .fromPath(params.groups, checkIfExists: true)
        .ifEmpty { exit 1, "Groups file not found: ${params.groups}" }
        .into { groups_1; groups_2 }
} else {
    groups = Channel.from(false).into{ groups_1; groups_2 }
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
summary['Pipeline Name']  = 'MethFlow - diff_meth'
summary['Pipeline Version'] = params.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Input dir'] = params.indir
summary['Comparisons'] = params.comparisons
if (params.groups) summary['Groups'] = params.groups
summary['Clusters'] = params.clusters ? 'Yes' : 'No'
summary['Minimal Diff Meth'] = params.minDiffMeth
summary['Q-value threshold'] = params.qval
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

process validateInput {    
    input:
    file comparisons from comparisons_1
    file groups from groups_1
    file indir from indir_1

    output:
    file "${comparisons.baseName}.valid" into valid_comparisons_1, valid_comparisons_2

    script:
    if (params.comprehensive) {
      comprehensive = "--comprehensive"
    } else {
      comprehensive = ""
    }
    if (params.groups) {
      """
      check_comparisons --ignore-diagonal $comprehensive \
      -c $comparisons -g $groups -i $indir -o ${comparisons.baseName}.valid
      """
    } else{
      """
      check_comparisons --ignore-diagonal $comprehensive \
      -c $comparisons -i $indir -o ${comparisons.baseName}.valid
      """
    }
}
