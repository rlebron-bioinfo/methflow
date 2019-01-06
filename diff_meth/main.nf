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
} else {
  exit 1, "No input directory specified!"
}

if ( params.comparisons ){
  comparisons = Channel
    .fromPath(params.comparisons, checkIfExists: true)
    .ifEmpty { exit 1, "Comparisons file not found: ${params.comparisons}" }
} else {
  exit 1, "No comparisons file specified!"
}

if ( params.groups ){
  groups = Channel
    .fromPath(params.groups, checkIfExists: true)
    .ifEmpty { exit 1, "Groups file not found: ${params.groups}" }
} else {
  exit 1, "No groups file specified!"
}

if ( params.clusters ){
  if ( params.fasta ){
      fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
  } else {
      exit 1, "No Fasta reference specified! This is required by GenomeCluster."
  }
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
summary['Groups'] = params.groups
summary['Flatten'] = params.flatten ? 'Yes' : 'No'
summary['Clusters'] = params.clusters ? 'Yes' : 'No'
if (params.clusters) summary['Fasta Ref'] = params.fasta
summary['Minimal Coverage'] = params.minCoverage
summary['Minimal Diff Meth'] = params.minDiffMeth
summary['Q-value threshold'] = params.qval
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveIntermediates ? 'Yes' : 'No'
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
 * STEP 0 - Flatten Input Directories
 */

if(params.flatten){
    process flattenInputDirectories {
        publishDir path: { params.saveIntermediates ? "${params.outdir}" : params.outdir },
          saveAs: { params.saveIntermediates ? it : null }, mode: 'copy'

        input:
        file indir from indir.collect()

        output:
        file "flat_input/*" into methylation_profiles

        script:
        if (params.comprehensive) {
          """
          meFlatten --indir $indir --outdir flat_input --comprehensive
          """
        } else {
          """
          meFlatten --indir $indir --outdir flat_input --no-comprehensive
          """
        }
    }
} else {
    process getInputDirectories {
        input:
        file indir from indir.collect()

        output:
        file "flat_input/*" into methylation_profiles

        script:
        """
        mkdir flat_input
        cp `find -L $indir -name \"*CG.output\"` flat_input/
        cp `find -L $indir -name \"*CHG.output\"` flat_input/
        cp `find -L $indir -name \"*CHH.output\"` flat_input/
        """
    }
}
