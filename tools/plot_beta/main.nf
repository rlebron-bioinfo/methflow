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
summary['Comparisons'] = params.comparisons
summary['Groups'] = params.groups
summary['Flatten'] = params.flatten ? 'Yes' : 'No'
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
}

/*
 * STEP 1 - Sort methylation profiles
 */

  process sortMethylationProfile {
      publishDir path: { params.saveIntermediates ? "${params.outdir}/sorted_methylation_profiles" : params.outdir },
        saveAs: { params.saveIntermediates ? it : null }, mode: 'copy'

      input:
      file infile from methylation_profiles.flatten()

      output:
      file "*.sorted" into sorted_methylation_profiles

      script:
      """
      meSort --input $infile --output ${infile}.sorted --file
      """
  }

/*
 * STEP 2 - Convert to methylKit
 */

  process convertToMethylKit {
      publishDir path: { params.saveIntermediates ? "${params.outdir}/methylKit_profiles" : params.outdir },
        saveAs: { params.saveIntermediates ? it : null }, mode: 'copy'

      input:
      file infile from sorted_methylation_profiles.flatten()

      output:
      file "*.mk" into methylkit_profiles

      script:
      """
      #!/usr/bin/env python3

      import os, shlex
      from subprocess import Popen, PIPE

      infile = \"${infile}\"
      outfile = infile.replace(\".output.sorted\", \".mk\")

      if infile.endswith(\"CG.output.sorted\"):
        cmd = f\"meToMethylKit --infile {infile} --outfile {outfile} --context CG --destrand\"
      elif infile.endswith(\"CHG.output.sorted\"):
        cmd = f\"meToMethylKit --infile {infile} --outfile {outfile} --context CHG\"
      elif infile.endswith(\"CHH.output.sorted\"):
        cmd = f\"meToMethylKit --infile {infile} --outfile {outfile} --context CHH\"
      else:
        raise Exception(\"Unknown context!\")

      p = Popen(shlex.split(cmd), stdin=PIPE, stdout=PIPE, stderr=PIPE)
      out, err = p.communicate()
      if p.returncode != 0:
        raise Exception(\"Conversion failed!\")
      """
  }

/*
 * STEP 3 - Plot beta
 */

  process plotBeta {
      publishDir "${params.outdir}/beta_plots", mode: 'copy'

      input:
      file profiles from methylkit_profiles.flatten()

      output:
      file "*.pdf" into plots

      script:
      """
      for F in *.mk
      do
        context=${F%.mk}
        context=${F##*.}
        mkBetaPlot $F $context
      done
      """ 
  }
