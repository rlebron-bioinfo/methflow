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

params.name = false
params.bismark_index = false
params.fasta = false

// Validate inputs
if( params.bismark_index ){
    bismark_index = Channel
        .fromPath(params.bismark_index, checkIfExists: true)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index}" }
}
if ( params.fasta ){
    fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
}
else {
    exit 1, "No Fasta reference specified! This is required by MethylExtract."
}

multiqc_config = file(params.multiqc_config)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Library prep presets
params.rrbs = false
params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.zymo = false
params.cegx = false
if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 9
} else if(params.single_cell){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
} else if(params.epignome){
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
} else if(params.accel || params.zymo){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
} else if(params.cegx){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
} else {
    params.clip_r1 = 0
    params.clip_r2 = 0
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
}

log.info """========================================================================
MethFlow : DNA Methylation (BS-Seq) Analysis Pipeline v${params.version}
========================================================================"""

def summary = [:]
summary['Pipeline Name']  = 'MethFlow'
summary['Pipeline Version'] = params.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
if(params.bismark_index) summary['Bismark Index'] = params.bismark_index
summary['Fasta Ref'] = params.fasta
if(params.rrbs) summary['RRBS Mode'] = 'On'
if(params.relaxMismatches) summary['Mismatch Func'] = "L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)"
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if(params.pbat)         summary['Trim Profile'] = 'PBAT'
if(params.single_cell)  summary['Trim Profile'] = 'Single Cell'
if(params.epignome)     summary['Trim Profile'] = 'TruSeq (EpiGnome)'
if(params.accel)        summary['Trim Profile'] = 'Accel-NGS (Swift)'
if(params.zymo)         summary['Trim Profile'] = 'Zymo Pico-Methyl'
if(params.cegx)         summary['Trim Profile'] = 'CEGX'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Trimming']  = params.notrim ? 'No' : 'Yes'
summary['Deduplication']  = params.nodedup || params.rrbs ? 'No' : 'Yes'
summary['Indel Realignment']  = params.norealign ? 'No' : 'Yes'
summary['M-plot Analysis']  = params.nombias ? 'No' : 'Yes'
summary['Directional Mode'] = params.single_cell || params.zymo || params.non_directional ? 'No' : 'Yes'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
if(params.mindepth) summary['Minimum Depth'] = params.mindepth
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Unmapped']  = params.unmapped ? 'Yes' : 'No'
summary['Save Ambiguous'] = params.ambiguous ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
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

fasta.into { fasta_1; fasta_2; fasta_3 }

/*
 * STEP 0 - Get FASTQ files
 */

if(params.reads){
    Channel
    .fromFilePairs( params.reads, checkIfExists: true, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }
} else {
    exit 1, "Cannot find any FASTQ file. Please use --reads argument."
}

/*
 * STEP 1 - Build Bismark index
 */

if(!params.bismark_index && params.fasta){
    process makeBismarkIndex {
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_1

        output:
        file "BismarkIndex" into bismark_index

        script:
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation BismarkIndex
        """
    }
}

bismark_index.into { bismark_index_1; bismark_index_2 }

/*
 * STEP 2 - FastQC
 */

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 3 - Trim Galore!
 */

if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

/*
 * STEP 4 - Align with Bismark
 */

process bismark_align {
    tag "$name"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
        saveAs: {filename ->
            if (params.unmapped && filename.indexOf("_unmapped_reads_") > 0) "unmapped/$filename"
            else if (params.ambiguous && filename.indexOf("_ambiguous_reads_") > 0) "ambiguous/$filename"
            else if (filename.indexOf(".fq.gz") == -1 && filename.indexOf(".bam") == -1) "logs/$filename"
            else params.saveAlignedIntermediates ? filename : null
        }

    input:
    set val(name), file(reads) from trimmed_reads
    file index from bismark_index_1.collect()

    output:
    file "*.bam" into bam_aligned_1, bam_aligned_2 
    file "*report.txt" into bismark_align_log_1, bismark_align_log_2, bismark_align_log_3  
    file "*unmapped_reads*" into bismark_unmapped
    file "*ambiguous_reads*" into bismark_ambiguous

    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    ambiguous = params.ambiguous ? "--ambiguous" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    multicore = ''
    if (task.cpus){
        // Numbers based on recommendation for human genome
        if(params.single_cell || params.zymo || params.non_directional){
            cpu_per_multicore = 5
            mem_per_multicore = (18.GB).toBytes()
        } else {
            cpu_per_multicore = 3
            mem_per_multicore = (13.GB).toBytes()
        }
        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore) as int
        // Check that we have enough memory, assuming 13GB memory per instance
        try {
            tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
            mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.debug "Not able to define bismark align multicore based on available memory"
        }
        if(ccore > 1){
            multicore = "--multicore $ccore"
        }
    }
    if (params.singleEnd) {
        """
        bismark \\
            --bam $pbat $non_directional $unmapped $ambiguous $mismatches $multicore \\
            --genome $index \\
            $reads
        """
    } else {
        """
        bismark \\
            --bam $pbat $non_directional $unmapped $ambiguous $mismatches $multicore \\
            --genome $index \\
            -1 ${reads[0]} \\
            -2 ${reads[1]}
        """
    }
}

/*
 * STEP 5 - Bismark Sample Report
 */

process bismark_report {
    tag "$name"
    publishDir "${params.outdir}/bismark_reports", mode: 'copy'

    input:
    file bismark_align_log_1

    output:
    file '*{html,txt}' into bismark_reports_results

    script:
    name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
    """
    bismark2report --alignment_report $bismark_align_log_1
    """
}

/*
 * STEP 6 - Bismark Summary Report
 */

process bismark_summary {
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'

    input:
    file ('*') from bam_aligned_1.collect()
    file ('*') from bismark_align_log_2.collect()

    output:
    file '*{html,txt}' into bismark_summary_results

    script:
    """
    bismark2summary
    """
}

/*
 * STEP 7 - Sort by coordinates with samtools
 */

process samtools_sort_by_coordinates {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/alignments_sorted_by_coordinates", mode: 'copy',
        saveAs: {filename ->
            if (params.saveAlignedIntermediates) filename
            else null
        }

    input:
    file bam from bam_aligned_2

    output:
    file "${bam.baseName}.sorted.bam" into bam_sorted
    file "${bam.baseName}.sorted.bam.bai" into bam_sorted_index
    file "${bam.baseName}_flagstat.txt" into flagstat_results
    file "${bam.baseName}_stats.txt" into samtools_stats_results

    script:
    """
    samtools sort \\
        -m ${task.memory.toBytes() / task.cpus} \\
        -@ ${task.cpus} \\
        $bam \\
        > ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    samtools flagstat ${bam.baseName}.sorted.bam > ${bam.baseName}_flagstat.txt
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}_stats.txt
    """
}

/*
 * STEP 8 - Merge with samtools
 */

process samtools_merge {
    if (params.norealign && (params.nodedup || params.rrbs)) {
        publishDir "${params.outdir}/merged_alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") > 0) "alignments/$filename"
                else null
            }
    } else {
        publishDir "${params.outdir}/merged_alignments", mode: 'copy',
            saveAs: {filename ->
                if (params.saveAlignedIntermediates) filename
                else null
            }
    }

    input:
    file bam from bam_sorted.collect()
    file bam_index from bam_sorted_index.collect()

    output:
    file "${params.name}.merged.bam" into bam_merged
    file "${params.name}.merged.bam.bai" into bam_merged_index

    script:
    """
    samtools merge -n ${params.name}.tmp.bam $bam
    samtools sort \\
        -m ${task.memory.toBytes() / task.cpus} \\
        -@ ${task.cpus} \\
        ${params.name}.tmp.bam \\
        > ${params.name}.merged.bam
    rm ${params.name}.tmp.bam
    samtools index ${params.name}.merged.bam
    """
}

/*
 * STEP 9 - Remove duplicates
 */

if (params.nodedup || params.rrbs) {
    bam_merged.set { bam_dedup }
    bam_merged_index.set { bam_dedup_index }
    dedup_metrics = Channel.from(false)
} else {
    process removeDuplicates {
        if (params.norealign) {
            publishDir "${params.outdir}/dedup_alignments", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else if (filename.indexOf(".bam") > 0) "alignments/$filename"
                    else params.saveAlignedIntermediates ? filename : null
                }
        } else {
            publishDir "${params.outdir}/dedup_alignments", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else params.saveAlignedIntermediates ? filename : null
                }
        }

        input:
        file bam from bam_merged
        file bam_index from bam_merged_index

        output:
        file "${bam.baseName}.dedup.bam" into bam_dedup
        file "${bam.baseName}.dedup.bam.bai" into bam_dedup_index
        file "${bam.baseName}.dedup_metrics.txt" into dedup_metrics

        script:
        """
        picard MarkDuplicates \\
            INPUT=$bam \\
            OUTPUT=${bam.baseName}.dedup.bam \\
            METRICS_FILE=${bam.baseName}.dedup_metrics.txt \\
            REMOVE_DUPLICATES=true \\
            ASSUME_SORTED=true \\
            PROGRAM_RECORD_ID='null' \\
            VALIDATION_STRINGENCY=LENIENT
        samtools index ${bam.baseName}.dedup.bam
        """
    }
}

/*
 * STEP 10 - Local indel realignment
 */

if (params.norealign) {
    bam_dedup.into { bam_realign_1; bam_realign_2 }
    bam_dedup_index.into { bam_realign_index_1; bam_realign_index_2 }
} else {
    process indelRealign {
        publishDir "${params.outdir}/realign_alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") > 0) "alignments/$filename"
                else null
            }

        input:
        file bam from bam_dedup
        file bam_index from bam_dedup_index
        file genome from bismark_index_2

        output:
        file "${bam.baseName}.realign.bam" into bam_realign_1, bam_realign_2
        file "${bam.baseName}.realign.bam.bai" into bam_realign_index_1, bam_realign_index_2

        script:
        """
        M-IndelRealigner \\
            -R $genome \\
            -I $bam \\
            -T temp \\
            -O ${bam.baseName}.realign.bam \\
            -C ${task.cpus} \\
        """
    }
}

/*
 * STEP 11 - Qualimap
 */

process qualimap {
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    file bam from bam_realign_1
    file bam_index from bam_realign_index_1

    output:
    file "${bam.baseName}_qualimap" into qualimap_results

    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    qualimap bamqc \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
 * STEP 12 - Sort by qname with samtools
 */

process samtools_sort_by_qname {
    publishDir "${params.outdir}/alignments_sorted_by_qname", mode: 'copy',
        saveAs: {filename ->
            if (params.saveAlignedIntermediates) filename
            else null
        }

    input:
    file bam from bam_realign_2
    file bam_index from bam_realign_index_2

    output:
    file "${bam.baseName}.sorted.bam" into bam_final_1, bam_final_2

    script:
    """
    samtools sort \\
        -n \\
        -m ${task.memory.toBytes() / task.cpus} \\
        -@ ${task.cpus} \\
        $bam \\
        > ${bam.baseName}.sorted.bam
    """
}

/*
 * STEP 13 - Bismark M-bias analysis
 */

if (params.nombias) {
    bismark_splitting_report = Channel.from(false)
    bismark_mbias = Channel.from(false)
    bismark_methXtract_results = Channel.from(false)
} else {
    process mbias_analysis {
        publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("splitting_report.txt") > 0) "logs/$filename"
                else if (filename.indexOf("M-bias") > 0) "m-bias/$filename"
                else if (filename.indexOf(".cov") > 0) "methylation_coverage/$filename"
                else if (filename.indexOf("bedGraph") > 0) "bedGraph/$filename"
                else "methylation_calls/$filename"
            }

        input:
        file bam from bam_final_1

        output:
        file "${bam.baseName}_splitting_report.txt" into bismark_splitting_report
        file "${bam.baseName}.M-bias.txt" into bismark_mbias
        file '*.{png,gz}' into bismark_methXtract_results

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        multicore = ''
        if (task.cpus){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 10) as int
            if(ccore > 1){
                multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if (task.memory){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // only set if we have more than 6GB available
            if(mbuffer.compareTo(4.GB) == 1){
                buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }
        if (params.singleEnd) {
            """
            bismark_methylation_extractor $comprehensive \\
                $multicore $buffer \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -s \\
                --report \\
                $bam
            """
        } else {
            """
            bismark_methylation_extractor $comprehensive \\
                $multicore $buffer \\
                --ignore_r2 2 \\
                --ignore_3prime_r2 2 \\
                --bedGraph \\
                --counts \\
                --gzip \\
                -p \\
                --no_overlap \\
                --report \\
                $bam
            """
        }
    }
}

/*
 * STEP 14 - MethylExtract
 */

process MethylExtract {
    publishDir "${params.outdir}/MethylExtract", mode: 'copy'

    input:
    file bam from bam_final_2
    file fasta from fasta_3

    output:
    file "*{output,vcf}" into methylextract_results
    file "*.log" into methylextract_log

    script:
    if (params.singleEnd) {
        flagW = '0'
        flagC = '16'
    } else {
        flagW = '99,147'
        flagC = '83,163'
    }
    memNumReads = 100000 * (task.memory.toGiga() / task.cpus).round()
    if (params.comprehensive) {
        context = 'ALL'
    } else {
        context = 'CG'
    }
    """
    MethylExtract.pl \\
        seq=$fasta \\
        inDir=. \\
        flagW=$flagW \\
        flagC=$flagC \\
        qscore=phred33-quals \\
        delDup=N \\
        peOverlap=N \\
        minDepthMeth=${params.mindepth} \\
        minDepthSNV=${params.mindepth} \\
        minQ=20 \\
        methNonCpGs=0.9 \\
        varFraction=0.1 \\
        maxPval=0.05 \\
        maxStrandBias=0.7 \\
        p=${task.cpus} \\
        memNumReads=$memNumReads \\
        context=$context \\
        outDir=. \\
        &> process.log || echo "ERROR" 
    """
}

/*
 * STEP 15 - Get software versions
 */

process get_software_versions {

    output:
    file 'software_versions_mqc.yml' into software_versions_yml

    script:
    """
    software_versions --module main &> software_versions_mqc.yml
    """
}

/*
 * STEP 16 - MultiQC
 */

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bismark/*') from bismark_align_log_3.toList()
    file ('bismark/*') from bismark_reports_results.toList()
    file ('bismark/*') from bismark_summary_results.toList()
    file ('samtools/*') from flagstat_results.flatten().toList()
    file ('samtools/*') from samtools_stats_results.flatten().toList()
    file ('picard/*') from dedup_metrics.flatten().toList()
    file ('qualimap/*') from qualimap_results.toList()
    file ('bismark/*') from bismark_splitting_report.toList()
    file ('bismark/*') from bismark_mbias.toList()
    file ('methylextract/*') from methylextract_log.toList()
    file ('software_versions/*') from software_versions_yml.toList()

    output:
    file "*_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr

    script:
    rtitle = custom_runName ? "--title \"MethFlow - $custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}
