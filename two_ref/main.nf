#!/usr/bin/env nextflow
/*
=========================================================================
                        METHFLOW PIPELINE - TWO REF
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
params.bismark_index_1 = false
params.fasta_1 = false
params.bismark_index_2 = false
params.fasta_2 = false

// First ref
if( params.bismark_index_1 ){
    bismark_index_1 = Channel
        .fromPath(params.bismark_index_1)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index_1}" }
}
if ( params.fasta_1 ){
    fasta_1 = Channel
        .fromPath(params.fasta_1)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta_1}" }
}
else {
    exit 1, "No First Fasta reference specified!"
}

// Second ref
if( params.bismark_index_2 ){
    bismark_index_2 = Channel
        .fromPath(params.bismark_index_2)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index_2}" }
}
if ( params.fasta_2 ){
    fasta_2 = Channel
        .fromPath(params.fasta_2)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta_2}" }
}
else {
    exit 1, "No Second Fasta reference specified!"
}
if (!params.use_unmapped && !params.use_ambiguous) {
    exit 1, "Please select unmapped (--use_unmapped), ambiguous (--use_ambiguous) or both!"
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
MethFlow.two_ref : DNA Methylation (BS-Seq) Analysis Pipeline v${params.version}
========================================================================"""

def summary = [:]
summary['Pipeline Name']  = 'MethFlow.two_ref'
summary['Pipeline Version'] = params.version
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
if(params.bismark_index_1) summary['First Bismark Index'] = params.bismark_index_1
else if(params.fasta_1)    summary['First Fasta Ref'] = params.fasta_1
if(params.bismark_index_2) summary['Second Bismark Index'] = params.bismark_index_2
else if(params.fasta_2)    summary['Second Fasta Ref'] = params.fasta_2
summary['Use Unmapped']  = params.use_unmapped ? 'Yes' : 'No'
summary['Use Ambiguous'] = params.use_ambiguous ? 'Yes' : 'No'
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
summary['Save First Reference'] = params.saveReference1 ? 'Yes' : 'No'
summary['Save Second Reference'] = params.saveReference2 ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
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

fasta_1.into { fasta_1A; fasta_1B; fasta_1C }
fasta_2.into { fasta_2A; fasta_2B; fasta_2C }

/*
 * STEP 0 - Get FASTQ files
 */

if(params.reads){
    Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }
} else {
    exit 1, "Cannot find any FASTQ file. Please use --reads argument."
}

/*
 * STEP 1A - Build Bismark index 1
 */

if(!params.bismark_index_1 && params.fasta_1){
    process makeBismarkIndex_1 {
        publishDir path: { params.saveReference1 ? "${params.outdir}/first_reference_genome" : params.outdir },
                   saveAs: { params.saveReference1 ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_1A

        output:
        file "BismarkIndex" into bismark_index_1

        script:
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation BismarkIndex
        """
    }
}

/*
 * STEP 1B - Build Bismark index 2
 */

if(!params.bismark_index_2 && params.fasta_2){
    process makeBismarkIndex_2 {
        publishDir path: { params.saveReference2 ? "${params.outdir}/second_reference_genome" : params.outdir },
                   saveAs: { params.saveReference2 ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_2A

        output:
        file "BismarkIndex" into bismark_index_2

        script:
        """
        mkdir BismarkIndex
        cp $fasta BismarkIndex/
        bismark_genome_preparation BismarkIndex
        """
    }
}

bismark_index_1.into { bismark_index_1A; bismark_index_1B }
bismark_index_2.into { bismark_index_2A; bismark_index_2B }

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
 * STEP 4A - First Align with Bismark
 */

process bismark_align_1 {
    tag "$name"
    publishDir "${params.outdir}/first_bismark_alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_unmapped_reads_") > 0) "unmapped/$filename"
            else if (filename.indexOf("_ambiguous_reads_") > 0) "ambiguous/$filename"
            else if (filename.indexOf(".fq.gz") == -1 && filename.indexOf(".bam") == -1) "logs/$filename"
            else params.saveAlignedIntermediates ? filename : null
        }

    input:
    set val(name), file(reads) from trimmed_reads
    file index from bismark_index_1A.collect()

    output:
    file "*.bam" into bam_aligned_1A, bam_aligned_1B 
    file "*report.txt" into bismark_align_log_1A, bismark_align_log_1B, bismark_align_log_1C  

    file "*unmapped_reads_{1,2}*" into bismark_unmapped
    file "*ambiguous_reads_{1,2}*" into bismark_ambiguous

    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = "--unmapped"
    ambiguous = "--ambiguous"
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

Channel
.fromFilePairs( bismark_unmapped, size: params.singleEnd ? 1 : 2 )
.set { bismark_unmapped_A }

Channel
.fromFilePairs( bismark_ambiguous, size: params.singleEnd ? 1 : 2 )
.set { bismark_ambiguous_A }

/*
 * STEP 4B - Second Align with Bismark
 */

process bismark_align_2 {
    tag "$name"
    publishDir "${params.outdir}/second_bismark_alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_unmapped_reads_") > 0) "unmapped/$filename"
            else if (filename.indexOf("_ambiguous_reads_") > 0) "ambiguous/$filename"
            else if (filename.indexOf(".fq.gz") == -1 && filename.indexOf(".bam") == -1) "logs/$filename"
            else params.saveAlignedIntermediates ? filename : null
        }

    input:
    set val(name), file(unmapped_reads) from bismark_unmapped_A
    set val(name), file(ambiguous_reads) from bismark_ambiguous_A
    file index from bismark_index_2A.collect()

    output:
    file "*.bam" into bam_aligned_2A, bam_aligned_2B 
    file "*report.txt" into bismark_align_log_2A, bismark_align_log_2B, bismark_align_log_2C  

    file "*unmapped_reads_{1,2}*" into bismark_unmapped_B
    file "*ambiguous_reads_{1,2}*" into bismark_ambiguous_B

    script:
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = "--unmapped"
    ambiguous = "--ambiguous"
    if (!params.use_unmapped) {
        unmapped_reads = Channel.from(false)
    }
    if (!params.use_ambiguous) {
        ambiguous_reads = Channel.from(false)
    }
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
            $unmapped_reads $ambiguous_reads
        """
    } else {
        """
        bismark \\
            --bam $pbat $non_directional $unmapped $ambiguous $mismatches $multicore \\
            --genome $index \\
            -1 ${unmapped_reads[0]} ${ambiguous_reads[0]} \\
            -2 ${unmapped_reads[1]} ${ambiguous_reads[0]}
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
    file bismark_align_log_1A
    file bismark_align_log_2A

    output:
    file '*{html,txt}' into bismark_reports_results

    script:
    name = bismark_align_log_1.toString() - ~/(_R1)?(_trimmed|_val_1).+$/
    """
    bismark2report --alignment_report $bismark_align_log_1A $bismark_align_log_2A
    """
}

/*
 * STEP 6 - Bismark Summary Report
 */

process bismark_summary {
    tag "${params.name}"
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'

    input:
    file ('*') from bam_aligned_1A.collect()
    file ('*') from bam_aligned_2A.collect()
    file ('*') from bismark_align_log_1B.collect()
    file ('*') from bismark_align_log_2B.collect()

    output:
    file '*{html,txt}' into bismark_summary_results

    script:
    """
    bismark2summary
    """
}

/*
 * STEP 7 - Sort with samtools
 */

process samtools_sort {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/sorted_alignments", mode: 'copy',
        saveAs: {filename ->
            if (params.saveAlignedIntermediates) filename
            else null
        }

    input:
    file bam1 from bam_aligned_1B
    file bam2 from bam_aligned_2B

    output:
    file "${bam.baseName}.sorted.bam" into bam_sorted
    file "${bam.baseName}.sorted.bam.bai" into bam_sorted_index
    file "${bam.baseName}_flagstat.txt" into flagstat_results
    file "${bam.baseName}_stats.txt" into samtools_stats_results

    script:
    """
    samtools sort \\
        -n \\
        -m ${task.memory.toBytes() / task.cpus} \\
        -@ ${task.cpus} \\
        $bam1 $bam2 \\
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
    tag "${params.name}"
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
    samtools merge -n ${params.name}.merged.bam $bam
    samtools index ${params.name}.merged.bam
    """
}

/*
 * STEP 9 - Remove duplicates
 */

if (params.nodedup || params.rrbs) {
    bam_merged.into { bam_dedup }
    bam_merged_index.set { bam_dedup_index }
    dedup_metrics = Channel.from(false)
} else {
    process removeDuplicates {
        tag "${bam.baseName}"
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
 * STEP 10A - Collapse FASTA files
 */

process collapseFasta {
    tag "$name"
    publishDir path: { params.saveReference2 ? "${params.outdir}/collapsed_reference_genome" : params.outdir },
        saveAs: { params.saveReference2 ? it : null }, mode: 'copy'

    input:
    file fasta1 from fasta_1C
    file fasta2 from fasta_2C

    output:
    file "${name}.fa" into collapsed_fasta
    file "BismarkIndex" into collapsed_fasta_index

    script:
    """
    cat $fasta1 $fasta2 | fasta_formatter | fasta_formatter -t \\
    | sort -k1,1V | uniq \\
    | sed 's/^/>/g' | sed -e 's/[\t]/\n/g' \\
    fasta_formatter -w 70 > collapsed.fa
    mkdir BismarkIndex
    cp collapsed.fa BismarkIndex/
    bismark_genome_preparation BismarkIndex
    """
}

/*
 * STEP 10B - Local indel realignment
 */

if (params.norealign) {
    bam_dedup.set { bam_final }
    bam_dedup_index.set { bam_final_index }
} else {
    process indelRealign {
        tag "${bam.baseName}"
        publishDir "${params.outdir}/realign_alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") > 0) "alignments/$filename"
                else null
            }

        input:
        file bam from bam_dedup
        file bam_index from bam_dedup_index
        file genome from collapsed_fasta_index

        output:
        file "${bam.baseName}.realign.bam" into bam_final
        file "${bam.baseName}.realign.bam.bai" into bam_final_index

        script:
        """
        M-IndelRealigner \\
            -R $genome \\
            -I $bam\\
            -T temp \\
            -O ${bam.baseName}.realign.bam\\
            -C ${task.cpus} \\
        """
    }
}

bam_final.into { bam_final_1; bam_final_2; bam_final_3 }
bam_final_index.into { bam_final_index_1; bam_final_index_2; bam_final_index_3 }

/*
 * STEP 11 - Qualimap
 */

process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    file bam from bam_final_1
    file bam_index from bam_final_index_1

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
 * STEP 12 - Bismark M-bias analysis
 */

if (params.nombias) {
    bismark_splitting_report = Channel.from(false)
    bismark_mbias = Channel.from(false)
    bismark_methXtract_results = Channel.from(false)
} else {
    process mbias_analysis {
        tag "${bam.baseName}"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("splitting_report.txt") > 0) "logs/$filename"
                else if (filename.indexOf("M-bias") > 0) "m-bias/$filename"
                else if (filename.indexOf(".cov") > 0) "methylation_coverage/$filename"
                else if (filename.indexOf("bedGraph") > 0) "bedGraph/$filename"
                else "methylation_calls/$filename"
            }

        input:
        file bam from bam_final_2
        file bam_index from bam_final_index_2

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
 * STEP 13 - MethylExtract ===============================>>>>>>>>>>>>> *************
 */

process MethylExtract {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/MethylExtract", mode: 'copy'

    input:
    file bam from bam_final_3
    file bam_index from bam_final_index_3
    file fasta from collapsed_fasta

    output:
    file "*{output,vcf}" into methylextract_results
    file "*{log}" into methylextract_log

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
        outDir=.
    """
}

/*
 * STEP 13 - Get software versions
 */

process get_software_versions {

    output:
    file 'software_versions_mqc.yml' into software_versions_yml

    script:
    """
    software_versions &> software_versions_mqc.yml
    """
}

/*
 * STEP 14 - MultiQC
 */

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bismark/*') from bismark_align_log_1C.toList()
    file ('bismark/*') from bismark_align_log_2C.toList()
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
    rtitle = custom_runName ? "--title \"MethFlow.two_ref - $custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}
