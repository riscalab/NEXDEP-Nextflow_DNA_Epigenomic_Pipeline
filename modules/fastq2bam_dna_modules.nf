
process fastp_SE_adapter_known {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir "${params.base_out_dir}/fastp_qc_single_end", mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir "${params.base_out_dir}/fastp_qc_single_end/html_reports", mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel
    val(adapter_seq)


    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --adapter_sequence "${adapter_seq}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}







// This next process will run some qc to look at the fastq files and trim adapters from the single end reads
process fastp_SE {
    // using this conda yml file that I created. Nextflow will now make its own conda environment from the dependencies found within.
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastp_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir "${params.base_out_dir}/fastp_qc_single_end", mode: 'copy', pattern:'*_fp_filt.fastq'
    publishDir "${params.base_out_dir}/fastp_qc_single_end/html_reports", mode: 'copy', pattern:'*.html'

    input:
    // input names dont have to be the exact same as what is seen in the workflow section of this script
    // as long and they are in the same order, you will have the correct input
    path(fastq_files) // the files
    val(fastq_names) // the names of the files as a value input channel



    output:
    path("${out_name}"), emit: filtered_fastqs
    path("${fastq_names}*.html"), emit: fastp_html_reports


    script:

    // getting the output name
    out_name = "${fastq_names}_fp_filt.fastq"


    """
    #!/usr/bin/env bash

    # ASK HERA ABOUT THE ADAPTERS USED AND FIND THEM FROM THE ILLUMINA KIT

    # remember this process is for the single end reads so only one read in and one filtered qc read out
    # adapter trimming is enabled by defualt for single end but use --detect_adapter_for_pe to enable it in paired end data
    # you can specify the adapter sequence if you have them. look at documentation here: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp
    
    # --adapter_sequence: string sequence that represents the adapter used. in the illumina NEBNext libraries adapter read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA. Adapter read 2 is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    # --dedup enables the deduplication to drop the duplicated reads
    # --dup_calc_accuracy: defualt is 3 but can go from 1~6. higher levels use more memory
    # --trim_poly_g: novaSeq data polyG can happen in read tails since G means no signal in the illumina two-color system. fastp can detect and trim them.
    # --qualified_quality_phred: the quality value that a base is qualified defualt is 15
    # --unqualified_percent_limit: how many percent of bases are allowed to be unqualified defualt is 40
    # --umi: not using this yet have to use --umi_loc also and i dont know that, fastp can extract the unique molecular identifiers and append them to the first part of the read names, so the umi's will be present in SAM/BAM records If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.
    # --overrepresentation_analysis and --overrepresentation_sampling see fastp documentation
    # --html: this is how the report file will be made

    # their fastq2bam single end pipeline used trim_galore tool for trimming adapters.
    # I want to implement a way to choose to let fastp do this by default or use a known sequence that the user can input.
    ###### can add or remove more options when needed #####
    
    fastp \
    --in1 "${fastq_files}" \
    --out1 "${out_name}" \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --html "${fastq_names}_fastp.html"






    """

}

// this next process is for fastqc tool
// dont forget about multiqc m
process fastqc_SE {
    // using the conda environment 
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'
    publishDir "${params.base_out_dir}/fastqc_htmls", mode: 'copy', pattern: '*.html'


    input:

    path(fastq_filt_files)
    val(fastq_filt_names)


    output:
    path("*.html"), emit: fastqc_htmls
    path("*.zip"), emit: fastqc_zip_files

    script:
    out_name = "${fastq_filt_files}"

    """
    #!/usr/bin/env bash

    # I need to add the adapter sequences later, ask hera.

    fastqc "${fastq_filt_files}"
    



    """
}

process multiqc_SE {
    // this yml file doesnt work
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/multiqc_rj'
    
    publishDir "${params.base_out_dir}/multiQC_collection", mode: 'copy', pattern: '*.html'

    input:
    path(fastp_filt_html)



    output:

    path("*.html"), emit: multiqc_html_collection

    script:


    """
    #!/usr/bin/env bash

    # I think i jsut have to pass the html files to multiqc for it to complie them
    multiqc . \
    --interactive \
    --profile-runtime \
    --title "Single-end QC"




    """

}


// Creating two processes that will index the reference genome

process bwa_index_genome {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    
    publishDir './genome_index_bwa', mode: 'copy', pattern: '*'


    input:
    path(ref_genome)
    


    output:

    path("*"), emit: genome_index_files


    script:
    // not getting the basename because bwa expects the index files will have the exact file name of the genome but with an .{ext} on the end. 
    // example genome.fa, genome.fa.amb, genome.fa.ann
    // not genome.fa, genome.amb, genome.ann
    //genome_file_name = "${ref_genome.baseName}"

    """
    #!/usr/bin/env bash

    ############### parameters used ###############
    # -p: a string representing the prefix of the output database [same as the db filename] so i think just the base name of the genome file
    # -a: a string, choosing the algorithm to construct the BWT index. either is or bwtsw. I'll use bwtsw


    bwa index \
    -p "${ref_genome}" \
    -a bwtsw \
    "${ref_genome}"


    """

}


// creating a process that will align the reads to the genome. i will take in the reference genome, the index files, the filtered fastq's and their names

process bwa_align_SE {
    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir "${params.base_out_dir}/bwa_outputs_singleEnd_SAM", mode: 'copy', pattern: '*.sam'
    publishDir "${params.base_out_dir}/sai_alignment_files", mode: 'copy', pattern: '*.sai'


    input:
    path(ref_genome)
    path(genome_index_files)
    path(fastq_filt_files)
    val(fastq_filt_names)


    output:

    path("*.sam"), emit: sam_se_files
    path("*.sai"), emit: sai_align_files


    script:
    sai_output_file = "${fastq_filt_names}_out.sai"
    sam_name = "${fastq_filt_names}.sam"

    """
    #!/usr/bin/env bash

    ############# Parameters used ############
    # first i need to get the sai file by using bwa aln. This gives the SA coordinates of the input reads
    # -t (nThrds) number of threads for multi threading mode. defualt is 1

    # using bwa samse : this will generate alignments in the SAM format given single-end reads
    # the two parameters may be used in the future.
    # -n: takes an integer. max number of alignments to output in the XA tag for reads paired properly
    # -r: takes a string. specify the read group
    #
    ##########################################

    ls .

    bwa aln \
    -t 20 \
    "${ref_genome}"  \
    "${fastq_filt_files}" \
    > "${sai_output_file}"


    bwa samse \
    "${ref_genome}" \
    "${sai_output_file}" \
    "${fastq_filt_files}" \
    > "${sam_name}"


    """
}




process samtools_sort {
    // using the conda yml file for samtools
    // it doesnt work
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'

    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj' // this was samtools version 1.3 which doesnt have samtools fixmate option -m

    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools-1.21_spec_env_rj.txt' // that is the explicit file but if that doesnt work try the yml file samtools-1.21_env_rj.yml; and if that doesnt work use the path to the 1.21 environment

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools-1.21_rj'


    
    if (!params.BL) {
        publishDir "${params.base_out_dir}/sorted_bam_files", mode: 'copy', pattern: '*_sorted.bam'
        publishDir "${params.base_out_dir}/sorted_bam_files", mode: 'copy', pattern: '*.{bai, csi}'

    }else {

        

    }
    // publishDir "${params.base_out_dir}/flag_stat_log", mode: 'copy', pattern: '*stats.log'
    // publishDir "${params.base_out_dir}/stats_tsv_files", mode: 'copy', pattern: '*stats.tsv'

    

    
    //publishDir './sorted_bam_files', mode: 'copy', pattern: '*_sorted.bam'
    //publishDir './indexed_bam_files', mode: 'copy', pattern: '*.{bai, csi}'
    //publishDir './flag_stat_log', mode: 'copy', pattern: '*stat.log'

    input:
    path(sam_files)


    output:

    path("*_sorted.bam"), emit: sorted_bams
    //tuple path("*.{bai,csi}"), emit: indexed_bams
    //tuple path("*.bai"), path("*.csi"), emit: indexed_bams
    path("*.bai"), emit: indexed_bams
    tuple path("*_sorted.bam"), path("*.bai"), emit: bam_index_tuple

    // path("*stats.log"), emit: flag_stats_log
    // path("*stats.txt"), emit: norm_stats_txt
    // path("*stats.tsv"), emit: tsv_SN_stats

    script:

    // i will start using baseName inside the process since its easier to keep track of different names an uses less inputs into a process
    out_bam_filt = "${sam_files.baseName}_bam_filt.bam"
    out_bam_name_sort = "${sam_files.baseName}_name_ordered.bam"
    out_bam_coor_sort = "${sam_files.baseName}_filt_coor_sorted.bam"
    out_bam_fixmate = "${sam_files.baseName}_fixmate.bam"
    // out_bam_final = "${sam_files.baseName}_markdup_filt_coor_sorted.bam"
    // flagstats_log = "${sam_files.baseName}_flag_stats.log"
    // samtools_stats_log = "${sam_files.baseName}_stats.txt"
    // tsv_file_with_stats = "${sam_files.baseName}_SN_stats.tsv"


    

    """
    #!/usr/bin/env bash

    ################# samtools parameters used ################
    # for samtools view
    # --min-MQ or -q : takes an INT and will skip alignments with a MAPQ smaller than INT
    # --bam or -b : output in the bam format
    # this version of bwa didnt recognize --bam or --min-MQ so i just used -b and -q respetively.

    # for samtools sort
    # -o : takes a file. it writes the final sorted output to file rather than standard output
    # -O : write the final output as sam, bam, or cram

    # samtools fixmate : preparing for finding the duplicates
    # -O : specify the format and i choose bam
    # -m : add ms(mate score) tags. these are used by markdup to select the best reads to keep 
    # other possible option to look at is -r : remove secondary and unmapped reads ?
    
    # for samtools markdup : can only be done on coordinate sorted bam files and run it through samtools fixmate first
    # -r : remove duplicate reads
    # --mode or -m : duplicate decision method for paired reads. values are "t" or "s". read documentation but i choose s becasue it tends to return more results. just incase i will remove this option in the pair end mode by adding an if else statement in this process.

    # now for samtools index to get index files
    # -b, --bai: create a bai index; this version of samtools does not support --bai --csi just use -b -c
    # -c, --csi: create a csi index
    # -o, --output: write the output index to a file specified  only when one alignment file is being indexed

    # using samtools flagstat: generate log files so i can use multiqc to get stats of all files into one html file
    # no parameters needed. just need to give the final bam file that went through all the processing
    ###########################################################

    # I should add a samtools filtering. looking to only get mapq scores higher than 30

    samtools view \
    -q 30 \
    -b \
    "${sam_files}" \
    > "${out_bam_filt}" 
    
    
    # first i have to name sort to use fixmate
    samtools sort \
    -o "${out_bam_name_sort}" \
    -n \
    -O bam \
    "${out_bam_filt}"


    samtools fixmate \
    -O bam \
    -m \
    "${out_bam_name_sort}"\
    "${out_bam_fixmate}"


    # now i will coordinate sort here 
    samtools sort \
    -o "${out_bam_coor_sort}" \
    -O bam \
    "${out_bam_fixmate}"


    # i might need to put coordinate sorted bam into markdup
    # this works but removing the duplicates results in the file being very small meaning too many reads were removed that were considered duplicates.
    # this results in the next process deeptools not being able to create a normalized bedgraph file

    #samtools markdup \
    #"\${out_bam_coor_sort}" \
    #"\${out_bam_final}"

    # so i will just use the out file from the coordinate sort samtools sort section instead of using out_bam_final
    samtools index \
    -b \
    "${out_bam_coor_sort}"

    #samtools flagstat \
    #"\${out_bam_coor_sort}" \
    #> "\${flagstats_log}"

    # adding another way to get stats from each bam file
    #samtools stats \
    #"\${out_bam_coor_sort}" \
    #> "\${samtools_stats_log}"

    # now only putting the stats into a tsv file
    #less "\${samtools_stats_log}" | grep ^SN | cut -f 2-3 >  "\${tsv_file_with_stats}"
    
    """





    
}

process bam_log_calc {

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools-1.21_rj'

    publishDir "${params.base_out_dir}/flag_stat_log", mode: 'copy', pattern: '*stats.log'
    publishDir "${params.base_out_dir}/stats_tsv_files", mode: 'copy', pattern: '*stats.tsv'


    input:
    tuple path(bam), path(index_bam)


    output:
    path("*stats.log"), emit: flag_stats_log
    path("*stats.txt"), emit: norm_stats_txt
    path("*stats.tsv"), emit: tsv_SN_stats



    script:
    flagstats_log = "${bam.baseName}_flag_stats.log"
    samtools_stats_log = "${bam.baseName}_stats.txt"
    tsv_file_with_stats = "${bam.baseName}_SN_stats.tsv"


    """
    #!/usr/bin/env bash

    # using samtools flagstat: generate log files so i can use multiqc to get stats of all files into one html file
    # no parameters needed. just need to give the final bam file that went through all the processing


    samtools flagstat \
    "${bam}" \
    > "${flagstats_log}"

    # adding another way to get stats from each bam file
    samtools stats \
    "${bam}" \
    > "${samtools_stats_log}"

    # now only putting the stats into a tsv file
    less "${samtools_stats_log}" | grep ^SN | cut -f 2-3 >  "${tsv_file_with_stats}"




    """
}



process deeptools_make_bed {
    // this conda env yml file didnt work. have to use the actual env
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/deeptools_rj_env.yml'
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    // this section is just a simple if else statement controlling the directories that are created and when the files end up
    // will copy and paste in other processes that need it.
    

    if (params.BL) {

        publishDir "${params.base_out_dir}/bl_filt_bed/bed_graphs_deeptools/", mode: 'copy', pattern: '*'

    }
    else {
        publishDir "${params.base_out_dir}/no_bl_filt/bed_graphs_deeptools/", mode: 'copy', pattern: '*'
    }
    
    


    input:

    tuple path(bams), path(index)


    output:
    path("${out_bed_name}*"), emit: bed_files_normalized

    script:

    out_bed_name="${bams.baseName}_normalized_cpm.bed"

    if (params.use_effectiveGenomeSize) {

        """

        ###### Using deeptools parameters ###############

        # first converting the bam file to a bed file using bamCoverage. I can also make a bigwig file if needed, it stores data better but is binary and cannot be opened in text editor
        # -b or --bam: takes the bam file that will be processed
        # -o or --outFileName: is the name you want the output file to have
        # -of or --outFileFormat: is the type of output file you want; either "bigwig" or "bedgraph"
        # --scaleFactor: the computed scaling factor (or 1, if not applicable) will be multiplied by this.
        # -bs or --binSize: are the size of the bins in bases, for output of bigwig or bedgraph. default is 50
        # -p or --numberOfProcessors: this is the number of processers you want to use. Not using this option yet but if needed I will use it.
        # --normalizeUsing: choose the type of normalization
        # bamCoverage offers normalization by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM), counts per million (CPM), bins per million mapped reads (BPM) and 1x depth (reads per genome coverage, RPGC).
        # --effectiveGenomeSize: choose the mappable genome size for your organism of choice used as reference. find length here: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
        # not using effectiveGenomeSize since multiple users will use this pipeline and might not be using the same organism.
        # actually i decided to use effectiveGenomeSize afterall since i can split this process into using it or not.

        # NOTE: since all the files will be processed using this tool and parameters, they will all be directly comparable in UCSC or IGV without needing to edit track heights.
        #################################################


        bamCoverage \
        --bam "${bams}" \
        --outFileName "${out_bed_name}" \
        --outFileFormat "bedgraph" \
        --scaleFactor 1 \
        --binSize 50 \
        --normalizeUsing CPM \
        --effectiveGenomeSize "${params.num_effectiveGenomeSize}"


        """


    }
    else {

        """

        bamCoverage \
        --bam "${bams}" \
        --outFileName "${out_bed_name}" \
        --outFileFormat "bedgraph" \
        --scaleFactor 1 \
        --binSize 50 \
        --normalizeUsing CPM

        """

    }

}

process bedtools_filt_blacklist {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    //publishDir './blacklist_filt_bam', mode: 'copy', pattern: '*.bam'

    
    // dont need to publish this since the samtools_bl_index will publish the sorted 2 version of the bam
    //publishDir "${params.base_out_dir}/blacklist_filt_bam", mode: 'copy', pattern: '*.bam'
    
    

    input:

    tuple path(bam), path(index)
    path(blacklist_file)

    output:
    path("${out_bl_filtered_bam}"), emit: bl_filtered_bams


    script:
    out_bl_filtered_bam = "${bam.baseName}_BL_filt.bam"

    """
    ############ using bedtools intersect parameters ##############
    # -a : takes a bam/bed/gff/vcf file A
    # -b : takes one or more of the files bam/bed/gff/vcf and calls it B
    # -v : if any regions in A do not have a overlap in B then keep only those regions. that will remove the blacklist regions found in the blacklist B file

    #####################################################

    bedtools intersect \
    -a "${bam}" \
    -b "${blacklist_file}" \
    -v \
    > "${out_bl_filtered_bam}"


    """

}

process samtools_bl_index {
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj'

    

    //publishDir './blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern:'*.bai'
   

    // if (params.ATAC) {

    //     publishDir "${params.base_out_dir}/ATAC_blacklist_filt_bam", mode: 'copy', pattern: '*.bai'
    //     publishDir "${params.base_out_dir}/ATAC_blacklist_filt_bam", mode: 'copy', pattern: '*_sort2.bam'

    // }
    // else {

    //     publishDir "${params.base_out_dir}/blacklist_filt_bam", mode: 'copy', pattern: '*.bai'
    //     publishDir "${params.base_out_dir}/blacklist_filt_bam", mode: 'copy', pattern: '*_sort2.bam'
    // }

    publishDir "${params.base_out_dir}/sorted_bam_files/blacklist_filt_bam", mode: 'copy', pattern: '*.bai'
    publishDir "${params.base_out_dir}/sorted_bam_files/blacklist_filt_bam", mode: 'copy', pattern: '*_sort2.bam'
    
    

    input:
    path(bl_filt_bam)


    output:

    tuple path("${out_bam_name_sort}"), path("*.bai"), emit: bl_filt_bam_index_tuple


    script:
    
    out_bam_name_sort = "${bl_filt_bam.baseName}_sort2.bam"

    """
    ####### parameters for indexing bam ######
    # -b : will create a bai file

    ##########################################

    # just do some sorting 

    samtools sort \
    -o "${out_bam_name_sort}" \
    -O bam \
    "${bl_filt_bam}"



    samtools index \
    -b \
    "${out_bam_name_sort}" 

    """
}

process fastp_PE {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/fastp_rj'

    publishDir "${params.base_out_dir}/fastp_pe_results/filt_fastqs", mode: 'copy', pattern: '*_filt_{R1,R2}*'

    publishDir "${params.base_out_dir}/fastp_pe_results/merged_filt_fastqs", mode: 'copy', pattern: '*_merged*'

    publishDir "${params.base_out_dir}/fastp_pe_results/failed_qc_reads", mode: 'copy', pattern: '*_failed_filter*'

    publishDir "${params.base_out_dir}/fastp_pe_results/htmls", mode: 'copy', pattern: '*fastp.html'


    input:

    tuple val(fastq_name), path(fastq)


    output:

    tuple val(fastq_name), path("${out_name_1}"), path("${out_name_2}"), emit: filt_PE_tuple
    //path("${merged_reads_file}"), emit: merged_filt_reads
    path("${html_file_name}"), emit: html_fastp_out
    path("${failed_reads_file}"), emit: failed_reads_out


    script:

    out_name_1 = "${fastq_name}_filt_R1_0.fastq"
    out_name_2 = "${fastq_name}_filt_R2_0.fastq"
    failed_reads_file = "${fastq_name}_failed_filter_reads.fastq"
    merged_reads_file = "${fastq_name}_merged_file_reads.fastq"
    html_file_name = "${fastq_name}_R1_R2_fastp.html"

    """
    #!/usr/bin/env bash

    # this below is just debugging
    echo "this is the file name: ${fastq_name}; this is the forward read (r1): ${fastq[0]}; this is the reverse read(r2): ${fastq[1]}"


    ############# Now the PE fastp parameters ################
    # --in1 or -i : takes the first pair end read
    # --in2 or -I : takes the second pair end read
    # --out1 or -o : is the output file name of the fastp first read
    # --out2 or -O : is the output file name of the fastp second read
    # --failed_out : specify the file to store reads that cannot pass the filters
    # --merge or -m : for paired end input, merge each pair of reads into a single read if they are overlapped. the merged reads are written to --merged_out, the unmerged reads will be written to the --out1 and --out2
    # --detect_adapter_for_pe : to enable auto detection for PE data. The auto detection for adapter is for SE only, so have to turn on for PE
    # --dedup : enable deduplication to drop the duplicated reads pairs
    # --dup_calc_accuracy : accuracy level to calculate duplication (1~6) higher level uses more memory
    # --trim_poly_g : force polyG tail trimming. polyG can happen if there is no signal in the illumina two-color systems
    # --trim_poly_x : enable polyX trimming in 3' ends. useful because polyX(polyA) can be found in the tails of mRNA seq reads. DONT NEED THIS IN THIS PIPELINE
    # --qualified_quality_phred : this is the quality value that a base is qualified. default 15 but i choose 20 meaning phred quality >= 20
    # --unqualified_percent_limit : how many percent of bases are allowed to be unqualified. defualt is 40 meaning 40%
    # --n_base_limit : if one read's number of N base is > n ase imit  , then this read/pair is discarded default is 5. see fastp documentation
    # --average_qual :  if one read's average quality score < avg ual , then this read/pair is discarded. Default 0 means no requirement
    # --correction or -c : enable base correction in overlapped regions (only for PE data), default is disabled. fastp performs overlap analysis for PE data, which try to find an overlap of each pair of reads. When this option is enabled, and if an proper overlap is found, it can correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base will be assigned to it so that they will share the same quality.
    # --overlap_len_require : the min length to detect oveerlapped region of PE reads. this will affect overlap analysis based PE merge, adapter trimming. default is 30 but i put 20 since the risca_lab snakemake pipeline looks for 20 
    # --overlap_diff_limit : the maximum number of mismatched bases to detect overlapped region of PE reads. it will affect any trimming or merge parameters. the default is 5 but i think the risc_lab snakemake pipeline allowed only for 1 mismatched bases. the snake make pipeline does only allow for 1 mismatches. But the author said "for now" so I assumed they will look to change it in the future. maybe i can change it back to the default which is 5
    # --overrepresentation_analysis : enable overrepresented sequence analysis
    # --overrepresentation_sampling : the number of reads computed for overrepresentation analysis (1~10000). default 20 i used 30
    # --html or -h : the html format report file name. default is fastp.html but i made my own name using the base name string (key) of the two fastq files that were input in this channel
    # --thread or -w : worker thread number default is 2; i used 15
    ###########################################################

    # NOTE: I will remove the merged reads and only keep all the reads that pass the filtering in their corresponding forward and reverse reads
    #--merge \
    #--merged_out "\${merged_reads_file}" 

    fastp \
    --in1 "${fastq[0]}" \
    --in2 "${fastq[1]}" \
    --out1 "${out_name_1}" \
    --out2 "${out_name_2}" \
    --failed_out "${failed_reads_file}" \
    --detect_adapter_for_pe \
    --dedup \
    --dup_calc_accuracy 5 \
    --trim_poly_g \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --average_qual 0 \
    --correction \
    --overlap_len_require 20 \
    --overlap_diff_limit 1 \
    --overrepresentation_analysis \
    --overrepresentation_sampling 30 \
    --html "${html_file_name}" \
    --thread 15

    """

}

process fastqc_PE {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/fastqc_rj_env.yml'

    publishDir "${params.base_out_dir}/fastqc_pe_files", mode: 'copy', pattern: '*'

    input:
   
    tuple val(fastq_name), path(filt_r1), path(filt_r2)


    output:

    path ("*.html"), emit: fastqc_htmls
    path("*.zip"), emit: fastqc_zip_files



    script:


    """
    #!/usr/bin/env bash

    ########### parameters for fastqc pe reads ###############
    # no real parameters 

    ###########################################################

    fastqc \
    "${filt_r1}" 

    fastqc \
    "${filt_r2}" 




    """


}


process multiqc_PE {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    publishDir "${params.base_out_dir}/multiqc_PE_output", mode: 'copy', pattern: '*'


    input:
    path(fastqc_zip_files)


    output:

    path("*"), emit: summary_of_PE_filt
    


    script:


    """
    #!/usr/bin/env bash

    ######### parameters for making the html file with all the fastqc output data #############

    
    # --ai : generate an AI summary of the report
    # --ai-summary-full : generate a detailed AI summary of the report
    # --interactive : use only interactive plots
    # --profile-runtime : add an analysis of how long multiqc takes to run to the report
    # --ai-provider : choose ai provider defualt sequra . openai or anthropic
    ############################################################################

    multiqc . \
    --ai-summary-full \
    --ai-provider openai \
    --interactive \
    --profile-runtime \
    --title "Pair end QC"


    """
}

process bwa_PE_aln {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/bwa_rj_env.yml'

    publishDir "${params.base_out_dir}/pe_bwa_files/pe_sam_files", mode: 'copy', pattern: '*.{sam, sai}'
    //publishDir "${params.base_out_dir}/pe_bwa_files/pe_sai_index_files", mode: 'copy', pattern: '*.sai'
    //cache false 


    input:
    tuple val(filt_fastq_name), path(fastq_r1), path(fastq_r2)
    path(genome)
    path(genome_index)


    output:

    path("*.sam"), emit: pe_sam_files
    path("*.sai"), emit: pe_sai_files



    script:

    sai_out_file_r1 = "${filt_fastq_name}_filt_r1.sai"
    sai_out_file_r2 = "${filt_fastq_name}_filt_r2.sai"

    out_sam_file = "${filt_fastq_name}_filt_r1_r2.sam"


    """
    #!/usr/bin/env bash

    ######## bwa aln parameters / bwa sampe params #########
    # -t : allows for the amout of threads you want this process to use

    #



    #########################################################


    bwa aln \
    -t 20 \
    "${genome}" \
    "${fastq_r1}" \
    > "${sai_out_file_r1}"

    bwa aln \
    -t 20 \
    "${genome}" \
    "${fastq_r2}" \
    > "${sai_out_file_r2}"


    bwa sampe \
    "${genome}" \
    "${sai_out_file_r1}" \
    "${sai_out_file_r2}" \
    "${fastq_r1}" \
    "${fastq_r2}" \
    > "${out_sam_file}"

    """
}

process multiqc_bam_stats {

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/multiqc_rj_env.yml'

    

    publishDir "${params.base_out_dir}/flag_stat_log/complete_log", mode: 'copy', pattern: '*.html'

   
    

    input:

    path(stats_log_files)
    path(norm_stats_files)


    output:
    path("*.html"), emit: bams_multiqc_html


    script:

    if (params.PE) {

        """

        #### parameters for multiqc ###

        # no real parameters
        ###############################
        
        multiqc . \
        --title "Pair end bams QC"



        """

    }
    else if (params.SE) {

        """

        #### parameters for multiqc ###

        # no real parameters
        ###############################
        
        multiqc . \
        --title "Single end bams QC"



        """

    }
    
}

process deeptools_aln_shift {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    // commenting this out since the samtools_index_sort will have the sorted bam published so dont need another
    /*if (params.PE) {

        publishDir './results_PE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }
    else {

        publishDir './results_SE/atac_shift_bam', mode: 'copy', pattern: '*shift.bam'
    }*/


    input:
    tuple path(bam), path(index)


    output:
    path("*shift.bam"), emit: atac_shifted_bam


    script:
    out_file = "${bam.baseName}_ATAC_shift.bam"

    """
    #!/usr/bin/env bash 

    ##### deeptools alignmentSieve params #####
    # --bam : this takes your input bam file
    # --numberofProcessors : number of threads 
    # --ATACshift : Shift the produced BAM file or BEDPE regions as commonly done for ATAC-seq
    ###########################################

    alignmentSieve \
    --bam "${bam}" \
    --numberOfProcessors 12 \
    --ATACshift \
    -o "${out_file}"




    """
}





/////////////////////////////////////////////////////////////////////////
// Above is where i copy and pasted the processes from the fastq2bam_nextflow_pipeline.nf main script
/////////////////////////////////////////////////////////////////////////











process samtools_index_sort {
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/samtools_rj'

    //publishDir './blacklist_filt_bam/bl_filt_index', mode: 'copy', pattern:'*.bai'
    

    if (params.ATAC && params.BL) {
        // I want to put both the bam and the index (bai) in the same channel
        publishDir "${params.base_out_dir}/sorted_bam_files/ATAC_blacklist_filt_bam", mode: 'copy', pattern: '*.bai'
        publishDir "${params.base_out_dir}/sorted_bam_files/ATAC_blacklist_filt_bam", mode: 'copy', pattern: '*_sort2.bam'

    }
    else if(params.ATAC) {

        publishDir "${params.base_out_dir}/sorted_bam_files/ATAC_filt_bam", mode: 'copy', pattern: '*.bai'
        publishDir "${params.base_out_dir}/sorted_bam_files/ATAC_filt_bam", mode: 'copy', pattern: '*_sort2.bam'
    }



          
    

    input:
    path(bam) // changed this to be just bam


    output:

    tuple path("${out_bam_name_sort}"), path("*.bai"), emit: bam_index_tuple
    
    

    script:
    
    out_bam_name_sort = "${bam.baseName}_sort2.bam"

    """
    ####### parameters for indexing bam ######
    # -b : will create a bai file

    ##########################################

    # just do some sorting 

    samtools sort \
    -o "${out_bam_name_sort}" \
    -O bam \
    "${bam}"



    samtools index \
    -b \
    "${out_bam_name_sort}" 

    """
}


process mk_break_points {
    // this is my attempt at creating the break density script, but i will try another process where I just call the break density wrapper.

    // all the wrapper does is make each peak file work with each bam file. so just make a process that takes the bam file and run all the break density scripts with each peak file. each bam will work with all the peak files in their own process instance.

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    publishDir "${params.base_out_dir}/break_point_bed", mode: 'copy', pattern: '*_breaks.bed'
    //publishDir './results_PE/'

    input:
    // take the bam files either bl filtered or not bl filtered from only the pair end path
    path(bams)
    // i think i need to have a narrowpeak file or a bed file that has peak information as input also for the next part,
    // have to call peaks and make a narrowpeak or bed peak file
    //path(bed_peaks)

    output:


    script:

    out_bampe_name = "${bams.baseName}_bampe.bed"
    break_point_name = "${bams.baseName}_breaks.bed"
    sorted_break_point = "${bams.baseName}_sorted.bed"

    // making sure the bed file is sorted so just doing it again
    sorted_bed_file = 

    """
    #!/usr/bin/env bash

    #### parameters bedtools ######
    # -bedpe : write bam alignments in bedpe format will have second field as start coordinates for forward read and 6th field as start coordinates for reverse reads
    # -i : the input bam file


    ###############################

    # this will create a bampe bed file so we can get the first field the second field and the sixth field

    bedtools bamtobed \
    -bedpe \
    -i "${bams}" \
    > "${out_bampe_name}"

    awk '{print \$1"\t"\$2"\t"\$2}' "${out_bampe_name}" > "${break_point_name}"
    awk '{print \$1"\t"\$6"\t"\$6}' "${out_bampe_name}" >> "${break_point_name}"

    # sorting the break bed files now
    bedtools sort \
    -i "${break_point_name}" \
    > "${sorted_break_point}"


    numBreaks=\$(wc -l "${break_point_name}")


    bedtools sort \
    -i "${bed}" \
    > 




    """
}




process breakDensityWrapper_process {

    publishDir "${params.base_out_dir}/break_density_calc", mode: 'copy', pattern: '*'
    

    input:
    // input has to be files that the breakDensityWrapper.sh script takes
    
    // it takes the bam files that have reads aligned to the reference genome. So I think it is best to collect all of the generated bams and ensure PLC is one of the bams.
    path(bams)
    
    // then it takes the peak files found in the directory  /lustre/fs4/home/ascortea/store/ascortea/beds
    path(peak_files)


    output:

    // output will be an AdjustedEnrichment.tsv
    path("adjustedEnrichment.tsv"), emit: adjusted_E_tsv
    path("Adjusted_Enrichment_of_*_Plot.pdf"), emit: break_plot_pdf
    path("densityCalculations.log"), emit: density_calc_log


    script:


    """
    #!/usr/bin/env bash

    # nextflow can find stuff in the bin dir but not recursively, so i have to specify the sub dir

    breakDensityWrapper.sh "${bams}" "${peak_files}"
    
    # this works but for some reason it is not seeing the output so it can be put in the published dir and also in the emit channels.




    """



}



process py_calc_stats_log {

    //debug true

    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/python_w_packages_rj_env.yml'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/python_w_packages_rj'

    //shell '/bin/python3'

    publishDir "${params.base_out_dir}/py_calc_stats_log", mode: 'copy', pattern: '*.tsv'


    input:
    tuple val(tsv_names), path(tsv_sn_stats)


    output:

    path("bam_*.tsv"), emit: pe_tsv_log
    
    


    script:

    

    name_of_file = "${tsv_names}"
    //pe_log_file_out = "pe_bam_stats_log.tsv"
    //se_log_file_out = "se_bam_stats_log.tsv"
    log_file_out = "bam_stats_log.tsv"

    """
    #!/usr/bin/env python
           
    # for testing purposes on the hpc use: 
    # list_tsv_files = ['results_PE/stats_tsv_files/B_HB1_0GyP_T1_S00_filt_r1_r2_SN_stats.tsv','results_PE/stats_tsv_files/I_HB2_PLC_T1_S00_filt_r1_r2_SN_stats.tsv']
    # list_of_names = ['B_HB1_0GyP_T1_S00_filt_r1_r2_SN_stats', 'I_HB2_PLC_T1_S00_filt_r1_r2_SN_stats']

    import pandas as pd

    print("before the split : ${tsv_sn_stats}")
    print("before the split: ${name_of_file}")
    

    # need to remove the brackets and split
    str_names = "${name_of_file}"  #.strip('[]').split(',')
    new_list_names = str_names.strip('[]').split(',')

    list_tsv_files = "${tsv_sn_stats}".split()
    #list_of_names = "\${name_of_file}".split()
    
    print("after the split:",  list_tsv_files)
    print("after the split:",  new_list_names)

    files_dict_df = {}

    for name, file in zip(new_list_names, list_tsv_files):
        files_dict_df[name] = pd.read_table(file, header=None, sep='\t').set_index(0).T

    # concat the dictionary of dataframes
    combined_stats_df = pd.concat( files_dict_df.values(), axis = 0)

    # now  putting sample names in the df
    combined_stats_df.insert(0, "sample_names", files_dict_df.keys())

    # use this code below to see the column names to choose from
    # combined_stats_df.columns

    final_log_stats_df = combined_stats_df.loc[:, ['sample_names','raw total sequences:','reads mapped:','reads mapped and paired:', 'reads duplicated:', 'reads MQ0:', 'percentage of properly paired reads (%):']]

    # calculating the percentage of reads mapped, and adding it to a column called percent_reads_mapped
    final_log_stats_df['percent_reads_mapped'] = (final_log_stats_df.iloc[:,2]/final_log_stats_df.iloc[:,1])*100

    # calculating the percentage of reads duplicated and adding that as a column
    final_log_stats_df['percent_reads_duplicated'] = (final_log_stats_df.iloc[:,4]/final_log_stats_df.iloc[:,1])*100


    # right before I make the log file I want to clean up the column names and replace and spaces with an underscore
    # then also remove the colons 
    final_log_stats_df.columns = final_log_stats_df.columns.str.replace(' ','_')
    final_log_stats_df.columns = final_log_stats_df.columns.str.replace(':','')

    # writing the tsv
    final_log_stats_df.to_csv("${log_file_out}", sep = '\t')

    
   
    
    """

}

