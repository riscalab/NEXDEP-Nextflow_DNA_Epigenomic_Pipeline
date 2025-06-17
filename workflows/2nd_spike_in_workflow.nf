// i need to put this before the modules are loaded 
// try this here


// I need access to the processes so I put them in the module file

// didn't need this here. was able to define this above the main workflow
//params.base_out_dir = params.PE ? './results_PE' : (params.SE ? './results_SE' : '')



params.t7_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/T7_Phage/genome/Sequence/WholeGenomeFasta/genome.fa')

t7_genome_tuple =Channel.value(params.t7_genome).map{file -> tuple '_t7', file}




params.yeast_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/S_pombe_EF2/genome/Sequence/WholeGenomeFasta/genome.fa')

yeast_genome_tuple = Channel.value(params.yeast_genome).map{file -> tuple '_yeast', file}




params.lambda_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Lambda_cI857ind_1_Sam_7/genome/Sequence/Bowtie2Index/genome.fa')


lambda_genome_tuple =Channel.value(params.lambda_genome).map{file -> tuple '_lambda', file }



// maybe this has to be outside of the workflow to work


// if i pair the genome with the spike in name I will be able to pass the proper spike in name along the process




include {
    fastp_SE_adapter_known_spike_in;
    fastp_SE_spike_in;
    fastqc_SE_spike_in;
    multiqc_SE_spike_in;
    bwa_index_genome_spike_in;
    bwa_align_SE_spike_in;
    samtools_sort_spike_in;
    deeptools_make_bed_spike_in;
    fastp_PE_spike_in;
    fastqc_PE_spike_in;
    multiqc_PE_spike_in;
    bwa_PE_aln_spike_in;
    samtools_index_sort_spike_in;
    deeptools_aln_shift_spike_in
    
} from '../modules/spike_in_modules.nf'


workflow pe_t7_spike_in_workflow {

    // take:
    // base_dir_results

    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch



    main:
    
    
    // testing something where I make a variable based on the spike in being ran
    //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

    // this will take the paired end reads and keep them together
    //params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

    //pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)

    if (params.test) {
            
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
        
    }else {

        
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    }

    //pe_fastqs_ch.view()

    fastp_PE_spike_in(pe_fastqs_ch,t7_genome_tuple)

    // checking the channels to see if everything works
    //fastp_PE_spike_in.out.filt_PE_tuple.view()

    //fastp_PE_spike_in.out.html_fastp_out.view()
    //fastp_PE_spike_in.out.failed_reads_out.view()
    //fastp_PE_spike_in.out.

    pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

    fastqc_PE_spike_in(pe_filt_tuple_ch,t7_genome_tuple)

    //fastqc_PE_spike_in.out.fastqc_zip_files.view()


    // then collect them so i can view in one file using multiqc.
    collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

    multiqc_PE_spike_in(collection_fastqc_ch,t7_genome_tuple)

    //multiqc_PE.out.summary_of_PE_filt.view()

    // first use the process to index the reference genome since the process exists already for the se

    bwa_index_genome_spike_in(t7_genome_tuple)

    //bwa_index_genome_spike_in.out.genome_index_files.view()

    genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

    // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

    //pe_filt_tuple_ch.view()
    //t7_phage_genome_ch.view()
    //genome_index_ch.view()

    bwa_PE_aln_spike_in(pe_filt_tuple_ch, genome_index_ch, t7_genome_tuple)

    // now check to see if the output channels are good
    //bwa_PE_aln_spike_in.out.pe_sam_files.view()
    
    
    // now i need to make the parameters for  if the bam file will be blacklist filtered or not
    
    // using this channel for both if Blacklist or not
    sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

    /*
    if (params.BL) {

        
        // i have to make a bam file to then use bedtools intersect to get the blacklist
        // using the same samtools sort process found in the SE part of the pipeline
        samtools_sort_spike_in(sam_files_pe_ch, t7_genome_tuple)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        // this will give a blacklist filtered bam but i need to index it again
        bedtools_filt_blacklist_spike_in(spike_in_bam_index_tuple_ch, blacklist_ch, t7_genome_tuple)
        //bedtools_filt_blacklist_spike_in.out.bl_filtered_bams.view()

        bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
        // so using the process to only index which means it will take the blacklist bam file
        
        samtools_bl_index(bl_filt_bams_ch) 

        spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

        

    }
    else {
        samtools_sort_spike_in(sam_files_pe_ch, t7_genome_tuple)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        
    }*/

    // not running blacklist filtering on spike ins
    samtools_sort_spike_in(sam_files_pe_ch, t7_genome_tuple)

    //samtools_sort_spike_in.out.bam_index_tuple.view()

    /*
    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
    


    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, t7_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }


    

}




workflow pe_lambda_spike_in_workflow {

    // take:
    // base_dir_results

    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch


    main:

    // testing something where I make a variable based on the spike in being ran
    //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

    // this will take the paired end reads and keep them together
    //params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

    //pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)

    if (params.test) {
            
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
    
    }else {

        
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    }


    //pe_fastqs_ch.view()

    fastp_PE_spike_in(pe_fastqs_ch,lambda_genome_tuple)

    // checking the channels to see if everything works
    //fastp_PE_spike_in.out.filt_PE_tuple.view()

    //fastp_PE_spike_in.out.html_fastp_out.view()
    //fastp_PE_spike_in.out.failed_reads_out.view()
    //fastp_PE_spike_in.out.

    pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

    fastqc_PE_spike_in(pe_filt_tuple_ch,lambda_genome_tuple)

    //fastqc_PE_spike_in.out.fastqc_zip_files.view()


    // then collect them so i can view in one file using multiqc.
    collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

    multiqc_PE_spike_in(collection_fastqc_ch,lambda_genome_tuple)

    //multiqc_PE.out.summary_of_PE_filt.view()

    // first use the process to index the reference genome since the process exists already for the se

    bwa_index_genome_spike_in(lambda_genome_tuple)

    //bwa_index_genome_spike_in.out.genome_index_files.view()

    genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

    // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

    //pe_filt_tuple_ch.view()
    //t7_phage_genome_ch.view()
    //genome_index_ch.view()

    bwa_PE_aln_spike_in(pe_filt_tuple_ch,  genome_index_ch, lambda_genome_tuple)

    // now check to see if the output channels are good
    //bwa_PE_aln_spike_in.out.pe_sam_files.view()
    
    
    // now i need to make the parameters for  if the bam file will be blacklist filtered or not
    
    // using this channel for both if Blacklist or not
    sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

    /*if (params.BL) {

        
        // i have to make a bam file to then use bedtools intersect to get the blacklist
        // using the same samtools sort process found in the SE part of the pipeline
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        // this will give a blacklist filtered bam but i need to index it again
        bedtools_filt_blacklist_spike_in(spike_in_bam_index_tuple_ch, blacklist_ch, lambda_genome_tuple)
        //bedtools_filt_blacklist_spike_in.out.bl_filtered_bams.view()

        bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
        // so using the process to only index which means it will take the blacklist bam file
        
        samtools_bl_index(bl_filt_bams_ch) 

        spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

        

    }
    else {
        samtools_sort_spike_in(sam_files_pe_ch, lambda_genome_tuple)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        
    }*/

    // not running blacklist filtering on spike ins
    //lambda_genome_tuple.view()
    samtools_sort_spike_in(sam_files_pe_ch, lambda_genome_tuple)

    //samtools_sort_spike_in.out.bam_index_tuple.view()
    /*
    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, lambda_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }



}

workflow pe_yeast_spike_in_workflow {

    
    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch


    main:

    // testing something where I make a variable based on the spike in being ran
    //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

    // this will take the paired end reads and keep them together
    //params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

    //pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)

    if (params.test) {
            
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
    
    }else {

        
        pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    }

    //pe_fastqs_ch.view()

    fastp_PE_spike_in(pe_fastqs_ch,yeast_genome_tuple)

    // checking the channels to see if everything works
    //fastp_PE_spike_in.out.filt_PE_tuple.view()

    //fastp_PE_spike_in.out.html_fastp_out.view()
    //fastp_PE_spike_in.out.failed_reads_out.view()
    //fastp_PE_spike_in.out.

    pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

    fastqc_PE_spike_in(pe_filt_tuple_ch,yeast_genome_tuple)

    //fastqc_PE_spike_in.out.fastqc_zip_files.view()


    // then collect them so i can view in one file using multiqc.
    collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

    multiqc_PE_spike_in(collection_fastqc_ch, yeast_genome_tuple)

    //multiqc_PE.out.summary_of_PE_filt.view()

    // first use the process to index the reference genome since the process exists already for the se

    bwa_index_genome_spike_in(yeast_genome_tuple)

    //bwa_index_genome_spike_in.out.genome_index_files.view()

    genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

    // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

    //pe_filt_tuple_ch.view()
    //t7_phage_genome_ch.view()
    //genome_index_ch.view()

    bwa_PE_aln_spike_in(pe_filt_tuple_ch,  genome_index_ch, yeast_genome_tuple)

    // now check to see if the output channels are good
    //bwa_PE_aln_spike_in.out.pe_sam_files.view()
    
    
    // now i need to make the parameters for  if the bam file will be blacklist filtered or not
    
    // using this channel for both if Blacklist or not
    sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

    /*if (params.BL) {

        
        // i have to make a bam file to then use bedtools intersect to get the blacklist
        // using the same samtools sort process found in the SE part of the pipeline
        samtools_sort_spike_in(sam_files_pe_ch)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        // this will give a blacklist filtered bam but i need to index it again
        bedtools_filt_blacklist_spike_in(spike_in_bam_index_tuple_ch, blacklist_ch, yeast_genome_tuple)
        //bedtools_filt_blacklist_spike_in.out.bl_filtered_bams.view()

        bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
        // so using the process to only index which means it will take the blacklist bam file
        
        samtools_bl_index(bl_filt_bams_ch) 

        spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

        

    }
    else {
        samtools_sort_spike_in(sam_files_pe_ch, yeast_genome_tuple)

        //samtools_sort_spike_in.out.bam_index_tuple.view()

        flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
        norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
        tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

        spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

        
    }*/

    // not running blacklist filtering on spike ins
    //yeast_genome_tuple.view()
    samtools_sort_spike_in(sam_files_pe_ch, yeast_genome_tuple)

    //samtools_sort_spike_in.out.bam_index_tuple.view()
    /*
    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple


    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, yeast_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }


}

workflow se_t7_spike_in_workflow {


    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch


    main:


    //params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    //se_reads_files = Channel.fromPath(params.single_end_reads).take(3)

    if (params.test) {
        se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
        
    
    }else {

        se_reads_files = Channel.fromPath(params.single_end_reads)
        

    }
    
    // now let's get the basename of the single end reads
    // removing both the .gz and the .fastq
    // I would normally use file.baseName here but it had the .gz and the .fastq
    se_reads_files.flatten()
                    .map{ file -> file.name.replace('.fastq.gz','')}
                    .set{se_reads_name}
    
    // let's view both the files and the names to make sure they match in order
    //se_reads_files.view()
    //se_reads_name.view()
    // this is where i send both the input file and their corresponding basenames to the fastp_SE process
    

    // if the adapter sequence is known then input it as a string if not dont use the parameter
    if ( params.ada_seq ) {

        //params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
        adapter_ch = Channel.value(params.adapter_seq_str)

        fastp_SE_adapter_known_spike_in(se_reads_files, se_reads_name, adapter_ch, t7_genome_tuple) // will have to make a new process for if the adapter sequence is known

        fastq_filts = fastp_SE_adapter_known_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_adapter_known_spike_in.out.fastp_html_reports

    }    
    else {

        fastp_SE_spike_in(se_reads_files, se_reads_name, t7_genome_tuple)

            // take all of the filtered fastq files and put them in a channel name
        // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
        fastq_filts = fastp_SE_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_spike_in.out.fastp_html_reports



    }
    

    //fastp_SE(se_reads_files, se_reads_name) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

    

    
    //fastp_SE.out.view()

    //fastq_filts.view()
    //fastp_filt_html.view()
    //fastp_filt_html.collect().view()


    // now creating a fastqc process
    fastqc_SE_spike_in(fastq_filts, fastq_filts_name, t7_genome_tuple)

    fastqc_html_files = fastqc_SE_spike_in.out.fastqc_htmls
    fastqc_zips = fastqc_SE_spike_in.out.fastqc_zip_files

    // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
    multiqc_SE_spike_in(fastqc_zips.collect(), t7_genome_tuple)

    // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
    bwa_index_genome_spike_in(t7_genome_tuple)

    // collecting the genome index files from the last process 
    // not sure if i should keep track of the order the files are in first
    // it looks like they are in the same order that they appeared in the published dir using ll
    //bwa_index_genome.out.genome_index_files.view()
    genome_index_files_ch = bwa_index_genome_spike_in.out.genome_index_files

    // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
    // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
    bwa_align_SE_spike_in(genome_index_files_ch, fastq_filts, fastq_filts_name, t7_genome_tuple )

    //bwa_align_SE.out.sam_se_files.view()

    // making a channel for the sam files generated
    sam_files = bwa_align_SE_spike_in.out.sam_se_files

    // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
    // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
    // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
    
    samtools_sort_spike_in(sam_files, t7_genome_tuple) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

    //samtools_sort_spike_in.out.sorted_bams.view()
    //samtools_sort_spike_in.out.indexed_bams.view()
    //samtools_sort_spike_in.out.bam_index_tuple.view()

    sorted_bams_ch = samtools_sort_spike_in.out.sorted_bams
    indexed_bams_ch = samtools_sort_spike_in.out.indexed_bams
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
    
    /*flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    // if you want to filter black list use the param --BL in the command line when calling nextflow
    
    // dont want to do bl yet
    // if ( params.BL ) {

    //     // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
    //     // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
    //     //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
        
    //     //blacklist_ch = Channel.value(params.blacklist_path)

    //     bedtools_filt_blacklist_spike_in(bam_index_tuple_ch, blacklist_ch)

    //     bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
    //     // i will need to index the black list filtered bam again so i have to create a different samtools process for this
    //     samtools_bl_index(bl_filt_bams_ch)

    //     bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

    //     /*if ( params.ATAC ) {


    //         // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

    //         deeptools_aln_shift(bam_index_tuple_ch)

    //         atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
    //         atac_shift_bam_ch.view()

    //         // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

    //         samtools_index_sort_spike_in(atac_shift_bam_ch )

    //         // so now name the tuple channel output appropriately 
    //         atac_shift_bam_index_ch = samtools_index_sort_spike_in.out.bl_filt_bam_index_tuple

    //         // now making the bed files for atac seq

    //         deeptools_make_bed(atac_shift_bam_index_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



    //     }
    //     else {



    //         // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

    //         deeptools_make_bed(bam_index_tuple_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


    //     }*/

    //     // then i need to pass the indexed_bl_bam and the bam to the deeptools process

    //     //deeptools_make_bed(bam_index_tuple_ch)
    //     //spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

    // }

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, t7_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, t7_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }
    


    // i will either use the bams or the bed files for any future processes depending on what tool needs what.



}


workflow se_lambda_spike_in_workflow {


    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch


    main:


    //params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    //se_reads_files = Channel.fromPath(params.single_end_reads).take(3)

    if (params.test) {
        se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
        
    
    }else {

        se_reads_files = Channel.fromPath(params.single_end_reads)
        

    }
    
    // now let's get the basename of the single end reads
    // removing both the .gz and the .fastq
    // I would normally use file.baseName here but it had the .gz and the .fastq
    se_reads_files.flatten()
                    .map{ file -> file.name.replace('.fastq.gz','')}
                    .set{se_reads_name}
    
    // let's view both the files and the names to make sure they match in order
    //se_reads_files.view()
    //se_reads_name.view()
    // this is where i send both the input file and their corresponding basenames to the fastp_SE process
    

    // if the adapter sequence is known then input it as a string if not dont use the parameter
    if ( params.ada_seq ) {

        //params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
        adapter_ch = Channel.value(params.adapter_seq_str)

        fastp_SE_adapter_known_spike_in(se_reads_files, se_reads_name, adapter_ch, lambda_genome_tuple) // will have to make a new process for if the adapter sequence is known

        fastq_filts = fastp_SE_adapter_known_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_adapter_known_spike_in.out.fastp_html_reports

    }    
    else {

        fastp_SE_spike_in(se_reads_files, se_reads_name, lambda_genome_tuple)

            // take all of the filtered fastq files and put them in a channel name
        // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
        fastq_filts = fastp_SE_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_spike_in.out.fastp_html_reports



    }
    

    //fastp_SE(se_reads_files, se_reads_name) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

    

    
    //fastp_SE.out.view()

    //fastq_filts.view()
    //fastp_filt_html.view()
    //fastp_filt_html.collect().view()


    // now creating a fastqc process
    fastqc_SE_spike_in(fastq_filts, fastq_filts_name, lambda_genome_tuple)

    fastqc_html_files = fastqc_SE_spike_in.out.fastqc_htmls
    fastqc_zips = fastqc_SE_spike_in.out.fastqc_zip_files

    // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
    multiqc_SE_spike_in(fastqc_zips.collect(), lambda_genome_tuple)

    // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
    bwa_index_genome_spike_in(lambda_genome_tuple)

    // collecting the genome index files from the last process 
    // not sure if i should keep track of the order the files are in first
    // it looks like they are in the same order that they appeared in the published dir using ll
    //bwa_index_genome.out.genome_index_files.view()
    genome_index_files_ch = bwa_index_genome_spike_in.out.genome_index_files

    // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
    // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
    bwa_align_SE_spike_in(genome_index_files_ch, fastq_filts, fastq_filts_name, lambda_genome_tuple )

    //bwa_align_SE.out.sam_se_files.view()

    // making a channel for the sam files generated
    sam_files = bwa_align_SE_spike_in.out.sam_se_files

    // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
    // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
    // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
    
    samtools_sort_spike_in(sam_files, lambda_genome_tuple) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

    //samtools_sort_spike_in.out.sorted_bams.view()
    //samtools_sort_spike_in.out.indexed_bams.view()
    //samtools_sort_spike_in.out.bam_index_tuple.view()

    sorted_bams_ch = samtools_sort_spike_in.out.sorted_bams
    indexed_bams_ch = samtools_sort_spike_in.out.indexed_bams
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
    
    /*
    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    // if you want to filter black list use the param --BL in the command line when calling nextflow
    
    // dont want to do bl yet
    // if ( params.BL ) {

    //     // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
    //     // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
    //     //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
        
    //     //blacklist_ch = Channel.value(params.blacklist_path)

    //     bedtools_filt_blacklist_spike_in(bam_index_tuple_ch, blacklist_ch)

    //     bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
    //     // i will need to index the black list filtered bam again so i have to create a different samtools process for this
    //     samtools_bl_index(bl_filt_bams_ch)

    //     bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

    //     /*if ( params.ATAC ) {


    //         // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

    //         deeptools_aln_shift(bam_index_tuple_ch)

    //         atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
    //         atac_shift_bam_ch.view()

    //         // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

    //         samtools_index_sort_spike_in(atac_shift_bam_ch )

    //         // so now name the tuple channel output appropriately 
    //         atac_shift_bam_index_ch = samtools_index_sort_spike_in.out.bl_filt_bam_index_tuple

    //         // now making the bed files for atac seq

    //         deeptools_make_bed(atac_shift_bam_index_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



    //     }
    //     else {



    //         // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

    //         deeptools_make_bed(bam_index_tuple_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


    //     }*/

    //     // then i need to pass the indexed_bl_bam and the bam to the deeptools process

    //     //deeptools_make_bed(bam_index_tuple_ch)
    //     //spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

    // }

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, lambda_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, lambda_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }
    


    // i will either use the bams or the bed files for any future processes depending on what tool needs what.



}

workflow se_yeast_spike_in_workflow {


    emit:
    spike_in_bam_index_tuple_ch
    spike_in_bed_files_norm_ch


    main:


    //params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    //se_reads_files = Channel.fromPath(params.single_end_reads).take(3)

    if (params.test) {
        se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
        
    
    }else {

        se_reads_files = Channel.fromPath(params.single_end_reads)
        

    }

    
    // now let's get the basename of the single end reads
    // removing both the .gz and the .fastq
    // I would normally use file.baseName here but it had the .gz and the .fastq
    se_reads_files.flatten()
                    .map{ file -> file.name.replace('.fastq.gz','')}
                    .set{se_reads_name}
    
    // let's view both the files and the names to make sure they match in order
    //se_reads_files.view()
    //se_reads_name.view()
    // this is where i send both the input file and their corresponding basenames to the fastp_SE process
    

    // if the adapter sequence is known then input it as a string if not dont use the parameter
    if ( params.ada_seq ) {

        //params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
        adapter_ch = Channel.value(params.adapter_seq_str)

        fastp_SE_adapter_known_spike_in(se_reads_files, se_reads_name, adapter_ch, yeast_genome_tuple) // will have to make a new process for if the adapter sequence is known

        fastq_filts = fastp_SE_adapter_known_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_adapter_known_spike_in.out.fastp_html_reports

    }    
    else {

        fastp_SE_spike_in(se_reads_files, se_reads_name, yeast_genome_tuple)

            // take all of the filtered fastq files and put them in a channel name
        // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
        fastq_filts = fastp_SE_spike_in.out.filtered_fastqs
        //fastp_SE.out.view()
        fastq_filts.map{file -> file.baseName}
                    .set{fastq_filts_name}

        // now getting the html files since i think fastqc combines them into one, that might be multiqc
        fastp_filt_html = fastp_SE_spike_in.out.fastp_html_reports



    }
    

    //fastp_SE(se_reads_files, se_reads_name) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

    

    
    //fastp_SE.out.view()

    //fastq_filts.view()
    //fastp_filt_html.view()
    //fastp_filt_html.collect().view()


    // now creating a fastqc process
    fastqc_SE_spike_in(fastq_filts, fastq_filts_name, yeast_genome_tuple)

    fastqc_html_files = fastqc_SE_spike_in.out.fastqc_htmls
    fastqc_zips = fastqc_SE_spike_in.out.fastqc_zip_files

    // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
    multiqc_SE_spike_in(fastqc_zips.collect(), yeast_genome_tuple)

    // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
    bwa_index_genome_spike_in(yeast_genome_tuple)

    // collecting the genome index files from the last process 
    // not sure if i should keep track of the order the files are in first
    // it looks like they are in the same order that they appeared in the published dir using ll
    //bwa_index_genome.out.genome_index_files.view()
    genome_index_files_ch = bwa_index_genome_spike_in.out.genome_index_files

    // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
    // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
    bwa_align_SE_spike_in(genome_index_files_ch, fastq_filts, fastq_filts_name, yeast_genome_tuple )

    //bwa_align_SE.out.sam_se_files.view()

    // making a channel for the sam files generated
    sam_files = bwa_align_SE_spike_in.out.sam_se_files

    // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
    // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
    // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
    
    samtools_sort_spike_in(sam_files, yeast_genome_tuple) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

    //samtools_sort_spike_in.out.sorted_bams.view()
    //samtools_sort_spike_in.out.indexed_bams.view()
    //samtools_sort_spike_in.out.bam_index_tuple.view()

    sorted_bams_ch = samtools_sort_spike_in.out.sorted_bams
    indexed_bams_ch = samtools_sort_spike_in.out.indexed_bams
    spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
    
    /*
    flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
    */
    // if you want to filter black list use the param --BL in the command line when calling nextflow
    
    // dont want to do bl yet
    // if ( params.BL ) {

    //     // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
    //     // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
    //     //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
        
    //     //blacklist_ch = Channel.value(params.blacklist_path)

    //     bedtools_filt_blacklist_spike_in(bam_index_tuple_ch, blacklist_ch)

    //     bl_filt_bams_ch = bedtools_filt_blacklist_spike_in.out.bl_filtered_bams
    //     // i will need to index the black list filtered bam again so i have to create a different samtools process for this
    //     samtools_bl_index(bl_filt_bams_ch)

    //     bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

    //     /*if ( params.ATAC ) {


    //         // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

    //         deeptools_aln_shift(bam_index_tuple_ch)

    //         atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
    //         atac_shift_bam_ch.view()

    //         // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

    //         samtools_index_sort_spike_in(atac_shift_bam_ch )

    //         // so now name the tuple channel output appropriately 
    //         atac_shift_bam_index_ch = samtools_index_sort_spike_in.out.bl_filt_bam_index_tuple

    //         // now making the bed files for atac seq

    //         deeptools_make_bed(atac_shift_bam_index_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



    //     }
    //     else {



    //         // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

    //         deeptools_make_bed(bam_index_tuple_ch)

    //         deeptools_make_bed.out.bed_files_normalized.view()

    //         spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


    //     }*/

    //     // then i need to pass the indexed_bl_bam and the bam to the deeptools process

    //     //deeptools_make_bed(bam_index_tuple_ch)
    //     //spike_in_bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

    // }

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort_spike_in(atac_shift_bam_ch, yeast_genome_tuple)

        // this will overwrite the spike_in_bam_index_tuple_ch if the data was atac data
        spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch, yeast_genome_tuple)

        //deeptools_make_bed.out.bed_files_normalized.view()

        spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized        


    }
    


    // i will either use the bams or the bed files for any future processes depending on what tool needs what.



}