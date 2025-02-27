
// i need to put this before the modules are loaded 
// try this here
// params.spike_name = params.gloe_seq ? '_t7' : (params.end_seq ? '_lambda' : (params.ricc_seq ? '_yeast': ''))

// // I need access to the processes so I put them in the module file

// include {
//     fastp_SE_adapter_known_spike_in;
//     fastp_SE_spike_in;
//     fastqc_SE_spike_in;
//     multiqc_SE_spike_in;
//     bwa_index_genome_spike_in;
//     bwa_align_SE_spike_in;
//     samtools_sort_spike_in;
//     deeptools_make_bed_spike_in;
//     fastp_PE_spike_in;
//     fastqc_PE_spike_in;
//     multiqc_PE_spike_in;
//     bwa_PE_aln_spike_in;
//     samtools_index_sort_spike_in;
//     deeptools_aln_shift_spike_in
    
// } from '../modules/spike_in_modules.nf'



// // used the t7 genome for the gloe_seq run of spike in
// params.t7_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/T7_Phage/genome/Sequence/WholeGenomeFasta/genome.fa')
// t7_phage_genome_ch = Channel.value(params.t7_genome)

// // next have to do ricc_seq runs with the yeast genome
// params.yeast_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/S_pombe_EF2/genome/Sequence/WholeGenomeFasta/genome.fa')
// yeast_genome_ch = Channel.value(params.yeast_genome)

// // now for the lambda with end_seq
// params.lambda_genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/Lambda_cI857ind_1_Sam_7/genome/Sequence/Bowtie2Index/genome.fa')
// lambda_genome_ch = Channel.value(params.lambda_genome)

// // maybe this has to be outside of the workflow to work
// //params.spike_name = params.gloe_seq ? '_t7' : (params.end_seq ? '_lambda' : (params.ricc_seq ? '_yeast': ''))

// workflow spike_in_runs_workflow {

//     emit:
//     spike_in_bam_index_tuple_ch
//     spike_in_bed_files_norm_ch



//     main:

//     // gloe_seq is pair end so if I copy and paste the pe path just say params.gloe_seq

//     // might have to initialize to a defualt
//     /*def spike_name = '' 
//     if (params.gloe_seq) {
//         spike_name =  '_t7' 
//     } else if (params.end_seq ) {
//         spike_name =  '_lambda' 
//     } else if (params.ricc_seq ) { 
//         spike_name = '_yeast'
//     }*/

//     //params.spike_name = params.gloe_seq ? '_t7' : (params.end_seq ? '_lambda' : (params.ricc_seq ? '_yeast': ''))

//     // if the parameters are gloe_seq and pair end then use the pair end workflow and reads to align to t7 and lambda genome
//     if (params.gloe_seq) {
        
//         // testing something where I make a variable based on the spike in being ran
//         //spike_name = '_t7' // this will be placed dynamically in the process as part of the output file names and dir

//         // this will take the paired end reads and keep them together
//         params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

//         pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

//         //pe_fastqs_ch.view()

//         fastp_PE_spike_in(pe_fastqs_ch.take(3))

//         // checking the channels to see if everything works
//         //fastp_PE_spike_in.out.filt_PE_tuple.view()

//         //fastp_PE_spike_in.out.html_fastp_out.view()
//         //fastp_PE_spike_in.out.failed_reads_out.view()
//         //fastp_PE_spike_in.out.

//         pe_filt_tuple_ch = fastp_PE_spike_in.out.filt_PE_tuple

//         fastqc_PE_spike_in(pe_filt_tuple_ch)

//         //fastqc_PE_spike_in.out.fastqc_zip_files.view()


//         // then collect them so i can view in one file using multiqc.
//         collection_fastqc_ch =fastqc_PE_spike_in.out.fastqc_zip_files.collect()

//         multiqc_PE_spike_in(collection_fastqc_ch)

//         //multiqc_PE.out.summary_of_PE_filt.view()

//         // first use the process to index the reference genome since the process exists already for the se

//         bwa_index_genome_spike_in(t7_phage_genome_ch)

//         //bwa_index_genome_spike_in.out.genome_index_files.view()

//         genome_index_ch = bwa_index_genome_spike_in.out.genome_index_files

//         // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

//         //pe_filt_tuple_ch.view()
//         //t7_phage_genome_ch.view()
//         //genome_index_ch.view()

//         bwa_PE_aln_spike_in(pe_filt_tuple_ch, t7_phage_genome_ch, genome_index_ch)

//         // now check to see if the output channels are good
//         //bwa_PE_aln_spike_in.out.pe_sam_files.view()
        
        
//         // now i need to make the parameters for  if the bam file will be blacklist filtered or not
        
//         // using this channel for both if Blacklist or not
//         sam_files_pe_ch = bwa_PE_aln_spike_in.out.pe_sam_files

//         /*if (params.BL) {

            
//             // i have to make a bam file to then use bedtools intersect to get the blacklist
//             // using the same samtools sort process found in the SE part of the pipeline
//             samtools_sort_spike_in(sam_files_pe_ch)

//             //samtools_sort_spike_in.out.bam_index_tuple.view()

//             spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
//             flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
//             norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
//             tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

//             // this will give a blacklist filtered bam but i need to index it again
//             bedtools_filt_blacklist(spike_in_bam_index_tuple_ch, blacklist_ch)
//             //bedtools_filt_blacklist.out.bl_filtered_bams.view()

//             bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
//             // so using the process to only index which means it will take the blacklist bam file
            
//             samtools_bl_index(bl_filt_bams_ch) 

//             spike_in_bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

            

//         }
//         else {
//             samtools_sort_spike_in(sam_files_pe_ch)

//             //samtools_sort_spike_in.out.bam_index_tuple.view()

//             flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
//             norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
//             tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

//             spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

            
//         }*/

//         // not running blacklist filtering on spike ins
//         samtools_sort_spike_in(sam_files_pe_ch)

//         //samtools_sort_spike_in.out.bam_index_tuple.view()

//         flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
//         norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
//         tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()

//         spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple

 
//     }


//     // since end_seq is single end, i will just use the single end path and call it params.end_seq with changes like no blacklist filtering

//     if ( params.end_seq ) {

//         // making the param dynamic option for end_seq also
//         //spike_name = '_lambda'


//         // lets get the channel for the single end reads first
//         // only use the single end read 1 data from the end seq which are already stored here: /rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs

    

//         params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
//         se_reads_files = Channel.fromPath(params.single_end_reads)
        
//         // now let's get the basename of the single end reads
//         // removing both the .gz and the .fastq
//         // I would normally use file.baseName here but it had the .gz and the .fastq
//         se_reads_files.flatten()
//                         .map{ file -> file.name.replace('.fastq.gz','')}
//                         .set{se_reads_name}
        
//         // let's view both the files and the names to make sure they match in order
//         //se_reads_files.view()
//         //se_reads_name.view()
//         // this is where i send both the input file and their corresponding basenames to the fastp_SE_spike_in process
        

//         // if the adapter sequence is known then input it as a string if not dont use the parameter
//         if ( params.ada_seq ) {

//             params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
//             adapter_ch = Channel.value(params.adapter_seq_str)

//             fastp_SE_adapter_known_spike_in(se_reads_files.take(3), se_reads_name.take(3), adapter_ch) // will have to make a new process for if the adapter sequence is known

//             fastq_filts = fastp_SE_adapter_known_spike_in.out.filtered_fastqs
//             //fastp_SE_spike_in.out.view()
//             fastq_filts.map{file -> file.baseName}
//                         .set{fastq_filts_name}

//             // now getting the html files since i think fastqc combines them into one, that might be multiqc
//             fastp_filt_html = fastp_SE_adapter_known_spike_in.out.fastp_html_reports

//         }    
//         else {

//             fastp_SE_spike_in(se_reads_files.take(3), se_reads_name.take(3))

//                 // take all of the filtered fastq files and put them in a channel name
//             // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
//             fastq_filts = fastp_SE_spike_in.out.filtered_fastqs
//             //fastp_SE_spike_in.out.view()
//             fastq_filts.map{file -> file.baseName}
//                         .set{fastq_filts_name}

//             // now getting the html files since i think fastqc combines them into one, that might be multiqc
//             fastp_filt_html = fastp_SE_spike_in.out.fastp_html_reports



//         }
        

//         //fastp_SE_spike_in(se_reads_files.take(3), se_reads_name.take(3)) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

        

        
//         //fastp_SE_spike_in.out.view()

//         //fastq_filts.view()
//         //fastp_filt_html.view()
//         //fastp_filt_html.collect().view()


//         // now creating a fastqc process
//         fastqc_SE_spike_in(fastq_filts, fastq_filts_name)

//         fastqc_html_files = fastqc_SE_spike_in.out.fastqc_htmls
//         fastqc_zips = fastqc_SE_spike_in.out.fastqc_zip_files

//         // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
//         multiqc_SE_spike_in(fastqc_zips.collect())

//         // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
//         bwa_index_genome_spike_in(lambda_genome_ch)

//         // collecting the genome index files from the last process 
//         // not sure if i should keep track of the order the files are in first
//         // it looks like they are in the same order that they appeared in the published dir using ll
//         //bwa_index_genome_spike_in.out.genome_index_files.view()
//         genome_index_files_ch = bwa_index_genome_spike_in.out.genome_index_files

//         // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
//         // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
//         bwa_align_SE_spike_in(lambda_genome_ch, genome_index_files_ch, fastq_filts, fastq_filts_name )

//         //bwa_align_SE_spike_in.out.sam_se_files.view()

//         // making a channel for the sam files generated
//         sam_files = bwa_align_SE_spike_in.out.sam_se_files

//         // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
//         // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
//         // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
        
//         samtools_sort_spike_in(sam_files) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

//         //samtools_sort_spike_in.out.sorted_bams.view()
//         //samtools_sort_spike_in.out.indexed_bams.view()
//         //samtools_sort_spike_in.out.bam_index_tuple.view()

//         sorted_bams_ch = samtools_sort_spike_in.out.sorted_bams
//         indexed_bams_ch = samtools_sort_spike_in.out.indexed_bams
//         spike_in_bam_index_tuple_ch = samtools_sort_spike_in.out.bam_index_tuple
//         flagstat_log_ch = samtools_sort_spike_in.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
//         norm_stats_txt_ch = samtools_sort_spike_in.out.norm_stats_txt.collect()
//         tsv_SN_stats_ch = samtools_sort_spike_in.out.tsv_SN_stats.collect()
//         // if you want to filter black list use the param --BL in the command line when calling nextflow
        
//     }






//     // what I am thinking is after gloe_seq, end_seq and ricc_seq runs I will get a spike_in_bam_index_tuple_ch from each (only one shoud be chosen per run ) and check if the data was atac_seq or not.
//     // so the atac_seq check will be here at the bottom

//     if (params.ATAC) {


//                 // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

//                 deeptools_aln_shift_spike_in(spike_in_bam_index_tuple_ch)

//                 atac_shift_bam_ch = deeptools_aln_shift_spike_in.out.atac_shifted_bam

//                 // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

//                 samtools_index_sort_spike_in(atac_shift_bam_ch)

//                 // so now name the tuple channel output appropriately 
//                 spike_in_bam_index_tuple_ch = samtools_index_sort_spike_in.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple. I also made it so the it doesnt sepcify its the atac-seq bam_index channel becuse having it the same name as the other spike_in_bam_index_tuple_ch means if the user specified params.ATAC there will only be the atac shifted bams in the bam channel. it gets overridden and that makes it easy for me to use channels with the same name.

//                 // now making the bed files for atac seq

//                 deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch)

//                 //deeptools_make_bed_spike_in.out.bed_files_normalized.view()

//                 spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized



//             }
//             else {



//                 // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

//                 deeptools_make_bed_spike_in(spike_in_bam_index_tuple_ch)

//                 //deeptools_make_bed_spike_in.out.bed_files_normalized.view()

//                 spike_in_bed_files_norm_ch = deeptools_make_bed_spike_in.out.bed_files_normalized          


//             }
    

    


// }