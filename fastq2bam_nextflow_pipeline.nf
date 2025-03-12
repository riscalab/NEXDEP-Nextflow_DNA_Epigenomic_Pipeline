// first i will activate dsl2
nextflow.enable.dsl=2


// params.pe_out_dir = './results_PE'
// params.se_out_dir = './results_SE'

params.base_out_dir = params.PE ? './results_PE' : (params.SE ? './results_SE' : '')
// finally using modules



include {
    fastp_SE_adapter_known;
    fastp_SE;
    fastqc_SE;
    multiqc_SE;
    bwa_index_genome;
    bwa_align_SE;
    samtools_sort;
    bam_log_calc;
    deeptools_make_bed;
    bedtools_filt_blacklist;
    samtools_bl_index;
    fastp_PE;
    fastqc_PE;
    multiqc_PE;
    bwa_PE_aln;
    multiqc_bam_stats;
    deeptools_aln_shift;
    samtools_index_sort;
    mk_break_points;
    breakDensityWrapper_process;
    py_calc_stats_log;





} from './modules/fastq2bam_dna_modules.nf'



//include {samtools_index_sort} from './modules/fastq2bam_dna_modules.nf'



// i want to include another workflow script

include {breakDensityWrapper_workflow} from './workflows/breakDensityWrapper_workflow.nf'

//include {py_calc_stats_log} from './modules/fastq2bam_dna_modules.nf'

//include {spike_in_runs_workflow} from './workflows/spike_in_workflow.nf'
include {
    pe_t7_spike_in_workflow;
    pe_lambda_spike_in_workflow;
    pe_yeast_spike_in_workflow;
    se_t7_spike_in_workflow;
    se_lambda_spike_in_workflow;
    se_yeast_spike_in_workflow  
            
            
} from './workflows/2nd_spike_in_workflow.nf'

include {
    align_depth_in_peaks_workflow

} from './workflows/align_depth_in_peaks_workflow.nf'


include {
    align_depth_in_peaks_spike_in_workflow
}from './workflows/align_depth_in_peaks_spike_in_workflow.nf'


workflow {

    // this is the end seq alignment steps first


    // i will use a path already in the hpc as the defualt human genome but the user can change the genome by using -genome parameter and putting the path to a new genome in the command line when calling nextflow run
    // this hg19 genome that i wanted to use did not have mitochondrial chromosome.
    //params.genome = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome.fa')

    // this genome version to use as default will have the mitochondrial genome and it was downloaded from ucsc.
    // i will create a process that will download all of the needed genomes and give the user acces to choose which one through the use of parameters. At a later time
    // I like it this way anyway, becasue the genome now has an actual name I can use if I need to reference it in the pipeline.
    //params.genome = file('/lustre/fs4/home/rjohnson/downloads/genomes/hg19/hg19.p13.plusMT.fa')

    // lets try this
    //params.genome = file('/lustre/fs4/home/rjohnson/downloads/genomes/hg19/hg19.p13.plusMT_M_both.fa')

    //next try this. it didnt work as the mitochondrial genome either
    //params.genome = file('/lustre/fs4/home/rjohnson/downloads/genomes/hg19/hg19.p13.plusMT_only2.fa')
    
    // trying the analysis set recommended by ucsc
    params.genome = file('/lustre/fs4/home/rjohnson/downloads/genomes/analysis_set_hg19/hg19.p13.plusMT.no_alt_analysis_set.fa')

    // this is the path to hg38 /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa
    // putting the human genome in a channel
    // keeping the human genome in a value channel so i can have other processes run more than once.
    genome_ch = Channel.value(params.genome)

    // hopefully uscs has a corresponding blacklist bed file I can use.
    // this is the only one i see on ucsc. I dont think theres a specific version for hg19 with mitochondrial
    params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
                
    blacklist_ch = Channel.value(params.blacklist_path)

    // i want to add an if then logic to the pipeline so i know which type of reads are comming in paired end or single end

    

    // params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    // params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        

    // // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
    // // or it will run in normal mode where you want all your data processed
    // if (params.test) {
    //     se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
    //     pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
    
    // }else {

    //     se_reads_files = Channel.fromPath(params.single_end_reads)
    //     pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

    // }
    
    
    if ( params.SE ) {

        
            

            

            // lets get the channel for the single end reads first
            // only use the single end read 1 data from the end seq which are already stored here: /rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs

            params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
    

        

            // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
            // or it will run in normal mode where you want all your data processed
            if (params.test) {
                se_reads_files = Channel.fromPath(params.single_end_reads).take(3)
                
            
            }else {

                se_reads_files = Channel.fromPath(params.single_end_reads)
               

            }



            //params.single_end_reads = file('/rugpfs/fs0/risc_lab/store/hcanaj/HC_ENDseq_Novaseq_010925/read1_fastqs/*_1.fastq.gz')
            //se_reads_files = Channel.fromPath(params.single_end_reads)
            
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

                params.adapter_seq_str = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' // this is just a place holder value for the adapter sequence
                adapter_ch = Channel.value(params.adapter_seq_str)

                fastp_SE_adapter_known(se_reads_files, se_reads_name, adapter_ch) // will have to make a new process for if the adapter sequence is known

                fastq_filts = fastp_SE_adapter_known.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE_adapter_known.out.fastp_html_reports

            }    
            else {

                fastp_SE(se_reads_files, se_reads_name)

                    // take all of the filtered fastq files and put them in a channel name
                // since the fastq files might be in a different order, if i need to get their base names I will have to do it from this new channel below
                fastq_filts = fastp_SE.out.filtered_fastqs
                //fastp_SE.out.view()
                fastq_filts.map{file -> file.baseName}
                            .set{fastq_filts_name}

                // now getting the html files since i think fastqc combines them into one, that might be multiqc
                fastp_filt_html = fastp_SE.out.fastp_html_reports



            }
            

            //fastp_SE(se_reads_files, se_reads_name) // REMEMBER TO REMOVE THIS TESTING FEATURE WHERE IT WILL ONLY TAKE THE FIRST 3

            

            
            //fastp_SE.out.view()

            //fastq_filts.view()
            //fastp_filt_html.view()
            //fastp_filt_html.collect().view()


            // now creating a fastqc process
            fastqc_SE(fastq_filts, fastq_filts_name)

            fastqc_html_files = fastqc_SE.out.fastqc_htmls
            fastqc_zips = fastqc_SE.out.fastqc_zip_files

            // now using multiqc to combine all of the se zip files. Multiqc takes the zip files generated by fastqc and puts them in a single html file
            multiqc_SE(fastqc_zips.collect())

            // first have a seprate process that indexes the reference genome using bwa or bwa mem so this part doesnt have to be done again and will be cached
            bwa_index_genome(genome_ch)

            // collecting the genome index files from the last process 
            // not sure if i should keep track of the order the files are in first
            // it looks like they are in the same order that they appeared in the published dir using ll
            //bwa_index_genome.out.genome_index_files.view()
            genome_index_files_ch = bwa_index_genome.out.genome_index_files

            // Now I need to pass the human genome file to the process to index the genome file. Also I will add the filtered fastq files from fastp into this process that will be aligned to the genome. 
            // each run of this only takes 20-30 min to run but since the hpc only is allowing 2-3 to run at one time it takes 3 hours
            bwa_align_SE(genome_ch, genome_index_files_ch, fastq_filts, fastq_filts_name )

            //bwa_align_SE.out.sam_se_files.view()

            // making a channel for the sam files generated
            sam_files = bwa_align_SE.out.sam_se_files

            // now I want to take any sam files generated by the bwa and use samtools to order them and convert them into bam files
            // I will hopefully be able to do this outside of the if else statement so the sam file from both conditions can be passed to the same samtools process
            // since i am just testing the pipeline i should find a way to do this on only a few files (about 3-4)
            
            samtools_sort(sam_files) // using take 3 should only take the first 3 files from the sam_files channel which should have 64 sam files. This is just for production and testing. will remove when running pipeline for real.

            //samtools_sort.out.sorted_bams.view()
            //samtools_sort.out.indexed_bams.view()
            //samtools_sort.out.bam_index_tuple.view()

            sorted_bams_ch = samtools_sort.out.sorted_bams
            indexed_bams_ch = samtools_sort.out.indexed_bams
            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple
            // flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            // norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            // tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()
            // if you want to filter black list use the param --BL in the command line when calling nextflow
            if ( params.BL ) {

                // using bedtools to filter black list but first giving the user an option to put the correct black list for an organism
                // by defualt it will use the hg19 v2 blacklist, but if you used a different organism or human genome use the appropriate blacklist
                //params.blacklist_path = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed')
                
                //blacklist_ch = Channel.value(params.blacklist_path)

                bedtools_filt_blacklist(bam_index_tuple_ch, blacklist_ch)

                bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
                // i will need to index the black list filtered bam again so i have to create a different samtools process for this
                samtools_bl_index(bl_filt_bams_ch)

                bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

                /*if ( params.ATAC ) {


                    // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                    deeptools_aln_shift(bam_index_tuple_ch)

                    atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
                    atac_shift_bam_ch.view()

                    // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                    samtools_index_sort(atac_shift_bam_ch )

                    // so now name the tuple channel output appropriately 
                    atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                    // now making the bed files for atac seq

                    deeptools_make_bed(atac_shift_bam_index_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



                }
                else {



                    // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                    deeptools_make_bed(bam_index_tuple_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }*/

                // then i need to pass the indexed_bl_bam and the bam to the deeptools process

                //deeptools_make_bed(bam_index_tuple_ch)
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

            }

            /*
            else {


                if (params.ATAC) {


                    // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                    deeptools_aln_shift(bam_index_tuple_ch)

                    atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

                    // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                    samtools_index_sort(atac_shift_bam_ch)

                    // so now name the tuple channel output appropriately 
                    atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                    // now making the bed files for atac seq

                    deeptools_make_bed(atac_shift_bam_index_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



                }
                else {



                    // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                    deeptools_make_bed(bam_index_tuple_ch)

                    deeptools_make_bed.out.bed_files_normalized.view()

                    bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


                }

                // Now i want to pass the tuple that has the bam and its corresponding index file into a process that will create a bigwig file for visulization, created from the bam file. This will show read coverage in the genome without looking for significant areas
                //deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()
                //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized
            }*/
            
    

            // i will either use the bams or the bed files for any future processes depending on what tool needs what.
    }

    // now I need to make the paired-end part of this pipeline.
    else if (params.PE) {

        // this will take the paired end reads and keep them together
        // params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

        // pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        
        params.paired_end_reads = '/rugpfs/fs0/risc_lab/store/hcanaj/HC_GLOEseq_Novaseq_010925/fastqs_read1_read2/*_{R1,R2}*'

            

        // this will give the user to run in test mode where the pipeline will only take 3 of the fastq files in the directory full of fastq files
        // or it will run in normal mode where you want all your data processed
        if (params.test) {
            
            pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads).take(3)
        
        }else {

            
            pe_fastqs_ch = Channel.fromFilePairs(params.paired_end_reads)

        }


        //pe_fastqs_ch.view()

        fastp_PE(pe_fastqs_ch)

        // checking the channels to see if everything works
        //fastp_PE.out.filt_PE_tuple.view()

        //fastp_PE.out.html_fastp_out.view()
        //fastp_PE.out.failed_reads_out.view()
        //fastp_PE.out.

        pe_filt_tuple_ch = fastp_PE.out.filt_PE_tuple

        fastqc_PE(pe_filt_tuple_ch)

        //fastqc_PE.out.fastqc_zip_files.view()


        // then collect them so i can view in one file using multiqc.
        collection_fastqc_ch =fastqc_PE.out.fastqc_zip_files.collect()

        multiqc_PE(collection_fastqc_ch)

        //multiqc_PE.out.summary_of_PE_filt.view()

        // first use the process to index the reference genome since the process exists already for the se

        bwa_index_genome(genome_ch)

        //bwa_index_genome.out.genome_index_files.view()

        genome_index_ch = bwa_index_genome.out.genome_index_files

        // now to use bwa aln and bwa sampe to align the filtered pair end reads to the reference genome and then to create the sam file respectively

        //pe_filt_tuple_ch.view()
        //genome_ch.view()
        //genome_index_ch.view()

        bwa_PE_aln(pe_filt_tuple_ch, genome_ch, genome_index_ch)

        // now check to see if the output channels are good
        //bwa_PE_aln.out.pe_sam_files.view()
        
        
        // now i need to make the parameters for  if the bam file will be blacklist filtered or not
        
        // using this channel for both if Blacklist or not
        sam_files_pe_ch = bwa_PE_aln.out.pe_sam_files

        if (params.BL) {

            
            // i have to make a bam file to then use bedtools intersect to get the blacklist
            // using the same samtools sort process found in the SE part of the pipeline
            samtools_sort(sam_files_pe_ch)

            //samtools_sort.out.bam_index_tuple.view()

            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple
            // flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            // norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            // tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()

            // this will give a blacklist filtered bam but i need to index it again
            bedtools_filt_blacklist(bam_index_tuple_ch, blacklist_ch)
            //bedtools_filt_blacklist.out.bl_filtered_bams.view()

            bl_filt_bams_ch = bedtools_filt_blacklist.out.bl_filtered_bams
            // so using the process to only index which means it will take the blacklist bam file
            
            samtools_bl_index(bl_filt_bams_ch) 

            bam_index_tuple_ch = samtools_bl_index.out.bl_filt_bam_index_tuple

            // so this would give a bam that is bl filtered and has an index

            /*if ( params.ATAC ) {


                // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                deeptools_aln_shift(bam_index_tuple_ch)

                atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam
                //atac_shift_bam_ch.view()

                // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                samtools_index_sort(atac_shift_bam_ch )

                // so now name the tuple channel output appropriately 
                atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                // now making the bed files for atac seq

                deeptools_make_bed(atac_shift_bam_index_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



            }
            else {



                // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


            }*/
            
            
            // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

            //deeptools_make_bed(bam_index_tuple_ch)

            //deeptools_make_bed.out.bed_files_normalized.view()

            //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized

            

        }
        else {
            samtools_sort(sam_files_pe_ch)

            //samtools_sort.out.bam_index_tuple.view()

            // flagstat_log_ch = samtools_sort.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
            // norm_stats_txt_ch = samtools_sort.out.norm_stats_txt.collect()
            // tsv_SN_stats_ch = samtools_sort.out.tsv_SN_stats.collect()

            bam_index_tuple_ch = samtools_sort.out.bam_index_tuple
            

            // add the if logic for ATAC here

            /*if (params.ATAC) {


                // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

                deeptools_aln_shift(bam_index_tuple_ch)

                atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

                // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

                samtools_index_sort(atac_shift_bam_ch)

                // so now name the tuple channel output appropriately 
                atac_shift_bam_index_ch = samtools_index_sort.out.bl_filt_bam_index_tuple

                // now making the bed files for atac seq

                deeptools_make_bed(atac_shift_bam_index_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



            }
            else {



                // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

                deeptools_make_bed(bam_index_tuple_ch)

                //deeptools_make_bed.out.bed_files_normalized.view()

                bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


            }*/






            // just using the original sorted and indexed bam
            
            // i will comment this out below since it was the original but the above if logic should work if there is ATAC data or if there is not

            //deeptools_make_bed(bam_index_tuple_ch)

            //deeptools_make_bed.out.bed_files_normalized.view()

            //bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized


        }

        

 
    }



    // I think i can just put the if ATAC script here separately and just grab all of the bam tuple channels

    if (params.ATAC) {


        // now if there is atac-seq data I need to take the bam and shift the alignment. I will do this using deeptools alignmentSieve in both pair end vs single end and bl vs no bl filter

        deeptools_aln_shift(bam_index_tuple_ch)

        atac_shift_bam_ch = deeptools_aln_shift.out.atac_shifted_bam

        // now i have to re index this new atac shifted bam. dispite the name of the process I can just pass any future created bam to this channel to be indexed

        samtools_index_sort(atac_shift_bam_ch)

        // so now name the tuple channel output appropriately 
        //atac_shift_bam_index_ch = samtools_index_sort.out.bam_index_tuple // i changed the emit ch to just be bam_index_tuple

        // just make it bam_index_tuple_ch
        bam_index_tuple_ch = samtools_index_sort.out.bam_index_tuple

        // now making the bed files for atac seq

        deeptools_make_bed(bam_index_tuple_ch)

        //deeptools_make_bed.out.bed_files_normalized.view()

        bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized



    }
    else {



        // now i want to take the bl filt bam files and pass them to deep tools to be converted into bed files

        deeptools_make_bed(bam_index_tuple_ch)

        //deeptools_make_bed.out.bed_files_normalized.view()

        bed_files_norm_ch = deeptools_make_bed.out.bed_files_normalized          


    }



    // making a multiqc process for the samtools flagstat log files. this should be able to take the flagstat_log_ch from any part of the choosen paths
    //multiqc_bam_stats(flagstat_log_ch, norm_stats_txt_ch)


    
    // // I want to have a parameter that takes peakfiles. The default will be the IMR90 narrowPeak files
    // params.peak_files_IMR90 = files('/lustre/fs4/home/ascortea/store/ascortea/beds/IMR90/*.narrowPeak')
    // //now making the channel for the files
    // peak_files_ch = Channel.fromPath(params.peak_files_IMR90)

    // // now i need to make channels with other cell lines with the narrow peaks but also get the nulls. first other peaks
    // peak_file_HUVEC = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/HUVEC/*.narrowPeak')
    // peak_file_k562 = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/*.narrowPeak')
    // peak_file_BJ = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/BJ/*.narrowPeak')

    // // wondering if i can do this
    // all_peaks_ch = peak_files_ch.concat(peak_file_HUVEC, peak_file_k562, peak_file_BJ)

    // some narrow peaks are under HUVEC/gkmsvm_null_regions/
    // more narrow peaks are under IMR90/scrambles/ , but these ones are chunked bed files.

    if (params.give_peak_files == null) {


        // I want to have a parameter that takes peakfiles. The default will be the IMR90 narrowPeak files
        params.peak_files_IMR90 = files('/lustre/fs4/home/ascortea/store/ascortea/beds/IMR90/*.narrowPeak')
        //now making the channel for the files
        peak_files_ch = Channel.fromPath(params.peak_files_IMR90)

        // now i need to make channels with other cell lines with the narrow peaks but also get the nulls. first other peaks
        peak_file_HUVEC = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/HUVEC/*.narrowPeak')
        peak_file_k562 = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/k562/*.narrowPeak')
        peak_file_BJ = Channel.fromPath('/lustre/fs4/home/ascortea/store/ascortea/beds/BJ/*minimal*.narrowPeak') // making sure to use only the minimal peaks.
        //E055-DNase.macs2.narrowPeak causing problems
        // wondering if i can do this
        all_peaks_ch = peak_files_ch.concat(peak_file_HUVEC, peak_file_k562, peak_file_BJ)


    }else {

        params.give_peak_files = files('/lustre/fs4/home/ascortea/store/ascortea/beds/IMR90/*.narrowPeak')

        all_peaks_ch = Channel.fromPath(params.give_peak_files)
    }

    // putting this lower since I need to use the bam_index that has all the files including the spike-ins
    // if (params.calc_break_density){
    // // i want to call the workflow breakDensityWrapper_workflow and pass the bam_index_tuple_ch as an input from either path where the user chose to do blacklist filter or not. Then also pass the peak files that already exists or are created later as input
    
    //     if (params.PE) {
            
    //         breakDensityWrapper_workflow(bam_index_tuple_ch, all_peaks_ch)

    //     }

    //     if (params.SE) {
         
    //         breakDensityWrapper_workflow(bam_index_tuple_ch, all_peaks_ch)
  
    //     }
    // }


    // first get a multi map of the different conditions (0gyr, PLC, 200(cells))

    // bed_files_norm_ch
    //     .multiMap{file ->

    //         zero_gy: file.name ==~ /.*0Gy.*\.bed/ ? file: null
    //         PLC: file.name ==~ /.*PLC.*\.bed/ ? file: null
    //         cells: file.name ==~ /.*Cells.*\.bed/ ? file: null
    //     }
    //     .set{cell_condition}
    

    // // cell_condition.zero_gy.view{file -> "0gy: $file"}
    // // cell_condition.PLC.view{file -> "PLC: $file"}
    // // cell_condition.cells.view{file -> "cells: $file"}
    // //cell_condition.zero_gy.view()
    // // cell_condition.PLC.view()
    // // cell_condition.cells.view()

    // // filtering the null out of the channels

    // zero_gy_ch = cell_condition.zero_gy.filter{it != null}
    // plc_ch = cell_condition.PLC.filter{it != null}
    // cells_ch = cell_condition.cells.filter{it != null}


    // try this
    // if (params.depth_intersection){

    //     bed_files_norm_ch
    //         .map { file -> tuple(file.baseName, file)}
    //         .map { name, file -> 
    //             tokens = name.tokenize("_") // there are 16 fields in the tokens now. I want the 3rd field 0Gy, cells, plc
    //             tuple(tokens, name, file)
    //         }
    //         .map {tokens, name, file ->
            
    //             ["${tokens[0]}_${tokens[1]}",tokens[2], name, file] // I needed to recreate the first two fields to get the files that share the same base name so i can group them and put them into the workflow.
            
    //         }
    //         .groupTuple(by:0) // the first element of the tuple is grouped by default. it is 0 based counting.
    //         .set { grouped_bed_ch}
    //         //.view()
            

    //     // trying to cross the channels
    //     // grouped_bed_ch.view()
    //     // all_peaks_ch.view()
    //     grouped_bed_ch.combine(all_peaks_ch).set{combined_bed_peak} // this has the grouped name, condition, basename, bed files, peak files. in that order

    //     align_depth_in_peaks(combined_bed_peak)

    // }
    // I want to make a log file with all the stats from using samtools stats on each bam file

    //tsv_SN_stats_ch.view{"These are the tsv files $it"}
    // tsv_SN_stats_ch
    //     .map{ file -> tuple(file.baseName, file)}
    //     .set{tsv_SN_stats_tuple_ch}
    //     //.view {"These are the files with their basenames in a tuple transposed: $it"} // i dont need this transposed, i want all the files and names
        
    // //tsv_SN_stats_tuple_ch.view{"These are the files with their basenames in a tuple: $it"}


    // if (params.PE) {

    //     py_calc_stats_log(tsv_SN_stats_tuple_ch)
    // }
    // if (params.SE) {

    //     py_calc_stats_log(tsv_SN_stats_tuple_ch)

    // }

    
    // looking to run the spike_in workflow
    

    if (params.spike_in) {

        // all this is doing is running the normal fastq2bam2 pipeline but with the specified genomes for spike in
       
        if (params.PE) {

            if (params.t7 ) {

                pe_t7_spike_in_workflow()
                // now getting the output channels from the workflows just incase i need to use them in a downstream analysis
                pe_t7_bam_index_tuple_ch = pe_t7_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                pe_t7_bed_files = pe_t7_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.lambda) {

                pe_lambda_spike_in_workflow()

                pe_lambda_bam_index_tuple_ch = pe_lambda_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                pe_lambda_bed_files = pe_lambda_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.yeast) {

                // so i would put the yeast spike in workflow here, for example.
                pe_yeast_spike_in_workflow()

                pe_yeast_bam_index_tuple_ch = pe_yeast_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                pe_yeast_bed_files = pe_yeast_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
        }
        if (params.SE) {

            //for this i need to make the single end workflow for spike ins
            if (params.t7) {

                se_t7_spike_in_workflow()

                se_t7_bam_index_tuple_ch = se_t7_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                se_t7_bed_files = se_t7_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }
            if (params.lambda) {

                se_lambda_spike_in_workflow()
                se_lambda_bam_index_tuple_ch = se_lambda_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                se_lambda_bed_files = se_lambda_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }           
            if (params.yeast) {

                // so i would put the yeast spike in workflow here, for example.
                se_yeast_spike_in_workflow()
                se_yeast_bam_index_tuple_ch = se_yeast_spike_in_workflow.out.spike_in_bam_index_tuple_ch 
                se_yeast_bed_files = se_yeast_spike_in_workflow.out.spike_in_bed_files_norm_ch
            }

        }

        // not sure if i would need to make this channel the normal 'bam_index_tuple_ch', or do what i did with atac-seq and give it a unique name calling it 'spike_in_bam_index_ch'
        // since spike in can occur at the same time as either pe or se, i need to give this channel a unique name. I dont think spike in bam files need to have break density calculated
        //spike_in_bam_index_ch = spike_in_runs_workflow.out.spike_in_bam_index_tuple_ch
        //spike_in_bed_files_norm_ch = spike_in_runs_workflow.out.spike_in_bed_files_norm_ch

    }


    if (params.PE) {
        // join the normal aligned reads bam with the spike in bam

        // first just initilize the channels to empty channels if their parameters were not set
        // doing this works
        if (!params.t7) { pe_t7_bam_index_tuple_ch = Channel.empty(); pe_t7_bed_files = Channel.empty() }
        if (!params.lambda) { pe_lambda_bam_index_tuple_ch = Channel.empty(); pe_lambda_bed_files = Channel.empty() }
        if (!params.yeast) { pe_yeast_bam_index_tuple_ch = Channel.empty(); pe_yeast_bed_files = Channel.empty() }
        
        all_spike_bam_index_tuple_ch = bam_index_tuple_ch.concat(pe_t7_bam_index_tuple_ch,
                                                pe_lambda_bam_index_tuple_ch, 
                                                pe_yeast_bam_index_tuple_ch
        )
        
        //all_spike_bam_index_tuple_ch.view()

        // now do the same for the bed files
        only_spike_beds = pe_t7_bed_files.concat(pe_lambda_bed_files,
                                            pe_yeast_bed_files
            
        )
    }else if (params.SE) {

        //doing the same for the single end path
        if (!params.t7) { se_t7_bam_index_tuple_ch = Channel.empty(); se_t7_bed_files = Channel.empty()}
        if (!params.lambda) { se_lambda_bam_index_tuple_ch = Channel.empty(); se_lambda_bed_files = Channel.empty()}
        if (!params.yeast) { se_yeast_bam_index_tuple_ch = Channel.empty(); se_yeast_bed_files = Channel.empty()}

        all_spike_bam_index_tuple_ch = bam_index_tuple_ch.concat(
                                                se_t7_bam_index_tuple_ch,
                                                se_lambda_bam_index_tuple_ch,
                                                se_yeast_bam_index_tuple_ch
        )
        
        //all_spike_bam_index_tuple_ch.view()

        // now do the same for the bed files
        only_spike_beds = se_t7_bed_files.concat(se_lambda_bed_files,
                                                se_yeast_bed_files
            
        )
    }

    // now that i have the bed file channel that contains all the normal bed files and the bed files for the spike ins if they are specified, i want to put the workflow to get the depth from intersection here

    if (params.depth_intersection){

        bed_files_norm_ch
            .map { file -> tuple(file.baseName, file)}
            .map { name, file -> 
                tokens = name.tokenize("_") // there are 16 fields in the tokens now. I want the 3rd field 0Gy, cells, plc
                tuple(tokens, name, file)
            }
            .map {tokens, name, file ->
            
                ["${tokens[0]}_${tokens[1]}",tokens[2], name, file] // I needed to recreate the first two fields to get the files that share the same base name so i can group them and put them into the workflow.
            
            }
            .groupTuple(by:0) // the first element of the tuple is grouped by default. it is 0 based counting.
            .set { grouped_bed_ch}
            
            //.view()
            

        // trying to cross the channels
        // grouped_bed_ch.view()
        // all_peaks_ch.view()
        grouped_bed_ch.combine(all_peaks_ch).set{combined_bed_peak} // this has the grouped name, condition, basename, bed files, peak files. in that order

        //combined_bed_peak.view()
        align_depth_in_peaks_workflow(combined_bed_peak)

        // here just say if spike_in has been specified i will do everything above but just with the spike ins only 

        if (params.spike_in) {

            // probably just start by combining the spike files with the peak files and work with it that way.

            
            only_spike_beds
                .map{ bed -> tuple(bed.baseName, bed)}
                .map{name, bed -> 
                    tokens = name.tokenize("_")
                    tuple(tokens, name, bed)

                }
                .map{tokens, name , bed ->

                    ["${tokens[0]}_${tokens[1]}_${tokens[5]}", tokens[2], tokens[5], name, bed]

                }
                .groupTuple(by:0)
                .combine(all_peaks_ch)
                
                
                
                
                
                .set{combined_bed_peak_base}

            //combined_bed_peak_base.view()
            align_depth_in_peaks_spike_in_workflow(combined_bed_peak_base)


        }

    }

    // make it so that if the user specifies spike-ins then use the bam_index tuple with spike-ins, if not then only use the normal bam_index tuple

    // fix this to params.spike_in later so i can get the full thing to run
    if (params.spike_in_) {

        if (params.calc_break_density){
            // i want to call the workflow breakDensityWrapper_workflow and pass the bam_index_tuple_ch as an input from either path where the user chose to do blacklist filter or not. Then also pass the peak files that already exists or are created later as input
            
            breakDensityWrapper_workflow(all_spike_bam_index_tuple_ch, all_peaks_ch)

        }


    } else {

        if (params.calc_break_density){
                // i want to call the workflow breakDensityWrapper_workflow and pass the bam_index_tuple_ch as an input from either path where the user chose to do blacklist filter or not. Then also pass the peak files that already exists or are created later as input
            
                breakDensityWrapper_workflow(bam_index_tuple_ch, all_peaks_ch)
       
        }


    }


    // Now I want to take the bam_index_tuple or the atac_shift_bam_index_ch  and get the stats. percent mitochondrial alignment, number of reads, percent reads aligned
    // just take a few of the stats from the py_calc_sats_log

    //bam_log_calc(bam_index_tuple_ch)
    bam_log_calc(all_spike_bam_index_tuple_ch)

    flagstat_log_ch = bam_log_calc.out.flag_stats_log.collect() // will make another process or send this to the multiqc process
    norm_stats_txt_ch = bam_log_calc.out.norm_stats_txt.collect()
    tsv_SN_stats_ch = bam_log_calc.out.tsv_SN_stats.collect()

    // making a multiqc process for the samtools flagstat log files. this should be able to take the flagstat_log_ch from any part of the choosen paths
    multiqc_bam_stats(flagstat_log_ch, norm_stats_txt_ch)

    tsv_SN_stats_ch
        .map{ file -> tuple(file.baseName, file)}
        .set{tsv_SN_stats_tuple_ch}
        //.view {"These are the files with their basenames in a tuple transposed: $it"} // i dont need this transposed, i want all the files and names

    //tsv_SN_stats_tuple_ch.view{"This is the channel tuple $it"}

    //tsv_SN_stats_tuple_ch.view{"These are the files with their basenames in a tuple: $it"}


    if (params.PE) {

        py_calc_stats_log(tsv_SN_stats_tuple_ch)
    }
    if (params.SE) {

        py_calc_stats_log(tsv_SN_stats_tuple_ch)

    }
}