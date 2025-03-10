
// making a parameter to put all the break density files in the process maybe


// note: this workflow is in the workflow dir so i need to go back one dir to get the modules dir
include {
    breakDensityWrapper_process;
    mk_break_points;
    break_concat_results




} from '../modules/fastq2bam_dna_modules.nf'

workflow breakDensityWrapper_workflow {

    // named workflows should receive their inputs explicitly though the take section
    take:
    bam_index_tuple 
    peak_files


    main:
    //bam_index_tuple.view() // if i collect everything the script should only look for the bam files anyway

    // this collects all of the elements in the channel, flattens them, and then gets only the bam files
    bam_index_tuple.collect()
                .flatten()
                .filter(~/.*bam/)
                .set{only_bams}

    //only_bams.view()

    // this workflow will take the files from the main workflow and pass them to the breakDensityWrapper process/module
    // if (params.test) {
    //     breakDensityWrapper_process(only_bams.collect(), peak_files.flatten().take(1)) // just for test, but i just want to take 1 peak file
    // }
    // else {
    //     breakDensityWrapper_process(only_bams.collect(), peak_files.flatten())

    // }

    // test putting all the bam index tuple in the process

    // bam_index_tuple.collect().flatten().view()
    // bam_index_tuple
    //     .multiMap{
    //         file -> 
    //         bams: file(file[0]).renameTo("${file[0].Parent}/${file[0].baseName}.sorted.bam") // this doesnt rename the file it makes a new file name "${file[0].Parent}/${file[0].baseName}.sorted.bam" //filter(~/.*bam/) not sure why filter doesnt do what i wanted
    //         index: file(file[1]).renameTo("${file[1].Parent}/${file[1].baseName}.sorted.bam") // this doesnt rename the file it makes a new file name "${file[1].Parent}/${file[1].baseName}.sorted.bam.bai" //filter(~/.*bai/)
    //     }
    //     .set{multi_bam_index_ch}
    // multi_bam_index_ch.bams.view{file -> "bams $file"}
    // multi_bam_index_ch.index.view{file -> "index $file"}
    
    
    
  
    bam_index_tuple
        .multiMap{
            file -> 
            bams:  file[0] //filter(~/.*bam/) not sure why filter doesnt do what i wanted
            index: file[1] //filter(~/.*bai/)
        }
        .set{multi_bam_index_ch}
    //multi_bam_index_ch.bams.view{file -> "bams $file"}
    //multi_bam_index_ch.index.view{file -> "index $file"}
    
    // downstream of this workflow I need to have a process that will sort the bam files by read group -t RG using samtools sort
    // // I might not need this process since I can just use  the -t RG option in the other sort process and it will first sort by the tag then the coordinate. so thats what i want, for it to be coordinate sorted but also have some kind of tag sort
    //rg_sort(multi_bam_index_ch.bams)


    // getting the break points sorted bed file first
    mk_break_points(multi_bam_index_ch.bams)

    //break_sorted_bed = mk_break_points.out.sorted_break_bed
    norm_break_files = mk_break_points.out.break_files
    // find a way to get all bams in one channel. probaly have to collect bams then collect index using multi map to have the same order
    // if using andrews scripts the break files still arent found so it still calculates its own.
    //peak_files.collect().view()

    if (params.test) {
        breakDensityWrapper_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files.take(1)) // just for test, but i just want to take 1 peak file
    }
    else {
        breakDensityWrapper_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files) // instead of collecting the peak files just have each process use all the bams and process on one peak file. then i can combine those output files

        // make a channel to collect the breakDensityWrapper_process outputs
        adj_enrich_tsv_files = breakDensityWrapper_process.out.adjusted_E_tsv
        break_plot_pdfs = breakDensityWrapper_process.out.break_plot_pdf
        density_calc_logs = breakDensityWrapper_process.out.density_calc_log

        // now using the process to get all the files produced and concatnate them
        // This process will concat the tsv files to each other and the log files to eachother

        break_concat_results(adj_enrich_tsv_files.collect(), density_calc_logs.collect())
    }

}



