
// making a parameter to put all the break density files in the process maybe


// note: this workflow is in the workflow dir so i need to go back one dir to get the modules dir
include {
    breakDensityWrapper_process;
    mk_break_points;
    break_concat_results;
    tally_break_density;
    cell_plc_tally_break_density;
    r_heatmap;
    get_scaling_numerator;
    mk_break_points as alias_mk_break_points




} from '../modules/fastq2bam_dna_modules.nf'

include {
    breakDensityWrapper_spike_in_process;
    mk_break_points_spike_in;
    break_concat_results_spike_in

}from '../modules/spike_in_modules.nf'

workflow breakDensityWrapper_workflow {

    // named workflows should receive their inputs explicitly though the take section
    take:
    bam_index_tuple
    spike_bam_index_tuple 
    peak_files

    emit:
    complete_adj_enrich_ch
    complete_density_calc_ch


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
            bams:  file[0] //filter(~/.*bam/) not sure why filter doesnt do what i wanted. also adding the string normal to be able to have dynamic out dir
            index: file[1] //filter(~/.*bai/)
        }
        .set{multi_bam_index_ch}
    //multi_bam_index_ch.bams.view{file -> "bams $file"}
    //multi_bam_index_ch.index.view{file -> "index $file"}

    // doing this for spike in bams also
    spike_bam_index_tuple
        .multiMap{
            file -> 
            bams:  file[0] //filter(~/.*bam/) not sure why filter doesnt do what i wanted. also adding the string normal to be able to have dynamic out dir
            index: file[1] //filter(~/.*bai/)
        }
        .set{spike_multi_bam_index_ch}
    
    

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
    sorted_break_files = mk_break_points.out.sorted_break_bed
    // i would need to tokenize and use group tuple to get the cells and plc together so i can divide the cells by the plc and get a better heatmap
    

    // now get the break points for spike in also
    alias_mk_break_points(spike_multi_bam_index_ch.bams)
    spike_break_files = alias_mk_break_points.out.break_files
    spike_sorted_break_files = alias_mk_break_points.out.sorted_break_bed

    // so now i need to get make a process that gets all the counts for each spike in file and get the scaling factor score that i want to use and put it in a channel
    // I need to just get the number that represents the max out of all spike ins. in this next process numerator
    // i will find the denominator for the equation in the actual tally process below

    //get_scaling_numerator(spike_sorted_break_files.collect())



    if (params.gloe_seq){
        sorted_break_files
            .filter(~/.*(?:Cell|PLC).*/)
            .map { file -> tuple(file.baseName, file.name, file)}
            .map{file_basename, file_name, file_path -> 
                tokens = file_basename.tokenize("_")
                tuple(tokens, file_basename, file_name, file_path)
            }
            .map {tokens, basename, filename, file ->
            
                ["${tokens[0]}_${tokens[1]}",tokens[2], basename, filename, file] // I needed to recreate the first two fields to get the files that share the same base name so i can group them and put them into the workflow.
            
            }
            .groupTuple(by:0, sort:true)
            .set{grouped_cells_plc}
        //grouped_cells_plc.view()
        // the first in gloe seq is cells and the second is plc. pass this to the tally_break_density process and divide cells by plc.
    }else if (params.end_seq){

    sorted_break_files
        .filter(~/.*(?:Cell|PLC).*/)
        .map { file -> tuple(file.baseName, file.name, file)}
        .map{file_basename, file_name, file_path -> 
            tokens = file_basename.tokenize("_")
            tuple(tokens, file_basename, file_name, file_path)
        }
        .map {tokens, basename, filename, file ->
        
            ["${tokens[0]}_${tokens[3]}",tokens[1], basename, filename, file] // I needed to recreate the first two fields to get the files that share the same base name so i can group them and put them into the workflow.
        
        }
        .groupTuple(by:0, sort:true)
        .set{grouped_cells_plc}
    //grouped_cells_plc.view()
    // the first in gloe seq is cells and the second is plc. pass this to the tally_break_density process and divide cells by plc.

    }

    cell_plc_tally_break_density(grouped_cells_plc)
    cell_plc_break_counts_chr_tsv_ch = cell_plc_tally_break_density.out.tsv_break_counts_files
    cell_plc_break_counts_chr_tsv_ch
        .collectFile(name: "complete_${params.pe_or_se_name}_${params.expr_type}_break_counts_cells_vs_plc_chr.tsv", keepHeader:true, sort: true, storeDir:"${params.base_out_dir}/break_point_bed/complete_break_chr_matrix_tsv" )
        .set {complete_break_counts_cell_plc_chr_ch}
    
    r_heatmap(complete_break_counts_cell_plc_chr_ch)


    // the code below is what i used to get the tsv files and heat map for all samples. but i switched to the above so i can divide, for each group, the cells vs the plc to make the heatmap more accurate at showing the biology
    //tally_break_density(sorted_break_files)
    //break_counts_chr_tsv_ch = tally_break_density.out.tsv_break_counts_files

    // break_counts_chr_tsv_ch
    //     .collectFile(name: "complete_${params.pe_or_se_name}_${params.expr_type}_break_counts_chr.tsv", keepHeader:true, sort: true, storeDir:"${params.base_out_dir}/break_point_bed/complete_break_chr_matrix_tsv" )
    //     .set {complete_break_counts_chr_ch}

    // now that i have the tsv of the break points for both single end and pair end also gloe and end seq
    // lets build the process to make the heatmap
    //r_heatmap(complete_break_counts_chr_ch)

    if (params.test) {
        breakDensityWrapper_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files.take(1)) // just for test, but i just want to take 1 peak file
    }
    else {
        breakDensityWrapper_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files) // instead of collecting the peak files just have each process use all the bams and process on one peak file. then i can combine those output files

        // make a channel to collect the breakDensityWrapper_process outputs
        adj_enrich_tsv_files = breakDensityWrapper_process.out.adjusted_E_tsv
        //break_plot_pdfs = breakDensityWrapper_process.out.break_plot_pdf
        density_calc_logs = breakDensityWrapper_process.out.density_calc_log

        // now using the process to get all the files produced and concatnate them
        // This process will concat the tsv files to each other and the log files to eachother

        break_concat_results(adj_enrich_tsv_files.collect(), density_calc_logs.collect())
        
        complete_adj_enrich_ch = break_concat_results.out.complete_adj_enrichment
        complete_density_calc_ch = break_concat_results.out.complete_density_calc
    }

}

workflow breakDensityWrapper_spike_in_workflow {

    // named workflows should receive their inputs explicitly though the take section
    take:
    bam_index_tuple 
    peak_files

    emit:
    complete_adj_enrich_ch
    complete_density_calc_ch


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
    mk_break_points_spike_in(multi_bam_index_ch.bams)

    //break_sorted_bed = mk_break_points.out.sorted_break_bed
    norm_break_files = mk_break_points_spike_in.out.break_files
    // find a way to get all bams in one channel. probaly have to collect bams then collect index using multi map to have the same order
    // if using andrews scripts the break files still arent found so it still calculates its own.
    //peak_files.collect().view()
    sorted_break_files = mk_break_points_spike_in.out.sorted_break_bed

    tally_break_density(sorted_break_files)

    break_counts_chr_tsv_ch = tally_break_density.out.tsv_break_counts_files

    break_counts_chr_tsv_ch
        .collectFile(name: 'complete_spike_in_break_counts_chr.tsv', keepHeader:true, sort: true, storeDir:"${params.base_out_dir}/break_point_bed/spike_in/complete_break_chr_matrix_tsv" )
        .set {complete_break_counts_chr_ch}



    if (params.test) {
        breakDensityWrapper_spike_in_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files.take(1)) // just for test, but i just want to take 1 peak file
    }
    else {
        breakDensityWrapper_spike_in_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files) // instead of collecting the peak files just have each process use all the bams and process on one peak file. then i can combine those output files

        // make a channel to collect the breakDensityWrapper_process outputs
        adj_enrich_tsv_files = breakDensityWrapper_spike_in_process.out.adjusted_E_tsv
        break_plot_pdfs = breakDensityWrapper_spike_in_process.out.break_plot_pdf
        density_calc_logs = breakDensityWrapper_spike_in_process.out.density_calc_log

        // now using the process to get all the files produced and concatnate them
        // This process will concat the tsv files to each other and the log files to eachother

        break_concat_results_spike_in(adj_enrich_tsv_files.collect(), density_calc_logs.collect())
        
        complete_adj_enrich_ch = break_concat_results_spike_in.out.complete_adj_enrichment
        complete_density_calc_ch = break_concat_results_spike_in.out.complete_density_calc
    }


    


}



