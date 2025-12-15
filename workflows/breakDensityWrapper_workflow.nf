
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
    mk_break_points as alias_mk_break_points;
    gloe_imr90_wrapper_process;
    gloe_imr90_wrapper_process as gloe_celltypes_wrapper_process;
    make_bigwig_gloe_process;
    filt_gloe_bam_samtools_process;
    get_ratio_cell_vs_plc_bigwig_process;
    gloe_celltypes_wrapper_scrm_process;
    endseq_celltypes_wrapper_process;
    endseq_celltypes_wrapper_scrm_process;
    make_bigwig_endseq_process;
    get_ratio_cell_vs_plc_bigwig_process as get_endseq_ratio_cell_vs_plc_bigwig_process;
    make_scrambled_peaks_process




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
    
    //comment out for now
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
    else {
        breakDensityWrapper_process(multi_bam_index_ch.bams.collect(),multi_bam_index_ch.index.collect(),norm_break_files.collect(), peak_files.take(1)) // instead of collecting the peak files just have each process use all the bams and process on one peak file. then i can combine those output files

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

// this will generate the scrambled peaks

workflow generate_scrambled_peaks_workflow {


    take:
    rpe_peaks
    blacklist_ch
    mappa_scrm_ch



    main:

    // now I want to make a process that will run andrews scrambled scripts
    make_scrambled_peaks_process(rpe_peaks.flatten(), blacklist_ch, mappa_scrm_ch)



    //emit:
}

workflow breakDensityWrapper_Gloe_workflow {

    take:
    bam_index_tuple_ch
    peak_file_imr90
    peak_file_k562
    peak_file_BJ
    peak_file_rpe1

    // scrambled peaks
    scrm_imr90_peaks
    scrm_k562_peaks
    scrm_bj_peaks
    scrm_rpe1_peaks


    main:

    // first separate the bam index tuple into a tuple that has all of one cell type

    bam_index_tuple_ch
        .map { bam, bai ->

        basename = bam.baseName
        bam_filename = bam.name

        tokens = basename.tokenize("_")

        cell_type = tokens[0]
        biorep_type = tokens[1]
        expr_type = tokens[2]

        tuple(cell_type, biorep_type, expr_type, bam_filename, bam, bai)


        }
        .groupTuple(by:0, sort:true)
        .filter { cell_type, biorep_type, expr_type, bam_filename, bam, bai ->

        cell_type in ['I', 'B', 'K', 'R']

        }
        .branch { cell_type, biorep_type, expr_type, bam_filename, bam, bai ->

        imr90: cell_type == 'I'
        bj: cell_type == 'B'
        k562: cell_type == 'K'
        rpe1: cell_type == 'R'

        }
        //.view() // reminder, we only have I(imr90), b(BJ), k(k562)
        // so i want to remove the R and undetermined
        .set{gloe_celltype}

        // gloe_celltype.imr90.view{it -> "this is imr90 branch: $it"}
        // gloe_celltype.bj.view{it -> "this is the bj cell type : $it"}
        // gloe_celltype.k562.view{it -> "this is the k562 cell type: $it"}

        // the above works with the branching

        // now to add the correct peeak set to each
        // something like using combine
        gloe_celltype.imr90
            .combine(peak_file_imr90.collect().toList())
            //.view{it -> "this is the imr90 bams and then combined with the imr90 peaks: $it"}
            .set{gloe_cell_imr90_w_peak_ch} // this is what i want. do it for the other cell types

        gloe_celltype.bj
            .combine(peak_file_BJ.collect().toList())
            .set{gloe_cell_bj_w_peak_ch}

        gloe_celltype.k562
            .combine(peak_file_k562.collect().toList())
            .set{gloe_cell_k562_w_peak_ch}

        gloe_celltype.rpe1
            .combine(peak_file_rpe1.collect().toList())
            .set{gloe_cell_rpe1_w_peak_ch}

        ///////////////////////////////////////////////////////////////////////////////
        // now doing the same but adding the scrambled peaks instead of the normal peaks
        gloe_celltype.imr90
            .combine(scrm_imr90_peaks.collect().toList())
            //.view{it -> "this is the imr90 bams and then combined with the imr90 peaks: $it"}
            .set{gloe_cell_imr90_w_scrm_peak_ch} // this is what i want. do it for the other cell types

        gloe_celltype.bj
            .combine(scrm_bj_peaks.collect().toList())
            .set{gloe_cell_bj_w_scrm_peak_ch}

        gloe_celltype.k562
            .combine(scrm_k562_peaks.collect().toList())
            .set{gloe_cell_k562_w_scrm_peak_ch}

        gloe_celltype.rpe1
            .combine(scrm_rpe1_peaks.collect().toList())
            .set{gloe_cell_rpe1_w_scrm_peak_ch}

        // then concate them like I did with the normal peaks below
        gloe_celltypes_w_scrm_peak_concat = gloe_cell_imr90_w_scrm_peak_ch.concat(gloe_cell_bj_w_scrm_peak_ch, gloe_cell_k562_w_scrm_peak_ch, gloe_cell_rpe1_w_scrm_peak_ch )

        ///////////////////////////////////////////////////////////////////////////////


        // so now I need three processes that will run each of these, or one process that will run a concatenated version of all three of these
        ////// using this or the next one for all celltypes later  ////////
        //gloe_imr90_wrapper_process(gloe_cell_imr90_w_peak_ch)
        ///////////////////////////////////////////////////////////////////
        // I will create a test process that will take all three channels in parallel
        // first concatenate them
        gloe_celltypes_w_peak_concat = gloe_cell_imr90_w_peak_ch.concat(gloe_cell_bj_w_peak_ch, gloe_cell_k562_w_peak_ch, gloe_cell_rpe1_w_peak_ch )
        // first check if it looks good
        gloe_celltypes_w_peak_concat.view(it -> "this is the concat channel with all three celltype tuples: $it")
        
        // I want to make a process that will generate the scrampled peaks and then output it into a channel like above but with scrambled peaks instead of normal peaks
        // not doing this yet

        
        // now the process as an alias
        //  WILL UNCOMMENT THIS BUT DONT NEED IT RIGHT NOW
        gloe_celltypes_wrapper_process(gloe_celltypes_w_peak_concat)

        // now run the process to get the break density for scrm gloe peaks
        // just copy the entire process but make it output with a change to file name from norm peaks to scrm peaks
        gloe_celltypes_wrapper_scrm_process(gloe_celltypes_w_scrm_peak_concat)


        // i would need to filter the bam using samtools to get only read1 for gloe seq
        filt_gloe_bam_samtools_process(bam_index_tuple_ch)
        r1_filt_gloe_bams = filt_gloe_bam_samtools_process.out.bam_filt_r1_gloe
        //r1_filt_gloe_bams.view{it -> "these are the gloe bams index tuple that are filt by read1: $it"}

        //get the bigwig files from each of the bam files
        //make_bigwig_gloe_process(bam_index_tuple_ch)
        make_bigwig_gloe_process(r1_filt_gloe_bams)

        // getting the chrmt bigwigs
        gloe_bigwig_chrmt_ch = make_bigwig_gloe_process.out.gloe_bigwig_mt

        // now to filter for only plc and cells
        gloe_bigwig_chrmt_ch
            .filter(~/.*(?:Cell|PLC).*/)
            .map{ file ->
            
            bigwig_basename = file.baseName
            bigwig_filename = file.name

            tokens = bigwig_basename.tokenize("_")

            biorep = tokens[1]
            expr_type = tokens[2]
            cell_type = tokens[0]

            tuple(cell_type, biorep, expr_type, bigwig_filename, file)

            
            }
            .groupTuple(by:[0,1], sort:true)
            .view{it -> "this is the chrmt bigwigs grouped by biorep: $it"}
            .set{gloe_bigwig_tuple_cell_plc_for_ratio_ch}

        
        // then with this output to find the ratio between cells over plc only filter for cells and plc, then group channel by the first field
        // gb1 should then have the cell and plc for that tuple
        get_ratio_cell_vs_plc_bigwig_process(gloe_bigwig_tuple_cell_plc_for_ratio_ch)

    



    //emit:
}

// now doing the same thing for end seq as i did for gloe seq above

workflow breakDensityWrapper_Endseq_workflow {

    take:
    bam_index_tuple_ch
    peak_file_imr90
    peak_file_k562
    peak_file_BJ
    peak_file_rpe1

    // scrambled peaks
    scrm_imr90_peaks
    scrm_k562_peaks
    scrm_bj_peaks
    scrm_rpe1_peaks


    main:

    // first separate the bam index tuple into a tuple that has all of one cell type

    bam_index_tuple_ch
        .map { bam, bai ->

        basename = bam.baseName
        bam_filename = bam.name

        tokens = basename.tokenize("_")

        // the biorep_type for end seq should be token 3 instead of token 1
        // and the expr_type for end seq should be token 1 instead of token 2
        cell_type = tokens[0]
        biorep_type = tokens[3]
        expr_type = tokens[1]

        tuple(cell_type, biorep_type, expr_type, bam_filename, bam, bai)


        }
        .groupTuple(by:0, sort:true)
        .filter { cell_type, biorep_type, expr_type, bam_filename, bam, bai ->

        cell_type in ['I', 'B', 'K', 'R']

        }
        .branch { cell_type, biorep_type, expr_type, bam_filename, bam, bai ->

        imr90: cell_type == 'I'
        bj: cell_type == 'B'
        k562: cell_type == 'K'
        rpe1: cell_type == 'R'

        }
        //.view() // reminder, we only have I(imr90), b(BJ), k(k562)
        // so i want to remove the R and undetermined
        .set{endseq_celltype}

        // gloe_celltype.imr90.view{it -> "this is imr90 branch: $it"}
        // gloe_celltype.bj.view{it -> "this is the bj cell type : $it"}
        // gloe_celltype.k562.view{it -> "this is the k562 cell type: $it"}

        // the above works with the branching

        // now to add the correct peeak set to each
        // something like using combine
        endseq_celltype.imr90
            .combine(peak_file_imr90.collect().toList())
            //.view{it -> "this is the imr90 bams and then combined with the imr90 peaks: $it"}
            .set{endseq_cell_imr90_w_peak_ch} // this is what i want. do it for the other cell types

        endseq_celltype.bj
            .combine(peak_file_BJ.collect().toList())
            .set{endseq_cell_bj_w_peak_ch}

        endseq_celltype.k562
            .combine(peak_file_k562.collect().toList())
            .set{endseq_cell_k562_w_peak_ch}

        endseq_celltype.rpe1
            .combine(peak_file_rpe1.collect().toList())
            .set{endseq_cell_rpe1_w_peak_ch}

        ///////////////////////////////////////////////////////////////////////////////
        // now doing the same but adding the scrambled peaks instead of the normal peaks
        endseq_celltype.imr90
            .combine(scrm_imr90_peaks.collect().toList())
            //.view{it -> "this is the imr90 bams and then combined with the imr90 peaks: $it"}
            .set{endseq_cell_imr90_w_scrm_peak_ch} // this is what i want. do it for the other cell types

        endseq_celltype.bj
            .combine(scrm_bj_peaks.collect().toList())
            .set{endseq_cell_bj_w_scrm_peak_ch}

        endseq_celltype.k562
            .combine(scrm_k562_peaks.collect().toList())
            .set{endseq_cell_k562_w_scrm_peak_ch}

        endseq_celltype.rpe1
            .combine(scrm_rpe1_peaks.collect().toList())
            .set{endseq_cell_rpe1_w_scrm_peak_ch}

        // then concate them like I did with the normal peaks below
        endseq_celltypes_w_scrm_peak_concat = endseq_cell_imr90_w_scrm_peak_ch.concat(endseq_cell_bj_w_scrm_peak_ch, endseq_cell_k562_w_scrm_peak_ch, endseq_cell_rpe1_w_scrm_peak_ch )

        ///////////////////////////////////////////////////////////////////////////////


        // so now I need three processes that will run each of these, or one process that will run a concatenated version of all three of these
        ////// using this or the next one for all celltypes later  ////////
        //gloe_imr90_wrapper_process(gloe_cell_imr90_w_peak_ch)
        ///////////////////////////////////////////////////////////////////
        // I will create a test process that will take all three channels in parallel
        // first concatenate them
        endseq_celltypes_w_peak_concat = endseq_cell_imr90_w_peak_ch.concat(endseq_cell_bj_w_peak_ch, endseq_cell_k562_w_peak_ch, endseq_cell_rpe1_w_peak_ch )
        // first check if it looks good
        endseq_celltypes_w_peak_concat.view(it -> "this is the concat channel with all three celltype tuples for endseq: $it")
        
        // I want to make a process that will generate the scrampled peaks and then output it into a channel like above but with scrambled peaks instead of normal peaks
        // not doing this yet

        
        // now the process as an alias
        
        endseq_celltypes_wrapper_process(endseq_celltypes_w_peak_concat)

        // now run the process to get the break density for scrm gloe peaks
        // just copy the entire process but make it output with a change to file name from norm peaks to scrm peaks
        endseq_celltypes_wrapper_scrm_process(endseq_celltypes_w_scrm_peak_concat)


        // i would need to filter the bam using samtools to get only read1 for gloe seq
        // HAVE TO FIND HOW TO DO THIS FOR END SEQ
        // I DO NOT NEED TO FILTER FOR READ 1 HERE, because it is pair end and is only using one read
        // filt_gloe_bam_samtools_process(bam_index_tuple_ch)
        // r1_filt_gloe_bams = filt_gloe_bam_samtools_process.out.bam_filt_r1_gloe
        // //r1_filt_gloe_bams.view{it -> "these are the gloe bams index tuple that are filt by read1: $it"}

        //get the bigwig files from each of the bam files
        //make_bigwig_gloe_process(bam_index_tuple_ch)
        make_bigwig_endseq_process(bam_index_tuple_ch)

        // getting the chrmt bigwigs
        endseq_bigwig_chrmt_ch = make_bigwig_endseq_process.out.endseq_bigwig_mt

        // I HAVE TO CHANGE TOKENS HERE
        // now to filter for only plc and cells
        endseq_bigwig_chrmt_ch
            .filter(~/.*(?:Cell|PLC).*/)
            .map{ file ->
            
            bigwig_basename = file.baseName
            bigwig_filename = file.name

            tokens = bigwig_basename.tokenize("_")
    
            biorep = tokens[3]  // was 1
            expr_type = tokens[1] // was 2
            cell_type = tokens[0] // was 0

            tuple(cell_type, biorep, expr_type, bigwig_filename, file)

            
            }
            .groupTuple(by:[0,1], sort:true)
            .view{it -> "this is the chrmt bigwigs grouped by biorep: $it"}
            .set{endseq_bigwig_tuple_cell_plc_for_ratio_ch}

        
        // then with this output to find the ratio between cells over plc only filter for cells and plc, then group channel by the first field
        // gb1 should then have the cell and plc for that tuple
        get_endseq_ratio_cell_vs_plc_bigwig_process(endseq_bigwig_tuple_cell_plc_for_ratio_ch)

    



    //emit:




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



