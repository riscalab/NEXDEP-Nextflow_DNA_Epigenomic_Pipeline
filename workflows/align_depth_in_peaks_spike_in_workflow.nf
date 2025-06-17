include {
    overlap_window_spike_in

} from '../modules/spike_in_modules.nf'

workflow align_depth_in_peaks_spike_in_workflow {


    take:

    combined_bed_peak   //dont have to name it the exact same but it helps for continuity
    // zero_gy
    // plc
    // cells
    // all_peaks

    // grouped_beds

    // all_peaks



    main:

    //zero_gy.view()

    // make a process that takes all of these files conditions and 

    // zero_gy.view{it -> "0gy: $it"}
    // PLC.view{it -> "plc: $it"}
    // cells.view{it -> "cells: $it"}
    //overlap_window(zero_gy, plc, cells, all_peaks)
    
    if (params.test) {

        overlap_window_spike_in(combined_bed_peak.take(10))

        multi_intersect_tsv_ch = overlap_window_spike_in.out.tsv_qc_files


    }else {
        
        
        overlap_window_spike_in(combined_bed_peak)

        multi_intersect_tsv_ch = overlap_window_spike_in.out.tsv_qc_files

        // t7_intersect_tsv_ch = overlap_window_spike_in.out.t7_qc_files
        // lambda_intersect_tsv_ch = overlap_window_spike_in.out.lambda_qc_files
        // yeast_intersect_tsv_ch = overlap_window_spike_in.out.yeast_qc_files
    
    
    }

    // multi_intersect_tsv_ch
    //     .collectFile(name: 'spike_in_alignmentReads_in_peaks_depth.tsv', keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth" )
    //     .set{spike_in_reads_in_peaks_depth_ch}



    multi_intersect_tsv_ch
        .collectFile(name: 'spike_in_alignmentReads_in_peaks_depth.tsv', keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth" )
        .map{file ->
            lines = file.text.readLines()
            header = lines[0]
            data = lines[1..-1]
            .collect {it.split('\t')}
            .sort{row -> row[1]}
            .collect{it.join("\t")}
            return ([header] + data).join("\n")
        
        
        }
        .subscribe { sorted_data ->
        
            file_name = file("${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth/spike_in_alignmentReads_in_peaks_depth_sorted2.tsv")
            file_name.text = sorted_data
        
        
        }
        .set{spike_in_reads_in_peaks_depth_sorted_ch}
    // multi_intersect_tsv_ch
    //     .collectFile(name: 'spike_in_alignmentReads_in_peaks_depth_sorted.tsv', sort: true, keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth" )
    //     .set{spike_in_reads_in_peaks_depth_sorted_ch}

    // this is me using collect files operator to get the 1071 files and save them to the new file specified in name
    // newLine ensures that every appended file's contents gets saved in a new line... hopefully. probably dont want an empty line every other line
    // keepHeader takes all the lines from the first collected file and including the first line as a header, but the other files will not have the first line appended. 
    // can use skip to choose the number of lines you want the resulting collected files to skip, but default is 0 so hopefully it skips the header still without using skip, if not then i have to use skip: 1
    // might have to use the storeDir option to put the resulting file in a dir of my choice
    // t7_intersect_tsv_ch
    //     .collectFile (name: 't7_alignmentReads_in_peaks_depth.tsv', keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth")
    //     .set {combined_t7_depth_intersect_ch}  

    // //combined_t7_depth_intersect_ch.view()

    // lambda_intersect_tsv_ch
    //     .collectFile (name: 'lambda_alignmentReads_in_peaks_depth.tsv', keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth")
    //     .set {combined_lambda_depth_intersect_ch} 

    // yeast_intersect_tsv_ch
    //     .collectFile (name: 'yeast_alignmentReads_in_peaks_depth.tsv', keepHeader: true, storeDir: "${params.base_out_dir}/alignment_peak_overlap_qc/spike_in/complete_intersection_depth")
    //     .set {combined_yeast_depth_intersect_ch} 

    // then I would have to make a process to do the calculation of percentages
    // awk 'NR>1 {print ($3/($7+1))*100}' alingmentReads_in_peaks_depth.tsv > percent_test.txt use someting like this but dont make a new file










}