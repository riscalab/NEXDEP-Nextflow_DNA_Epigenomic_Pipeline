

include {
    overlap_window

} from '../modules/fastq2bam_dna_modules.nf'

workflow align_depth_in_peaks {


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

        overlap_window(combined_bed_peak.take(10))


    }else {
        
        
        overlap_window(combined_bed_peak)
    
    
    }

    











}