

include{
    make_bigwig;
    make_heatmap;
    get_insert_size_metrics_gatk_process
}from '../modules/downstream_analysis_modules.nf'

include{
    merge_bam_lanes_process

}from '../modules/fastq2bam_dna_modules.nf'

workflow down_stream_workflow {



    take:
    bam_ch
    index_ch


    emit:




    main:

    // if i plot each of the groups then combine them after


    //bam_ch.view()

    make_bigwig(bam_ch, index_ch)

    big_wigs = make_bigwig.out.big_wig_for_plotting

    make_heatmap(big_wigs.collect())


}

workflow gatk_analysis_workflow {


    take:
    bam_index_tuple_ch


    main:

    // do some channel manipulation here
    // check how the bam_index_tuple_ch channel looks
    bam_index_tuple_ch.view{it -> "this is the bam index channel for the current data in gatk workflow: $it"}

    // run gatk picard for collecting insert size metrics in a process
    get_insert_size_metrics_gatk_process(bam_index_tuple_ch)




    
}

workflow merge_by_lane_or_techrep_workflow {


    take:
    bam_index_tuple_ch


    main:

    if (params.cadc_grouping_key != false) {


        bam_index_tuple_ch 
            .map{ bam_file, bai_file -> 

            bam_name = bam_file.baseName

            file_tokens = bam_name.tokenize("_")


            
            
            // what i actually need to do is recreate the file basename without the lane number and group based on that key.
            // that will get all the files that have the same name but different lanes
            condition_name = file_tokens[params.condition_type_field_num]
            experiment_name = file_tokens[params.experiment_type_field_num]
            replicate_name = file_tokens[params.replicate_type_field_num]
            lane_name = file_tokens[params.lane_type_field_num]

            // for the cadc data I need to just use the first field as a way to group the data. for now, and call that the grouping field
            
            grouping_key = file_tokens[params.cadc_grouping_key]
            // in the future put the replicate number here so that information is recorded, instead of adding it to the grouping key abovefeil
            tuple(grouping_key, bam_file)

        }
        .groupTuple(by:0)
        // .view {it -> "these are the bams to merge based on the lane number: $it"}
        .set{bams_to_merge_ch}


    } else if (!params.cadc_grouping_key) {

        bam_index_tuple_ch 
            .map{ bam_file, bai_file -> 

            bam_name = bam_file.baseName

            file_tokens = bam_name.tokenize("_")


            
            
            // what i actually need to do is recreate the file basename without the lane number and group based on that key.
            // that will get all the files that have the same name but different lanes
            condition_name = file_tokens[params.condition_type_field_num]
            experiment_name = file_tokens[params.experiment_type_field_num]
            replicate_name = file_tokens[params.replicate_type_field_num]
            lane_name = file_tokens[params.lane_type_field_num]

            // for the cadc data I need to just use the first field as a way to group the data. for now, and call that the grouping field

            grouping_key = "${condition_name}_${experiment_name}_${replicate_name}"
            

            // in the future put the replicate number here so that information is recorded, instead of adding it to the grouping key abovefeil
            tuple(grouping_key, bam_file)

        }
        .groupTuple(by:0)
        // .view {it -> "these are the bams to merge based on the lane number: $it"}
        .set{bams_to_merge_ch}
        
    }

    bams_to_merge_ch.view{it -> "these are the bams to merge based on the lane number: $it"}
    

    // then here make a process to merge bam files, then generate the index files and emit it to a bam index tuple channel and save it into that name also.

    merge_bam_lanes_process(bams_to_merge_ch)

    bam_index_tuple_ch_out = merge_bam_lanes_process.out.merged_bam_index_tuple_no_dup

    // now here i need to put a channel that has the bam before the dedup
    // so the process above will dedup but i will save that to bam_index_tuple_ch_out
    // while the bams generated after merging and indexing without dedup first will go into a name called non_dedup_stats_ch or bam_index_tuple_for_stats_ch

    bam_index_tuple_for_stats_ch = merge_bam_lanes_process.out.merged_bam_index_tuple_dup

    emit:
    bam_index_tuple_ch_out
    bam_index_tuple_for_stats_ch



}