

include{
    make_bigwig;
    make_heatmap;
    get_insert_size_metrics_gatk_process
}from '../modules/downstream_analysis_modules.nf'



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