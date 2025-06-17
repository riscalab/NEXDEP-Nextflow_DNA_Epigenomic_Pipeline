

include{
    make_bigwig;
    make_heatmap
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