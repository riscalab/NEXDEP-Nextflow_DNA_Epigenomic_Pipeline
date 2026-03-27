


process make_bigwig {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    label 'normal_small_resources'

    input:
    
    path(bam_file)
    path(index_files)



    output:

    path("${output_bw}" ), emit: big_wig_for_plotting

    script:

    output_bw = "${bam_file.baseName}_bw"

    """
    #!/usr/bin/env bash

    bamCoverage \
    -b ${bam_file} \
    -o "${output_bw}" \
    --normalizeUsing RPKM

    # put this in another process
    #plotHeatmap \
    #-m "${output_bw}" \

    """


}

process make_heatmap {
    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    publishDir "./heatmap_from_bigwigs", mode:'copy', pattern: '*'

    input:

    path(all_bigwigs)


    output:

    path("heatmap.png"), emit: bigwig_heatmaps


    script:


    """
    #!/usr/bin/env bash

    multiBigwigSummary bins \
    -b ${all_bigwigs} \
    --binSize 10000 \
    -out results.npz \
    --outRawCounts results.tab

    plotHeatmap \
    -m results.npz \
    -out heatmap.png



    """



}


process get_insert_size_metrics_gatk_process {

    label 'normal_big_resources'
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/gatk_4.6.2.0_rj'
    
    // cant have two environments, just need to add R in the gatk environment
    //conda '/ru-auth/local/home/rjohnson/miniconda3/envs/R_lan_2_rj'

    publishDir "${params.base_out_dir}/GATK_results", mode: 'copy', pattern: '*'

    input:
    // this will be a bam index tuple
    tuple path(bam), path(index)


    output:

    path("${output_txt_filename}"), emit: insert_size_metrics
    path("${output_histogram_filename}"), emit: insert_size_histogram



    script:

    name_split = "${bam.baseName}".split("_")
    // i should get the second and third field also for the files that are named properly
    main_name = name_split[0]

    output_txt_filename = "${main_name}_insert_size_metrics.txt"

    output_histogram_filename = "${main_name}_insert_size_histogram.pdf"


    """
    #!/usr/bin/env bash

    gatk CollectInsertSizeMetrics \
    --INPUT ${bam} \
    --OUTPUT ${output_txt_filename} \
    --Histogram_FILE ${output_histogram_filename} \
    --MINIMUM_PCT 0.05



    """
}