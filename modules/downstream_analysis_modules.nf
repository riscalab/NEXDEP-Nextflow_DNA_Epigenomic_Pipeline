


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