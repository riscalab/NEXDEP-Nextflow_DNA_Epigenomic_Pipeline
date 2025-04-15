
nextflow.enable.dsl=2

// first get the parameters above everything else

params.sorted_break_beds_se = '/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/hera_pipeline/results_SE/break_point_bed/*sort2.breaks.sorted*'
se_sorted_break_beds = Channel.fromPath(params.sorted_break_beds_se)

params.sorted_break_beds_pe = '/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/hera_pipeline/results_PE/break_point_bed/*sort2.breaks.sorted*'
pe_sorted_break_beds = Channel.fromPath(params.sorted_break_beds_pe)


//getting all the spike bedfiles and finding what the largest number is so I can use a numerator
params.spike_sorted_break_beds_se = '/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/hera_pipeline/results_SE/break_point_bed/*{t7,lambda}*sorted.breaks.sorted*'
se_spike_beds = Channel.fromPath(params.spike_sorted_break_beds_se)

params.spike_sorted_break_beds_pe = '/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/hera_pipeline/results_PE/break_point_bed/*{t7,lambda}*sorted.breaks.sorted*'
pe_spike_beds = Channel.fromPath(params.spike_sorted_break_beds_pe)


include {
    get_scaling_numerator


}from './modules/fastq2bam_dna_modules.nf'


workflow {


    //se_spike_beds.view()
    
    // i used stdout in this process to get the value of the spike in file that contains the highest number of reads.
    // this value will be the numerator i use in the formula to normalize by spike in.
    get_scaling_numerator(se_spike_beds.collect())

    max_spike_count = get_scaling_numerator.out

    max_spike_count.view()



    // now I need to get the grouped_cells_plc for all the break files including spike in i think

    if (params.gloe_seq){
        pe_sorted_break_beds
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

    se_sorted_break_beds
        .combine(se_spike_beds)
        .flatten()
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
    }

    
    
    
    grouped_cells_plc.view()
    //se_sorted_break_beds.view()


}