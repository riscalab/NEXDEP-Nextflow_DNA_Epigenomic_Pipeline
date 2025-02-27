


// note: this workflow is in the workflow dir so i need to go back one dir to get the modules dir
include {breakDensityWrapper_process} from '../modules/fastq2bam_dna_modules.nf'

workflow breakDensityWrapper_workflow {

    // named workflows should receive their inputs explicitly though the take section
    take:
    bam_index_tuple 
    peak_files


    main:
    //bam_index_tuple.view() // if i collect everything the script should only look for the bam files anyway

    // this collects all of the elements in the channel, flattens them, and then gets only the bam files
    bam_index_tuple.collect()
                .flatten()
                .filter(~/.*bam/)
                .set{only_bams}

    only_bams.view()

    // this workflow will take the files from the main workflow and pass them to the breakDensityWrapper process/module
    if (params.test) {
        breakDensityWrapper_process(only_bams.collect(), peak_files.flatten().take(1)) // just for test, but i just want to take 1 peak file
    }
    else {
        breakDensityWrapper_process(only_bams.collect(), peak_files.flatten())

    }

}



