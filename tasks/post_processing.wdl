version development

task AggregateInteractionResults {
    input {
        Array[File] listOfFiles # to figure out (how to get these files)
        Float FDR_threshold
    }

    command {
        python summarise.py pathResults --geneFiles {join(" ", listOfFiles)} 

        # inputs:
        #   ../inputs/1231239/geneName.txt
        #   ../inputs/5126313/geneName2.txt

        # summarise.py --geneFiles file1 file2 file3
    }
    
    output {
        File all_results
        File significant_results_only
    }
}

task AggregateBetasResults {
    input {
        Array[File] listOfFiles
    }

    command {
        python summarise_betas.py --geneFiles {join(" ", listOfFiles)}
    }

    output {
        File all_betaG_s
        File all_betaGxC_s
    }
}