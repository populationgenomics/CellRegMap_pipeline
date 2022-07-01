version development # what is this again?

task GetGeneChrPairs {

    input {
        File featureVariantFile
        String outputName
    }

    command {
        # bash environment
        python get_scatter.py ${featureVariantFile} --outputName "${outputName}"
    }

    runtime {
        memory: "20Gb" # I have no sense how much memory would be needed for this - I think very little? Just needs to open a not-very-long txt file and make a new one
    }

    output {
        File thisIsMyOutputFile = "GeneChromosomePairs.txt"
    }
}

task RunInteraction {
    input {
        Int chrom # how do I go from the GeneChromosomePairs file to getting a specific chrom geneName pair for this task?
        Float geneName
        File sampleMappingFile
        File genotypeFile
        File phenotypeFile
        File contextFile
        File kinshipFile
        File featureVariantFile
        Int nContexts = 10
    }


    command {
        conda activate cellregmap_notebook
        python run_interaction.py chrom geneName sampleMappingFile genotypeFile phenotypeFile contextFile kinshipFile featureVariantFile nContexts --outputFile ${geneName + ".csv"}
    }

    output {
        File geneOutput = geneName + ".csv"
    }

    runtime {
        # static
        memory: "400Gb"
        
        # # calculated
        ## Unable to allocate 9.78 GiB for an array with shape (133716, 9820)
        ## shape is number of cells X (numbers of donors X number of contexts)
        ## here, 133,716 cells across 982 individuals, 10 contexts
        # memory: size(inputFile) + size(interval)
        
        # # passed in
        # memory: memory # from input
    }
}

task AggregateInteractionResults {
    input {
        Array[File] listOfFiles
        Float FDR_threshold
    }

    command {
        python summarise.py pathResults --geneFiles {join(" ", listOfFiles)} 
    }
    
    output {
        File all_results
        File significant_results_only
    }
}

# here I need GetScatter again but with a new fvf? significant_results_only from above, can i reuse task from above?
task EstimateBetas {
    input {
        Int chrom # how do I go from the GeneChromosomePairs file to getting a specific chrom geneName pair for this task?
        Float geneName
        File sampleMappingFile
        File genotypeFile
        File phenotypeFile
        File contextFile
        File kinshipFile
        File significant_results_only
        Int nContexts = 10
    }


    command {
        conda activate cellregmap_notebook
        python estimateBetas.py chrom geneName sampleMappingFile genotypeFile phenotypeFile contextFile kinshipFile featureVariantFile nContexts --outputFile ${geneName + ".csv"}
    }

    output {
        File geneOutput2 = geneName + "_betaG.csv"
        File geneOutput2 = geneName + "_betaGxC.csv"
    }

    runtime {
        # static
        memory: "400Gb"
        
        # # calculated
        # memory: size(inputFile) + size(interval)
        
        # # passed in
        # memory: memory # from input
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



# { WorkflowName.inputName: "value" }
# { RunCellRegMap.contextFile: "" }
workflow RunCellRegMap { # I am confused, does this just effectively call the tasks in the right order?
    input {
        File [genotypeFile] # again need to sort out how to deal with this (one file per chromsome)
        File [phenotypeFile]
        File contextFile
        File covariateFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
        Array[File] outputFiles


        # how to determine which genes to run
        Array[String] genes


        Array[File] intervals
    }

    # if not is_defined(genes) {
    #     call GetGeneList...
    # }

    # call UseGeneList {
    #     geneList=select_first([genes, GetGeneList.genes])
    # }

    call GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile
    }

    # call FeatureDoes {
    #     input:
    #         file=featureVariantFile
    # }


    scatter ((chr, gene) in GetGeneChrPairs.outputPairs) {

        call RunInteraction {
            input:
                # smallFile is a couple of genes
                inputFile=inputFile,
                chr=chr,
                gene=gene,
                featureVariantFile,
                sampleMappingFile,
                phenotypeFile,
                genotypeFile,

        }
        }

    }

    scatter (chrom_gene in ChromGenePairs) {
        call RunInteraction {
            input: 
                Int chrom=GetScatter.out[0]
                Int i=GetScatter.out[1]
            conda activate my_conda_env
            python run_interaction.py chrom i
        }
    }

    call AggregateInteractionResults {
        input:
        # geneOutput is implicitly a Array[File]
            listOfFiles=RunInteraction.geneOutput

    }

    call CreateBetaFvf {
        input 
            File interaction_results.txt
        ouput 
            File beta_fvf
    }
    call EstimateBetas {
        input {
            File beta_fvf
        }
    }

    runtime {
        memory: memory_requirements
        cpus: 4
        container: "limix/cellregmap:v1.0.0"
    }

    output {
        File out_interaction = AggregateIntegrationResults.out
        File out_betas = AggregateBetaResults.out
    }

}