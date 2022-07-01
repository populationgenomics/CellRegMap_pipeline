version development # tells WDL to use most recent version

task GetGeneChrPairs {

    input {
        File featureVariantFile
        String outputName
    }

    command {
        # bash environment
        python get_scatter.py ${featureVariantFile} --outputName "${outputName}"

        # for output_pairs
        echo 'chr1\tGeneName\nchr2\tGeneName2' > outputPairs.tsv
    }

    runtime {
        memory: "2Gb" # I have no sense how much memory would be needed for this - I think very little? Just needs to open a not-very-long txt file and make a new one
    }

    output {
        # File thisIsMyOutputFile = "GeneChromosomePairs.txt"
        # [["chr1", "GeneName"], ["chr2", "GeneName2"]]
        Array[Array[String]] output_pairs = read_tsv("./outputPairs.tsv")
    }
}

task RunInteraction {
    input {
        String chrom # how do I go from the GeneChromosomePairs file to getting a specific chrom geneName pair for this task?
        String geneName
        File sampleMappingFile
        File genotypeFile
        File phenotypeFile
        File contextFile
        File kinshipFile
        File featureVariantFile # still need this to select specific SNPs
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
        File geneOutput1 = geneName + "_betaG.csv"
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
workflow RunCellRegMap {
    input {
        Map[String, File] genotypeFiles # one file per chromsome
        Map[String, File] phenotypeFiles
        File contextFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
        Array[File] outputFiles # what is this? do I need one for betas results too?
    }

    call GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile
    }

    scatter (outputPair in GetGeneChrPairs.outputPairs) {

        call RunInteraction {
            input:
                inputFile=inputFile,
                chr=outputPair[0],
                gene=outputPair[1],
                featureVariantFile,
                sampleMappingFile,
                phenotypeFile=phenotypeFiles[outputPair[0]],
                genotypeFile=genotypeFiles[outputPair[0]],
                kinshipFile,
                contextFile,

        }
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