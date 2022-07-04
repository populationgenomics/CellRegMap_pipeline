version development

task RunInteraction {
    input {
        String chrom 
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

task EstimateBetas {
    input {
        Int chrom # how do I go from the GeneChromosomePairs file to getting a specific chrom geneName pair for this task?
        Float geneName
        File sampleMappingFile
        File genotypeFile
        File phenotypeFile
        File contextFile
        File kinshipFile
        File betaFeatureVariantFile
        Int nContexts = 10
        File mafFile
    }


    command {
        conda activate cellregmap_notebook
        python estimateBetas.py chrom geneName sampleMappingFile genotypeFile phenotypeFile contextFile kinshipFile betaFeatureVariantFile nContexts mafFile --outputFile ${geneName}
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