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
        Int nContexts = 10 # add additional flag for cases when the contexts tested and those in the background aren't the same

    }

    String outputFilename = geneName + ".csv"

    command <<<
        # for now, use conda, but when we're closer,
        # remove this in favor of 'container' in the runtime section
        conda activate cellregmap_notebook 

        python /Shared/annacuomo/scripts/run_interaction.py \
            --chrom ~{chrom} \
            --gene-name ~{geneName} \
            --sample-mapping-file ~{sampleMappingFile} \
            --genotype-file ~{genotypeFile} \ 
            --phenotype-file ~{phenotypeFile} \
            --context-file ~{contextFile} \
            --kinship-file ~{kinshipFile} \
            --feature-variant-file ~{featureVariantFile} \
            --n-context ~{nContexts} \
            --outputFile ~{outputFilename}
    >>>

    output {
        File geneOutput = outputFilename
    }

    runtime {
        # container: "annacuomo/limix:dev"

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
        Int chrom 
        String geneName
        File sampleMappingFile
        File genotypeFile
        File phenotypeFile
        File contextFile
        File kinshipFile
        File betaFeatureVariantFile
        Int nContexts = 10
        File mafFile # if MAFs are known / pre-calculated they should be given as inputs
    }


    command <<<
        # leave this until containers are used
        conda activate cellregmap_notebook

        python /Shared/annacuomo/scripts/estimate_betas.py \
            --chrom ~{chrom} \
            --gene-name ~{geneName} \
            --sample-mapping-file ~{sampleMappingFile} \
            --genotype-file ~{genotypeFile} \ 
            --phenotype-file ~{phenotypeFile} \
            --context-file ~{contextFile} \
            --kinship-file ~{kinshipFile} \
            --beta-feature-variant-file ~{betaFeatureVariantFile} \
            --n-context ~{nContexts} \
            --maf-file ~{mafFile} \
            --outputFile ${geneName}
    >>>

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
