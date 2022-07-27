version development

# Tasks to intersect all input files and perform CellRegMap 
# functions for a specific chromosome-gene pair


task RunInteraction { # CellRegMap's run_interaction()
    input {
        String chrom 
        String geneName
        File sampleMappingFile
        File genotypeFile
        File genotypeFilesBim
        File genotypeFilesFam
        File phenotypeFile
        File contextFile
        File kinshipFile
        File featureVariantFile 
        Int nContexts = 10 
    }

    command <<<
        # make sure secondaries are paired together
        cp -f '~{genotypeFilesBim}' "$(echo '~{genotypeFile}' | sed 's/\.[^.]*$//').bim"
        cp -f '~{genotypeFilesFam}' "$(echo '~{genotypeFile}' | sed 's/\.[^.]*$//').fam"

        # for now, use conda, but when we're closer,
        # remove this in favor of 'container' in the runtime section
        # eval "$(conda shell.bash hook)" 
        # conda activate cellregmap_notebook 

        python /app/run_interaction.py \
            --chrom ~{chrom} \
            --gene-name ~{geneName} \
            --sample-mapping-file ~{sampleMappingFile} \
            --genotype-file ~{genotypeFile} \
            --phenotype-file ~{phenotypeFile} \
            --context-file ~{contextFile} \
            --kinship-file ~{kinshipFile} \
            --feature-variant-file ~{featureVariantFile} \
            --n-contexts ~{nContexts} 
    >>>

    output {
        File geneOutputPvalues = geneName + ".csv"
    }

    runtime {
        docker: "annasecuomo/cellregmap_pipeline:dev"

        memory: "400Gb" # static
    }
}

task EstimateBetas { # CellRegMap's estimate_betas()
    input {
        Int chrom 
        String geneName
        File sampleMappingFile
        File genotypeFile
        File genotypeFilesBim
        File genotypeFilesFam
        File phenotypeFile
        File contextFile
        File kinshipFile
        File betaFeatureVariantFile
        Int nContexts = 10
        File? mafFile # if MAFs are known / pre-calculated they should be given as inputs
    }

    command <<<
        # make sure secondaries are paired together
        cp -f '~{genotypeFilesBim}' "$(echo '~{genotypeFile}' | sed 's/\.[^.]*$//').bim"
        cp -f '~{genotypeFilesFam}' "$(echo '~{genotypeFile}' | sed 's/\.[^.]*$//').fam"
        
        # leave this until containers are used
        # eval "$(conda shell.bash hook)" 
        # conda activate cellregmap_notebook

        python /app/estimate_betas.py \
            --chrom ~{chrom} \
            --gene-name ~{geneName} \
            --sample-mapping-file ~{sampleMappingFile} \
            --genotype-file ~{genotypeFile} \
            --phenotype-file ~{phenotypeFile} \
            --context-file ~{contextFile} \
            --kinship-file ~{kinshipFile} \
            --beta-feature-variant-file ~{betaFeatureVariantFile} \
            --n-contexts ~{nContexts} \
            ~{if defined(mafFile) then "--maf-file ~{mafFile}" else ""}
    >>>

    output {
        File geneOutputBetaG = geneName + "_betaG.csv"
        File geneOutputBetaGxC = geneName + "_betaGxC.csv"
    }

    runtime {
        docker: "annasecuomo/cellregmap_pipeline:dev"

        # static
        memory: "400Gb"
    }
}
