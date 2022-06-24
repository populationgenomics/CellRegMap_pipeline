version development

task GetScatter {

    input {
        File featureVariantFile
        String outputName
    }

    command {
        # bash environment
        python get_scatter.py ${featureVariantFile} --outputName "${outputName}"
    }

    runtime {
        memory: ""
    }

    output {
        File thisIsMyOutputFile = "this is the output name.txt"
    }
}

task RunInteraction {
    input {
        Int chrom
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
        conda activate my_conda_env
        python run_interaction.py chrom geneName sampleMappingFile genotypeFile phenotypeFile contextFile kinshipFile featureVariantFile nContexts --outputFile ${geneName + ".tsv"}
    }

    output {
        File geneOutput = geneName + ".tsv"
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

task AggregateInteractionResults {
    input {
        Array[File] listOfFiles
        Float FDR_threshold
    }

    command {
        python summarise.py --geneFiles {join(" ", listOfFiles)}
    }

    output {
        File all_results
        File significant_results_only
    }
}


task EstimateBetas {
    input {}

    command {}

    output {}

    runtime {

    }
}

task AggregateBetasResults {
    input {
        Array[File] listOfFiles
    }

    command {
        python summarise_betas.py --geneFiles {join(" ", listOfFiles)}
    }
}

task GetGeneChrPairs {
    input {
        File featureVariantsFile
    }

    command {
        # does a thing
        # write a TSV
        echo 'chr1  BRCA1' > output.tsv
    }

    output {
        Array[Array[String, String]] outputPairs = read_tsv("output.tsv")
    }

}


# { WorkflowName.inputName: "value" }
# { RunCellRegMap.contextFile: "" }
workflow RunCellRegMap {
    input {
        File [genotypeFile]
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

    if not is_defined(genes) {
        call GetGeneList...
    }

    call UseGeneList {
        geneList=select_first([genes, GetGeneList.genes])
    }

    call GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile
    }

    call FeatureDoes {
        input:
            file=featureVariantFile
    }


    scatter ((chr, gene) in GetGeneChrPairs.outputPairs) {

        call RunInteraction {
            input:
                # smallFile is a couple of genes
                inputFile=inputFile,
                chr=chr,
                gene=gene,
                featureVariantFile

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