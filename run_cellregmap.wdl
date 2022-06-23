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
        File interval
        File inputFile

        String memory
        String geneGame
    }


    command {
        conda activate my_conda_env
        python run_interaction.py chrom i --outputFile ${geneName + ".txt"}
    }

    output {
        File geneOutput = geneName + ".txt"
    }

    runtime {
        # static
        memory: "4Gb"
        # calculated
        memory: size(inputFile) + size(interval)
        # passed in
        memory: memory # from input
    }
}

task AggregateIntegrationResults {
    input {
        Array[File] listOfFiles
    }

    command {
        python summarise.py --geneFiles {join(" ", listOfFiles)}
    }
}

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


        Array[File] intervals
    }

    call GetScatter {
        input:
            featureVariantFile=featureVariantFile
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

    call AggregateIntegrationResults {
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