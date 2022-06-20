version development

workflow RunCellRegMap {
    input {
        File genotypeFile
        File phenotypeFile
        File contextFile
        File covariateFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
    }

    scatter (chrom_gene in ChromGenePairs) {
        call RunInteraction {
            input: 
                Int chrom
                Int i
            conda activate my_conda_env
            python run_interaction.py chrom i
        }
    }

    call AggregateIntegrationResults {
        input:
            listOfFiles=RunInteraction.out

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
        File out = stdout()
    }

}