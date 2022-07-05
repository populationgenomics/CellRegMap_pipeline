version development # tells WDL to use most recent version

import "tasks/CellRegMap.wdl" as C
import "tasks/utils.wdl" as u
import "tasks/post_processing.wdl" as pp


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

    call u.GetGeneChrPairs as GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile
    }

    scatter (outputPair in GetGeneChrPairs.outputPairs) {

        call C.RunInteraction as RunInteraction {
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

    call pp.AggregateInteractionResults as AggregateIntegrationResults{
        input:
        # geneOutput is implicitly a Array[File]
            listOfFiles=RunInteraction.geneOutput

    }

    # call u.CreateBetaFvf as CreateBetaFvf { # not needed hopefully
    #     input 
    #         File interaction_results.txt
    #     ouput 
    #         File beta_fvf
    # }

    call C.EstimateBetas as EstimateBetas{
        input {
            Map[String, File] genotypeFiles # one file per chromsome
            Map[String, File] phenotypeFiles
            File contextFile
            File kinshipFile
            File sampleMappingFile
            File betaFeatureVariantFile = AggregateIntegrationResults.out[1] #syntax ok?
            Array[File] betaOutputFiles
            Map[String, File] mafFiles
        }
    }

    # runtime { # do we need to specify a runtime for the entire workflow?
    #     memory: memory_requirements
    #     cpus: 4
    #     container: "limix/cellregmap:v1.0.0"
    # }

    output {
        File out_interaction = AggregateIntegrationResults.out
        File out_betas = AggregateBetaResults.out
    }

}