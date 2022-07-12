# tells WDL to use most recent version
version development

import "tasks/CellRegMap.wdl" as C
import "tasks/utils.wdl" as u
import "tasks/post_processing.wdl" as pp


workflow RunCellRegMap {
    input {
        Map[String, File] genotypeFiles # one file per chromsome
        Map[String, File] phenotypeFiles
        Map[String, File] mafFiles
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

    scatter (outputPair in GetGeneChrPairs.output_pairs) {

        call C.RunInteraction as RunInteraction {
            input:
                inputFile=inputFile, # what is this??
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

    call u.GetGeneChrPairs as GetGeneChrPairsBetas {
        input:
            betaFeatureVariantFile=AggregateInteractionResults.out[1]
    }

    scatter (outputPair in GetGeneChrPairsBetas.output_pairs) { 

        call C.EstimateBetas as EstimateBetas {
            input:
                genotypeFile=genotypeFiles[outputPair[0]],
                phenotypeFile=phenotypeFiles[outputPair[0]],
                contextFile,
                kinshipFile,
                betaFeatureVariantFile=AggregateInteractionResults.out[1],
                chr=outputPair[0],
                gene=outputPair[1],
                mafFile=mafFiles[outputPair[0]],

            # Map[String, File] genotypeFiles # one file per chromsome
            # Map[String, File] phenotypeFiles
            # File contextFile
            # File kinshipFile
            # File sampleMappingFile
            # File betaFeatureVariantFile = AggregateIntegrationResults.out[1] #syntax ok?
            # Array[File] betaOutputFiles
            # Map[String, File] mafFiles
        }
    }

    call pp.AggregateBetaResults as AggregateBetaResults{
        input:
            listOfFiles1=EstimateBetas.geneOutput1,
            listOfFiles2=EstimateBetas.geneOutput2

    }

    output {
        File out_interaction_all_results = AggregateIntegrationResults.all_results
        File out_interaction_significant_results = AggregateIntegrationResults.significant_results
        File out_betaG = AggregateBetaResults.all_betaG
        File out_betaGxC = AggregateBetaResults.all_betaGxC
    }

}
