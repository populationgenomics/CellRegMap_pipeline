# tells WDL to use most recent version
version development

import "tasks/CellRegMap.wdl" as C
import "tasks/utils.wdl" as u
import "tasks/post_processing.wdl" as pp


workflow RunCellRegMap {
    input {
        Map[String, File] genotypeFiles # one file per chromosome
        Map[String, File] phenotypeFiles
        Map[String, File] mafFiles
        File contextFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
        Int nContexts
        Int FDR_threshold
        # Array[File] outputFiles # what is this? do I need one for betas results too?
    }

    call u.GetGeneChrPairs as GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile,
    }

    scatter (outputPair in GetGeneChrPairs.output_pairs) {

        call C.RunInteraction as RunInteraction {
            input:
                # inputFile=inputFile, # what is this??
                chrom=outputPair[0],
                geneName=outputPair[1],
                sampleMappingFile,
                genotypeFile=genotypeFiles[outputPair[0]],
                phenotypeFile=phenotypeFiles[outputPair[0]],
                contextFile,
                kinshipFile,
                featureVariantFile,
                nContexts,
        }
    }

    call pp.AggregateInteractionResults as AggregateInteractionResults{
        input:
            # geneOutput is implicitly a Array[File]
            listOfFiles=RunInteraction.geneOutput,
            FDR_threshold=FDR_threshold,

    }

    call u.GetGeneChrPairs as GetGeneChrPairsBetas {
        input:
            featureVariantFile=AggregateInteractionResults.significant_results,

    }

    scatter (outputPair in GetGeneChrPairsBetas.output_pairs) { 

        call C.EstimateBetas as EstimateBetas {
            input:
                chrom=outputPair[0],
                geneName=outputPair[1],
                sampleMappingFile,
                genotypeFile=genotypeFiles[outputPair[0]],
                phenotypeFile=phenotypeFiles[outputPair[0]],
                contextFile,
                kinshipFile,
                betaFeatureVariantFile=AggregateInteractionResults.significant_results,
                nContexts,
                mafFile=mafFiles[outputPair[0]],

            # Map[String, File] genotypeFiles # one file per chromosome
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
        File out_interaction_all_results = AggregateInteractionResults.all_results
        File out_interaction_significant_results = AggregateInteractionResults.significant_results
        File out_betaG = AggregateBetaResults.all_betaG
        File out_betaGxC = AggregateBetaResults.all_betaGxC
    }

}
