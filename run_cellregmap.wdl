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
        Int nContexts=10
        Float FDR_threshold=0.05
        # Array[File] outputFiles # what is this? do I need one for betas results too?
    }

    call u.GetGeneChrPairs as GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile,
    }

    scatter (outputPair in GetGeneChrPairs.output_pairs) {

        String RunInteractionChr = outputPair[1]

        call C.RunInteraction as RunInteraction {
            input:
                # inputFile=inputFile, # what is this??
                chrom=RunInteractionChr,
                geneName=outputPair[0],
                sampleMappingFile=sampleMappingFile,
                genotypeFile=genotypeFiles[RunInteractionChr],
                phenotypeFile=phenotypeFiles[RunInteractionChr],
                contextFile=contextFile,
                kinshipFile=kinshipFile,
                featureVariantFile=featureVariantFile,
                nContexts=nContexts,
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

        String EstimateBetaChr = outputPair[1]

        call C.EstimateBetas as EstimateBetas {
            input:
                chrom=EstimateBetaChr,
                geneName=outputPair[0],
                sampleMappingFile=sampleMappingFile,
                genotypeFile=genotypeFiles[EstimateBetaChr],
                phenotypeFile=phenotypeFiles[EstimateBetaChr],
                contextFile=contextFile,
                kinshipFile=kinshipFile,
                betaFeatureVariantFile=AggregateInteractionResults.significant_results,
                nContexts=nContexts,
                mafFile=mafFiles[EstimateBetaChr],

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
