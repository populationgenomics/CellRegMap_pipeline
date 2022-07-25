# tells WDL to use most recent version
version development

import "tasks/CellRegMap.wdl" as C
import "tasks/utils.wdl" as u
import "tasks/post_processing.wdl" as pp


workflow RunCellRegMap {
    input {
        Map[String, File] genotypeFiles # one file per chromosome
        Map[String, File] genotypeFilesBims
        Map[String, File] genotypeFilesFams
        Map[String, File] phenotypeFiles
        Map[String, File] mafFiles
        File contextFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
        Int nContexts=10
        Float FDR_threshold=1
        # Array[File] outputFiles # what is this? do I need one for betas results too?
    }

    call u.GetGeneChrPairs as GetGeneChrPairs {
        input:
            featureVariantFile=featureVariantFile,
            columnsToSelect=["chrom", "gene"],
    }

    scatter (outputPair in GetGeneChrPairs.output_pairs) {

        String RunInteractionChr = outputPair[0]

        call C.RunInteraction as RunInteraction {
            input:
                # inputFile=inputFile, # what is this??
                chrom=RunInteractionChr,
                geneName=outputPair[1],
                sampleMappingFile=sampleMappingFile,
                genotypeFile=genotypeFiles[RunInteractionChr],
                genotypeFilesBim=genotypeFilesBims[RunInteractionChr],
                genotypeFilesFam=genotypeFilesFams[RunInteractionChr],
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
            columnsToSelect=["chrom", "gene"],
    }

    scatter (outputPair in GetGeneChrPairsBetas.output_pairs) { 

        String EstimateBetaChr = outputPair[0]

        call C.EstimateBetas as EstimateBetas {
            input:
                chrom=EstimateBetaChr,
                geneName=outputPair[1],
                sampleMappingFile=sampleMappingFile,
                genotypeFile=genotypeFiles[EstimateBetaChr],
                genotypeFilesBim=genotypeFilesBims[EstimateBetaChr],
                genotypeFilesFam=genotypeFilesFams[EstimateBetaChr],
                phenotypeFile=phenotypeFiles[EstimateBetaChr],
                contextFile=contextFile,
                kinshipFile=kinshipFile,
                betaFeatureVariantFile=AggregateInteractionResults.significant_results,
                nContexts=nContexts,
                mafFile=mafFiles[EstimateBetaChr],
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
