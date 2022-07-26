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
        Float fdrThreshold=1
    }

    call u.CsvPairExtractor as GetGeneChrPairs { # returns [chrom, gene] pairs
        input:
            csvFile=featureVariantFile,
            columnsToSelect=["chrom", "gene"],
    }

    scatter (outputPair in GetGeneChrPairs.outputPairs) {

        String RunInteractionChr = outputPair[0]

        call C.RunInteraction as RunInteraction {
            input:
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

    call pp.AggregateInteractionResults as AggregateInteractionResults {
        input:
            listOfFiles=RunInteraction.geneOutputPvalues,
            fdrThreshold=fdrThreshold,

    }

    call u.CsvPairExtractor as GetGeneChrPairsBetas {
        input:
            csvFile=AggregateInteractionResults.significantResults,
            columnsToSelect=["chrom", "gene"],
    }

    scatter (outputPair in GetGeneChrPairsBetas.outputPairs) { 

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
                betaFeatureVariantFile=AggregateInteractionResults.significantResults,
                nContexts=nContexts,
                mafFile=mafFiles[EstimateBetaChr],
        }
    }

    call pp.AggregateBetaResults as AggregateBetaResults {
        input:
            listOfFilesBetaG=EstimateBetas.geneOutputBetaG,
            listOfFilesBetaGxC=EstimateBetas.geneOutputBetaGxC

    }

    output {
        File outInteractionAllResults = AggregateInteractionResults.allResults
        File outInteractionSignificantResults = AggregateInteractionResults.significantResults
        File outBetaG = AggregateBetaResults.allResultsBetaG
        File outBetaGxC = AggregateBetaResults.allResultsBetaGxC
    }

}
