version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
        String outputName
    }

    command {
        # bash environment
        python get_scatter.py ${featureVariantFile} --outputName "${outputName}"

        # for output_pairs
        echo 'chr1\tGeneName\nchr2\tGeneName2' > outputPairs.tsv
    }

    runtime {
        memory: "2Gb" 
    }

    output {
        # File thisIsMyOutputFile = "GeneChromosomePairs.txt"
        # [["chr1", "GeneName"], ["chr2", "GeneName2"]]
        Array[Array[String]] output_pairs = read_tsv("./outputPairs.tsv")
    }
}

# task CreateBetaFvf{ # this may not actually be needed - just use the "significant only results" output from aggregate interactions

#     input {
#         File AggregateInteractionResults.out
#     }
# }