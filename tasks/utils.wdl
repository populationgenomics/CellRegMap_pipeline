version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
        String outputName
    }

    command <<<
        # bash environment
        # delimiter=comma, select 2nd and 4th columnd, omit first line
        cut -d "," -f 2,4 ${featureVariantFile | tail -n +2 | sort | unique > outputPairs.csv
        # python get_scatter.py ${featureVariantFile} --outputName "${outputName}"

        # for output_pairs
        # echo 'chr1\tGeneName\nchr2\tGeneName2' > outputPairs.tsv
    >>>

    runtime {
        memory: "2Gb" 
    }

    output {
        # File thisIsMyOutputFile = "GeneChromosomePairs.txt"
        # [["chr1", "GeneName"], ["chr2", "GeneName2"]]
        Array[Array[String]] output_pairs = read_tsv("./outputPairs.csv")
    }
}

# task CreateBetaFvf{ # this may not actually be needed - just use the "significant only results" output from aggregate interactions

#     input {
#         File AggregateInteractionResults.out
#     }
# }
