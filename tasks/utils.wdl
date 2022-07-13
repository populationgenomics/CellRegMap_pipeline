version development

task GetGeneChrPairs {

    input {
        File featureVariantFile
        # String outputName # do I need this??
    }

    command <<<
        # bash environment
        # delimiter=comma, select 2nd and 4th columnd, omit first line
        cut -d "," -f 2,4 ~{featureVariantFile | tail -n +2 | sort | uniq > outputPairs.csv
        # python get_scatter.py ${featureVariantFile} --outputName "${outputName}"
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
