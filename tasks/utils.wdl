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
        memory: "2Gb" # I have no sense how much memory would be needed for this - I think very little? Just needs to open a not-very-long txt file and make a new one
    }

    output {
        # File thisIsMyOutputFile = "GeneChromosomePairs.txt"
        # [["chr1", "GeneName"], ["chr2", "GeneName2"]]
        Array[Array[String]] output_pairs = read_tsv("./outputPairs.tsv")
    }
}

task CreateBetaFvf{
    
}