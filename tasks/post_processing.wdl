version development

task AggregateInteractionResults {
    input {
        Array[File] listOfFiles # to figure out (how to get these files)
        Float FDR_threshold
    }

    command <<<
        # michael: to work out how to more cleanly pass these files in
        cat << EOF >> path-results.txt
~{sep("\n", listOfFiles)}
EOF

        python summarise.py \
            --fileWithFilenames path-results.txt
    >>>
    
    output {
        File all_results = "file.txt"
        File significant_results = "significant_results.csv"
    }
}

task AggregateBetasResults {
    input {
        Array[File] listOfFiles
    }

    command <<<
        python summarise_betas.py --geneFiles ~{sep(" ", listOfFiles)}
    >>>

    output {
        File all_betaG = "summary_betaG.csv"
        File all_betaGxC = "summary_betaGxC.csv"
    }
}
