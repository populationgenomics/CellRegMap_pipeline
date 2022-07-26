version development

''' 
Tasks to aggregate results from multiple scatters 
(chromosome-gene pairs) each running CellRegMap functions
'''

# takes as input one file per scatter, returns both 
# all results and all significant results (at a given FDR)

task AggregateInteractionResults { # results from RunInteraction
    input {
        Array[File] listOfFiles 
        Float FDR_threshold # false discovery rate
    }

    command <<<
        # michael: to work out how to more cleanly pass these files in
cat << EOF >> path-results.txt
~{sep("\n", listOfFiles)}
EOF

    eval "$(conda shell.bash hook)" 
    conda activate cellregmap_notebook

    python /share/ScratchGeneral/anncuo/github_repos/CellRegMap_pipeline/images/cellregmap/summarise.py \
        --fdr-threshold ~{FDR_threshold} \
        --file-with-filenames path-results.txt
    >>>
    
    output {
        File all_results = "summary.csv"
        File significant_results = "significant_results.csv"
    }
}

# takes as input one file per scatter, returns both 
# two sets of parameters (betaG and betaGxC) across all scatters

task AggregateBetaResults { # results from EstimateBetas
    input {
        Array[File] listOfFiles1
        Array[File] listOfFiles2
    }

    command <<<

cat << EOF >> path-results-betaG.txt
~{sep("\n", listOfFiles1)}
EOF
cat << EOF >> path-results-betaGxC.txt
~{sep("\n", listOfFiles2)}
EOF

    eval "$(conda shell.bash hook)" 
    conda activate cellregmap_notebook

    python /share/ScratchGeneral/anncuo/github_repos/CellRegMap_pipeline/images/cellregmap/summarise_betas.py \
        --file-with-filenames-1 path-results-betaG.txt \
        --file-with-filenames-2 path-results-betaGxC.txt
    >>>

    output {
        File all_betaG = "summary_betaG.csv"
        File all_betaGxC = "summary_betaGxC.csv"
    }
}
