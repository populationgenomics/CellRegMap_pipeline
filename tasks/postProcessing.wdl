version development

# Tasks to aggregate results from multiple scatters
# (chromosome-gene pairs) each running CellRegMap functions

# takes as input one file per scatter, returns both
# all results and all significant results (at a given FDR)

task AggregateInteractionResults { # results from RunInteraction
    input {
        Array[File] listOfFiles
        Float fdrThreshold # false discovery rate
    }

    command <<<
    # pass files to the task, there can be a lot of files which would break the command line
cat << EOF >> path-results.txt
~{sep("\n", listOfFiles)}
EOF

    # eval "$(conda shell.bash hook)"
    # conda activate cellregmap_notebook

    python /app/summarise.py \
        --fdr-threshold ~{fdrThreshold} \
        --file-with-filenames path-results.txt
    >>>

    output {
        File allResults = "summary.csv"
        File significantResults = "significant_results.csv"
    }

    runtime {
        docker: "annasecuomo/cellregmap_pipeline:dev"
        memory: "8GB"
    }
}

# takes as input one file per scatter, returns both
# two sets of parameters (betaG and betaGxC) across all scatters

task AggregateBetaResults { # results from EstimateBetas
    input {
        Array[File] listOfFilesBetaG
        Array[File] listOfFilesBetaGxC
    }

    command <<<

cat << EOF >> path-results-betaG.txt
~{sep("\n", listOfFilesBetaG)}
EOF
cat << EOF >> path-results-betaGxC.txt
~{sep("\n", listOfFilesBetaGxC)}
EOF

    # eval "$(conda shell.bash hook)"
    # conda activate cellregmap_notebook

    python /app/summarise_betas.py \
        --file-with-filenames-1 path-results-betaG.txt \
        --file-with-filenames-2 path-results-betaGxC.txt
    >>>

    output {
        File allResultsBetaG = "summary_betaG.csv"
        File allResultsBetaGxC = "summary_betaGxC.csv"
    }

    runtime {
        docker: "annasecuomo/cellregmap_pipeline:dev"
        memory: "8GB"
    }
}
