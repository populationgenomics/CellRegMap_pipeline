task R_Task {

    input {
        File myInput
    }

    command <<<
        Rscript <myscriptInsideTheContainer>.R_Task \
            --input ~{myInput} \
            --outputPath "file.txt"

    >>>

    runtime {
        docker: "annacuomo/cell-reg-map:"
    }

    output {
        File out = "file.txt"
    }
}