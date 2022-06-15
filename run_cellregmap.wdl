version development

workflow RunCellRegMap {
    input {
        File genotypeFile
        File phenotypeFile
        File contextFile
        File covariateFile
        File kinshipFile
        File sampleMappingFile
        File featureVariantFile
    }
    command {
        activate cellregmap_notebook
        python ~{pythonScript}
    }
    output {
        File results.txt
    }
}