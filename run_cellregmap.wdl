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
    call RunInteraction {
        conda activate my_conda_env
        python run_interaction.py
    }

    output {
        File interaction_results.txt
    }
}