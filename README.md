# CellRegMap_pipeline
WDL workflow to facilitate running [CellRegMap](https://github.com/limix/CellRegMap).

A container to run the workflow is available [on Dockerhub](https://hub.docker.com/repository/docker/annasecuomo/cellregmap_pipeline).

## To run on a High Performance Computing (HPC) system
```
java \
    -Dconfig.file=qsub.conf \
    -jar /cromwell/path/cromwell-56.jar run \
    runCellRegMap.wdl \
    --imports tasks.zip \
    --inputs inputs.json
```

## To run on a the Google Cloud Platform (GCP) using analysis-runner
```
analysis-runner cromwell submit  \
    --dataset "tob-wgs" \
    --description "Run CellRegMap WDL pipeline" \
    --access-level "full" \
    --output-dir "scrna-seq/CellRegMap_input_files/"
    -i input_gcs.json \
    -p tasks/ \
    runCellRegMap.wdl
```
