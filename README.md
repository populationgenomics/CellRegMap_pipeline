# CellRegMap_pipeline
WDL workflow to facilitate running [CellRegMap](https://github.com/limix/CellRegMap).

A container to run the workflow is available [on Dockerhub](https://hub.docker.com/repository/docker/annasecuomo/cellregmap_pipeline).

## To run on a High Performance Computing (HPC) system

### PBS (Portable Batch System, qsub)
```
java \
    -Dconfig.file=qsub.conf \
    -jar /cromwell/path/cromwell-56.jar run \
    runCellRegMap.wdl \
    --imports tasks.zip \
    --inputs inputs.json
```
### TO DO: SLURM (sbatch) / LFS (platform Load Sharing Facility, bsub)

## To run on a the Google Cloud Platform (GCP) using analysis-runner
```
analysis-runner cromwell submit  \
    --dataset "tob-wgs" \
    --description "Run CellRegMap WDL pipeline" \
    --access-level "full" \
    --output-dir "scrna-seq/CellRegMap_input_files/2022-08-10_output/"
    -i input_gcs.json \
    -p tasks/ \
    runCellRegMap.wdl
```
