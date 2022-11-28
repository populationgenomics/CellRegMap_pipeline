# CellRegMap pipeline

## Usage overview

In its original version, this is a WDL workflow to facilitate running [CellRegMap](https://github.com/limix/CellRegMap).

[CellRegMap v2](https://github.com/annacuomo/CellRegMap) (new name to come) now comes with new functionalities, particularly around the implementation of tests to map effects of rare genetic variants (called from whole-genome sequencing, WGS).
To assess these new methods, this pipeline needs to be run on both real and simulated data.

For simulations, the pipeline also needs to run alternative methods (including the [SKAT](https://github.com/leelabsg/SKAT)- and [ACAT](https://github.com/yaowuliu/ACAT)-derived methods) for comparison purposes, which are implemented in R.

For real data, the pipeline also needs to use Hail Query to query genetic variants and annotations stored as Hail Matrix Tables and Hail Tables.

### Structure of repo:
* WDL workflow 
  * this does not contain input files generation (does not require Hail Query) and takes pre-generated input files
  * should work for both simulated and real data (see below)
  * contains both python and R scripts as tasks (and respective Docker images)
  * returns p-values and other summary statistics for all results (original CellRegMap, new CellRegMap, R-implemented tests) 
* Hail Batch workflow
  * single-script end-to-end (Python) workflow from data collection / input file generation (using Hail Query) to association testing to returning summary stats (for new CellRegMap specifically, to run internally)
* (Python) scripts to generate both real and simulated input files to feed to the WDL workflow

## Hail Batch workflow
Single [python script](batch.py), key steps (implemented as distinct functions):
* preselection of full WGS hail matrix table to includ only biallelic (n alleles=2) rare (freq<5%) SNVs
* given one gene and the prefiltered WGS hail matrix table, selects relevant variants and export as plink files
* plink files are read in and turned into genotype input files, other input files (phenotype, kinship) are processed
* genes are parallelised into chunks
* CellRegMap-RV various tests are run (each chunk independently)
* resulting p-values are aggregated across jobs and saved

To run:
```
analysis-runner \
    --dataset tob-wgs \
    --access-level test \
    --output-dir "tob_wgs_rv/pseudobulk_rv_association" \
    --description "CellRegMap batch job" \
    python3 batch.py \
      --expression-file-prefix scrna-seq/grch38_association_files \
      --sample-mapping-file scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv \
      -- genes VPREB3 \
      --chromosomes 22 \
      --cell-types B_intermediate
```

## CellRegMap pipeline v1

A WDL workflow to facilitate running [CellRegMap](https://github.com/limix/CellRegMap).
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
<!-- ### TO DO: SLURM (sbatch) / LFS (platform Load Sharing Facility, bsub) -->

## To run on a the Google Cloud Platform (GCP) using analysis-runner
```
analysis-runner cromwell submit  \
    --dataset "tob-wgs" \
    --description "Run CellRegMap WDL pipeline" \
    --access-level "full" \
    --output-dir "scrna-seq/CellRegMap_input_files/2022-08-10_output/" \
    -i input_gcs.json \
    -p tasks/ \
    runCellRegMap.wdl
```
