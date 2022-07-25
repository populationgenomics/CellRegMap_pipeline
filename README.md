# CellRegMap_pipeline
Pipeline / workflow to facilitate running [CellRegMap](https://github.com/limix/CellRegMap).


```git pull && zip -r tasks.zip tasks && java -Dconfig.file=qsub.conf -jar /share/ScratchGeneral/anncuo/cromwell/cromwell-56.jar run run_cellregmap.wdl --imports tasks.zip --inputs inputs.json```
