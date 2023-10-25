# StarFusionFFPEpipe

## Requirements

1) Docker 
2) python >= 3.*
3) Java 11
4) Minimum 40 GB of RAM 


## Getting started

Installing Docker:
https://docs.docker.com/engine/install/
  
To run containers in rootless mode:
https://docs.docker.com/engine/install/
(you may also need to increase number of subordinate UIDs/GIDs up to 1 million)


1) Clone STAR-Fusion docker image
```
docker pull docker.io/trinityctat/starfusion:1.12.0
```

2) Download, extract The CTAT genome library and move it into project directory:
```
git clone
cd starfusionffpe
wget https://nextcloud.ispras.ru/index.php/s/rLcdrHTZ7CB3LFY/download/GRCh38_gencode_v43_CTAT_lib_Jul012023.plug-n-play.tar.gz
tar -xvf GRCh38_gencode_v43_CTAT_lib_Jul012023.plug-n-play.tar.gz
```
For short reads (50bp, single-end) use [this](https://nextcloud.ispras.ru/index.php/s/5tix8ABytZcg5zW/download/GRCh38_gencode_v43_CTAT_lib_Jul012023.50bp.plug-n-play.tar.gz) genome library.

3) Download Cromwell jar file and move it into project directory:
```
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
```

## Running pipelines

The simpliest way to start pipeline is using run_pipe.py:
```
python run_pipe.py -fq input_fastq_file
```

The another way is to add path for the input fastq file and The CTAT genome library in input_config_template_[PE|SE].json and run:
```
java -Xmx2g -Dconfig.file=docker.conf -jar cromwell-85.jar run StarFusionFFPE.wdl --options options.json --inputs input_config_template_[PE|SE].json
```

## Results

When pipeline successfully finished result file stored in results directory.
Intermediate files stored in cromwell-executions directory.


## Possible issues

Issues may occur if some files are missing or placed in wrong directory (not in project directory).\
You may move your files into project directory or specify paths for them.
For more info use:
```
python run_pipe.py --help
```
