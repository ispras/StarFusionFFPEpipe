# Run the Docker daemon as a non-root user (Rootless mode)
(from https://docs.docker.com/engine/security/rootless/)

## Requirements

1) python >= 3.*
2) Java 11
3) Minimum 40 GB of RAM 

## Getting started

1) Clone STAR-Fusion docker image
```
docker pull docker.io/trinityctat/starfusion:1.12.0
```

2) Download, extract The CTAT genome library and move it into project directory:
```
git clone
cd starfusionffpe
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
tar -xvf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
mv GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir ctat_genome_lib_build_dir
```

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
