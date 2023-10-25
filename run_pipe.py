from argparse import ArgumentParser
import os
import json
import subprocess


def get_arguments():
    parser = ArgumentParser()
    parser.add_argument('-fq', '--fastq',
                        help='fastq file for single-end and left fastq file for paired-end',
                        dest="fq",
                        metavar="PATH",
                        required=True)
    parser.add_argument('-fq2', '--fastq2',
                        help='right fastq file for paired-end',
                        dest="fq2",
                        metavar="PATH",
                        default=None,
                        required=False)
    parser.add_argument('-cl', '--ctat-lib',
                        help='Directory of CTAT genome library. (default ctat_genome_lib_build_dir)',
                        dest="ctat",
                        metavar="PATH",
                        default="ctat_genome_lib_build_dir",
                        required=False)
    parser.add_argument('-k', '--kinases',
                        help='File with kinases genes. (default kinase_domain_transcripts.tsv)',
                        dest="kinases",
                        metavar="PATH",
                        default="kinase_domain_transcripts.tsv",
                        required=False)
    parser.add_argument('-e', '--exons',
                        help='File with exon numbers. (default exons.json)',
                        dest="exons",
                        metavar="PATH",
                        default="exons.json",
                        required=False)
    parser.add_argument('-l', '--lncrna',
                        help="File with lncRNA transcript names. If no file was given pipeline wouldn't search lncRNA.",
                        dest="lncrna",
                        metavar="PATH",
                        default=None,
                        required=False)
    parser.add_argument('-t', '--threads',
                        help='Number of threads. (16)',
                        dest="threads",
                        type=int,
                        default=16,
                        required=False)
    parser.add_argument('-cm', '--cromwell-memory',
                        help='Maximum gigabytes used by Cromwell. (2)',
                        dest="cromwell_mem",
                        type=int,
                        default=2,
                        required=False)
    parser.add_argument('-sm', '--star-memory',
                        help='Maximum megabytes used by STAR-Fusion. (46080)',
                        dest="star_mem",
                        type=int,
                        default=46080,
                        required=False)
    parser.add_argument('-O', '--out-name',
                        help='Name of result config file (default input_config.json)',
                        dest="out",
                        metavar="String",
                        default="input_config.json",
                        required=False)
    parser.add_argument('-cr', '--cromwell',
                        help='Cromwell jar file. (default cromwell-85.jar)',
                        dest="cromwell",
                        metavar="PATH",
                        default="cromwell-85.jar",
                        required=False)
    parser.add_argument('-wdl', '--wdl',
                        help='WDL file of StarFusionFFPE pipe. (default StarFusionFFPE.wdl)',
                        dest="wdl",
                        metavar="PATH",
                        default="StarFusionFFPE.wdl",
                        required=False)
    parser.add_argument('-conf', '--config',
                        help='Cromwell config file. (default docker.conf)',
                        dest="config",
                        metavar="PATH",
                        default="docker.conf",
                        required=False)
    parser.add_argument('-options', '--options',
                        help='Cromwell options config file. (default options.json)',
                        dest="options",
                        metavar="PATH",
                        default="options.json",
                        required=False)
    parser.add_argument('--no_kinase_info',
                        help='Use this flag to not analyze kinase exons and domains',
                        action='store_true')
    return parser.parse_args()


def check_restrictions(_args):
    if not os.path.exists(_args.fq):
        print(f"ERROR! fastq file {_args.fq} doesn't exists.\nExiting...")
        raise SystemExit
    if _args.fq2 is not None and not os.path.exists(_args.fq2):
        print(f"ERROR! fastq file {_args.fq2} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.ctat):
        print(f"ERROR! Compressed (.zip) file of CTAT genome {_args.ctat} library doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.kinases):
        print(f"ERROR! File with kinases genes {_args.kinases} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.exons):
        print(f"ERROR! File with exon numbers {_args.exons} doesn't exists.\nExiting...")
        raise SystemExit
    if _args.lncrna is not None and not os.path.exists(_args.lncrna):
        print(f"ERROR! File with lncRNA transcript names {_args.lncrna} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.cromwell):
        print(f"ERROR! Cromwell jar file {_args.cromwell} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.wdl):
        print(f"ERROR! WDL file of StarFusionFFPE pipe {_args.wdl} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.config):
        print(f"ERROR! Cromwell config file {_args.config} doesn't exists.\nExiting...")
        raise SystemExit
    if not os.path.exists(_args.options):
        print(f"ERROR! Cromwell options config file {_args.options} doesn't exists.\nExiting...")
        raise SystemExit


if __name__ == '__main__':
    args = get_arguments()
    check_restrictions(args)

    input_json = {"StarFusionFFPE.input_fastq": args.fq}
    if args.fq2 is not None:
        input_json["StarFusionFFPE.input_fastq2"] = args.fq2
    input_json["StarFusionFFPE.ctat_lib"] = os.path.abspath(args.ctat)
    input_json["StarFusionFFPE.kinases_domains"] = os.path.abspath(args.kinases)
    input_json["StarFusionFFPE.exons"] = os.path.abspath(args.exons)
    if args.lncrna is not None:
        input_json["StarFusionFFPE.lncrna_list"] = os.path.abspath(args.lncrna)
    input_json["StarFusionFFPE.threads"] = args.threads
    input_json["StarFusionFFPE.star_memory"] = args.star_mem
    if args.no_kinase_info:
        input_json["StarFusionFFPE.with_kinase_info"] = 'false'
    input_json["StarFusionFFPE.docker"] = {"StarFusion": "docker.io/trinityctat/starfusion:1.12.0"}

    with open(args.out, 'w') as f:
        json.dump(input_json, f, indent=4)

    p = subprocess.Popen(
        ["java", f"-Xmx{args.cromwell_mem}g", f"-Dconfig.file={args.config}", "-jar", f"{args.cromwell}",
         "run", f"{args.wdl}", "--options", f"{args.options}", "--inputs", f"{args.out}"])
    p.wait()
