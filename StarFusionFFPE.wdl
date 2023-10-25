version 1.0

workflow StarFusionFFPE {

    input {
        File input_fastq
        File? input_fastq2
        String ctat_lib
        File kinases_domains
        File exons
        File? lncrna_list
        Int threads
        Int star_memory
        Boolean with_kinase_info = true
        Map[String, String] docker
    }
    String res_name = basename(basename(basename(input_fastq, ".gz"), ".fastq"), ".fq")
    Boolean paired_end = defined(input_fastq2)

    call StarFusion {
        input:
            left_fastq = input_fastq,
            right_fastq = input_fastq2,
            ctat_lib = ctat_lib,
            threads = threads,
            star_memory = star_memory,
            docker = docker["StarFusion"]
    }

    call FilterStarFusion as FilterStarFusionMain{
        input:
            nofilter_tsv = StarFusion.tsv,
            res_name = res_name,
            docker = docker["StarFusion"]
    }

    if (defined(lncrna_list)) {
        call FilterStarFusion  as FilterStarFusionLncRNA {
            input:
                nofilter_tsv = StarFusion.tsv,
                lncrna_list = lncrna_list,
                res_name = res_name,
                docker = docker["StarFusion"]
        }
    }

    call FusionInspector {
        input:
            left_fastq = input_fastq,
            right_fastq = input_fastq2,
            ctat_lib = ctat_lib,
            fusions = FilterStarFusionMain.tsv,
            fusions_lncrna = FilterStarFusionLncRNA.tsv,
            threads = threads,
            star_memory = star_memory,
            docker = docker["StarFusion"]
    }

    call BamToSamRmDup as BamToSamRmDupJunction {
        input:
            bam = FusionInspector.bam_j,
            res_name = res_name,
            junction = true,
            docker = docker["StarFusion"]
    }

    if (paired_end) {
        call BamToSamRmDup as BamToSamRmDupSpanning {
            input:
                bam = FusionInspector.bam_s,
                res_name = res_name,
                junction = false,
                docker = docker["StarFusion"]
        }
    }

    call FilterFusionInspector {
        input:
            paired_end = paired_end,
            nofilter_tsv = FusionInspector.tsv,
            nofilter_sam_j = BamToSamRmDupJunction.sam,
            nofilter_sam_s = BamToSamRmDupSpanning.sam,
            lncrna_list = lncrna_list,
            res_name = res_name,
            docker = docker["StarFusion"]
    }

    if (with_kinase_info) {
        call KinasesAndExons {
            input:
                fusions = FilterFusionInspector.tsv,
                star_fusion_data = false,
                exons = exons,
                kinases_domains = kinases_domains,
                res_name = res_name,
                docker = docker["StarFusion"]
        }

        call AddReads {
            input:
                noreads_tsv = KinasesAndExons.tsv,
                sam = FilterFusionInspector.sam_j,
                res_name = res_name,
                docker = docker["StarFusion"]
        }
    }

    call SamToBam as SamToBamJunction {
        input:
            sam = FilterFusionInspector.sam_j,
            res_name = res_name,
            junction = true,
            docker = docker["StarFusion"]
    }

    if (paired_end) {
        File sam_s = select_first([FilterFusionInspector.sam_s,])
        call SamToBam as SamToBamSpanning {
            input:
                sam = sam_s,
                res_name = res_name,
                junction = false,
                docker = docker["StarFusion"]
        }
    }

    output {
        File fusions_tsv = select_first([AddReads.tsv, FilterFusionInspector.tsv])
        File fusions_j_bam = SamToBamJunction.bam
        File fusions_j_bai = SamToBamJunction.bai
        File? fusions_s_bam = SamToBamSpanning.bam
        File? fusions_s_bai = SamToBamSpanning.bai
    }
}

task StarFusion {
    input {
        File left_fastq
        File? right_fastq
        String ctat_lib
        Int threads
        Int star_memory
        String docker
    }

    command <<<
        STAR-Fusion --left_fq ~{left_fastq} \
        ~{if defined(right_fastq) then "--right_fq " + right_fastq else ""} \
        --genome_lib_dir /~{basename(ctat_lib)} \
        --CPU ~{threads} --require_LDAS 0 \
        --skip_FFPM --min_sum_frags 1 \
        --examine_coding_effect
    >>>

    runtime {
        ctat_dir: ctat_lib
        cpu: threads
        memory: star_memory
        docker: docker
    }

    output {
        File tsv = "STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    }
}

task FilterStarFusion {
    input {
        File? lncrna_list
        File nofilter_tsv
        String res_name
        String docker
    }

    command <<<
        python3 <<CODE
        with open("~{nofilter_tsv}", "r") as f_in:
            data = f_in.read().strip().split("\n")

        ~{if defined(lncrna_list) then
            'with open("~{lncrna_list}", "r") as f_in: lncrna = f_in.read().strip().split("\\n")'
        else
            ''
        }

        header = data[0]
        data = [x for x in data if x.split("\t")[5] == "ONLY_REF_SPLICE" and \
        ~{if defined(lncrna_list) then
            '(x.split("\\t")[0].replace("--", "|").split("|")[0] in lncrna or ' +
            'x.split("\\t")[0].replace("--", "|").split("|")[1] in lncrna)'
        else
            'x.split("\\t")[20] == "INFRAME"'
        }]

        with open("~{res_name}~{if defined(lncrna_list) then '.lnc' else ''}.sf.tsv", "w") as f_out:
            f_out.write(header + "\n")
            if len(data): f_out.write("\n".join(data) + "\n")
        CODE
    >>>

    runtime {
        cpu: 1
        memory: 4096
        docker: docker
    }

    output {
        File tsv = "~{res_name}~{if defined(lncrna_list) then '.lnc' else ''}.sf.tsv"
    }
}

task FusionInspector {
    input {
        File left_fastq
        File? right_fastq
        String ctat_lib
        File fusions
        File? fusions_lncrna
        Int threads
        Int star_memory
        String docker
    }

    command <<<
        FusionInspector --fusions ~{fusions}~{if defined(fusions_lncrna) then ',' + fusions_lncrna  else ''} \
        --left_fq ~{left_fastq} \
        ~{if defined(right_fastq) then "--right_fq " + right_fastq else ""} \
        --min_junction_reads 1  --min_novel_junction_support 3 \
        --min_spanning_frags_only 5 --min_sum_frags 1 \
        --max_mate_dist 100000 --max_promiscuity 10 \
        --genome_lib_dir /~{basename(ctat_lib)} --CPU ~{threads} \
        --require_LDAS 0 --no_FFPM --annotate --vis --cleanup \
        --examine_coding_effect --fusion_contigs_only
    >>>

    runtime {
        ctat_dir: ctat_lib
        cpu: threads
        memory: star_memory
        docker: docker
    }

    output {
        File tsv = "FI/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect"
        File bam_j = "FI/finspector.junction_reads.bam"
        File bam_s = "FI/finspector.spanning_reads.bam"
    }
}

task BamToSamRmDup {
    input {
        File bam
        String res_name
        Boolean junction
        String docker
    }
    String prefix = if junction then '.junction' else '.spanning'

    command <<<
        samtools rmdup -s ~{bam} - | samtools view -h > \
        ~{res_name}~{prefix}.nofilter.sam
    >>>

    runtime {
        continueOnReturnCode: true
        cpu: 1
        memory: 1024
        docker: docker
    }

    output {
        File sam = "~{res_name}~{prefix}.nofilter.sam"
    }
}

task FilterFusionInspector {
    input {
        Boolean paired_end
        File nofilter_tsv
        File nofilter_sam_j
        File? nofilter_sam_s
        File? lncrna_list
        String res_name
        String docker
    }

    command <<<
        python3 <<CODE
        ~{if defined(lncrna_list) then
            'with open("~{lncrna_list}", "r") as f_in: lncrna = f_in.read().strip().split("\\n")'
        else
            ''
        }

        with open("~{nofilter_tsv}", "r") as f_in:
            data = f_in.read().strip().split("\n")
        header = data[0]
        data = [x.split("\t") for x in data[1:]]
        data = [x for x in data if x[11] == "ONLY_REF_SPLICE" and (\
        ~{if defined(lncrna_list) then
            'x[0].replace("--", "|").split("|")[0] in lncrna or ' +
            'x[0].replace("--", "|").split("|")[1] in lncrna or '
        else
            ''
        } x[28] == "INFRAME")]

        with open("~{nofilter_sam_j}", "r") as f_in:
            reads_j = f_in.read().strip().split("\n")
        header_j = [x for x in reads_j if x[0] == "@"]

        paired_end = ~{if paired_end then 'True' else 'False'}
        res_reads_j = []
        for i in range(0, len(data)):
            fusion = data[i]
            fusion_name = fusion[0]
            fusion_pos = [str(int(fusion[6]) + 1), str(int(fusion[9]) - 1)]

            fusion_reads = [x for x in reads_j[len(header_j):] if fusion_name == x.split("\t")[2] and
                            fusion_pos == x.split("\t")[18].split(",")[-2:]]

            fusion[1] = str(len(fusion_reads))
            data[i] = fusion

            if len(fusion_reads) > 1:
                res_reads_j = res_reads_j + fusion_reads

        res_reads_j = list(set(res_reads_j))
        data = [x for x in data if int(x[1]) > 1]
        data = sorted(data, key=lambda x: int(x[1]), reverse=True)

        res_fusion_names = [x[0] for x in data]
        header_j = [header_j[0]] + [x for x in header_j[1:] if x.split('\t')[1].split(':')[1] in res_fusion_names]

        if paired_end:
            with open("~{nofilter_sam_s}", "r") as f_in:
                reads_s = f_in.read().strip().split("\n")

            header_s = [x for x in reads_s if x[0] == "@"]
            reads_s = [x for x in reads_s[len(header_s):] if x.split('\t')[2] in res_fusion_names]
            header_s = [header_s[0]] + [x for x in header_s[1:] if x.split('\t')[1].split(':')[1] in res_fusion_names]

            with open("~{res_name}.spanning.sam", "w") as f_out:
                f_out.write("\n".join(header_s) + "\n")
                if len(reads_s): f_out.write("\n".join(reads_s) + "\n")

        data = ["\t".join(x) for x in data]
        with open("~{res_name}.nokinase.tsv", "w") as f_out:
            f_out.write(header + "\n")
            if len(data): f_out.write("\n".join(data) + "\n")

        with open("~{res_name}.junction.sam", "w") as f_out:
            f_out.write("\n".join(header_j) + "\n")
            if len(res_reads_j): f_out.write("\n".join(res_reads_j) + "\n")
        CODE
    >>>

    runtime {
        cpu: 1
        memory: 4096
        docker: docker
    }

    output {
        File tsv = "~{res_name}.nokinase.tsv"
        File sam_j = "~{res_name}.junction.sam"
        File? sam_s = "~{res_name}.spanning.sam"
    }
}

task KinasesAndExons {
    input {
        File fusions
        Boolean star_fusion_data
        File exons
        File kinases_domains
        String res_name
        String docker
    }
    String prefix = if star_fusion_data then '' else 'fi.'

    command <<<
        python3 <<CODE
        def parse_exons(fm, tr):
            en = []
            for elem in fm:
                start, end = [int(x) for x in elem[3:-3].split("-")]
                en += [x["number"] for x in tr if start == x["start"] or end == x["end"] or
                       max(0, min(end, x["end"]) - max(start, x["start"]))]
            return en


        def post_process_exons(exons, transcript=None, kinase=False):
            if len(exons) > 1:
                a, b = min(exons[0], exons[-1]), max(exons[0], exons[-1]) + 1
                concat = True if (list(range(a, b)) ==
                                  list(range(min(exons), max(exons) + 1))) else False
                if concat:
                    exons = f"{exons[0]}-{exons[-1]}"
                else:
                    exons = ",".join([str(x) for x in exons])
            elif len(exons) == 1:
                exons = str(exons[0])
            elif kinase:
                exons = "."
            elif transcript is not None and len(transcript) == 1:
                exons = f"{transcript[0]['number']}"
            else:
                exons = "Undefined_behaviour!"
            return exons


        def parse_uniprot(transcript_id, cds_position, exons, fm, tr_ex, kinase, strand):
            kin_en = []
            kin_en_extra = []
            uniprot = []
            cds_start, cds_end = cds_position.split("-")
            cds_start, cds_end = int(cds_start), int(cds_end)
            kinase_positions = fusion_transcripts[transcript_id][1].split("|")
            to_reverse = False

            for i in range(len(kinase_positions)):
                kinase_start, kinase_end = kinase_positions[i].split("-")
                kinase_start, kinase_end = int(kinase_start), int(kinase_end)

                if max(0, min(kinase_end * 3, cds_end) - max(kinase_start * 3, cds_start)):
                    kinase = True
                    if fusion_transcripts[transcript_id][0] not in uniprot:
                        uniprot.append(fusion_transcripts[transcript_id][0])

                    if kinase_start * 3 < cds_start:
                        cds_current = cds_start
                        exons_smaller = list(filter(lambda x: x['number'] < min(exons), tr_ex))
                        for elem in reversed(exons_smaller):
                            start, end = elem["start"], elem["end"]
                            diff = abs(end - start) + 1
                            if max(0, min(kinase_end * 3, cds_current) - max(kinase_start * 3, cds_current - diff)):
                                kin_en_extra.append(elem["number"])
                            else:
                                break
                            cds_current -= diff

                    exon_shift = 0
                    cds_current = cds_start
                    if strand == "-" and len(exons) > 0 and exons[0] > exons[-1]:
                        to_reverse = True
                        exons = list(reversed(exons))
                    for elem in fm:
                        start, end = elem[3:-3].split("-")
                        diff = abs(int(end) - int(start)) + 1
                        if max(0, min(kinase_end * 3, cds_current + diff) - max(kinase_start * 3, cds_current)):
                            kin_en.append(exons[exon_shift])
                        cds_current += diff
                        exon_shift += 1

                    if kinase_end * 3 > cds_end:
                        exons_bigger = list(filter(lambda x: x['number'] > max(exons), tr_ex))
                        for elem in exons_bigger:
                            start, end = elem["start"], elem["end"]
                            diff = abs(end - start) + 1
                            if max(0, min(kinase_end * 3, cds_current + diff) - max(kinase_start * 3, cds_current)):
                                kin_en_extra.append(elem["number"])
                            else:
                                break
                            cds_current += diff

            if len(kin_en) > 0:
                if kin_en[0] > kin_en[-1]:
                    kin_en = sorted(kin_en + kin_en_extra, reverse=True)
                else:
                    kin_en = sorted(kin_en + kin_en_extra, reverse=False)
            if to_reverse:
                kin_en = list(reversed(kin_en))
            return kinase, uniprot, kin_en


        def post_process_uniprot(is_kinase_, uniprot_):
            if is_kinase_:
                is_kinase_ = "TRUE"
                if len(uniprot_) == 2:
                    uniprot_ = "|".join(uniprot_)
                elif len(uniprot_) == 1:
                    uniprot_ = uniprot_[0]
                elif len(uniprot_) == 0:
                    uniprot_ = "."
                else:
                    uniprot_ = "Undefined_behaviour!"
            else:
                is_kinase_, uniprot_ = "FALSE", "."
            return is_kinase_, uniprot_


        with open("~{exons}", "r") as f_in:
            from json import load
            all_exons = load(f_in)

        fusion_transcripts = {}
        with open("~{kinases_domains}", "r") as f_in:
            data = f_in.read().strip().split("\n")
            for elem in data[1:]:
                s = elem.strip().split("\t")
                fusion_transcripts[s[1]] = (s[3], s[4])

        column_fm = ~{if star_fusion_data then '21' else '29'}
        column_id_left = ~{if star_fusion_data then '16' else '24'}
        column_id_right = ~{if star_fusion_data then '18' else '26'}
        column_pfam_left = ~{if star_fusion_data then '24' else '32'}
        column_pfam_right = ~{if star_fusion_data then '25' else '33'}

        with open("~{res_name}.noreads.~{prefix}tsv", "w") as f_out:
            with open("~{fusions}", "r") as f_in:
                s = f_in.readline().strip().split("\t")
                s.append("IS_KINASE")
                s.append("UNIPROT")
                s.append("EXONS_LEFT")
                s.append("EXONS_RIGHT")
                s.append("KINASE_EXONS_LEFT")
                s.append("KINASE_EXONS_RIGHT")
                s = "\t".join(s)
                f_out.write(f"{s}\n")

                s = f_in.readline()
                while len(s) > 0:
                    s = s.strip().split("\t")
                    # Add exons
                    fm_l, fm_r = None, None
                    en_l, en_r = None, None
                    if s[column_fm] == ".":
                        en_lpp = "."
                        en_rpp = "."
                    else:
                        fm_l, fm_r = s[column_fm].split("<==>")
                        fm_l_strand, fm_r_strand = fm_l.split("|")[1], fm_r.split("|")[1]
                        fm_l, fm_r = fm_l.split("|")[2:], fm_r.split("|")[2:]
                        tr_l, tr_r = all_exons[s[column_id_left]], all_exons[s[column_id_right]]
                        en_l = parse_exons(fm_l, tr_l)
                        en_r = parse_exons(fm_r, tr_r)
                        en_lpp = post_process_exons(en_l, tr_l)
                        en_rpp = post_process_exons(en_r, tr_r)
                    # Add kinase info
                    is_kinase = False
                    uniprots = []
                    kin_en_l, kin_en_r = [], []
                    uniprot_l, uniprot_r = [], []

                    if s[column_id_left] != "." and s[column_id_left].split(".")[0] in fusion_transcripts or \
                            s[column_id_right] != "." and s[column_id_right].split(".")[0] in fusion_transcripts:
                        # Left
                        if s[column_id_left].split(".")[0] in fusion_transcripts:
                            is_kinase, uniprot_l, kin_en_l = parse_uniprot(transcript_id=s[column_id_left].split(".")[0],
                                                                           cds_position=s[column_id_left + 1], exons=en_l,
                                                                           fm=fm_l, tr_ex=tr_l, kinase=is_kinase, strand=fm_l_strand)
                        # Right
                        elif s[column_id_right].split(".")[0] in fusion_transcripts:
                            is_kinase, uniprot_r, kin_en_r = parse_uniprot(transcript_id=s[column_id_right].split(".")[0],
                                                                           cds_position=s[column_id_right + 1], exons=en_r,
                                                                           fm=fm_r, tr_ex=tr_r, kinase=is_kinase, strand=fm_r_strand)
                        uniprots = uniprot_l + uniprot_r
                    if "kinase" in s[column_pfam_left] or "kinase" in s[column_pfam_right] or \
                       "Kinase" in s[column_pfam_left] or "Kinase" in s[column_pfam_right]:
                        is_kinase = True

                    is_kinase, uniprots = post_process_uniprot(is_kinase, uniprots)
                    s.append(is_kinase)
                    s.append(uniprots)
                    s.append(en_lpp)
                    s.append(en_rpp)
                    kin_en_lpp = post_process_exons(kin_en_l, kinase=True)
                    kin_en_rpp = post_process_exons(kin_en_r, kinase=True)
                    s.append(kin_en_lpp)
                    s.append(kin_en_rpp)
                    s = "\t".join(s) + "\n"
                    f_out.write(s)
                    s = f_in.readline()
        CODE
    >>>

    runtime {
        cpu: 1
        memory: 2048
        docker: docker
    }

    output {
        File tsv = "~{res_name}.noreads.~{prefix}tsv"
    }
}

task AddReads {
    input {
        File noreads_tsv
        File sam
        String res_name
        String docker
    }

    command <<<
        python3 <<CODE
        with open("~{noreads_tsv}", "r") as f_in:
            data = f_in.read().strip().split("\n")
        data[0] = data[0].strip() + "\tJUNCTION_READS"

        with open("~{sam}", "r") as f_in:
            reads = f_in.read().strip().split("\n")
        reads = [x for x in reads if x[0] != "@"]

        for i in range(1, len(data)):
            fusion = data[i].split("\t")
            fusion_name = fusion[0]
            fusion_pos = [str(int(fusion[6]) + 1), str(int(fusion[9]) - 1)]
            fusion_reads = [x.split("\t")[9] for x in reads if fusion_name == x.split("\t")[2] and
                            fusion_pos == x.split("\t")[18].split(",")[-2:]]
            fusion.append(",".join(fusion_reads))
            data[i] = fusion
        
        data = [data[0]] + sorted(data[1:], key=lambda x: int(x[1]), reverse=True)
        data = [data[0]] + ["\t".join(x) for x in data[1:]]

        with open("~{res_name}.tsv", "w") as f_out:
            f_out.write("\n".join(data) + "\n")
        CODE
    >>>

    runtime {
        cpu: 1
        memory: 2048
        docker: docker
    }

    output {
        File tsv = "~{res_name}.tsv"
    }
}

task SamToBam {
    input {
        File sam
        String res_name
        Boolean junction
        String docker
    }
    String prefix = if junction then '.junction' else '.spanning'

    command <<<
        touch ~{res_name}~{prefix}.bam; \
        touch ~{res_name}~{prefix}.bam.bai; \
        samtools sort ~{sam} | \
        samtools view -bh - > ~{res_name}~{prefix}.bam; \
        samtools index -b ~{res_name}~{prefix}.bam > \
        ~{res_name}~{prefix}.bam.bai;
    >>>

    runtime {
        continueOnReturnCode: true
        cpu: 1
        memory: 1024
        docker: docker
    }

    output {
        File bam = "~{res_name}~{prefix}.bam"
        File bai = "~{res_name}~{prefix}.bam.bai"
    }
}
