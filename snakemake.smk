import os

# 1. 加载配置
configfile: "config.yml"

DIRS = config["directories"]
PATTERN = config["file_pattern"]
PROJ = config["project"]

# 2. 初始化：自动创建所有输出目录
for d_path in DIRS.values():
    os.makedirs(d_path, exist_ok=True)

# 3. 动态识别样本名
# 假设文件名为 SampleA_L1_1.fq.gz，则 {sample} 会匹配到 SampleA
SAMPLES, = glob_wildcards(DIRS["raw"] + "/{sample}" + PATTERN["r1_suffix"])

# 4. 目标任务汇总
rule all:
    input:
        expand(DIRS["qc"] + "/{sample}_L1_1_fastqc.html", sample=SAMPLES),
        expand(DIRS["qc"] + "/{sample}_1.clean.fq.gz", sample=SAMPLES),
        expand(DIRS["mapping"] + "/{sample}.sorted.bam", sample=SAMPLES),
        expand(DIRS["novel"] + "/{sample}.gtf", sample=SAMPLES),
        expand(DIRS["expression"] + "/{sample}.counts", sample=SAMPLES),
        DIRS["expression"] + "/gene_expression_matrix.csv",
        DIRS["expression"] + "/gene_expression_FPKM.csv",
        DIRS["expression"] + "/gene_expression_TPM.csv",
        DIRS["qc"] + "/multiqc_report.html"



# --- 分析规则 ---

# 1. 原始数据质控
rule fastqc:
    input:
        r1 = DIRS["raw"] + "/{sample}" + PATTERN["r1_suffix"],
        r2 = DIRS["raw"] + "/{sample}" + PATTERN["r2_suffix"]
    output:
        # 必须匹配 FastQC 默认生成的名称：原文件名去掉 .fq.gz 加上 _fastqc.html
        html1 = DIRS["qc"] + "/{sample}_L1_1_fastqc.html",
        html2 = DIRS["qc"] + "/{sample}_L1_2_fastqc.html",
        zip1  = DIRS["qc"] + "/{sample}_L1_1_fastqc.zip",
        zip2  = DIRS["qc"] + "/{sample}_L1_2_fastqc.zip"
    threads: 4
    shell:
        "fastqc -t {threads} -o {DIRS[qc]} {input.r1} {input.r2}"

# 2.1 去接头与质量过滤
rule cutadapt:
    input:
        r1 = DIRS["raw"] + "/{sample}" + PATTERN["r1_suffix"],
        r2 = DIRS["raw"] + "/{sample}" + PATTERN["r2_suffix"]
    output:
        c1 = DIRS["qc"] + "/{sample}_1.clean.fq.gz",
        c2 = DIRS["qc"] + "/{sample}_2.clean.fq.gz",
        log = DIRS["qc"] + "/{sample}.cutadapt.log"
    threads: 8
    shell:
        """
        cutadapt -a {PROJ[adapter_r1]} -A {PROJ[adapter_r2]} \
        -q 20,20 --trim-n -m 75 \
        -o {output.c1} -p {output.c2} {input.r1} {input.r2} \
        > {output.log} 2>&1
        """

# 2.2 生成MultiQC 汇总报告

rule multiqc:
    input:
        # 这里列出需要被 MultiQC 扫描的目录
        # 注意：MultiQC 会自动递归扫描这些目录下的所有支持的工具日志
        qc_dir = DIRS["qc"],
        mapping_dir = DIRS["mapping"]
    output:
        # 报告输出路径
        report = DIRS["qc"] + "/multiqc_report.html",
        # 可选：如果你希望明确控制数据和图片目录的输出位置（通常 MultiQC 会自动创建）
        # 如果不需要在 rule all 中显式依赖这两个文件夹，可以注释掉下面两行
        data_dir = DIRS["qc"] + "/multiqc_data", 
        plots_dir = DIRS["qc"] + "/multiqc_plots"
    params:
        title = PROJ.get("name", "Bulk RNA-Seq Project"), # 增加 .get() 防止 config 中缺少 name 报错
        outdir = DIRS["qc"]
    threads: 4
    shell:
        """
        multiqc {input.qc_dir} {input.mapping_dir} \
            --title "{params.title}" \
            --outdir {params.outdir} \
            --force \
            --filename multiqc_report.html
        """

# 3. 比对并排序
rule hisat2:
    input:
        c1 = rules.cutadapt.output.c1,
        c2 = rules.cutadapt.output.c2
    output:
        bam = DIRS["mapping"] + "/{sample}.sorted.bam",
        bai = DIRS["mapping"] + "/{sample}.sorted.bam.bai"
    threads: PROJ["threads"]
    shell:
        """
        set -o pipefail
        hisat2 -p {threads} --dta -x {PROJ[ref_genome]} \
        -1 {input.c1} -2 {input.c2} | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """

# 4. 转录本组装（用于发现新转录本）
rule stringtie:
    input:
        bam = rules.hisat2.output.bam,
        gtf = PROJ["gtf"]
    output:
        gtf = DIRS["novel"] + "/{sample}.gtf"
    threads: 8
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gtf} -o {output.gtf} -l {wildcards.sample}"

# 5. 基因定量
rule htseq_count:
    input:
        bam = rules.hisat2.output.bam,
        gtf = PROJ["gtf"]
    output:
        counts = DIRS["expression"] + "/{sample}.counts"
    threads: 4 
    shell:
        # 注意：HTSeq 在处理大文件时较慢，增加 -n 参数可以提高多核效率（如果版本支持）
        "htseq-count -f bam -r pos -s no -t exon -i gene_id {input.bam} {input.gtf} > {output.counts}"

# 6. 生成基因差异表达矩阵
rule merge_counts:
    input:
        counts = expand(DIRS["expression"] + "/{sample}.counts", sample=SAMPLES)
    output:
        matrix = DIRS["expression"] + "/gene_expression_matrix.csv"
    params:
        script = "get_expr_matrix.py"
    shell:
        """
        python {params.script} {output.matrix} {input.counts}
        """

# 7. 计算 FPKM 和 TPM
rule calculate_normalization:
    input:
        matrix = rules.merge_counts.output.matrix,
        gtf = PROJ["gtf"],
        bams = expand(DIRS["mapping"] + "/{sample}.sorted.bam", sample=SAMPLES)
    output:
        fpkm = DIRS["expression"] + "/gene_expression_FPKM.csv",
        tpm = DIRS["expression"] + "/gene_expression_TPM.csv"
    params:
        script = "get_expr_normal_matrix.py"
    shell:
        """
        python {params.script} {input.matrix} {input.gtf} {output.fpkm} {output.tpm} {input.bams}
        """
       
