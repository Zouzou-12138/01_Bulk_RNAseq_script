#!/usr/bin/env python3
"""
get_expr_normal_matrix.py
功能：基于原始计数矩阵、GTF 注释和 BAM 文件，计算 FPKM 和 TPM 标准化矩阵。
用法：python get_expr_normal_matrix.py <input_matrix_csv> <gtf_file> <output_fpkm_csv> <output_tpm_csv> <bam_file_1> <bam_file_2> ...
注意：BAM 文件的顺序不需要与矩阵列严格对应，脚本会通过文件名自动匹配样本。
      假设 BAM 文件名为 {sample}.sorted.bam，矩阵列名为 {sample}。
"""
# !/usr/bin/env python3
"""
get_expr_normal_matrix.py
功能：基于原始计数矩阵、GTF 注释和 BAM 文件，计算 FPKM 和 TPM 标准化矩阵。
用法：python get_expr_normal_matrix.py <input_matrix_csv> <gtf_file> <output_fpkm_csv> <output_tpm_csv> <bam_file_1> <bam_file_2> ...
"""

import sys
import os
import re
import subprocess
import pandas as pd
import numpy as np


def calculate_union_length(intervals):
    """合并重叠的区间并计算总长度。"""
    if not intervals:
        return 0
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            last[1] = max(last[1], current[1])
        else:
            merged.append(current)
    return sum(iv[1] - iv[0] + 1 for iv in merged)


def parse_gtf_lengths(gtf_path):
    """解析 GTF 文件，计算每个基因的非重叠外显子（Union Exons）总长度。"""
    print(f"[*] Parsing GTF for gene lengths: {gtf_path}")
    gene_intervals = {}
    gene_id_pattern = re.compile(r'gene_id\s+"?([^";\s]+)"?')

    try:
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'exon': continue

                start, end = int(parts[3]), int(parts[4])
                match = gene_id_pattern.search(parts[8])
                if match:
                    gene_id = match.group(1)
                    if gene_id not in gene_intervals:
                        gene_intervals[gene_id] = []
                    gene_intervals[gene_id].append([start, end])
    except Exception as e:
        print(f"[-] Error parsing GTF: {e}")
        sys.exit(1)

    gene_lengths = {gid: calculate_union_length(ivs) for gid, ivs in gene_intervals.items()}
    print(f"[*] Found lengths for {len(gene_lengths)} genes.")
    return gene_lengths


def get_mapped_reads(bam_path):
    """使用 samtools flagstat 获取比对上的 Reads 总数。"""
    try:
        # 使用 -@ 4 加速 samtools
        result = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.split('\n'):
            # 匹配具体的 mapped 行，排除 secondary 和 supplementary
            if 'mapped (' in line and 'secondary' not in line and 'supplementary' not in line:
                return int(line.split()[0])
        return 0
    except Exception as e:
        print(f"[-] Error running samtools on {bam_path}: {e}")
        return 0


def main():
    if len(sys.argv) < 6:
        print("Usage: python get_expr_normal_matrix.py <input_matrix> <gtf> <out_fpkm> <out_tpm> <bam1> [bam2 ...]")
        sys.exit(1)

    input_matrix, gtf_file, out_fpkm, out_tpm = sys.argv[1:5]
    bam_files = sys.argv[5:]

    # 1. 获取基因长度
    gene_lengths = parse_gtf_lengths(gtf_file)

    # 2. 读取计数矩阵
    print(f"[*] Reading count matrix: {input_matrix}")
    df_counts = pd.read_csv(input_matrix, index_col=0)

    # 3. 匹配 BAM 文件与样本
    sample_bam_map = {}
    for bam in bam_files:
        basename = os.path.basename(bam)
        # 提取核心样本名 (例如 sample1.sorted.bam -> sample1)
        s_name = re.sub(r'(\.sorted)?\.bam$', '', basename)
        sample_bam_map[s_name] = bam

    # 4. 获取每个样本的库大小 (Library Size)
    print("[*] Calculating mapped reads from BAMs...")
    mapped_reads = {}
    for sample in df_counts.columns:
        bam_path = sample_bam_map.get(sample)
        if not bam_path:
            # 模糊匹配尝试
            for s_key, b_path in sample_bam_map.items():
                if sample in s_key or s_key in sample:
                    bam_path = b_path
                    break

        if bam_path:
            count = get_mapped_reads(bam_path)
            mapped_reads[sample] = count if count > 0 else 1
        else:
            print(f"  [!] Warning: No BAM found for {sample}. Using 1.")
            mapped_reads[sample] = 1

    # 5. 向量化计算 FPKM 和 TPM
    # 准备长度向量 (只保留矩阵中存在的基因)
    df_counts = df_counts.loc[df_counts.index.isin(gene_lengths.keys())]
    L = np.array([gene_lengths[g] for g in df_counts.index])
    N = np.array([mapped_reads[s] for s in df_counts.columns])

    print("[*] Calculating FPKM and TPM via vectorization...")

    # FPKM = (C * 10^9) / (L * N)
    # 使用 Pandas 的广播机制
    df_fpkm = (df_counts.mul(1e9)).div(L, axis=0).div(N, axis=1)

    # TPM = (FPKM / sum(FPKM)) * 10^6
    df_tpm = df_fpkm.div(df_fpkm.sum(axis=0), axis=1).mul(1e6)

    # 6. 保存
    print(f"[*] Saving: \n  - {out_fpkm}\n  - {out_tpm}")
    df_fpkm.to_csv(out_fpkm)
    df_tpm.to_csv(out_tpm)
    print("[*] Done!")


if __name__ == "__main__":
    main()
