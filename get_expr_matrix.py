#!/usr/bin/env python3
"""
get_expr_matrix.py
功能：合并多个样本的 HTSeq count 文件为一个基因表达矩阵 (CSV)。
用法：python get_expr_matrix.py <output_csv> <count_file_1> <count_file_2> ...
"""

import sys
import os
import pandas as pd


def main():
    if len(sys.argv) < 3:
        print("Usage: python get_expr_matrix.py <output_csv> <count_file_1> [count_file_2 ...]")
        sys.exit(1)

    output_file = sys.argv[1]
    count_files = sys.argv[2:]

    if not count_files:
        print("Error: No count files provided.")
        sys.exit(1)

    df_list = []

    print(f"Processing {len(count_files)} count files...")

    for f_path in count_files:
        if not os.path.exists(f_path):
            print(f"Warning: File not found: {f_path}, skipping.")
            continue

        # 从文件名提取样本名 (去掉路径和 .counts 后缀)
        sample_name = os.path.basename(f_path).replace(".counts", "")

        try:
            # 读取 HTSeq 输出 (无表头，两列: GeneID, Count)
            df = pd.read_csv(f_path, sep="\t", header=None, names=["gene_id", sample_name])

            # 过滤掉 HTSeq 的统计行 (以 __ 开头)
            df = df[~df["gene_id"].str.startswith("__")]

            # 设置索引
            df.set_index("gene_id", inplace=True)
            df_list.append(df)
            print(f"  - Loaded {sample_name}: {len(df)} genes")

        except Exception as e:
            print(f"Error reading {f_path}: {e}")
            sys.exit(1)

    if not df_list:
        print("Error: No valid data frames to merge.")
        sys.exit(1)

    # 横向合并 (Outer Join)，缺失值填 0
    print("Merging matrices...")
    merged_df = pd.concat(df_list, axis=1, join="outer").fillna(0)

    # 确保数据类型为整数
    merged_df = merged_df.astype(int)

    # 保存
    print(f"Saving to {output_file}...")
    merged_df.to_csv(output_file)
    print("Done.")


if __name__ == "__main__":
    main()