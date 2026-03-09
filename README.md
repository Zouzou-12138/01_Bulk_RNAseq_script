# 标准RNA-Seq分析流程

## 1. 项目概述

本项目是一个基于Snakemake构建的自动化Bulk RNA-seq分析流程，适用于小鼠（Mus musculus，GRCm39）转录组数据处理，遵循行业标准。

**核心特点：**

* **声明式工作流**：由Snakemake管理，支持目标驱动执行和错误恢复
* **标准化质量控制**：集成FastQC和Cutadapt过滤（Q20， 长度>75bp）
* **新转录本预测**：使用StringTie进行从头转录本组装
* **结构化输出**：自动生成00-13编号目录系统，符合专业生物信息学报告标准

---

## 2. 目录结构

执行流程时，脚本会自动创建以下目录结构：

* **00_Rawdata**：原始双端fastq文件输入目录（.fq.gz格式）
* **01_QC**：FastQC质控报告和Cutadapt修剪后数据
* **02_Mapping**：HISAT2比对结果（排序后的BAM文件及索引）
* **03_AS**：可变剪接分析占位目录
* **04_Novel**：StringTie转录本组装和新异构体GTF文件
* **05_SNV**：SNV/Indel分析占位目录
* **06_GeneExpression**：HTSeq-count生成的原始counts矩阵
* **07-13**：下游分析预创建目录（DE、GO、KEGG、PCA、GSEA等）

---

## 3. 脚本位置

**所有分析脚本均位于 `01_Bulk_RNAseq_script` 目录下**

```bash
# 进入脚本目录
cd /path/to/01_Bulk_RNAseq_script
```

---

## 4. 环境配置

建议使用提供的environment.yml文件创建一致的运行环境：

```bash
# 创建环境
conda env create -f environment.yml

# 激活环境
conda activate Bulk_RNA_Seq
```

---

## 5. 数据准备

### 5.1 复制原始数据

**将原始测序数据复制到脚本目录下的 `00_Rawdata` 文件夹**（如目录不存在，脚本会自动创建）：

```bash
# 确保在脚本目录下
cd /path/to/01_Bulk_RNAseq_script

# 复制原始数据到00_Rawdata目录
cp /path/to/your/rawdata/*.fq.gz ./00_Rawdata/
```

### 5.2 文件命名规范

原始fastq文件必须遵循以下命名规则：
- 样本名称必须与后续分析对应
- 格式：`{样本ID}_L1_1.fq.gz` 和 `{样本ID}_L1_2.fq.gz`

示例：
- `WT-testis-1_L1_1.fq.gz`
- `WT-testis-1_L1_2.fq.gz`
- `C2-KO-testis-1_L1_1.fq.gz`
- `C2-KO-testis-1_L1_2.fq.gz`

---

## 6. 配置文件（config.yml）

运行前，请根据实际情况修改 `config.yml` 中的绝对路径：

* **threads**：CPU核心分配（默认：36）
* **ref_genome**：HISAT2索引前缀路径
* **gtf**：参考注释文件路径（如 GRCm39.110.gtf）
* **adapter**：R1和R2接头序列

---

## 7. 运行流程

### 7.1 方法一：使用Python包装脚本（推荐）

```bash
# 确保在脚本目录下
cd /path/to/01_Bulk_RNAseq_script

# 后台运行完整流程（会自动初始化目录结构）
nohup python run_pipeline.py > pipeline.log 2>&1 &
```

### 7.2 方法二：直接运行Snakemake（更灵活）

```bash
# 确保在脚本目录下
cd /path/to/01_Bulk_RNAseq_script

# 试运行，检查逻辑
snakemake -n

# 正式运行
snakemake --cores 36 --printshellcmds
```

### 7.3 监控运行状态

```bash
# 查看运行日志
tail -f pipeline.log

# 查看进程
ps -ef | grep snakemake
```

---

## 8. 分析方法

1. **质量控制**：使用FastQC进行初始评估，Cutadapt修剪接头和低质量碱基（Q20），修剪后长度<75bp的reads被丢弃
2. **比对**：使用HISAT2将reads比对到GRCm39基因组，添加`--dta`标记便于下游组装
3. **转录本组装**：StringTie以参考GTF为指导重构转录本，识别新异构体
4. **定量**：HTSeq-count基于外显子特征计算基因水平的原始counts，用于下游差异表达分析

---

## 9. 下游分析

### 9.1 创建样本分组信息

**必须创建样本分组文件 `sample_group.csv`，确保样本名称与原始数据完全一致：**

```csv
sample,group
WT-testis-1,WT
WT-testis-2,WT
C2-KO-testis-1,C2_KO
C2-KO-testis-2,C2_KO
C6-KO-testis-1,C6_KO
C6-KO-testis-2,C6_KO
C26-DKO-testis-1,C26_DKO
C26-DKO-testis-2,C26_DKO
```

### 9.2 运行下游分析脚本

```bash
# 确保在脚本目录下
cd /path/to/01_Bulk_RNAseq_script

# 使用默认R环境运行
nohup Rscript downstream_analysis.R > downstream.log 2>&1 &

# 或指定特定R环境运行（例如使用/opt/R/4.4.1）
nohup /opt/R/4.4.1/bin/Rscript downstream_analysis.R > downstream.log 2>&1 &
```

### 9.3 下游分析输出目录

运行完成后，将在脚本目录下生成以下结果：

* **07_AdvancedQC**：相关性热图、基因表达热图、PCA图
* **08_DE**：差异表达基因列表（CSV格式）和火山图
* **09_GO**：GO富集分析结果（表格和可视化图）
* **10_KEGG**：KEGG通路富集分析结果（表格和多类型可视化图）
* **12_PCA**：PCA分析图

---

## 10. 注意事项

1. **路径问题**：所有路径建议使用绝对路径，避免相对路径导致的错误
2. **样本命名**：样本名称一旦确定，在整个分析流程中保持一致
3. **内存需求**：HISAT2比对和StringTie组装需要较大内存，建议服务器配置≥64GB
4. **检查点**：每个步骤完成后检查日志文件，确保无错误
5. **分组文件**：`sample_group.csv` 必须与`00_Rawdata`中的样本名称完全一致

---

## 11. 常见问题排查

| 问题 | 可能原因 | 解决方案 |
|------|----------|----------|
| 找不到输入文件 | 路径错误或文件命名不规范 | 检查`00_Rawdata`目录和文件命名 |
| 比对率低 | 接头未正确去除或基因组索引问题 | 检查config.yml中的adapter序列 |
| 下游分析报错 | 分组文件与样本名不匹配 | 检查`sample_group.csv` |
| 内存不足 | 数据量过大 | 减少并行任务数或增加内存 |

---

## 12. 引用和致谢

如使用本流程，请引用：
- Snakemake: Mölder et al. (2021)
- HISAT2: Kim et al. (2019)
- StringTie: Pertea et al. (2015)
- HTSeq: Anders et al. (2015)