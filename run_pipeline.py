import os
import subprocess
import sys
import yaml
import logging

# 日志配置 (Python 终端显示)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(config_path):
    """确保所有定义的目录无论是否使用都已创建"""
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)
    dirs = cfg.get('directories', {})
    for name, path in dirs.items():
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            logging.info(f"创建目录: {path}")


def run_snakemake():
    """执行 Snakemake 流程并记录详细日志"""
    logging.info(">>> 启动 C2C6_Project 自动化分析 (Snakemake 模式) <<<")

    # 构建 Snakemake 命令
    # 注意：移除了旧版本可能不支持的 --reason
    # 增加了 --rerun-incomplete (如果上次运行中断，自动处理未完成文件)
    cmd = [
        "snakemake",
        "-s", "snakemake.smk",
        "--cores", "36",
        "--printshellcmds",
        "--keep-going",
        "--rerun-incomplete"
    ]

    # 定义日志文件路径
    log_file = "workflow_snakemake.log"
    logging.info(f"详细运行日志将写入: {log_file}")

    try:
        # 使用 subprocess 将输出同时记录到文件
        with open(log_file, "w") as f:
            # stdbuf -oL 确保实时写入
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            # 同时在终端打印并写入文件
            for line in process.stdout:
                sys.stdout.write(line)  # 终端显示
                f.write(line)  # 写入文件

            process.wait()

            if process.returncode == 0:
                logging.info(">>> 分析任务圆满完成！ <<<")
            else:
                logging.error(f">>> 流程中断，Snakemake 返回码: {process.returncode}，请检查 {log_file} <<<")
                sys.exit(1)

    except Exception as e:
        logging.error(f">>> 脚本执行异常: {e} <<<")
        sys.exit(1)


if __name__ == "__main__":
    if not os.path.exists("config.yml"):
        logging.error("错误: 找不到 config.yml")
        sys.exit(1)

    setup_directories("config.yml")
    run_snakemake()
