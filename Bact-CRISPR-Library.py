#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Bact-CRISPR-Library.py

这是一个通用调用脚本，用于整合并调用不同的Cas9基因编辑文库设计模式，专为细菌（如大肠杆菌）设计。
它不会修改原始脚本，而是通过子进程运行相应的设计脚本，并根据模式传递参数。

支持的模式:
- Knockdown_Cas9: Cas9-mediated gene knockdown design.
- PromoterChange_Cas9: Promoter replacement design (Cas9).
- Cfusion_Cas9: C-terminal fusion design (Cas9).
- Knockout_Cas9: Gene knockout design (Cas9).
- Knockout_CASTs: Gene knockout design (CASTs).
- PromoterChange_CASTs: Promoter replacement design (CASTs).

使用示例:
python Bact-CRISPR-Library.py Knockdown_Cas9 --input_fna genome.fna --input_gff annotation.gff --output results.csv --synthesis_template template.txt [其他参数...]

运行 'python Bact-CRISPR-Library.py -h' 显示此帮助信息。
运行 'python Bact-CRISPR-Library.py <mode> -h' 显示指定模式的详细参数。
"""

import argparse
import subprocess
import sys
import logging

# --- 日志设置 ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def build_cmd(args, script_file, sub_mode=None):
    cmd = ["python3", script_file]
    arg_dict = vars(args)
    if sub_mode:
        cmd.extend(["--mode", sub_mode])
    for key, value in arg_dict.items():
        if key == 'mode':
            continue
        if value is not None:
            if isinstance(value, bool) and value:
                cmd.append(f"--{key}")
            elif isinstance(value, list):
                cmd.append(f"--{key}")
                cmd.extend(value)
            else:
                cmd.append(f"--{key}")
                cmd.append(str(value))
    return cmd

def main():
    parser = argparse.ArgumentParser(description="Bact-CRISPR-Library: A tool for designing CRISPR libraries in bacteria using Cas9 or CASTs.", 
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    # 共同参数 (使用 add_help=False 以避免重复帮助)
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--input_fna", required=True, help="输入的基因组FASTA文件。")
    common_parser.add_argument("--input_gff", required=True, help="输入的GFF3注释文件。")
    common_parser.add_argument("--output", required=True, help="输出的CSV文件路径。")
    common_parser.add_argument("--genome_type", type=str, default="linear", choices=['linear', 'circle'], help="基因组类型 (默认: linear)。")
    common_parser.add_argument("--sgRNA_num", type=int, default=1, help="每个基因的设计数量 (默认: 1)。")
    common_parser.add_argument("--barcode_len", type=int, default=10, help="条形码(barcode)的长度 (默认: 10)。")
    common_parser.add_argument("--restriction_site", type=str, nargs='+', help="需要避免的限制性酶切位点。")
    common_parser.add_argument("--synthesis_template", type=str, required=True, help="定义最终oligo结构的模板文件路径。")
    common_parser.add_argument("--cloning_site", type=str, help="将被豁免检查的克隆位点序列。")
    common_parser.add_argument("--relax", action="store_true", help="不考虑是否破坏临近基因的CDS区。")
    
    subparsers = parser.add_subparsers(dest="mode", required=True, help="设计模式: Knockdown_Cas9 (基因下调), PromoterChange_Cas9 (启动子替换), Cfusion_Cas9 (C端融合), Knockout_Cas9 (基因敲除) (Cas9-based), Knockout_CASTs, PromoterChange_CASTs (CASTs-based)。")
    
    # Knockdown_Cas9 子命令
    knockdown_parser = subparsers.add_parser("Knockdown_Cas9", parents=[common_parser], help="基因下调 (Knockdown) 模式 (Cas9)。")
    knockdown_parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    knockdown_parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度范围 '首选:最小'。")
    knockdown_parser.add_argument("--sgrna_upstream_range", type=str, default="20:60", help="sgRNA距离起始密码子上游的搜索范围 '最小:最大' bp (默认: 20:60)。")
    knockdown_parser.add_argument("--rha_min_dist_to_atg", type=int, default=15, help="插入点与ATG的最小距离 (默认: 15)。")
    
    # PromoterChange_Cas9 子命令
    promoter_parser = subparsers.add_parser("PromoterChange_Cas9", parents=[common_parser], help="启动子替换 (Promoter Change) 模式 (Cas9)。")
    promoter_parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    promoter_parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度范围 '首选:最小'。")
    promoter_parser.add_argument("--promoter_search_size", type=int, default=150, help="在ATG上游搜索sgRNA的区域大小 (默认: 150)。")
    
    # Cfusion_Cas9 子命令
    cfusion_parser = subparsers.add_parser("Cfusion_Cas9", parents=[common_parser], help="C端融合 (C-Fusion) 模式 (Cas9)。")
    cfusion_parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    cfusion_parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度范围 '首选:最小'。")
    cfusion_parser.add_argument("--sgrna_search_range", type=int, default=50, help="在终止密码子下游搜索sgRNA的初始范围 (默认: 50)。")
    cfusion_parser.add_argument("--max_search_expansion", type=int, default=300, help="sgRNA搜索的最大下游扩展距离 (默认: 300)。")
    cfusion_parser.add_argument("--strict", action="store_true", help="开启严格模式，避免破坏邻近基因。")
    
    # Knockout_Cas9 子命令
    knockout_parser = subparsers.add_parser("Knockout_Cas9", parents=[common_parser], help="基因敲除 (Knockout) 模式 (Cas9)。")
    knockout_parser.add_argument("--pam", type=str, default="NGG", help="PAM序列 (默认: 'NGG')。")
    knockout_parser.add_argument("--HR_len", type=str, default="50:40", help="同源臂长度范围 '首选:最小'。")
    knockout_parser.add_argument("--del_length", type=str, default="50:100", help="删除长度范围 '最小:最大' bp (默认: 50:100)。")
    knockout_parser.add_argument("--ko_search_range", type=str, default="5:80", help="在CDS中搜索sgRNA的百分比范围 (默认: 5:80)。")
    knockout_parser.add_argument("--promoter_region", type=str, default="50:150", help="在严格模式下，邻近基因启动子受保护的区域，格式'最小:最大'bp (默认: 50:150)。")
    knockout_parser.add_argument("--strict", action="store_true", help="开启严格模式，避免破坏邻近基因。")
    
    # Knockout_CASTs 子命令
    casts_ko_parser = subparsers.add_parser("Knockout_CASTs", parents=[common_parser], help="基因敲除 (Knockout) 模式 (CASTs)。")
    casts_ko_parser.add_argument("--target_cds_range", type=str, default="5:80", help="在CDS中搜索sgRNA的百分比范围, 格式 '最小:最大' (默认: 5:80)。")
    
    # PromoterChange_CASTs 子命令
    casts_pc_parser = subparsers.add_parser("PromoterChange_CASTs", parents=[common_parser], help="启动子替换 (Promoter Change) 模式 (CASTs)。")
    casts_pc_parser.add_argument("--insertion_range_promoter", type=str, default="25:50", help="插入点距离起始密码子上游的距离范围, 格式 '最小:最大' bp (默认: 25:50)。")
    
    args = parser.parse_args()
    
    # 根据模式选择脚本文件
    script_map = {
        "Knockdown_Cas9": "Cas9_knockdown_designer_v1.py",
        "PromoterChange_Cas9": "Cas9_PromoterChange_designer_v2.py",
        "Cfusion_Cas9": "Cas9_Cfusion_designer_v1.py",
        "Knockout_Cas9": "Cas9_knockout_designer_v1.py",
        "Knockout_CASTs": "CASTs_designer_v3.py",
        "PromoterChange_CASTs": "CASTs_designer_v3.py"
    }
    sub_mode_map = {
        "Knockout_CASTs": "Knockout_CASTs",
        "PromoterChange_CASTs": "PromoterChange_CASTs"
    }
    
    script_file = script_map[args.mode]
    sub_mode = sub_mode_map.get(args.mode)
    
    # 构建并运行命令
    cmd = build_cmd(args, script_file, sub_mode)
    logging.info(f"正在调用模式: {args.mode}，脚本: {script_file}")
    logging.info(f"命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            logging.error(result.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(f"调用失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
