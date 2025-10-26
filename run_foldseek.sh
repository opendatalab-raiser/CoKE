#!/bin/bash
# Foldseek 便捷运行脚本
# 自动设置环境和库路径

# 激活 conda 环境
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bioanalysis

# 设置数据库路径
# export FOLDSEEK_DB="/zhuangkai2/projects/lost_in_tokenization/foldseek_db/sp"

# 使用 LD_PRELOAD 修复库依赖问题
export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6.0.34"

# 运行 foldseek 命令
foldseek "$@"

