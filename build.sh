#!/usr/bin/env bash
# 一键生成: 遍历 templates/ 生成 sections.typ, 再编译出 xcpc.pdf
# 用法: ./build.sh          # 单次编译
#       ./build.sh --watch  # 监听模式(改模板文件即时重编译)
set -e
cd "$(dirname "$0")"

python3 gen.py
echo '--- typst ---'
if [[ "$1" == "--watch" ]]; then
  typst watch xcpc.typ
else
  typst compile xcpc.typ "$@"
  echo "OK -> $(pwd)/xcpc.pdf"
fi
