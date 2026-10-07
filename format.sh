#!/usr/bin/env bash
# 默认只整理模板本体与 base header;自测可显式传路径,不触碰参考代码。
set -euo pipefail
cd "$(dirname "$0")"

check=0
case "${1:-}" in
    -h|--help)
        printf '用法: ./format.sh [--check] [源码路径 ...]\n默认处理全部模板本体与 _check_base.hpp;需要 clang-format 14+。\n'
        exit 0 ;;
    --check) check=1; shift ;;
esac
command -v clang-format >/dev/null || { echo '找不到 clang-format(需要 14+)' >&2; exit 1; }
if [[ $# -gt 0 ]]; then
    sources=("$@")
else
    mapfile -d '' -t sources < <(find templates -type f \( -name '*.cpp' -o -name '*.hpp' \) ! -name '*.check.cpp' -print0 | LC_ALL=C sort -z)
fi
if [[ $check -eq 1 ]]; then
    clang-format --style=file --dry-run --Werror "${sources[@]}"
else
    clang-format --style=file -i "${sources[@]}"
fi
