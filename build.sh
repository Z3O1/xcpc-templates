#!/usr/bin/env bash
# 一键生成 sections.typ 并编译 xcpc.pdf;--watch 监听已有模板的内容变化。
set -euo pipefail
cd "$(dirname "$0")"

watch=0
case "${1:-}" in
    -h|--help)
        printf '用法: ./build.sh [--watch] [Typst 选项 ...]\n加删模板或改 manifest 后需重新启动构建。TYPST 可指定编译器。\n'
        exit 0 ;;
    --watch) watch=1; shift ;;
esac

# /snap/bin/typst 是 snapd 启动器:在 PID 隔离沙箱中无法经 host D-Bus 建 scope。
# 优先普通二进制,遇到 snap 启动器则直接用载荷;不要退回会失败的启动器。
if [[ -n "${TYPST:-}" ]]; then
    command -v "$TYPST" >/dev/null || { echo "找不到指定的 typst: $TYPST" >&2; exit 1; }
elif candidate=$(command -v typst) && [[ "$(readlink -f "$candidate")" != /usr/bin/snap ]]; then
    TYPST=$candidate
elif [[ -x /snap/typst/current/bin/typst ]]; then
    TYPST=/snap/typst/current/bin/typst
else
    echo '找不到 typst 二进制(PATH 或 /snap/typst/current/bin/typst)' >&2
    exit 1
fi

# snap 载荷需要可写运行目录且看不到 /tmp,调试与编译产物都放项目 tmp/。
mkdir -p tmp/xdg
export XDG_RUNTIME_DIR="$PWD/tmp/xdg"
python3 gen.py
echo '--- typst ---'
echo "using $TYPST ($("$TYPST" --version))"
if [[ $watch -eq 1 ]]; then
    "$TYPST" watch xcpc.typ xcpc.pdf "$@"
else
    "$TYPST" compile xcpc.typ xcpc.pdf "$@"
    echo "OK -> $PWD/xcpc.pdf"
fi
