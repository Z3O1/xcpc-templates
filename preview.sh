#!/usr/bin/env bash
# 只编译样式小样,不生成章节、不覆盖 xcpc.pdf; 默认监听,产物统一放 tmp/。
set -euo pipefail
cd "$(dirname "$0")"

mode=watch
case "${1:-}" in
    -h|--help)
        printf '用法: ./preview.sh [--once] [Typst 选项 ...]\n默认监听 style.typ 和小样,输出 tmp/style-preview.pdf。\n--once 只编译一次; --open 打开 PDF; TYPST 可指定编译器。\n'
        exit 0 ;;
    --once) mode=compile; shift ;;
esac

# 与 build.sh 一致:优先普通二进制,避开沙箱内不可用的 snap 启动器。
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

mkdir -p tmp/xdg
export XDG_RUNTIME_DIR="$PWD/tmp/xdg"
echo "using $TYPST ($("$TYPST" --version))"
echo "预览输出 -> $PWD/tmp/style-preview.pdf"
exec "$TYPST" "$mode" --root . tests/style-preview.typ tmp/style-preview.pdf "$@"
