#!/usr/bin/env bash
# 一键生成: 遍历 templates/ 生成 sections.typ, 再编译出 xcpc.pdf
# 用法: ./build.sh          # 单次编译
#       ./build.sh --watch  # 监听模式(改模板文件即时重编译)
set -e
cd "$(dirname "$0")"

python3 gen.py
echo '--- typst ---'

# 本机 typst 是 snap 版。注意 /snap/bin/typst 只是 snapd 启动器(→ /usr/bin/snap),
# 它在 DSH 沙箱里必然启动失败: snapd 要先经 D-Bus 向 systemd 申请 transient scope,
# 而沙箱的 PID 1 是 `bwrap --unshare-pid`, host 上的 systemd/dbus 看不到命名空间里的
# PID → `cannot create transient scope: ... UnixProcessIdUnknown`(PID 1 为 bwrap 时
# 报错可能变成 "job ... finished with result failed")。这是 PID 命名空间限制,
# 跟文件权限无关, 放宽沙箱也修不好; 直接跑 snap 载荷里的真二进制即可绕过 snapd。
# 顺带: typst 读不到 /tmp(snap 私有 tmpfs), 调试小样要写在项目目录里。
if command -v typst >/dev/null 2>&1; then
  TYPST=typst
elif [[ -x /snap/typst/current/bin/typst ]]; then
  TYPST=/snap/typst/current/bin/typst
elif [[ -x /snap/bin/typst ]]; then
  TYPST=/snap/bin/typst
else
  echo "找不到 typst(PATH 里没有, 也没有 /snap/typst/current/bin/typst)" >&2
  exit 1
fi

# snap 版 typst 需要一个可写的 XDG_RUNTIME_DIR, 否则报
# `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`
mkdir -p tmp/xdg
export XDG_RUNTIME_DIR="$PWD/tmp/xdg"
echo "using $TYPST ($("$TYPST" --version 2>/dev/null | head -1))"

if [[ "$1" == "--watch" ]]; then
  "$TYPST" watch xcpc.typ
else
  "$TYPST" compile xcpc.typ "$@"
  echo "OK -> $(pwd)/xcpc.pdf"
fi
