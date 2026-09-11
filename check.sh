#!/usr/bin/env bash
# 编译并运行 templates/ 下所有 *.check.cpp 自测(类似 skip2004 的 check)
# 用法: ./check.sh            # 跑全部
#       ./check.sh 数论        # 只跑名字里含"数论"的
#       ./check.sh -v 数论     # 顺便打印每个 check 的输出
#
# 注意:templates/计算几何/geo.check.cpp 目前**故意失败** —— 它复现了 geo.cpp:33-34
# 的两个真 bug(`-=` 与 `/=` 都写成了 `x = x + y`)。修好那两行它就自动转绿。
set -u
cd "$(dirname "$0")"

timeout_s="${XCPC_CHECK_TIMEOUT:-60}"   # 单个 check 的墙钟上限(卡死=FAIL,不再挂住整轮)
verbose=0
[[ "${1:-}" == "-v" ]] && verbose=1 && shift
filter="${1:-}"

mapfile -t checks < <(find templates -name '*.check.cpp' | sort)
[[ -n "$filter" ]] && mapfile -t checks < <(printf '%s\n' "${checks[@]}" | grep -- "$filter")

if [[ ${#checks[@]} -eq 0 ]]; then
    echo "没有匹配的 *.check.cpp(check 文件尚未创建?)"
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

pass=0; fail=0; failed=()
for src in "${checks[@]}"; do
    bin="$tmpdir/$(echo "$src" | tr '/.' '__')"
    if ! err=$(g++ -std=c++17 -O2 -o "$bin" "$src" 2>&1); then
        printf 'COMPILE FAIL  %s\n' "$src"
        printf '%s\n' "$err" | head -15
        failed+=("$src"); ((fail++)); continue
    fi
    if out=$(cd "$(dirname "$src")" && timeout "$timeout_s" "$bin" 2>&1); then
        printf 'ok    %s\n' "$src"
        ((pass++))
        [[ $verbose -eq 1 ]] && printf '%s\n' "$out" | sed 's/^/      /'
    else
        rc=$?
        if [[ $rc -eq 124 ]]; then
            printf 'FAIL  %s  (超时 %ss —— 很可能是 check 触发了模板的死循环用例)\n' "$src" "$timeout_s"
        else
            printf 'FAIL  %s\n' "$src"
        fi
        printf '%s\n' "$out" | tail -20 | sed 's/^/      /'
        failed+=("$src"); ((fail++))
    fi
done

echo
echo "== $pass passed, $fail failed =="
for f in "${failed[@]:-}"; do [[ -n "$f" ]] && echo "   $f"; done
exit $((fail > 0))
