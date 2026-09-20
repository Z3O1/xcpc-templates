#!/usr/bin/env bash
# 编译并运行 templates/ 下的 *.check.cpp 自测(类似 skip2004 的 check)
#
# 用法:
#   ./check.sh                          # 跑全部(慢:66 个 check 要几分钟,只在收尾/推送前跑)
#   ./check.sh 高斯整数                  # 只跑路径里含该子串的 check(可给多个子串,OR 匹配)
#   ./check.sh -x 高斯整数               # 精确匹配单个 check(给模板名/文件名/路径都行,忽略后缀)
#   ./check.sh -f templates/数论/高斯整数.check.cpp   # 直接指定 check 文件(可给多个)
#   ./check.sh -l 数论                   # 只列出会跑哪些,不编译不运行
#   ./check.sh -v -x 高斯整数            # 顺便打印这个 check 的断言输出
#   ./check.sh -j 4                     # 最多 4 个 check 并行(默认串行;并行只影响耗时,输出顺序不变)
#   ./check.sh -h                       # 帮助
#
# 改一个模板时只跑它那一个:./check.sh -x <模板名> —— 别每次都跑全局,全量编译+运行太慢。
set -u
cd "$(dirname "$0")"

timeout_s="${XCPC_CHECK_TIMEOUT:-60}"   # 单个 check 的墙钟上限(卡死=FAIL,不再挂住整轮)

usage() {
    cat <<'EOF'
用法:
  ./check.sh                          # 跑全部(慢:全部 check 要几分钟,只在收尾/推送前跑)
  ./check.sh 高斯整数                  # 只跑路径里含该子串的 check(可给多个子串,OR 匹配)
  ./check.sh -x 高斯整数               # 精确匹配单个 check(给模板名/文件名/路径都行,忽略后缀)
  ./check.sh -f templates/数论/高斯整数.check.cpp   # 直接指定 check 文件(可给多个)
  ./check.sh -l 数论                   # 只列出会跑哪些,不编译不运行
  ./check.sh -v -x 高斯整数            # 顺便打印这个 check 的断言输出
  ./check.sh -j 4                     # 最多 4 个 check 并行(默认串行;并行只影响耗时,输出顺序不变)
  ./check.sh -h                       # 本帮助

改一个模板时只跑它那一个:./check.sh -x <模板名> —— 别每次都跑全局,全量编译+运行太慢。
环境变量 XCPC_CHECK_TIMEOUT 可改单个 check 的墙钟上限(默认 60s,超时按 FAIL 报)。
EOF
}

jobs=1
verbose=0
list_only=0
exact=0
filters=()   # 子串过滤(-x 时是精确名字)
files=()     # -f 直接给的 check 文件

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        -v|--verbose) verbose=1; shift ;;
        -l|--list) list_only=1; shift ;;
        -x|--exact) exact=1; shift ;;
        -j|--jobs)
            jobs="${2:-}"
            [[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "-j 需要一个正整数(收到:'${jobs}')" >&2; exit 2; }
            shift 2 ;;
        -f|--file)
            [[ $# -ge 2 ]] || { echo "-f 后面要跟 check 文件" >&2; exit 2; }
            files+=("$2"); shift 2 ;;
        -*) echo "未知选项: $1" >&2; echo; usage >&2; exit 2 ;;
        *) filters+=("$1"); shift ;;
    esac
done

mapfile -t all < <(find templates -name '*.check.cpp' | sort)
if [[ ${#all[@]} -eq 0 ]]; then
    echo "没找到任何 *.check.cpp(check 文件尚未创建?)"
    exit 1
fi

# 按名字精确解析成 check 路径:接受 模板名 / xxx.check.cpp / 章/模板名 等写法
resolve_exact() {
    local name=$1 base sb src
    if [[ -f "$name" ]]; then printf '%s\n' "$name"; return 0; fi
    base=$(basename "$name"); base=${base%.check.cpp}; base=${base%.cpp}
    local hit=0
    for src in "${all[@]}"; do
        sb=$(basename "$src"); sb=${sb%.check.cpp}
        if [[ "$sb" == "$base" ]]; then printf '%s\n' "$src"; hit=1; fi
    done
    [[ $hit -eq 1 ]]
}

# —— 选定要跑的 check 列表(保持 find 的顺序,去重) ——
declare -A seen=()
checks=()
add_check() {
    [[ -n "${seen[$1]:-}" ]] && return 0
    seen[$1]=1; checks+=("$1")
}

if [[ ${#files[@]} -gt 0 ]]; then
    for f in "${files[@]}"; do
        matched=$(resolve_exact "$f") || matched=""
        if [[ -z "$matched" ]]; then echo "找不到 check:'$f'(用 -l 列出现有的)" >&2; exit 1; fi
        while IFS= read -r line; do add_check "$line"; done <<< "$matched"
    done
elif [[ ${#filters[@]} -gt 0 && $exact -eq 1 ]]; then
    for f in "${filters[@]}"; do
        matched=$(resolve_exact "$f") || matched=""
        if [[ -z "$matched" ]]; then echo "找不到 check:'$f'(用 -l 列出现有的)" >&2; exit 1; fi
        while IFS= read -r line; do add_check "$line"; done <<< "$matched"
    done
elif [[ ${#filters[@]} -gt 0 ]]; then
    for f in "${filters[@]}"; do
        while IFS= read -r line; do add_check "$line"; done < <(printf '%s\n' "${all[@]}" | grep -F -- "$f" || true)
    done
else
    for src in "${all[@]}"; do add_check "$src"; done
fi

if [[ ${#checks[@]} -eq 0 ]]; then
    echo "没有匹配的 *.check.cpp(用 -l 列出现有的)"
    exit 1
fi

if [[ $list_only -eq 1 ]]; then
    printf '%s\n' "${checks[@]}"
    echo "共 ${#checks[@]} 个(全部 ${#all[@]} 个)"
    exit 0
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

now() { if [[ -n "${EPOCHREALTIME:-}" ]]; then echo "$EPOCHREALTIME"; else date +%s.%N; fi; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.1f", b - a}'; }

# run_check <src> <log>:编译并运行一个 check,把结果写进 <log>
# 退出码:0 = ok,1 = FAIL(含超时),2 = 编译失败
run_check() {
    local src=$1 log=$2 bin t0 t1 rc r line body=''
    bin="$tmpdir/$(echo "$src" | tr '/.' '__')"
    t0=$(now)
    if ! err=$(g++ -std=c++17 -O2 -o "$bin" "$src" 2>&1); then
        t1=$(now)
        line="COMPILE FAIL  $src"
        body=$(printf '%s\n' "$err" | head -15)
        rc=2
    elif out=$(cd "$(dirname "$src")" && timeout "$timeout_s" "$bin" 2>&1); then
        t1=$(now)
        line="ok    $src"
        rc=0
        [[ $verbose -eq 1 ]] && body="$out"
    else
        r=$?
        t1=$(now)
        if [[ $r -eq 124 ]]; then
            line="FAIL  $src  (超时 ${timeout_s}s —— 很可能是 check 触发了模板的死循环用例)"
        else
            line="FAIL  $src"
        fi
        body=$(printf '%s\n' "$out" | tail -20)
        rc=1
    fi
    {
        printf '%s  (%ss)\n' "$line" "$(elapsed "$t0" "$t1")"
        [[ -n "$body" ]] && printf '%s\n' "$body" | sed 's/^/      /'
    } > "$log"
    return $rc
}

declare -a status=()
n=${#checks[@]}
echo "跑 ${n} 个 check(共 ${#all[@]} 个;单个上限 ${timeout_s}s;并行度 ${jobs})"
if [[ $jobs -le 1 ]]; then
    for ((k = 0; k < n; ++k)); do
        if run_check "${checks[$k]}" "$tmpdir/log_$k"; then status[k]=0; else status[k]=$?; fi
        cat "$tmpdir/log_$k"          # 串行时跑完一个就打印一个(长任务能看进度)
    done
else
    for ((i = 0; i < n; i += jobs)); do
        pids=()
        for ((k = i; k < i + jobs && k < n; ++k)); do
            run_check "${checks[$k]}" "$tmpdir/log_$k" &
            pids+=("$!")
        done
        for ((k = i, j = 0; k < i + jobs && k < n; ++k, ++j)); do
            if wait "${pids[$j]}"; then status[k]=0; else status[k]=$?; fi
        done
    done
    for ((k = 0; k < n; ++k)); do cat "$tmpdir/log_$k"; done
fi

pass=0; fail=0; failed=()
for ((k = 0; k < n; ++k)); do
    if [[ ${status[k]:-1} -eq 0 ]]; then ((pass++)); else ((fail++)); failed+=("${checks[$k]}"); fi
done

echo
echo "== $pass passed, $fail failed =="
for f in "${failed[@]:-}"; do [[ -n "$f" ]] && echo "   $f"; done
[[ $n -lt ${#all[@]} ]] && echo "(只跑了 $n / ${#all[@]} 个;跑全部:./check.sh)"
exit $((fail > 0))
