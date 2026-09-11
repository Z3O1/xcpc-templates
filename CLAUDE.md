# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

算法竞赛模板集:XCPC 板子库,`templates/` 下按章目录放代码模板,一键生成并编译成 `xcpc.pdf`(横向 A4、双栏正文、右上+右下大页码)。

## 常用命令

```bash
./build.sh                   # python3 gen.py && typst compile xcpc.typ → xcpc.pdf
./build.sh --watch           # 改 templates/ 后 typst 自动重编译
python3 gen.py               # 只看生成、不编译
```

- 改模板后跑 `./build.sh`,Typst 直接报 `error: ... templates/xxx.typ:行号` 加文件路径,照报错修即可
- **本环境 typst 是 snap 版且不在 PATH**:用 `/snap/bin/typst`,并且要先建可写的运行目录(否则报 `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`):`mkdir -p tmp/xdg && XDG_RUNTIME_DIR=$PWD/tmp/xdg /snap/bin/typst compile xcpc.typ`;`./build.sh` 会因找不到 `typst` 在最后一步失败(gen.py 部分已跑完)
- **标题必须 Typst 安全**:标题会原样进 `sections.typ` 的 `== 标题`,带 `(`/`)` 会被当函数调用报 `unclosed delimiter`,带 `*`/反引号会被当标记——标题写短名词,细节放注释
- **改了首行注释但标题没变**:`gen.py` 只给新条目从注释取标题,已有条目以 manifest 为准 → 想换标题就先从 `.manifest.json` 里摘掉该条再 `python3 gen.py`
- 调试单点 Typst 小样:写个 `_t.typ` 放到**项目目录**(snap 版 typst 读不到 /tmp!),`typst compile _t.typ`,用完删;要用 `#panic("...")` 把值打到 stderr——`pdftotext` 的页面折行会伪造换行字符,不能用来调试字符串
- 没有任何测试/CI;编译通过 + 抽查 PDF 页(如 `pdftoppm -png -f N -l N xcpc.pdf /tmp/x`)就是验证

## 架构:生成链路

```
templates/<章目录>/<模板>.cpp|sh|py  ─┐
templates/<章目录>/<模板>.typ(介绍)  ─┼─→ gen.py ─→ sections.typ ─→ #include 进 xcpc.typ ─→ xcpc.pdf
```

- Typst 的 `read()` 不能遍历目录,所以"扫描"全部由 `gen.py` 承担:扫磁盘 → 与 `.manifest.json` 对账(记忆顺序/标题,增/删/恢复都是增量) → 输出 `sections.typ` → 确保 `xcpc.typ` 正文为 `#include "sections.typ"` 模式
- `xcpc.typ.bak` 是迁移前原始内联版(39499 B,勿删):首次无 manifest 时用来还原章节顺序与标题
- `sections.typ` 由 gen.py 生成,**不要手改**
- **外部依赖**:gen.py 输出的 `sections.typ` 首行 `#import "@preview/zebraw:0.6.3"`(代码高亮包,需要该 Typst 包已装好——联网拉取或本地缓存;编译报包找不到就先检查 `~/.cache/typst`)

## 模板约定(模板文件上需要遵守的规则)

- **章节**:`templates/` 下每个目录是一章,顺序由 manifest 决定(新目录追加到最后)
- **标题**:取代码文件首行注释 `// xxx(): 章节 · 标题 (括号内的备注会被去掉)` 的「· 标题」部分;无注释则用文件名。标题存 manifest,与首行注释解耦
- **隐藏**:文件头注释区含 `// hide` 或 `// 隐藏` → 整块不进 PDF,去掉标记即恢复(`sam.cpp` 当前处于此状态)
- **介绍(推荐)**:同名 `.typ` 文件(如 `ntt.cpp` 配 `ntt.typ`)渲染在该代码前面,不进 manifest、不占条目;**判断规则是"存在同基底名非 `.typ` 文件"**——没有同名代码的 `.typ` 是独立可渲染模板
- **代码文件头注释规范**:第一空行之前的内容(头部注释)不渲染,`readcode` 会跳过 `//` 开头的行与空行,其余全部当代码——无注释无空行(纯代码文件)也能正常渲染
- 介绍 `.typ` 里不要再写与 `== 标题` 同级的标题;可用 `*文本*`/`_文本_` 加粗、``` `代码` ``` 等 Typst 语法

## 代码风格(house 宏,写模板必看)

模板是**竞赛代码片段**,不是可单独编译的独立单元:文件内不定义宏/类型,而是沿用竞赛 base header 提供的一套命名。书里的代码在 PDF 中是"参考实现",阅读者平时直接背这套宏,所以新写模板必须沿用同一套,别自造风格:

- `For(i, l, r)` / `rFor(i, r, l)` / `ForD` 代替 `for`;`vect<T>` 代替 `vector<T>`
- `ll` = `long long`、`db` = 浮点(常用 `long double`)、`mint` = 模数类(带 `.inv()`)、`poly` = 多项式/`vector<mint>`、`ksm` = 快速幂
- **`ksm` 的模数限制**:base header 的 `ksm` 内部按 `ll` 相乘,模数超过 2^32 时 `a * a` 会溢出算错 → 需要大模数(如 `u64` 素性检验、`p > 4e9` 的离散对数)时**自带模乘**,不要复用 `ksm`(见 `templates/数论/Miller-Rabin.cpp` 的 `MR::mul`)
- **`vect<T>` 没有 `push_back`**:它是带 `+=` 的自定义容器(`p += x` 才是追加);标准容器用 `std::vector`。写模板要按"读者自己 header 里那套"来,别混用
- 计算几何统一 `struct p2`(点),线段/直线用 `seg`;`eps`/`cmp`/`sgn`/`cross`/`det` 常规
- 模板顶部的 `// 标题 · 标题 (备注)` 首行注释是选标题用的,不用改内容

这类宏的具体定义不在本仓库,而在用户的 contbase header 里;拿不准时参考 `skip2004-ICPC-Templates/`(其 `contents/*/all.cpp` 是完整的 base header)与 `~/0/Code/` 里既有模板的写法。

## 移植来源(补模板时优先搬运,别从零写)

给 `templates/` 补新模板时,先找现成轮子,按 house 风格改写后放进去,再 `./build.sh` 验证:

- **`~/0/Code/`** —— 用户自己的板子(接近 house 风格),含 `一般图最大匹配.cpp`/`上下界.cpp`/`Hopcroft–Karp.cpp`/`FHQ.cpp`/`FFT.cpp`/`poly 1..3.cpp`/`Bostan-Mori.cpp`/`杜教筛.cpp`/`Gause.cpp`/`det.cpp`/`Matrix.cpp`/`BigInt.cpp`/`Hash.cpp`/`Geo.cpp` 等
- **`~/0/Lib/atcoder/`** —— 官方 atcoder 库(现成 `maxflow`/`mincostflow`/`scc`/`twosat`/`segtree`/`lazysegtree`/`convolution`/`math`/`string` 等,可作为"该算法该怎么做"的权威参考)
- **`skip2004-ICPC-Templates/contents/<分类>/`** —— 第三方**完整**模板集,覆盖图论/数学/数论/数据结构/字符串/几何/杂项 173 个文件(含 `2-sat`/`AC自动机`/`BM`/`CRT`/`D-MST`/`Delaunay`/`FFT`/`Fiduccia`/`GF`/`HK`/`KM`/`LP`/`bostan_mori`/`cactus`/`convex` 等)。**改动只写进本仓库 `templates/`**,该目录本身保持原样不动

## Typst 陷阱(本环境实测)

## Typst 陷阱(本环境实测)

- `$...$` 数学模式里**不要用带反斜杠的多字母符号**:`\log`、`\sqrt`、`\alpha` 等一律报 `unknown variable`;多字母标识符(`len`、`mcf`、`que`)也会被拆成变量序列。替代方案:复杂度等写反引号文本(`` `O(n log n)` ``)或 Unicode 文本;数学里只留单字母和基本运算
- 强调用单 `*` 或 `_`,**双星 `**` 无效**
- 反引号必须成对(不成对报 `unclosed raw text`,报错跨度可能跨行,很迷惑);`**`/`_` 同理
- Typst 0.15:`().join("\n")` 返回 `none`(空数组 join),需要显式 `if ... else { "" }`
- `#set page(header: ...)`/`footer` 里用 `#place(right + top/bottom, ...)` 放页码,`#text(size: 12pt)` 放大

## 已知的模板原始 bug(读代码时注意,未修)

- `templates/图论/dinic.cpp`:`bfs()` 首行 `d[i] = -(i == s)` **不是 bug**(本条原先误记为"最大流恒 0",已推翻并订正):`d[s] = -1` 只当"源点已访问"的哨兵,其余为 0;判空条件 `!~d[v]`(即 `d[v] != -1`)只挡住源点,BFS 从 `d[s]+1 = 0` 开始给邻居定层、逐层递增,`d[v] == d[u] + 1` 照常成立。`-1` 在这份实现里统一表示"已访问/废弃",与 `dfs()` 末尾 `if (!out) d[u] = -1` 的剪枝一致
- `templates/数学/类欧.cpp`:结构版 `f(n,a,b,c)` 与单值版 `f(a,b,c,n)` 签名同为 `(int,int,int,int)`,C++ 不允许只按返回类型重载,**不能同时编译**,按需只留一个(已实测:报 `ambiguating new declaration of 'mint f(int, int, int, int)'`)
- `templates/计算几何/geo.cpp`:`operator-=` 与 `operator/=` 都实现成了 `x = x + y`

## 其它

- 目录已是 git 仓库(远端 `git@github.com:Z3O1/xcpc-templates.git`,private);`tmp/` 是本地编译夹具与调试产物,已进 `.gitignore`;`a.cpp` 等为无关文件,不动;`skip2004-ICPC-Templates/` 是第三方仓库,只当参考来源(见「移植来源」),目录本身不要改、也不入库
