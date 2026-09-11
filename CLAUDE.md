# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

算法竞赛模板集:XCPC 板子库,`templates/` 下按章目录放代码模板,一键生成并编译成 `xcpc.pdf`(横向 A4、双栏正文、右上+右下大页码)。

## 常用命令

```bash
./build.sh                   # python3 gen.py && typst compile xcpc.typ → xcpc.pdf
./build.sh --watch           # 改 templates/ 后 typst 自动重编译
python3 gen.py               # 只看生成、不编译
./check.sh                   # 编译并运行 templates/**/*.check.cpp 自测(每个模板一份)
./check.sh -v 数论           # 只跑某章并打印每条断言的输出
```

- 改模板后跑 `./build.sh`,Typst 直接报 `error: ... templates/xxx.typ:行号` 加文件路径,照报错修即可
- **本环境 typst 是 snap 版且不在 PATH**:直接跑 `./build.sh` 就行(已自动选对二进制并设好 `XDG_RUNTIME_DIR`)。要点与坑:
  - **不要用 `/snap/bin/typst`**——它只是 snapd 启动器(→ `/usr/bin/snap`),在 DSH 沙箱里**必然启动失败**:snapd 要先经 D-Bus 向 systemd 申请 transient scope,而沙箱 PID 1 是 `bwrap --unshare-pid`,host 的 systemd/dbus 看不见命名空间里的 PID → `cannot create transient scope: DBus error ... UnixProcessIdUnknown`(PID 1 为 bwrap 时也可能报 `job .../job/N finished with result failed`),退出码 46。**这是 PID 命名空间隔离,与文件权限无关,放宽沙箱(danger-full-access)也修不好**;失败发生在 snapd 申请 scope 阶段,真二进制根本没被执行
  - 绕法:直接跑 snap 载荷里的真二进制 `/snap/typst/current/bin/typst`(typst 0.15.1,classic confinement,能读仓库目录)。手动编译:`mkdir -p tmp/xdg && XDG_RUNTIME_DIR=$PWD/tmp/xdg /snap/typst/current/bin/typst compile xcpc.typ`;不设 `XDG_RUNTIME_DIR` 会报 `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`
  - 沙箱里也没有网络(`curl https://github.com` 报 `SSL_ERROR_SYSCALL`),所以别指望现下 non-snap 版 typst;主机终端里 `/snap/bin/typst` 是正常的
- **标题必须 Typst 安全**:标题会原样进 `sections.typ` 的 `== 标题`,带 `(`/`)` 会被当函数调用报 `unclosed delimiter`,带 `*`/反引号会被当标记——标题写短名词,细节放注释
- **改了首行注释但标题没变**:`gen.py` 只给新条目从注释取标题,已有条目以 manifest 为准 → 想换标题就先从 `.manifest.json` 里摘掉该条再 `python3 gen.py`
- 调试单点 Typst 小样:写个 `_t.typ` 放到**项目目录**(snap 版 typst 读不到 /tmp!),`typst compile _t.typ`,用完删;要用 `#panic("...")` 把值打到 stderr——`pdftotext` 的页面折行会伪造换行字符,不能用来调试字符串
- **模板自测**:每个模板配一份 `X.check.cpp`(仿 skip2004 的 check),`./check.sh` 批量编译运行
  - 写法:`#include "../_check_base.hpp"`(提供 house 宏 + `CHECK`/`PASSED`/`rnd`)再 `#include "X.cpp"`(被测模板本体),main 里做暴力对照/断言,末尾 `PASSED("名字")`
  - `gen.py` 会跳过 `*.check.cpp`,所以 check 不会进 PDF
  - 单个 check 有 60s 墙钟上限(`XCPC_CHECK_TIMEOUT` 可调);卡死会报成 FAIL
  - `templates/计算几何/geo.check.cpp` 目前**故意 FAIL**(复现 geo.cpp 的 `-=`/`/=` bug),修好即转绿
- 视觉验证:抽查 PDF 页(`pdftoppm -png -f N -l N xcpc.pdf /tmp/x`)

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
- **`ksm` 的模数限制**:base header 的 `ksm` 内部按 `ll` 相乘,模数超过 2^32 时 `a * a` 会溢出算错 → 需要大模数(如 `u64` 素性检验、`p > 4e9` 的离散对数)时**自带模乘**,不要复用 `ksm`(见 `templates/数论/Miller-Rabin.cpp` 的 `MR::mul`、`templates/数论/二次剩余.cpp` 的 `QR::mul`)
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

## 已知的模板原始 bug

> 状态:下面除**带花树**一条外都已在 2024 这轮修复(修法与实测证据见各条描述);每条都有对应的
> `*.check.cpp` 守着,回归就会红。带花树那条只做了诊断,待重写缩花流程。

- `templates/图论/dinic.cpp`:`bfs()` 首行 `d[i] = -(i == s)` **是 bug(最大流恒 0)**:只有 s 得 -1、其余为 0,而判空条件 `!~d[v]` 只认 -1,邻居 `d[v] = 0` 被当"已访问"、BFS 一个点都进不去。实测(模板本体):`f.add(1,2,3); f.add(2,3,4); f.add(1,3,1); f.solve(1,3,3)` 返回 0(应 4)。**修法是 `d[i] = (i == s) - 1`(即 d[s] = 0、其余 -1),不是 `d[i] = -1`** —— 后者会把源点也当"未访问",残余图有边回 s 时 s 被重新定层入队,`while (bfs())` 恒真而 dfs 推不进流 → solve() 死循环(实测:6 条容量 2 的边 1-2/1-3/2-3/2-4/3-5/4-5、solve(1,5,5) 应 4 却挂死)。`dinic.typ` 的原记载("最大流恒 0")是对的,本条此前被我两次记反,现以实测为准
- `templates/图论/dinic.cpp`:`cut(auto *f)` 的 `auto` 形参是 C++20 语法,C++17 下只是 GCC 扩展(告警),`-Werror`/clang 会挂
- `templates/数据结构/全局平衡二叉树.cpp`:**模板是完整程序(自带 main)**,而 main 末尾原来没有 `return 0;` —— 这样单独编译运行没问题(C++ 规定 main 掉出末尾等价 return 0),但**一旦把 main 改名后 include 再调用**(写 check 时的常见做法),GCC 会把"非 main 的 non-void 函数掉出末尾"当 `__builtin_unreachable`,在 -O2 下**连循环退出判断和 ret 一起删掉** → 变成死循环。已补 `return 0;`(一行修复)。教训:驱动这类"自带 main 的模板"时,要么补 return,要么只 include 不调用它的 main。另外它依赖 base header 的 `using namespace std`(裸用 `max`),单独编译会报 `max was not declared`。其 5 参数版 `que(L,R,k,l,r)` 用 `(l+r)/2` 分界而 `build`/`upd` 用加权中位数 `k>>1`,权值不均时子区间查询会走到 `tr[0]`(模板自身调用点都传整链,所以当前不可达)
- `templates/数学/类欧.cpp`:结构版 `f(n,a,b,c)` 与单值版 `f(a,b,c,n)` 签名同为 `(int,int,int,int)`,C++ 不允许只按返回类型重载,**不能同时编译**,按需只留一个(已实测:报 `ambiguating new declaration of 'mint f(int, int, int, int)'`)- `templates/多项式/多项式复合逆.cpp`:`comp_inv(F, n)` 只在 `F.c.size() == n` 时正确 —— `pw_pj` 只写前 n 项却不 `resize(n)`,接着 `comp_inv` 又 `reverse` 整个向量,于是 `[n, |F|)` 的高次项被倒到低次参与运算,**静默算错**(不崩,返回长度仍是 n)。最小复现:n=2、F={0,1,0} → g={0,0}(应 {0,1});n=4、F={0,1,1,0,0} → 全 0。修法:`pw_pj` 末尾加 `F.c.resize(n)`。模板自带 main 会先 `resize(n)`,所以契约内用法没暴露
- `templates/多项式/多项式复合.cpp` 与 `多项式复合逆.cpp` 有**未文档化的越界读前置条件**(ASAN 实测 heap-buffer-overflow):`comp` 要求 `G.c.size() >= n`;`comp_inv` 要求 `F.c.size() >= n`,且 `n == 1` 时也要给到 2 个系数(内部读 `F.c[1]`)
- `templates/多项式/ntt.cpp`:变换长度 `l = 2^ceil(lg(n+m-1))` 必须 ≤ **2^20**(三模里 M1 = 1004535809 = 479·2^21+1 最小);超限结果错,但 P 恰为 998244353 / 469762049 时 CRT 退化成单模、反而"看着对",只拿这两个试上限会得出错误结论。另外**负系数不在契约内**(`mul({-1},{1},998244353)` 返回 -1,应传 `P-1`);非负但 ≥ P 的未归一系数反而是对的
- `templates/计算几何/geo.cpp`:`operator-=` 与 `operator/=` 都实现成了 `x = x + y`

## 其它

- 目录已是 git 仓库(远端 `git@github.com:Z3O1/xcpc-templates.git`,private);`tmp/` 是本地编译夹具与调试产物,已进 `.gitignore`;`a.cpp` 等为无关文件,不动;`skip2004-ICPC-Templates/` 是第三方仓库,只当参考来源(见「移植来源」),目录本身不要改、也不入库

### 一般图最大匹配.cpp 的进一步诊断(2024 实测,未修)
原始症结:`for(auto v : t[u])` 里"增广"那一支(判据仅 `!mat[v]`)写在最前面,遇到**搜索根
rt 自己**(已访问且未匹配,`vis[rt] = 2`、`mat[rt] = 0`)也会走增广:此时 `pr[rt]` 从未被
赋值,回溯 `while(v) mat[v] = pr[v], swap(mat[pr[v]], v)` 里 `mat[pr[rt]] = mat[0]`,
于是 `v` 在 `rt` 与 `mat[0]` 之间来回换 → 死循环,并把 `mat[0]` 写脏(最小复现:n=6,
边 `1-2 1-6 2-6 3-4 3-5 4-5`,标准驱动卡在 `getmat(6)`)。
已试过三种修法都**不行**(实测,均已回退):
1. 把 `if(vis[v] == 2)` 提到最前 → 死循环消失,但 3 元奇环上匹配变小(两个三角形只得 2,应 3);
2. 改成"未访问且已匹配才扩展 / 未访问且未匹配才增广" → `lca()` 依赖的并查集在根自身未匹配时没建好,同样错;
3. 增广时加 `v != rt` 守卫 + 回溯到根即停 + 同花不再缩 → 不挂不脏,但两个三角形仍只得 2。
结论:这不是"调换两个 if"能修好的,这份实现的缩花流程在**根进入搜索时未匹配**的语义下不完整
(参考实现如 skip2004 的带花树,根进搜索前是已匹配状态)。要修得重写缩花/回溯,工作量不小;
它的 check 保持"响亮打印 [BUG] + 用例硬断言"的形态,等以后重写。
