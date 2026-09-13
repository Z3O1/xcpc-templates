# CLAUDE.md —— XCPC 模板集 · 交接文档

**这是什么**:一套**算法竞赛纸质模板**(XCPC 板子书)。`templates/` 下按章放代码片段,
`gen.py` 扫描目录生成 `sections.typ`,再用 Typst 编译成 `xcpc.pdf`(A4 横向、正文双栏),
比赛时打印出来翻。每个模板都配一份**自测**(`X.check.cpp`)守着正确性。

**本文档是完整的交接文档**:仓库状态、怎么构建/验证、写模板的规矩、每个模板的契约与已知坑、
历史上修过的 bug、上一任踩过的坑、环境限制,全部都在这个文件里 —— 接手只需要读这一份。

## 0. 现状

| 项目 | 值 |
|---|---|
| 章节 | **8 章**:字符串 / 数据结构 / 图论 / 数学 / 数论 / 多项式 / 计算几何 / 通用(顺序由 `.manifest.json` 定) |
| 代码模板 | **66** 个 `.cpp`/`.sh`(其中 2 个标了 `// hide`,只留档不进 PDF) |
| PDF 代码块 | **65** 块(64 个模板 + `通用/常数速查表.typ` 这一个"独立 typ 模板") |
| 自测 | **65** 份 `X.check.cpp`, `./check.sh` → 65 passed, 0 failed |
| `.typ` 介绍 | **17** 份(渲染在对应代码前面,不占条目) |
| 产物 | `xcpc.pdf`,**35 页**(含 3 栏目录页),页眉右上 + 页脚右下都有页码 |
| 未完成 | §14 交接清单(2 个待补模板、3 个未修缺陷、3 个未移植模板) |

## 1. 五分钟上手

```bash
git clone git@github.com:Z3O1/xcpc-templates.git && cd xcpc-templates   # private 仓库,要权限
./build.sh                   # gen.py + typst compile → xcpc.pdf
./check.sh                   # 编译并运行全部自测(全绿 = 65 passed, 0 failed)
```

其它命令:

```bash
./build.sh --watch           # 改 templates/ 后自动重编译
python3 gen.py               # 只生成 sections.typ、不编译(会打印 "8 章 / N 个代码块")
./check.sh 数论              # 只跑路径里含"数论"的 check(参数是 grep 子串)
./check.sh -v 数论           # 顺便打印每个 check 的断言输出(平时只看 ok/FAIL)
pdftoppm -png -r 150 -f N -l N xcpc.pdf tmp/x    # 出图判断排版(别用 pdftotext,见 §10)
```

- **加/删模板的完整流程**:把片段丢进 `templates/<章>/`(写法见 §4)→ `./build.sh` → `./check.sh <章>`。
  `gen.py` 会自动把它追加到该章末尾、从首行注释取标题;顺序/标题想改就编辑 `.manifest.json`(见 §3)。
- `check.sh` 里每个 check 的墙钟上限 60s(`XCPC_CHECK_TIMEOUT` 可改),**超时按 FAIL 报**,不会挂住整轮。
  check 的运行目录是它所在的那一章目录。
- 循环卡死的模板(带花树、费用流、LCT、GBST、SA 这 5 个)在 check 里用看门狗(`alarm`+`sigsetjmp`,
  LCT 用 `fork`+`waitpid`),让死循环变成 FAIL 而不是挂住。
- **改完模板必须跑 `./build.sh`**,Typst 报错会直接给 `templates/xxx.typ:行号`,照着修。

## 2. 完整清单(每个模板的标题 + check 覆盖什么)

> ⚠ 标记 = 该模板不进 PDF(`// hide`);"有 .typ 介绍"表示 PDF 里代码前面有一段文字说明。
> 标题即 PDF 里的 `== 标题`(存在 `.manifest.json`)。

### 字符串(6)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `manacher.cpp` | Manacher | 每个中心的奇/偶回文半径与 O(n²) 暴力逐中心对照,覆盖全同串/单字符/交替串/哨兵/1e6 极端 |
| `sa.cpp` | SA | 与 sort 暴力对照 + sa 是合法排列 + rk 与 sa 互逆 + `lcp(x,y)` 与暴力对照(含 1e6 极端);**有已知缺陷,见 §8.1** |
| `zfunc.cpp` | Z 函数 | z 数组与暴力 LCP 逐位对照,覆盖全同/交替/周期/Fibonacci 串 + 1e6 极端;**有契约坑,见 §8.1** |
| `最小表示法.cpp` | 最小表示法 | 返回起点循环移位后确实最小,与"枚举全部循环移位取最小"暴力对照 |
| `sam.cpp` | SAM | ⚠ 不进 PDF(用户不用 SAM)。check 只提供 `N` 再 include 本体,靠"不同子串数"性质反查正确性 |
| `pam.cpp` | 广义 PAM | ⚠ 不进 PDF(`pam.typ` 介绍也一并被跳过)。check 对照多串在线构造的节点集合 / 回文个数 / fail 链 |

### 数据结构(7)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `lct.cpp` | LCT | 与暴力(邻接表 + BFS)对照 link/cut/upd/路径异或查询(fork 看门狗) |
| `rmq.cpp` | O(1) RMQ | 与逐段扫描对照(重复值 / 单元素 / 块边界 / `dt=0/1` 两种下标约定) |
| `wqs.cpp` | wqs 二分构造方案 | 与朴素 O(n²k) DP 对照最优值,并校验返回的划分本身 |
| `二分栈.cpp` | 二分栈 | 与朴素 O(n²) DP 对照最优值(值 + 段数,`nd` 的 `(x,c)` 字典序比较) |
| `全局平衡二叉树.cpp` | 全局平衡二叉树 | 与朴素「树上带权最大独立集」树 DP 对拍 |
| `广义串并联图.cpp` | 广义串并联图 | 与「同一套约简规则的独立参考实现」对照(缩点后的核 + `add` 的记账) |
| `李超树.cpp` | 李超树 | 与「把已插入直线在该点逐条求值取 min」暴力对照(值域 [1,1e6]) |

### 图论(4)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `dinic.cpp` | 网络流 | 小图枚举割 / 随机图与暴力最大流对照 + `cut()` 割集自洽性(割容量 == 最大流) |
| `mcmf.cpp` | 费用流 | `mcmf()` 与 check 自带的 SPFA 版 SSP 费用流对拍(穷举 n=3 全边集 + 随机 n≤8、含负边);`mcmf2` 覆盖上下界可行性流 |
| `一般图最大匹配.cpp` | 一般图最大匹配 | 与两种独立暴力对拍(① 枚举边子集 ② 子集 DP `dp[mask]`)+ 看门狗抓死循环 |
| `支配树.cpp` | 支配树 | 与**按定义**的暴力对照(逐点删点看可达性,再取最近支配点) |

### 数学(5)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `barrett.cpp` | Barret 约简 | 与内建 `%` 逐值对照(小模数穷举 + u64 全域随机 + `mod_without_chk` 契约) |
| `lagrange.cpp` | 拉格朗日插值 | `lag()/lag_i()/lag_p()` 与「已知系数的多项式」和朴素养值公式对照 |
| `min25.cpp` | min_25 | 与「乘性延拓的线性筛前缀和」对照(1..3000 全枚举 + 随机到 2e6 + 1e10 量级差分) |
| `pollard-rho.cpp` | Pollard Rho | `chkp()` 与试除/素数表对照,`fact()` 与试除分解对照 |
| `类欧.cpp` | 类欧 | 结构版 `nd f(n,a,b,c)` 与单值版 `mint f(a,b,c,n)` 都与暴力求和对照;**两版不能同时编译,见 §8.1** |

### 数论(13)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `线性筛.cpp` | 线性筛 | 五个函数(素数 / φ / μ / d / 最小质因子)与暴力、独立筛法对照 |
| `杜教筛.cpp` | 杜教筛 | 与线性筛前缀和对照 + 大 n 已知值 |
| `exgcd.cpp` | 扩展欧几里得 | 恒等式 `a*x + b*y == gcd(a,b)`,覆盖 0 / 负数 / 大数 / 极值 |
| `CRT.cpp` | 中国剩余定理 | 与暴力枚举对照(含无解、非互质模、极值) |
| `BSGS.cpp` | 离散对数 | 与暴力对照(契约:`gcd(a,p)=1` 且 p 为素数) |
| `原根.cpp` | 原根与阶 | 原根最小性与阶,均与暴力对照 |
| `Miller-Rabin.cpp` | Miller-Rabin 素性检验 | 与试除/筛法对照,**重点覆盖大模数溢出**(`ksm` 换不来的地方) |
| `Pollard-Rho.cpp` | Pollard-Rho 分解质因数 | 与试除对照,覆盖平方数 / 半素数 / 大数 |
| `二次剩余.cpp` | 二次剩余 | 欧拉判别与暴力对照,开根结果验证 `r² ≡ x`(含大素数) |
| `fgcd.cpp` | O(V)-O(1) GCD | 与 `std::gcd` 全量 + 随机对照 |
| `minmod.cpp` | 最小模线性值 | 与暴力对照 `min_{0<=i<n} (a*i+b) mod m` |
| `下取整和.cpp` | 下取整和 | 与暴力求和对照,覆盖各分支与大值 |
| `分数还原.cpp` | 分数还原 | `approx(p,q,A)` 求 `x/a ≡ q (mod p)`,`|x| <= A` 且 `|a|` 最小 |

### 多项式(13)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `ntt.cpp` | 任意模数 NTT | 三模 CRT vs O(n·m) 朴素卷积(`__int128` 精确累加)大规模对拍 |
| `FFT.cpp` | 复数 FFT 卷积 | 与 O(n·m) 暴力卷积对拍 + 精度边界 |
| `多项式求逆.cpp` | 多项式求逆 | 与"按定义解系数递推"的 O(n²) 暴力 + `f*g ≡ 1` 性质与边界 |
| `多项式ln.cpp` | 多项式对数 | 两套独立暴力(`f'/f` 再积分、`log(1+u)` 级数展开)+ 代数性质 |
| `多项式exp.cpp` | 多项式指数 | 两套独立暴力(`b' = f'·b` 系数递推、`Σ f^t/t!` 展开)+ 代数性质 |
| `多项式开根.cpp` | 多项式开根 | 与"逐项比较系数"的 O(n²) 暴力 + `g² ≡ f` 与完全平方还原 |
| `多项式除法.cpp` | 带余除法 | 与朴素长除法对照,并验证 `f = q*g + r` 且 `deg r < deg g` |
| `Berlekamp-Massey.cpp` | 最短线性递推 | 与"已知线性递推序列"对拍 + 高斯消元暴力求最短阶数 + 各种退化序列 |
| `Bostan-Mori.cpp` | 线性递推第 n 项 | 与"逐项递推第 k 项"对拍 |
| `多点求值.cpp` | 多点求值 | 与逐点 Horner 暴力 + 拉格朗日插值回代 + 退化输入 |
| `快速插值.cpp` | 快速插值 | 与 O(m²) 拉格朗日插值 + "求值后再插值应还原 f" + 退化输入 |
| `多项式复合.cpp` | 多项式复合 | `F(G(x)) mod x^n` 与 O(n³) 朴素参考大规模对拍 + 代数性质 |
| `多项式复合逆.cpp` | 多项式复合逆 | `g = comp_inv(f,n)` 必须满足 `f(g(x)) ≡ x (mod x^n)`(朴素复合验证) |

### 计算几何(17)

| 文件 | 标题 | check 覆盖 |
|---|---|---|
| `geo.cpp` | 基础几何库 | p2 运算 / 比较函数 / 距离 / 线段位置关系 / 凸包 / 圆相关(强化版) |
| `半平面交.cpp` | 半平面交 | `hpi`(有向直线左侧求交)与「逐半平面暴力裁剪多边形」对拍 + 性质断言(8 万组随机才复现的退化反例已钉进 check) |
| `上凸壳.cpp` | 上凸壳 | 与「翻转 y 后求标准下凸壳」的独立参考对拍 + 性质断言 |
| `凸包内点判定.cpp` | 凸包内点判定 | `in_convex`(O(log n))与 O(n) 半平面参考 + `geo.cpp` 的 `contain` 三方对照 |
| `多边形包含.cpp` | 多边形包含 | `wn_contain` 与「精确整数缠绕数」「有向角积分」「射线法奇偶」四方对照 |
| `多边形重心.cpp` | 多边形重心 | 与「逐边积分」独立参考 + 蒙特卡洛数值重心 + 已知形状对拍 |
| `最近点对.cpp` | 最近点对 | 与 O(n²) 全点对暴力 + 退化/极端/大规模用例 |
| `最小圆覆盖.cpp` | 最小圆覆盖 | 随机增量法与「枚举两点直径圆 / 三点外接圆」暴力对照 |
| `图形交.cpp` | 圆相关求交与面积交 | 圆与直线/线段/圆求交、两圆交面积、圆与多边形交面积(蒙特卡洛 + 条带积分 + 200 万边形裁剪独立核对) |
| `简单多边形三角剖分.cpp` | 简单多边形三角剖分 | 随机凸/凹/带共线点多边形上做性质级断言 + 独立参考判定全部校验 |
| `Delaunay.cpp` | Delaunay 三角剖分 | 空圆性质(O(n²) 暴力检查)/ 剖分合法性(覆盖凸包、边恰好被两个三角形共享) |
| `Voronoi.cpp` | Voronoi 图 | 与「暴力最近站点」逐点对照(格点 + 随机采样)/ 面积和 = 包围盒面积 |
| `三维向量.cpp` | 三维向量 | p3 运算与手算/独立参考对照(叉积、混合积、模长、共线共面) |
| `三维直线.cpp` | 三维直线与线段 | 点线距 / 线线距 / 线段最近点 / 相交判定 |
| `三维平面.cpp` | 三维平面 | 点面距 / 平面交线 / 线面交点 / 两平面夹角 |
| `三维凸包.cpp` | 三维凸包 | 与「枚举所有三点成面、判断是否所有点在同侧」暴力对拍 |
| `三维旋转.cpp` | 三维旋转 | 罗德里格斯公式 / 旋转矩阵 / 坐标基底 |

### 通用(2)

| 文件 | 标题 | 说明 |
|---|---|---|
| `常数表.cpp` | 常用素数与大整数 | NTT 三模(998244353 / 1004535809 / 469762049)、1e9+7 / 1e9+9、skip2004 的 10 个大素数、哈希底数、`INF_ADD`/`INF_MUL` 哨兵。**没有 check**(纯常数清单,无法对拍) |
| `常数速查表.typ` | 常数速查表 | **独立 typ 模板**(没有同名 `.cpp`):表 A = `n!` / `C(n,n/2)` / `LCM(1..n)` / `P_n` + `n ≤ 10^k` 的 `max ω(n)` / `max d(n)` / `π(n)`;表 B = 69 个高合成数(约数最多),左右两半并排 6 列 |

## 3. 生成链路与状态文件

```
templates/<章>/<模板>.cpp|sh|py  ─┐
templates/<章>/<模板>.typ(介绍)  ─┼─→ gen.py ─→ sections.typ ─→ #include 进 xcpc.typ ─→ xcpc.pdf
templates/<章>/<模板>.check.cpp  ─┘(gen.py 直接跳过,只给 check.sh 用)
```

- Typst 的 `read()` 不能遍历目录,所以"扫描"全部由 `gen.py` 承担:扫磁盘 → 与 `.manifest.json`
  对账(记忆章节顺序 + 标题,增删改都是增量) → 输出 `sections.typ` → 确保 `xcpc.typ` 正文是
  `#include "sections.typ"`。
- **`sections.typ` 是生成物,不要手改**;`.manifest.json` 也是状态文件,**只手工改它**来调顺序/标题。
- `.manifest.json` 结构:`{"sections": [章名…], "entries": [{key, title, missing, hidden}, …]}`(现有 79 条)。
  - `missing: true` = 曾经存在、后来删掉/改名的条目(**记忆保留**,将来同名文件回来会原位复原)。
    现在有 12 条是 missing(改名或合并了):`通用/快读.cpp`、`通用/随机.cpp`、`通用/对拍.sh`、
    `通用/typ演示.typ`、`通用/高合成数表.typ`(并入 `常数速查表.typ`)、`数学/modint.cpp`、
    `数学/ntt.cpp`(→ `多项式/ntt.cpp`)、`数学/fgcd.cpp`(→ `数论/fgcd.cpp`)、
    `数学/多项式复合.cpp`、`数学/多项式复合逆.cpp`(→ `多项式/`)、`数学/类欧测试.cpp`、
    `数据结构/线段树.cpp`。**它们不出现在 PDF 里,不用管**。
    (注意 `数学/pollard-rho.cpp` 与 `数论/Pollard-Rho.cpp` 是**两个都还在**的独立模板,不是同一份)
  - `hidden: true` = 文件头有 `// hide`(当前:SAM、广义 PAM)。
- **改标题**:`gen.py` 只给*新*条目从首行注释取标题,已有条目一律以 manifest 为准 →
  想换标题就先从 `.manifest.json` 里删掉该条,再跑 `python3 gen.py`(它会从注释重取)。
- 外部依赖:`sections.typ` 首行 `#import "@preview/zebraw:0.6.3": zebraw`(代码高亮)。
  报"包找不到"先看 `~/.cache/typst`;离线环境可把它换成空壳 `#let zebraw(..a) = ..`。
- `old_versions/xcpc.typ.bak`(39499 B)是**迁移前的原始内联版**,别删:manifest 丢失时 `gen.py`
  靠它按代码内容反查、还原章节顺序与标题。
- `xcpc.typ`(主文件,非生成物):A4 横向、页边距 `0.8cm/0.6cm`、`#set text(size: 9pt)`、
  标题编号 `1.1`、目录用 3 栏 `#outline()` 后 `#pagebreak()`,正文切 `columns: 2` 再
  `#include "sections.typ"`;页码同时放页眉右上和页脚右下。

## 4. 模板写法约定(写新模板必须遵守)

- **章节**:`templates/` 下每个目录是一章,顺序由 manifest 决定(新目录追加到最后)。
- **首行注释就是标题来源**:`// 接口名(): 章节 · 标题 (备注)`。
  - 冒号前**可以为空**(允许 `// poly_inv(): 多项式 · 求逆` 这种写法)
  - 尾部括号备注会被**配平剥离**(用括号配平,能处理嵌套括号)
  - 无注释则回退文件名;标题存在 manifest 里,与注释解耦
  - **标题必须 Typst 安全**:会原样进 `== 标题`,带 `(`/`)` 会被当函数调用、带 `*`/反引号会被当标记
    → 标题写短名词(如 `网络流`,不写 `Dinic 最大流 (含割集)`)
- **头注释区**:第一个空行之前的内容不渲染。
  - **头注释后必须留一个空行**,否则 `codelines()` 取不到正文、PDF 里渲染成**空代码块**
    (实测踩过 11 个文件)
  - 正文里的 `//` 行与空行都会被 `readcode` 跳过(它跳过开头连续的空行/注释行),正文中间可随意写注释
- **不写 `main`**:模板是给读者抄的片段。自带 main 的模板会把输入输出调度整块渲染进 PDF
  (已有 3 个这样的模板被抽象掉了)
- **不写 `#include` / `using namespace std;`**:片段默认读者有 base header
  (`多项式/多项式复合.cpp` 历史上带着这两个,已清)
- **隐藏**:头注释区含 `// hide` 或 `// 隐藏` → 整块不进 PDF(去掉标记即恢复)。用于"留档但不上书"。
- **介绍 `.typ`(推荐写)**:同名 `.typ` 渲染在该代码前面,不进 manifest、不占条目。
  判断规则是"**存在同基底名的非 `.typ` 文件**" —— 没有同名代码的 `.typ` 才是独立渲染模板
  (`通用/常数速查表.typ` 就是这种)。介绍文件首行习惯写 `// 介绍: <标题>`(Typst 注释,不渲染);
  介绍里**不要再写与 `== 标题` 同级的标题**;可用 `*文本*` / `_文本_` 加粗、`` `代码` ``,
  注意 Typst **不认 `**双星粗体**`**(会渲染成"星号吞掉")。

## 5. 代码风格(house 宏,写模板必看)

模板是**竞赛代码片段**,不是能单独编译的单元:文件内不定义宏/类型,沿用读者 base header 的那套命名。
拿不准就参考 `~/0/Code/`(用户自己的板子,最接近 house 风格)与 `skip2004-ICPC-Templates/`。

- `For(i, l, r)` / `rFor(i, r, l)` / `ForD(i, l, r)` 代替 `for`;`vect<T>` 代替 `vector<T>`
- `ll` = `long long`、`db` = 浮点(常用 `long double`)、`mint` = 模数类(带 `.inv()`)、
  `poly` = 多项式 / `vector<mint>`、`ksm` = 快速幂
- `vect<T>` 是带 `+=` 的自定义容器(**没有 `push_back`**,追加用 `p += x`);要标准语义就用 `std::vector`
- `pii` 是 `array<int,2>`(支持 `f[i] = {dfn[i]}`、`que(x)[1]`),不是 `std::pair`
- 计算几何统一 `struct p2`(点)、`seg`(线段)、`line`;`eps`/`cmp`/`sgn`/`cross`/`det` 常规
- **`ksm` 的模数限制**:它内部按 `ll` 相乘,模数 > 2^32 时 `a * a` 溢出算错 → 大模数场合**自带模乘**
  (见 `数论/Miller-Rabin.cpp` 的 `MR::mul`、`数论/二次剩余.cpp` 的 `QR::mul`;
  `__int128` **必须用 unsigned** —— signed 上限 1.7e38,两个 u64 相乘到 3.4e38 会静默溢出)
- **补形状用 `resize` 不要用 `assign`**:`assign` 会清空调用方传进来的数据(踩过,见 §10)
- 多测/多次调用的模板要在入口重置自己维护的数组(踩过:支配树、费用流 `clear()`)
- 头注释里**要写清前置条件/契约**(很多模板都这么写,PDF 里读者看得到,见 §8)

## 6. 自测(check)约定

每个模板配一份同目录 `X.check.cpp`,骨架固定:

```cpp
// X 自测:<测什么>(首行注释随意写,不渲染)
#include "../_check_base.hpp"   // house 宏 + CHECK/PASSED/ok/rnd/rng
#include "X.cpp"                // 直接 include 被测模板本体
int main() { /* 断言 + 暴力对照 */; PASSED("X"); }
```

- 跑法:`./check.sh [子串] [-v]`;编译用 `g++ -std=c++17 -O2`,运行目录 = 该章目录。
  `check.sh` 的逻辑就是 `find templates -name '*.check.cpp'` + 逐个编译运行,没有别的魔法。
- `templates/_check_base.hpp` 提供:`ll/u64/s64/i128/u128/db`、`vect<T>`(薄派生自 `vector`,
  带 `+=`/`substr`/`operator+`)、`pii = array<int,2>`、`all/cmin/cmax/For/rFor/ForD`、
  `ksm`、`mint`(成员 `int v`,还有 union 别名 `x`;`.inv()`)、`rng`/`rnd(l,r)`、
  `ok()`、`CHECK(cond, name)`、`PASSED(name)`。注意它的 `ksm` 同样有 §5 的模数限制。
- **标准**:「内置**独立**参考实现 + 大规模随机对拍 + 退化/极端用例 + 性质断言」,失败要打印具体输入。
  别写成"跑一遍不崩":写完用**变异测试**验证 check 有牙(故意把模板改错一处,看它是否变红)。
- 会死循环的模板:用 `alarm` + `sigsetjmp/siglongjmp`(或 `fork`+`waitpid`)加看门狗,让它报 FAIL
  而不是挂住整轮(带花树就是这么发现那个 bug 的:**穷举 n≤5 的 1099 个图都不触发,400 组随机
  n=6..12 里才卡死 4 组** → 随机 + 看门狗是必需的)。
- 内存越界类问题用 ASAN 抓:
  `g++ -std=c++17 -O1 -g -fsanitize=address,undefined -o tmp/asan_x X.check.cpp && ./tmp/asan_x`
  (产物写 `tmp/`,别写 `/tmp`,原因见 §13)。
- **模板本体出 bug 时先报告再改**;改完必须跑该章 check 与 `./build.sh`。

## 7. 移植来源(补模板时优先搬,别从零写)

- **`~/0/Code/`** —— 用户自己的板子(接近 house 风格):`一般图最大匹配.cpp`、`上下界.cpp`、
  `Hopcroft–Karp.cpp`、`FHQ.cpp`、`FFT.cpp`、`poly 1..3.cpp`、`Bostan-Mori.cpp`、`杜教筛.cpp`、
  `Gause.cpp`、`det.cpp`、`Matrix.cpp`、`BigInt.cpp`、`Hash.cpp`、`Geo.cpp` 等
- **`~/0/Lib/atcoder/`** —— 官方 atcoder 库(`maxflow`/`mincostflow`/`scc`/`twosat`/`segtree`/
  `lazysegtree`/`convolution`/`math`/`string`),当"该算法该怎么做"的权威参考
- **`skip2004-ICPC-Templates/contents/<分类>/`** —— 第三方完整模板集(173 个文件;`numbertheory`/
  `geometry`/`math`/`string` 等,还带 `tables/Addon.tex`、`tables/dn2.md` 这类常数表)。
  **只当参考来源:目录本身不要改、也不要入库**(已 gitignore)
- 移植时**只把改动写进本仓库的 `templates/`**;`~/0/**` 与 `skip2004/**` 一律不动

## 8. 契约与已知缺陷(用模板前必须知道的前提)

### 8.1 未修缺陷(模板本身的问题,check 里以"契约"形式绕开了)

| 模板 | 问题 | 现在怎么办 |
|---|---|---|
| `字符串/sa.cpp` | ① `n == 1` 段错误(重排名循环不跑,`rk[1]` 停在字符值,LCP 段自比自身无限增长);② `lcp()` 依赖调用方 `a[0]`;③ 字符值必须 ≥ 1(0-based 映射会坏) | **未修**。check 只在"n=1 且字符值 ≥ 1"的契约内测。要修就从这三处下手 |
| `字符串/zfunc.cpp` | 注释说哨兵是 `s[n]`,实际读 **`s[n+1]`**;且要求 z 数组**零初始化** | **未修**(check 里记录了注释)。补文档或改实现 |
| `数学/类欧.cpp` | 结构版 `f(n,a,b,c)` 与单值版 `f(a,b,c,n)` 签名都是 `(int,int,int,int)` → **不能同时编译**(`ambiguating new declaration of 'mint f(int,int,int,int)'`) | **未修**,按需只留一版;要两版就改名 |
| `数据结构/全局平衡二叉树.cpp` | 5 参 `que()` 在 `(l+r)/2` 处切分,而 `build/upd` 用加权中位数 `k>>1` —— 目前没有调用路径所以不可达,但语义不一致 | 潜在坑,改的时候留意 |
| `计算几何/geo.cpp` | `ons()` 对零长线段恒真 → `contain()` 在多边形有相邻重复顶点时会误判;`eps` 是**绝对**的,坐标量级 < 1e-5 会塌缩 | 潜在坑;可照 `图形交.cpp` 里"按半径缩放判据"的思路修 |

### 8.2 前置条件(头注释里写了,调用方必须满足)

- `多项式/ntt.cpp`:变换长度 `l = 2^ceil(lg(n+m-1))` **必须 ≤ 2^20**(三模里 `1004535809 = 479·2^21+1` 最小);
  负系数不在契约内
- `多项式/多项式复合.cpp`、`多项式/多项式复合逆.cpp`:**必须先调一次 `prep(22)`**(初始化 NTT 的
  `rev/omg`,长度上界按 2^22 准备);复合逆还要求 `F[1] != 0`、`F.c.size() >= n`(n=1 时也要 2 个系数)
- `多项式/多项式ln.cpp`、`多项式开根.cpp`:`f[0] == 1`;`多项式exp.cpp`:`f[0] == 0` 且 `n >= 1`
- `多项式/多项式求逆.cpp`:`f[0] != 0`;`多项式/多项式除法.cpp`:除数 `g.back() != 0`;
  `Bostan-Mori.cpp`:`Q[0] != 0`;`Berlekamp-Massey.cpp`:模 998244353 下的序列;
  `快速插值.cpp`:`xs` 两两不同(否则 `M'(x_i) = 0` 除零)
- `数论/BSGS.cpp`:`gcd(a,p) = 1` 且 `p` 素数;`二次剩余.cpp`:奇素数模、`0 <= a < p`;
  `原根.cpp`:`gcd(a,p) = 1`;`minmod.cpp`:`m > 0`;`分数还原.cpp`:首参是模数(与 skip2004 一致);
  `下取整和.cpp`:返回值/累加器是 u64,答案不能超过 2^64−1
- `计算几何/geo.cpp::convex_hull`:**输出缓冲必须 ≥ 2n+2**(实测 `p2 b[5]` + 5 个点会 stack-buffer-overflow);
  `n <= 1` 时返回 1 但**不写点**;会**原地排序输入数组**,且输出数组 `b` 与输入 `a` 不能是同一块内存
- `计算几何/上凸壳.cpp`:同样会原地排序输入;`最近点对.cpp`:要求 `n >= 2`;
  `三维向量.cpp`:`unit()` 要求非零向量;`三维平面.cpp`:`islp()` 要求不平行;
  `三维旋转.cpp`:轴必须非零;`图形交.cpp`:`arg2(u,v) ∈ (-pi, pi]` **不能**归一到 `[0, 2pi)`;
  `半平面交.cpp`:各直线方向两两不「几乎平行」(注意是**无向**平行),`vs` 非空;
  `简单多边形三角剖分.cpp`:要求**逆时针**简单多边形

## 9. 历史 bug 表(全部已修,每条的 check 都守着,回归会变红)

| 模板 | bug | 修法 | 怎么发现的 |
|---|---|---|---|
| `图论/dinic.cpp` | `bfs()` 首行 `d[i] = -(i == s)` → 邻居拿到 `0`,而判空 `!~d[v]` 只认 `-1` → BFS 一个点都进不去,**最大流恒 0** | 改成 `d[i] = -1` 再 `d[s] = 0`(只写 `d[i] = -1` 不设源点会让源点被当未访问重新入队 → `solve()` 死循环) | 最小复现 `add(1,2,3) add(2,3,4) add(1,3,1); solve(1,3,3)` 返回 0(应 4) |
| `图论/dinic.cpp` | `cut(auto *f)` 的 `auto` 形参是 C++20 语法,C++17 只有 GCC 扩展 | 改 `int *f` | `-Werror`/clang 会挂 |
| `图论/mcmf.cpp` | `mcmf2` 内层两次 `mcmf(...)` 漏传第 4 参 `_n` → 内层 n=0 → **任何输入死循环** | 补 `_n+2` | `mf.add(1,2,1,-5); mf.mcmf2(1,2,0,2);` 挂死 |
| `图论/mcmf.cpp` | `a1 += e[tot].w` 写在 `mcmf(...)` 之前(那时 `_t→_s` 上流量还是 0,且增广会追加边、下标会变) | 存 `fake = tot`,增广后再读 `e[fake].w` | 可行性阶段的循环流永远算不进流量 |
| `图论/mcmf.cpp` | `clear()` 只清到成员 `n`,而 `n` 只在 `mcmf()` 里赋值 → 首次 `clear()` 是空操作,残留 `hd[]` 让 dij 取到奇下标前驱 → **`mcmf()` 也死循环** | 新增 `mx`(用过的最大点数),`clear()` 清到 `mx` | 3/3 组子进程被看门狗强杀 |
| `图论/一般图最大匹配.cpp` | 遍历 `for(auto v : t[u])` **少一行 `if(fd(u) == fd(v)) continue;`**:花缩起来后花内相邻点还会互相走分支 → 死循环 + `mat[]` 写脏 | 循环体第一行加这句 | 400 组随机 n=6..12 里 4 组卡死;**穷举 n≤5 的 1099 个图不触发** |
| `图论/支配树.cpp` | 处理不可达前驱时没跳过 `dfn[x] == 0` → `sd[u]` 被拉成 0、`dm[u]` 算错 | `for(auto x : t2[u]) if(dfn[x]) cmin(...)` | 最小复现 n=3、边 `1→3, 2→3` → `dm[3]=0` 应 1;19408 用例里触发 2798 次 |
| `图论/支配树.cpp` | `getdom()` 不自重置 `dfn/pos/dt/q1`(只重置 `fa/f`)→ 同进程二次调用算错 | 入口重置算法自身状态(**不动 `t[]`/`t2[]`**,图是调用方建的) | 同一进程连跑两张图结果错 |
| `计算几何/geo.cpp` | `operator-=` 与 `operator/=` 都写成 `x = x + y`;且 `/=` 只有 `(p2&, p2)` 重载,`a /= 2.0` 编译不过 | `x - y` / 补 `(p2&, db)` 版 `x / y` | check 里用 SFINAE 探测重载存在性 |
| `计算几何/半平面交.cpp` | ① 两条近同向直线(极角差 1 ulp)去重时保留**较弱**的一条;② 真交集为空时双端队列交出**顺时针假三角形**(面积 −1.4e17) | 显式比强弱;带号面积 ≤ 0 一律当退化返回空 | 随机对拍 8 万组才碰到一次,check 里钉了 8 直线可复现反例 |
| `计算几何/上凸壳.cpp` | 同 x 判据用了 eps,圆上 15° 间隔的两点被当同 x → 丢真正的最高点 | 同 x 用**精确相等**,近似同 x 交给栈的叉积处理 | 半径 1000 圆上 24 点只出 9 个顶点 |
| `计算几何/凸包内点判定.cpp` | 先算对角线分支再算外边 → 顶点本身被判成"内部"(返回 2 而非 1) | 先算外边半平面,`c==0` 直接返回 1 | 与 `contain` 三方对照 |
| `计算几何/图形交.cpp` | `cir_poly_area` 在"圆完全在多边形外"时返回**整块圆面积**;精度只有 1e-7;相切/圆内判据的绝对 eps 未按半径缩放(半径 1e-6 时失效) | 按 skip2004 思路重写成解析版(逐边有向扇形 + 弦三角形、有向角不归一),O(n);判据按半径缩放 | 4 处 check 期望值也算错了,用蒙特卡洛/条带积分/200 万边形裁剪独立核对 |
| `计算几何/geo.cpp`(`convex_hull`) | 头注释没写"输出缓冲要 ≥ 2n+2"(实测 `p2 b[5]` + 点 `(i,i²)` 写 `b[5]` 时 stack-buffer-overflow);`n<=1` 返回 1 但不写点;会原地排序输入 | 契约写进头注释 | ASAN 实测 |
| `数论/二次剩余.cpp` | 复用 base header 的 `ksm`,模数 > 2^32 时溢出 → `quadres/sqrtp` 在大素数上全错 | 自带 `QR::mul/QR::pw`(u128);Cipolla 的 `F_p[√w]` 乘法一并改用;返回类型 `int → ll` | `ksm(4000000000,2,4294967291)` 得负数;2000 组"y² 构造"的二次剩余全被判 −1 |
| `数论/分数还原.cpp` | 移植时把连分数的 `x %= y` 抄成 `x *= y` → 对几乎所有输入溢出返回垃圾 | 改回 `%=`;参数顺序对齐 skip2004 的 `(p, q, A)` | `approx(5,2,1)` 曾返回 −1000 万级垃圾值 |
| `数论/Miller-Rabin.cpp` | 复用 `ksm` 做模乘,`p > 2^32` 时溢出 | 自带 `MR::mul`(用 **unsigned** `__int128`) | `4294967291`、`2^64-59` 被误判合数 |
| `多项式/多项式复合逆.cpp` | `pw_pj` 只写前 n 项却不 `resize(n)`,随后 `comp_inv` `reverse` 整个向量 → `[n,|F|)` 的高次项倒到低次,**静默算错** | `pw_pj` 末尾 `F.c.resize(n)`;`comp_inv` 入口补 F 的系数 | n=2、F={0,1,0} → g={0,0}(应 {0,1}) |
| `多项式/多项式复合.cpp` | `Gp[1][i] = mod - G.c[i]` 越界读(调用方给的 G 不足 n 个系数) | 入口把 `G`/`F` 补到 n 个系数 | ASAN 实测 heap-buffer-overflow |
| `多项式/多项式复合逆.cpp` | 同上,且 `n == 1` 时也要至少 2 个系数(内部读 `F.c[1]`) | 同上 | ASAN 实测 |
| `字符串/manacher.cpp` | 签名 `void manacher(int n, char *s, int d)` **少个 `*`**,函数体里 `d[i]` → 编译不过 | 改 `int *d` | 字符串章 check 实测 |
| `数据结构/全局平衡二叉树.cpp` | **模板曾是完整程序(自带 main)**,而 main 末尾没有 `return 0;` —— 改名 include 后 GCC 把"non-void 函数掉出末尾"当 `__builtin_unreachable`,`-O2` 下**连循环退出判断和 `ret` 一起删掉** → 死循环 | 补 `return 0;`,并按新约定去掉了整个 main | 反汇编 + 插桩(循环计数被优化掉) |
| `多项式/ntt.cpp` | 变换长度 `l` 必须 ≤ 2^20 但没写进注释;负系数不在契约内 | 契约写进头注释 | `mul({-1},{1},998244353)` 返回 −1 |

## 10. 上一任踩过的坑(不是模板的错,但会重复踩)

- **把"两个不相交三角形"的最大匹配算成 3**(正确答案 **2**:每个三角形内部只能取 1 条边,
  取两条必有公共端点)。因为这个错的期望值,把"实现正确"当成"漏增广路",试了三种错方向改法,白费一轮。
  → 教训:**期望值要用暴力枚举独立算,别口算**。
- **把 `d[i] = -1` 当成 dinic 的修法推荐出去**:那会让源点被当未访问重新入队,`solve()` 死循环。
  → 教训:改最小复现之后必须实测,别只做纸面推理。
- **`assign` 清空调用方数据**:写 charpoly 时用 `a.assign(n, poly(n))` 补形状,把输入矩阵元素全抹了
  (对角阵还"看起来对",随机对拍才发现)。→ 教训:补形状用 `resize`。
- **空代码块**:头注释后没留空行,`codelines()` 取不到正文,PDF 里渲染成空块(一次踩了 11 个文件)。
- **`pdftotext` 不能用来判断排版**:页面折行会伪造换行,上标会被压平(`10^6` → `106`)。
  要判断排版必须 `pdftoppm` 出图看。
- **`gen.py` 的标题正则曾不认空接口名**(`// poly_inv(): ...`)→ 标题回退成文件名;尾部括号曾用
  `[^)]*` 剥离,遇到嵌套括号会剥错(已改成配平剥离)。

## 11. 未采用的尝试(留档,避免重踩)

- **特征多项式(charpoly)**:按 skip2004 风格写过 Hessenberg 版,对角阵对、但 Jordan 块给不出
  `(x−7)³`、随机 n≤4 矩阵与 Leibniz 展开 **300 组里 221 组不一致** → 已删除,没进仓库。
  要写请照 `skip2004/contents/math/charpoly.cpp` 逐行核对,并自备 Leibniz 对拍。
- **带花树的"标准结构重排"**:把 `if(vis[v]==2)` 提前、"未访问且已匹配才扩展"等三种改法都会让匹配
  变小或仍挂(真正需要的只是那行 `if(fd(u) == fd(v)) continue;`)。
- **常数表的怪符号**:`⇝`、`~=>~` 这类符号在 Typst 里渲染很怪,已去掉;高合成数表原本是单列 69 行,
  太高,改成**左右两半并排 6 列**才压到 35 行高(这是"行间距太大"那次反馈的解法)。

## 12. Typst 排版陷阱(本环境实测)

- `$...$` 数学模式里**不要用带反斜杠的多字母符号**:`\log`、`\sqrt`、`\alpha` 报 `unknown variable`;
  多字母标识符(`len`、`mcf`)会被拆成变量序列 → 复杂度等写成反引号文本或 Unicode 文本,
  数学里只留单字母与基本运算
- **LaTeX 源不能直接抄**:`\times`、`\cdot`、`^{0}`、`\multicolumn` 这些要手工转成 Typst 写法
- 强调用**单** `*` 或 `_`;`**` 无效(编译 warning "no text within stars",排版乱)
- 反引号必须成对,不成对报 `unclosed raw text`(报错跨度可能跨行,很迷惑);
  内容块 `[...]` 里的 `_`、`^`、`#`、`$`、`[`、`]` 要转义
- **表格排版**(两栏很窄):`#table(inset: 1.5pt, row-gutter: 0pt, column-gutter: 3pt)` +
  `#set text(size: 6.5pt)` + `#set par(leading: 0.35em)` 才够紧;
  想让整张表留在同一页用 `#block(breakable: false)[...]`(会整体推到下一页 —— 注意别把标题和表拆散,
  把介绍文字和表包在一起);
  长表可改**左右两半并排**(6 列布局)把行数减半
- 标题要 Typst 安全(见 §4)
- Typst 0.15:`().join("\n")` 返回 `none`(空数组 join),需要显式 `if ... else { "" }`
- 调试单点小样:写 `_t.typ` 到**项目目录**(snap 版 typst 读不到 /tmp),用 `#panic("...")`
  把值打到 stderr

## 13. 环境与运维

- **typst 是 snap 版**(本机 0.15.1)。`./build.sh` 会自动选二进制并设 `XDG_RUNTIME_DIR`。要点:
  - **不要直接用 `/snap/bin/typst`** —— 它只是 snapd 启动器,在沙箱里必然失败
    (`cannot create transient scope: DBus error ... UnixProcessIdUnknown`):snapd 要经 D-Bus 让 host
    systemd 申请 scope,而沙箱 PID 1 是 `bwrap --unshare-pid`,host 看不见。
    **这是 PID 命名空间隔离,放宽沙箱也修不好**
  - 绕法:跑 snap 载荷里的真二进制 `/snap/typst/current/bin/typst`,并给一个可写的
    `XDG_RUNTIME_DIR`(`mkdir -p tmp/xdg && XDG_RUNTIME_DIR=$PWD/tmp/xdg …`),
    否则报 `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`
  - snap 版 typst **读不到 /tmp**,调试小样要写在项目目录里
- **沙箱里没有网络**:`curl https://github.com` 报 `SSL_ERROR_SYSCALL`。git 推送若报
  `Connection closed by ... port 22`,是网络层断 SSH,改用 gh token 走 HTTPS:
  ```bash
  TOK=$(grep oauth_token ~/.config/gh/hosts.yml | awk '{print $2}')
  git push https://x-access-token:$TOK@github.com/Z3O1/xcpc-templates.git main:main
  ```
- 远端是 `git@github.com:Z3O1/xcpc-templates.git`(**private**);SSH 拉取失败时本地 `origin/main`
  可能是旧的,要用 GitHub API 或重新 fetch 确认真实状态。
- 早期会话里 `git ls-tree` 对非 ASCII 路径会加引号转义,用 `git ls-tree -z … | tr '\0' '\n'` 才 grep 得到。

## 14. 交接清单(还没做的事)

1. **两个待补模板**(放 `templates/数据结构/`,写法照 §4/§6,再配 check + `./build.sh`):
   - **线性序列并查集** —— 序列上"删点后跳过"的并查集(`fd(i)` 跳到下一个未删位置;常配离线倒序加边/删点)
   - **线性树上并查集** —— 树上版本的跳过并查集(每个点往上跳时跳过已处理的祖先)
2. **三个未修缺陷**:`字符串/sa.cpp`(n=1 段错误 + `lcp` 依赖 `a[0]` + 字符值 ≥ 1)、
   `字符串/zfunc.cpp`(哨兵实际是 `s[n+1]` + 需零初始化)、`数学/类欧.cpp`(两版重载冲突)。详见 §8.1
3. **未移植的 skip2004 模板**:`seg_in_polygon`、在线卷积、幂投影(都在 `skip2004-ICPC-Templates/contents/` 下)
4. **`通用/常数表.cpp` 没有 check**(纯常数清单,没法对拍)—— 想要覆盖率的话可以加一个轻量 check:
   "能编译 + NTT 三模的原根是 3 + 那几个大素数过 Miller-Rabin"
5. 想调整 PDF:**标题/条目顺序**改 `.manifest.json`(§3);**章节顺序**也在 manifest 的 `sections`;
   **页码/栏数/字号**改 `xcpc.typ`;改完一律 `./build.sh` 并**肉眼过一遍 PDF**(`pdftoppm` 看关键页)

## 15. 杂项

- `tmp/` 是编译夹具与调试产物,**已 gitignore**(里面有 `tmp/数论/` 早期测试驱动、`tmp/xdg/` typst 运行目录、
  一堆 ASAN 二进制与 PNG)—— 从 GitHub 克隆下来不会有它
- `a.cpp`、`in2.txt` 是无关残留文件(早期手测输入),**不动也不删**
- `skip2004-ICPC-Templates/` 是参考用的第三方模板集,**已 gitignore,不进仓库**
- 仓库内容:`templates/`(模板 + check + 介绍)、`gen.py`(生成器)、`build.sh`、`check.sh`、
  `xcpc.typ`(主文件)、`sections.typ`(生成物)、`.manifest.json`(顺序/标题状态)、`xcpc.pdf`(产物)、
  `old_versions/xcpc.typ.bak`(迁移前原始版,别删)、`CLAUDE.md`/`BUGS.md`/`TODO.md`(文档)

> `BUGS.md` 与 `TODO.md` 的内容**已全部并入本文档**(§2、§9、§10、§11、§14);它们现在只是占位指针,
> 想删可以直接删,不会丢信息。
