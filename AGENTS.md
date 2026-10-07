# AGENTS.md —— XCPC 模板集 · 交接文档

**这是什么**:一套**算法竞赛纸质模板**(XCPC 板子书)。`templates/` 下按章放代码片段,
`gen.py` 扫描目录生成 `sections.typ`,再用 Typst 编译成 `xcpc.pdf`(A4 横向、正文双栏),
比赛时打印出来翻。每个模板配一份自测(`X.check.cpp`)守着正确性。

**本文档是完整交接文档**:仓库状态、怎么构建/验证、写模板的规矩、契约与已知坑、环境限制都在这一份里。

**改本文档的规矩**:**改动越小越好** —— 只做手术式修补(改个数字、补一行契约/坑),别整段重写、别扩写;
命令 `--help`、代码头注释、`.check.cpp` 首行、`git log` 里已有的东西不要抄进来,新信息并进现有小节。

## 0. 现状

| 项目 | 值 |
|---|---|
| 章节 | **8 章**:字符串 / 数据结构 / 图论 / 数学 / 数论 / 多项式 / 计算几何 / 通用(顺序见 `.manifest.json`) |
| 模板 | **69** 个代码文件,其中 2 个标了 `// hide`(只留档不进 PDF) |
| PDF | `xcpc.pdf` **40 页**(含 3 栏目录),71 个代码块 = 67 个模板 + 4 个独立 typ 页(`通用/大质数表.typ`、`通用/常数速查表.typ`、`数据结构/四边形不等式/四边形不等式.typ`、`数论/积性函数/积性函数.typ`);另有 17 份 `.typ` 介绍(渲染在对应代码前,不占条目);第 2 行起的头注释也进 PDF(见 §4) |
| 自测 | **69** 份 `X.check.cpp`,全跑 = 69 passed / 0 failed;改一个模板只跑 `./check.sh -x <模板名>` |
| 未完成 | 见 §11 → `TODO.md` |

## 1. 五分钟上手

```bash
git clone git@github.com:Z3O1/xcpc-templates.git && cd xcpc-templates   # private 仓库,要权限
./build.sh                   # gen.py + typst compile → xcpc.pdf
./check.sh -x 原根           # 只跑一个 check(改完某个模板默认就跑这个)
./check.sh -j 4              # 全量、4 路并行(全量要几分钟,只在收尾/推送前跑)
./check.sh -h                # 其余选项:-f 指定 check 文件 / -l 只列不跑 / -v 打印断言 / 子串过滤 / -j N
```

其它命令:

```bash
./build.sh --watch           # 改 templates/ 后自动重编译
python3 gen.py               # 只生成 sections.typ、不编译(打印 "8 章 / N 个代码块")
pdftoppm -png -r 150 -f N -l N xcpc.pdf tmp/x    # 出图判断排版(别用 pdftotext:折行会伪造换行、上标被压平)
```

- **加/删模板**:片段丢进 `templates/<章>/`(写法见 §4)→ `./build.sh` → `./check.sh -x <新模板名>`。
  `gen.py` 会自动追加到该章末尾、从首行注释取标题;顺序/标题想改就编辑 `.manifest.json`(见 §3)。
- check 的运行目录是它所在那一章;单个 check 墙钟上限 60s(`XCPC_CHECK_TIMEOUT` 可改),**超时按 FAIL 报**,不会挂住整轮。
- 可能死循环的模板(带花树、费用流、LCT、GBST、SA)在 check 里带看门狗(`alarm`+`sigsetjmp`,LCT 用 `fork`+`waitpid`),把死循环变成 FAIL。
- **改完模板必须跑 `./build.sh`**,Typst 报错会直接给 `templates/xxx.typ:行号`,照着修。

## 2. 完整清单

> 「这个 check 测了什么」写在各 `.check.cpp` 的**首行注释**里(本文档不再重复);标题与顺序见 `.manifest.json`。
> **⚠ = 头注释带 `// hide`,不进 PDF。**

| 章 | 模板 |
|---|---|
| 字符串(6) | `manacher` `sa` `zfunc` `最小表示法` `sam`⚠ `pam`⚠ |
| 数据结构(11) | `lct` `rmq` `wqs` `四边形不等式`(独立页) `决策单调性分治` `SMAWK` `Wilber` `二分栈` `全局平衡二叉树` `广义串并联图` `李超树` |
| 图论(4) | `dinic` `mcmf` `一般图最大匹配` `支配树` |
| 数学(4) | `barrett` `lagrange` `pollard-rho` `类欧` |
| 数论(15) | `积性函数`(小节:`线性筛` `杜教筛` `min_25`)+ `exgcd` `CRT` `BSGS` `原根` `Miller-Rabin` `Pollard-Rho` `二次剩余` `fgcd` `minmod` `下取整和` `分数还原` `高斯整数`(自配介绍页) |
| 多项式(13) | `ntt` `FFT` `多项式求逆` `多项式ln` `多项式exp` `多项式开根` `多项式除法` `Berlekamp-Massey` `Bostan-Mori` `多点求值` `快速插值` `多项式复合` `多项式复合逆` |
| 计算几何(17) | `geo` `半平面交` `上凸壳` `凸包内点判定` `多边形包含` `多边形重心` `最近点对` `最小圆覆盖` `图形交` `简单多边形三角剖分` `Delaunay` `Voronoi` `三维向量` `三维直线` `三维平面` `三维凸包` `三维旋转` |
| 通用(2) | 独立页 `大质数表.typ`($10^k$ 以上的前 10 个素数,写成 $10^k + x$,数据来自 cnblogs/ljxtt/p/13514346)+ 独立页 `常数速查表.typ` |

- `sam`(用户不用 SAM)与 `pam` 不进 PDF,它们的介绍 `.typ` 一并被跳过。
- 仍需留意的前提:`全局平衡二叉树` / `geo` → 见 §8;SA/Z/Min_25 的旧缺陷已修复并加回归断言。
- `高斯整数`(数论)配了介绍页 `高斯整数.typ`:那份资料页现在渲染在 `高斯整数.cpp` 前面,行文口味见 §4。

## 3. 生成链路与状态文件

```
templates/<章>/<模板>.cpp|sh|py       ─┐
templates/<章>/<模板>.typ(介绍)       ─┼─→ gen.py ─→ sections.typ ─→ #include 进 xcpc.typ ─→ xcpc.pdf
templates/<章>/<小节>/<模板>.cpp      ─┤(子目录 = 二级小节,里面的文件是它的三级子条)
templates/<章>/<模板>.check.cpp       ─┘(gen.py 直接跳过,只给 check.sh 用)
```

- Typst 的 `read()` 不能遍历目录 → "扫描"全由 `gen.py` 承担:扫磁盘 → 与 `.manifest.json` 对账
  (记忆章节顺序 + 标题,增删改都是增量) → 输出 `sections.typ` → 确保 `xcpc.typ` 正文是 `#include "sections.typ"`。
- **`sections.typ` 是生成物,不要手改**;`.manifest.json` 是状态文件,调顺序/标题**只手工改它**。
- `.manifest.json` 结构:`{"sections": [章名…], "entries": [{key,title,missing,hidden}, …]}`(现 90 条)。
  - `missing: true` = 曾存在、后来删除/改名的条目(**记忆保留**,同名文件回来会原位复原),现有 17 条
    (如 `通用/常数表.cpp`→`通用/大质数表.typ`、`数学/ntt.cpp`→`多项式/`、`数学/min25.cpp`→`数论/积性函数/`、`数论/两平方和.typ`→`高斯整数.typ`、`数据结构/线段树.cpp`)。**不进 PDF,别管也别清**。
  - `hidden: true` = 文件头有 `// hide`(当前 SAM、广义 PAM)。
  - 注意 `数学/pollard-rho.cpp` 与 `数论/Pollard-Rho.cpp` 是**两个都还在**的独立模板,不是同一份。
- **改标题**:`gen.py` 只给*新*条目从首行注释取标题,已有条目一律以 manifest 为准 →
  想换标题先从 `.manifest.json` 删掉该条,再 `python3 gen.py` 重取。
- **三级标题(靠目录,不用标记)**:章内的子目录 = 一个二级小节,子目录里的文件 = 它的三级子条(编号如 `2.5.1`);
  子目录里与目录同名的 `.typ`(如 `数据结构/四边形不等式/四边形不等式.typ`)是小节正文,不算子条。
  当前 `数据结构/四边形不等式/` 里挂 2.5.1 分治 / 2.5.2 SMAWK / 2.5.3 Wilber / 2.5.4 二分栈。
- 外部依赖:`sections.typ` 首行 `#import "@preview/zebraw:0.6.3": zebraw`。
  报"包找不到"先看 `~/.cache/typst`;离线环境可换成空壳 `#let zebraw(..a) = ..`。
- `old_versions/xcpc.typ.bak` 是**迁移前的原始内联版,别删**:manifest 丢失时 `gen.py`
  靠它按代码 token(忽略格式/注释)反查、还原章节顺序与标题。
- `xcpc.typ`(主文件,非生成物):A4 横向、页边距 `0.8cm/0.6cm`、9pt、标题编号 `1.1`、
  3 栏 `#outline()` 后 `#pagebreak()`、正文 `columns: 2` 再 `#include "sections.typ"`;页码在页眉右上 + 页脚右下。
  字体固定为 Libertinus Serif(正文西文)/DejaVu Sans Mono(代码西文)+Noto Sans CJK JP(中文),不能依赖系统 fallback(会选到楷体/新宋体)。

## 4. 模板写法约定(写新模板必须遵守)

- **章节** = `templates/` 下的目录,顺序由 manifest 定(新目录追加到最后)。
- **首行注释就是标题来源**:`// 接口名(): 章节 · 标题 (备注)` —— 冒号前可空;尾部括号按**配平**剥离(能处理嵌套);
  无注释回退文件名;标题存在 manifest 里,与注释解耦。
  - **标题必须 Typst 安全**(会原样进 `== 标题`):带 `(`/`)` 会被当函数调用、带 `*`/反引号会被当标记
    → 标题写短名词(如 `网络流`,不写 `Dinic 最大流 (含割集)`)。
- **头注释只有第 1 行是元数据**:`readcode()` 只丢掉第 1 个非空行(标题行,`gen.py` 同时用它取标题),
  **第 2 行起连同说明注释一律照原样渲染** —— 契约/接口写进头注释就会印在 PDF 上,`.typ` 介绍页只留长篇与公式。
  **第 1 行后要留一个空行**,否则 `codelines()`(按第一个空行切)取不到正文,manifest 丢失时反查不到该文件。
- **不写 `main`**:模板是给读者抄的片段,自带 main 会把输入输出调度整块渲染进 PDF(已有 3 个被抽象掉)。
  **不写 `#include` / `using namespace std;`**:片段默认读者有 base header。
- **隐藏**:开头连续行注释区(允许空行)含 `// hide` 或 `// 隐藏` → 整块不进 PDF(去掉标记即恢复),用于"留档但不上书"。
- **介绍 `.typ`(推荐写)**:同名 `.typ` 渲染在该代码前面,**不进 manifest、不占条目**。
  判断规则 = "**存在同基底名的非 `.typ` 文件**";没有同名代码的 `.typ` 才是独立渲染页(`常数速查表.typ`)。首行习惯 `// 介绍: <标题>`(Typst 注释,不渲染);
  介绍里**不要再写与 `== 标题` 同级的标题**;加粗用单 `*`/`_`,Typst **不认 `**双星粗体**`**。
- **资料页/介绍的行文(用户口味)**:同余条件别写成 `$a equiv b mod m$`,说人话 ——
  "4k+1 型素数""模 4 余 1 的因子个数";别用"有序/无序解数"这类术语(用户原话"根本不是人话"),
  直接说清什么算同一个;小标题用 `*加粗*` 行(见 §9);公式密度照 `数学/lagrange.typ`(前提 → 公式 → 怎么做)。
  样板 = `数论/高斯整数.typ` —— **用户会自己改这份文件,改前先读盘上真实内容,别照记忆覆盖**。

## 5. 代码风格(house 宏)

模板是**竞赛代码片段**,不是能单独编译的单元:文件内不定义宏/类型,沿用读者 base header 的那套命名。
拿不准就参考 `~/0/Code/`(最接近 house 风格)与 `skip2004-ICPC-Templates/`。

- `For(i,l,r)` / `rFor` / `ForD` 代替 `for`;`vect<T>` 代替 `vector<T>`;`ll` / `db` / `mint`(带 `.inv()`)/ `poly` / `ksm`。
- `vect<T>` 是带 `+=` 的自定义容器,**没有 `push_back`**(追加用 `p += x`);要标准语义就用 `std::vector`。
- `pii` 是 `array<int,2>`(支持 `f[i] = {dfn[i]}`、`que(x)[1]`),不是 `std::pair`。
- 计算几何统一 `struct p2`(点)/`seg`/`line`;`eps`/`cmp`/`sgn`/`cross`/`det` 常规。
- **`ksm` 的模数限制**:内部按 `ll` 相乘,模数 > 2^32 时 `a * a` 溢出算错 → 大模数场合**自带模乘**
  (见 `数论/Miller-Rabin.cpp` 的 `MR::mul`、`数论/二次剩余.cpp` 的 `QR::mul`);
  `__int128` **必须用 unsigned** —— signed 上限 1.7e38,两个 u64 相乘到 3.4e38 会静默溢出。
- **补形状用 `resize` 不要用 `assign`**(`assign` 会清空调用方传进来的数据)。
- 多测/多次调用的模板要在入口重置自己维护的数组(踩过:支配树、费用流 `clear()`)。
- 头注释里**要写清前置条件/契约**(很多模板都这么写,PDF 里读者看得到,见 §8)。
- 格式由 `.clang-format` + `format.sh` 管理(默认不改 check);生成/构建/选测工具回归用 `python3 -m unittest discover -s tests`。

## 6. 自测(check)约定

每个模板配一份同目录 `X.check.cpp`,骨架固定(模板在小节子目录里时,base 头写 `#include "../../_check_base.hpp"`):

```cpp
// X 自测:<测什么>            (首行注释写清测什么;不渲染进 PDF,但 §2 靠它)
#include "../_check_base.hpp"   // house 宏 + CHECK/PASSED/ok/rnd/rng
#include "X.cpp"                // 直接 include 被测模板本体
int main() { /* 断言 + 暴力对照 */; PASSED("X"); }
```

- 跑法:**改完一个模板只跑它那一个** —— `./check.sh -x <模板名>`(全量要几分钟;选项见 `./check.sh -h`)。
  编译 `g++ -std=c++17 -O2`,运行目录 = 该章目录,单个 check 60s 超时按 FAIL(`XCPC_CHECK_TIMEOUT` 可改)。
  逻辑就是 `find templates -name '*.check.cpp'` 逐个编译运行,没有别的魔法。
- `templates/_check_base.hpp` 提供:`ll/u64/s64/i128/u128/db`、`vect<T>`(薄派生自 `vector`,带 `+=`/`substr`/`operator+`)、
  `pii = array<int,2>`、`all/cmin/cmax/For/rFor/ForD`、`ksm`、`mint`(成员 `int v`,还有 union 别名 `x`;`.inv()`)、
  `rng`/`rnd(l,r)`、`ok()`、`CHECK(cond,name)`、`PASSED(name)`;它的 `ksm` 同样有 §5 的模数限制。
- **标准**:「内置**独立**参考实现 + 大规模随机对拍 + 退化/极端用例 + 性质断言」,失败要打印具体输入。
  别写成"跑一遍不崩":写完用**变异测试**验证 check 有牙(故意把模板改错一处,看它是否变红)。
- 会死循环的模板:用 `alarm` + `sigsetjmp/siglongjmp`(或 `fork`+`waitpid`)加看门狗,让它报 FAIL 而不是挂住整轮 ——
  带花树那个 bug 就是 400 组随机 n=6..12 里卡死 4 组才发现的(穷举 n≤5 的 1099 个图都不触发)。
- 内存越界类问题用 ASAN 抓:
  `g++ -std=c++17 -O1 -g -fsanitize=address,undefined -o tmp/asan_x X.check.cpp && ./tmp/asan_x`
  (产物写 `tmp/`,别写 `/tmp`,原因见 §10)。
- **模板本体出 bug 时先报告再改**;改完必须跑该章 check 与 `./build.sh`。

## 7. 移植来源(补模板时优先搬,别从零写)

- **`~/0/Code/`** —— 用户自己的板子(接近 house 风格):`一般图最大匹配.cpp`、`上下界.cpp`、
  `Hopcroft–Karp.cpp`、`FHQ.cpp`、`FFT.cpp`、`poly 1..3.cpp`、`Bostan-Mori.cpp`、`杜教筛.cpp`、
  `Gause.cpp`、`det.cpp`、`Matrix.cpp`、`BigInt.cpp`、`Hash.cpp`、`Geo.cpp` 等。
- **`~/0/Lib/atcoder/`** —— 官方 atcoder 库(`maxflow`/`mincostflow`/`scc`/`twosat`/`segtree`/
  `lazysegtree`/`convolution`/`math`/`string`),当"该算法该怎么做"的权威参考。
- **`skip2004-ICPC-Templates/contents/<分类>/`** —— 第三方完整模板集(173 个文件,还带 `tables/Addon.tex`、
  `tables/dn2.md` 这类常数表)。**只当参考来源:目录本身不要改、也不要入库**(已 gitignore)。
- 移植时**只把改动写进本仓库的 `templates/`**;`~/0/**` 与 `skip2004/**` 一律不动。

## 8. 前置条件与剩余限制

| 模板 | 限制 | 现在怎么办 |
|---|---|---|
| `数据结构/全局平衡二叉树.cpp` | `dfs1()` 的静态计数不重置,单进程只能建一次;只验证非负点权,中间和须在 int 内;大树递归有栈限制 | 按头注释初始化,多组另开进程;加权分界的子区间查询已修并对拍 |
| `计算几何/geo.cpp` | `eps` 是**绝对**的,坐标量级 < 1e-5 会塌缩;`proj/reflect` 要求非零长直线 | 先缩放坐标;`ons/nearest/disss` 的零长线段与重复多边形顶点误判已修 |

## 9. Typst 排版陷阱(本环境实测)

- `$...$` 数学模式里**不要用带反斜杠的多字母符号**:`\log`、`\sqrt`、`\alpha` 报 `unknown variable`;
  多字母标识符(`len`、`mcf`)会被拆成变量序列 → 复杂度等写成反引号文本或 Unicode 文本,数学里只留单字母与基本运算。
- **相邻变量必须留空格**:`$N(a+bi)$` 报 `unknown variable: bi`,要写 `$N(a + b i)$`;`$2^k u v^2$` 里的 `u v` 同理。
  `\pm` 不存在(`$\pm 1$` 排版成"pm1"),正负号写 `$plus.minus$`(`$plus.minus 1$` → ±1),别写文本 `±`。
- **数学函数名要带括号或装进 `op(...)`**:`$O(sqrt n)$` 排版成斜体"sqrt 𝑛",必须写 `$O(sqrt(n))$`;
  `gcd` 这类同理写 `$op("gcd")(a,b)$`。
- 共轭用 `macron(pi)`(渲染成 π̄);**`bar(pi)` 是错的** —— 那是模长的竖线,会渲染成 |π|。
- **LaTeX 源不能直接抄**:`\times`、`\cdot`、`^{0}`、`\multicolumn` 这些要手工转成 Typst 写法;
  `^{}`/`_{}` 的**花括号会原样印进 PDF**(`$10^{-5}$` 出 `10{−5}`),分组一律用圆括号:`$10^(-5)$`、`$F^(-k)$`。
  `\equiv` 也报 `unknown variable: quiv` —— Typst 里是**无斜杠**的 `equiv`(`$r^2 equiv -1 mod p$`),要 `≡` 就直接写文本。
- 强调用**单** `*` 或 `_`;`**` 无效(编译 warning "no text within stars",排版乱)。
- 反引号必须成对,不成对报 `unclosed raw text`(报错跨度可能跨行,很迷惑);
  内容块 `[...]` 里的 `_`、`^`、`#`、`$`、`[`、`]` 要转义。
- **表格排版**(两栏很窄):统一样式见 `templates/_table.typ`,用 `inset: (x: 2pt, y: 1pt)`、`column-gutter: 0pt`(连续网格) +
  `#set text(size: 6.5pt)` + `#set par(leading: 0.35em)` 才够紧;想让整张表留在同一页用
  `#block(breakable: false)[...]`(会整体推到下一页 —— 别把标题和表拆散,把介绍文字和表包在一起);
  长表改**左右两半并排**(6 列)可把行数减半。
- 标题要 Typst 安全(见 §4)。
- **内容页里的小标题用 `*加粗*` 单独一行,不要用 `===`**:`===` 是真标题,会进 3 栏目录、还参与 `1.1` 式编号
  (独立 typ 页会长出 `5.14.1` 这种),`数学/lagrange.typ` 与 `高斯整数.typ` 都只加粗。
- Typst 0.15:`().join("\n")` 返回 `none`(空数组 join),需要显式 `if ... else { "" }`。
- 调试单点小样:写 `_t.typ` 到**项目目录**(snap 版 typst 读不到 /tmp),用 `#panic("...")` 把值打到 stderr。

## 10. 环境与运维

- **typst 是 snap 版**(本机 0.15.1),`./build.sh` 会自动选二进制并设 `XDG_RUNTIME_DIR`:
  - **不要直接用 `/snap/bin/typst`** —— 它只是 snapd 启动器,在沙箱里必然失败
    (`cannot create transient scope: DBus error ... UnixProcessIdUnknown`):snapd 要经 D-Bus 让 host
    systemd 申请 scope,而沙箱 PID 1 是 `bwrap --unshare-pid`,host 看不见。**这是 PID 命名空间隔离,放宽沙箱也修不好**。
  - 绕法:跑 snap 载荷里的真二进制 `/snap/typst/current/bin/typst`,并给一个可写的 `XDG_RUNTIME_DIR`
    (`mkdir -p tmp/xdg && XDG_RUNTIME_DIR=$PWD/tmp/xdg …`),否则报 `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`。
  - snap 版 typst **读不到 /tmp**,调试小样要写在项目目录里。
- **不是"沙箱没网",是 22 端口不稳**:`curl https://github.com` 与 `gh api rate_limit` 都通(200、配额满额)。
  SSH 报 `Connection closed by <ip> port 22` 是 TCP 建得起来、客户端 banner 发完约 5s 后被断(收不到服务端 banner),
  **与权限/仓库无关**(git 那句 "check access rights" 是误导);实测时好时坏,别依赖。走 **443** 可以继续用 SSH key、命令不变:
  ```bash
  git ls-remote ssh://git@ssh.github.com:443/Z3O1/xcpc-templates.git HEAD   # 首次报 Host key verification failed
  ssh-keyscan -p 443 ssh.github.com >> ~/.ssh/known_hosts                    # 先把 443 的 host key 加进去
  # 长期写法:~/.ssh/config 里 Host github.com / HostName ssh.github.com / Port 443 / User git
  ```
- **或 gh token 走 HTTPS**(`head -1` 不能省:hosts.yml 里 `oauth_token` 有两行,不截会报
  `url contains a newline in its password component`):
  ```bash
  TOK=$(grep oauth_token ~/.config/gh/hosts.yml | awk '{print $2}' | head -1)
  git push https://x-access-token:$TOK@github.com/Z3O1/xcpc-templates.git main:main
  ```
- 远端是 `git@github.com:Z3O1/xcpc-templates.git`(**private**);SSH 拉取失败时本地 `origin/main`
  可能是旧的,要用 GitHub API 或重新 fetch 确认真实状态。
- 早期会话里 `git ls-tree` 对非 ASCII 路径会加引号转义,用 `git ls-tree -z … | tr '\0' '\n'` 才 grep 得到。

## 11. 待办

**待办不放这里,全部在 `TODO.md`**(待补模板、待改章节、遗留缺陷)。

## 12. 杂项

- `tmp/` 是编译夹具与调试产物,**已 gitignore**(`tmp/数论/` 早期测试驱动、`tmp/xdg/` typst 运行目录、一堆 ASAN 二进制与 PNG)—— 克隆下来没有它。
- `a.cpp` 是无关残留文件(早期手测输入),**不动也不删**(同期的 `in2.txt` 已删)。
- `skip2004-ICPC-Templates/` 是参考用的第三方模板集,**已 gitignore,不进仓库**。
- 仓库内容:`templates/`(模板 + check + 介绍)、`gen.py`、`build.sh`、`check.sh`、`format.sh`、`.clang-format`、`tests/`、`xcpc.typ`(主文件)、
  `sections.typ`(生成物)、`.manifest.json`(顺序/标题状态)、`xcpc.pdf`(产物)、`old_versions/xcpc.typ.bak`(别删)、
  `AGENTS.md`(交接文档)、`TODO.md`(待办)。
