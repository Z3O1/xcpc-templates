# STYLE —— XCPC 模板代码风格

> 本页是模板源码的统一规则。操作命令见 [DEVELOPING.md](DEVELOPING.md)，模板组织与生成约束见 [AGENTS.md](AGENTS.md)。
> 规则分两类：**新写与被触碰的代码必须遵守** §1–§4；**既有遗留**见 §6，已豁免但不得扩散。
> 排版以根目录 `.clang-format` 为唯一机器基准；既有文件不要求一次性重排，不做全仓重排。
> 校验：`python3 lint_style.py`（只读扫描，非零退出即有违规），回归在 `tests/test_style_lint.py`。

## 1. 类型别名与公共设施

- 使用 `templates/通用/头文件.cpp` 的别名，不写裸内置类型：
  `ll`（long long）、`u64`（unsigned long long）、`i128`/`u128`、`db`（long double）、`pii`/`pll`。
  不再写 `__int128`、`unsigned __int128`、`unsigned long long`；冗余别名（`i64`/`ull`/`u32`）不用于新代码。
- 遍历统一用宏 `For(i, l, r)` / `rFor(i, r, l)`，闭区间 `[l, r]`，端点缓存见头文件。
- 普通多项式统一 `poly = vector<int>`、模数 `M = 998244353`、系数在 `[0, M)`；乘法升 `ll` 逐步取模。
- 模幂统一 `qpow(a, b, p)`（原 `ksm` 已改名），定义在 `templates/_check_base.hpp`；非零逆元写 `qpow(x, M - 2, M)`。
- 不重复定义局部 `typedef vector<...> vi/vv`，不另起 `mod` 别名代替 `M`；算法自有结构体保留。

## 2. 命名与封装

- 变量、函数、成员用 lowercase snake_case；宏除外。
- 以下既有的大写缩写名是**契约名，不改写**：`NTT`、`FFT`、`SA`、`SAM`、`BSGS`、`CRT`、`BM`、
  `DJS_mu`、`DJS_phi`、`Wilber`、`SMAWK`、`chkp`，以及宏 `For`、`rFor`。
  `chkp` 与 `数论/Miller-Rabin.cpp` 绑定，是 `数论/Pollard-Rho.cpp` 的共用入口（AGENTS.md §1），不要改成 `is_prime`。
- 有状态的“工具型”结构体保留 `_t` 后缀：`Dinic_t`、`MCMF_t`、`rmq_t`、`sieve_t`、`fgcd_t` 等；
  纯数据类型不加后缀：`p2`、`plane`、`nd`、`mat`、`poly`，既有 `SPQR_tree`、`Tri` 保留。
- **替换拼音与不透明缩写**（限新代码）：`ksm → qpow` 已完成；新名字不用 `bd`/`ins` 这类残词。
  `que` 是四边形不等式章节的回调名契约（代码与介绍页都这么写），保留。
- 命名空间：既有的 `Geo`、`IGeo`、`DOM`、`Mul`、`LCT`、`DS`、`Probe`、`Delaunay`、`Voronoi` 保留；
  新命名空间用 lowercase snake_case。
- 单字母与惯用短循环变量保留：`n`、`m`、`i`、`j`、`k`、`l`、`r`、`x`、`y` 等；
  `多项式复合逆.cpp` 的局部 `pA`/`pB`/`pC`/`dM` 保留。
- 相互关联的自由函数收进小写命名空间（如 `mul`、`geo`、`ds`）；独立板子允许全局数组 + 自由函数，
  这是纸质查阅的默认形态（全仓 87 个源码都是这样），不强制包 namespace。
- 模板片段按单编译单元被包含：自由函数默认不写 `inline`；`_check_base.hpp` 的断言/随机工具与自测内的文件级辅助用 `static`。
- **自测 `.check.cpp` 一律 lowercase snake_case**：辅助函数、局部结构体、变量都写 `ref_area`、`cut_area` 这种形式，
  参考实现的结构体也不加 PascalCase。

## 3. 排版

以 `.clang-format` 为准，要点：

- 缩进 4 空格，**不使用 tab**；命名空间内容不缩进。
- 花括号 Attach（`){` / `) {`）在同行。
- `if(`、`for(`、`while(`、`switch(` 前置括号**不加空格**；模板调用 `For(`、`rFor(` 同理。
- 指针、引用右靠：`int *a`、`const vector<int> &b`。
- `{...}` 列表内不补空格：`{v, hd[u], w}`。
- 每行一个声明；尾随注释前留两个空格。
- 不强制列宽；长行可接受，但不得为凑行宽制造怪异换行。
- 规则外写法（`if (`、tab 缩进等）须改正；`clang-format` 只用于新代码与显式整理，**不做全仓一次性重排**。

## 4. 文件头与注释

- 模板源码（`.cpp`）：首行 `// 接口名(): 章节 · 短标题 (备注)`；需要隐藏时**第二行紧跟** `// hide`；
  再空一行写共用契约。标记写 `// hide`（有空格），不写 `//hide`。
- 介绍页（`.typ`）：首行 `// 介绍: <标题>`，隐藏标记同样紧跟第二行；`templates/_table.typ` 等内部工具不受此约束。
- 自测（`.check.cpp`）：首行 `// X 自测:<覆盖内容>`，见 AGENTS.md §5。
- 函数前按需写 0～2 行中文注释，交代调用、输入/数学前提、下标/区间、返回/副作用、初始化/重复调用、资源/特殊复杂度、溢出/精度。

## 5. 生效与验证

- 触碰某模板后：`./check.sh -x <模板名>`；同章改名波及多文件时 `./check.sh <章>`。
- 改模板代码后必跑 `./build.sh`；收尾跑 `./check.sh -j 4`。
- 改完跑 `python3 lint_style.py`；改 lint 本身时跑 `python3 -m unittest discover -s tests`。
- 重命名须同步 `.check.cpp` 与同名 `.typ` 中的引用，保证仓内无旧名残留（`grep` 为 0）。
- 提交按维度独立进行，提交前 `git diff --check`，只暂存本次改动；推送仅手动。

## 6. 既有遗留（已豁免，新代码不新增）

以下写法在本轮风格统一中按用户决定保留，**新代码不得仿照**；lint 把它们列为白名单而不是违规。

| 项 | 位置 | 数量 |
| --- | --- | --- |
| 裸 `long double` | `数论/Pollard-Rho.cpp`（`sqrtl` 参数）+ 三维/几何自测的 `(long double)` 转换 | 18 |
| `typedef` | `图论/mcmf.check.cpp`、`计算几何/平面剖分/Delaunay.cpp` 及 `.check.cpp` | 4 |
| `inline` 自由函数 | `多项式/多项式复合`、`多项式复合逆`、`数学/barrett`、`数据结构/lct`、`计算几何/基础/直线与线段` | 11 |
| `static` 自由函数/全局量 | `多项式/多点求值`、`多项式ln`、`快速插值`、`数论/Pollard-Rho`、`计算几何/三维/三维平面`、`圆/最小圆覆盖` | 8 |
| `#define ls/rs` 局部宏 | `数据结构/lct`、`数据结构/全局平衡二叉树` | 4 |

命名方面保留的既有写法见 §2（`NTT`/`SA`/`chkp`/PascalCase 命名空间/`pA` 等）。
