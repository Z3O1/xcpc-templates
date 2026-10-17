# STYLE —— XCPC 模板代码风格

> 本页是模板源码的统一规则。操作命令见 [DEVELOPING.md](DEVELOPING.md)，模板组织与生成约束见 [AGENTS.md](AGENTS.md)。
> 排版以根目录 `.clang-format` 为唯一机器基准；既有文件不要求一次性重排，**新写或被触碰的文件须符合本页**。
> 文件头标题、隐藏标记、注释格式仍按 `AGENTS.md` §4，不属于本页。

## 1. 类型别名与公共设施

- 使用 `templates/通用/头文件.cpp` 的别名，不写裸内置类型：
  `ll`（long long）、`u64`（unsigned long long）、`i128`/`u128`、`db`（long double）、`pii`/`pll`。
  不再写 `__int128`、`unsigned __int128`、`unsigned long long`；冗余别名（`i64`/`ull`/`u32`）不用于新代码。
- 遍历统一用宏 `For(i, l, r)` / `rFor(i, r, l)`，闭区间 `[l, r]`，端点缓存见头文件。
- 普通多项式统一 `poly = vector<int>`、模数 `M = 998244353`、系数在 `[0, M)`；乘法升 `ll` 逐步取模。
- 模幂统一 `qpow(a, b, p)`（原 `ksm` 已改名），定义在 `templates/_check_base.hpp`；非零逆元写 `qpow(x, M - 2, M)`。
- 不重复定义局部 `typedef vector<...> vi/vv`，不另起 `mod` 别名代替 `M`；算法自有结构体保留。

`qpow` 与标准缩略名白名单见 §2。

## 2. 命名与封装

- **一律 lowercase snake_case**：命名空间、结构体/类、函数、变量、成员。不写 PascalCase、camelCase。
- 有状态的“工具型”结构体保留 `_t` 后缀：`dinic_t`、`mcmf_t`、`rmq_t`、`sieve_t`、`fgcd_t` 等；
  纯数据类型不加后缀：`p2`、`plane`、`nd`、`mat`、`poly`。
- **替换拼音与不透明缩写**：`ksm → qpow`、`chkp → is_prime`、`bd → build`、`ins → insert`、`que → query` 等，
  以语义命名取代。
- **保留标准算法缩略名**（白名单，不改写为全称）：`ntt`、`fft`、`gcd`、`exgcd`、`crt`、`sa`、`sam`、
  `bsgs`、`lct`、`dinic`、`mcmf`、`wqs`、`spqr`、`smawk`、`min25`、`djs` 等。
- 单字母与惯用短循环变量保留：`n`、`m`、`i`、`j`、`k`、`l`、`r`、`x`、`y` 等。
- 相关联的自由函数收进小写命名空间（如 `mul`、`geo`、`ds`），不散落在全局。
- 模板片段按单编译单元被包含：自由函数默认不写 `inline`；仅头文件或确有 ODR 需求时保留。
  `_check_base.hpp` 的断言/随机工具与 check 内的文件级辅助函数用 `static`。

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

## 4. 生效与验证

- 触碰某模板后：`./check.sh -x <模板名>`；同章改名波及多文件时 `./check.sh <章>`。
- 改模板代码后必跑 `./build.sh`；收尾跑 `./check.sh -j 4`。
- 重命名须同步 `.check.cpp` 与同名 `.typ` 中的引用，保证仓内无旧名残留（`grep` 为 0）。
- 提交按维度独立进行，提交前 `git diff --check`，只暂存本次改动；推送仅手动。
