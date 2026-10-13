# AGENTS.md — XCPC 模板集

算法竞赛纸质板子书：`templates/` 放代码与介绍，`gen.py` 生成章节，Typst 编译 `xcpc.pdf`（A4 横向，正文双栏、目录三栏）；每份源码配同目录的 `.check.cpp`。

本文件只维护执行规则和必要环境信息。模板清单/顺序看 `.manifest.json`，待办看 `TODO.md`，函数约束看源码注释，测试覆盖看各 `.check.cpp` 首行；数量、页数、验证结果以实际检查为准，不维护易过期的快照。通常只更新变动处，不复制命令 help、算法详解或提交历史；用户明确要求重整时可重写。

## 1. 优先遵守

- **“删/去掉板子”默认隐藏，不物理删除**：在源码头部加 `// hide`，保留源码、自测及介绍；只有明确要求物理删除才删文件。同名 `.typ` 也可单独隐藏，代码继续显示。
- **`sections.typ` 不手改**；标题与顺序改 `.manifest.json`，历史 `missing` 条目不要清理。保留 `old_versions/xcpc.typ.bak`，用于状态丢失时恢复。
- **模板本体发现 bug，先报告再修**；按 §2 验证，不能把编译成功当作正确性验证。
- **来源目录只读**：`~/0/**` 与 `skip2004-ICPC-Templates/**` 只作参考；移植写入本仓库，第三方参考目录不入库。无关残留 `a.cpp` 不动、不删。
- **不重新引入 ModInt/mint**；普通模运算直接用整数，任意模数卷积保留 `Mul::mul(a, b, P)`。
- 介绍不写“坑:”说明段。`templates/数论/高斯整数.typ` 由用户自行编辑，改前必须读取盘上真实内容，不能照记忆覆盖。

## 2. 工作流程与验证

在仓库根目录运行：

```bash
./check.sh -x 原根                    # 单模板自测（隐藏模板也要测）
./check.sh 图论                       # 算法修复还要测受影响的同章模板
./build.sh                           # gen.py + Typst → sections.typ / xcpc.pdf
./check.sh -j 4                       # 收尾/推送前全量自测
python3 -m unittest discover -s tests # 改生成、构建、选测工具时必跑
```

- **一个独立改动一个提交**：单个模板的修复/新增、一次生成或版式调整、一次工具或文档改动分别提交；完成该步并通过相应验证后立即 `git commit`，不攒大提交、不提交半成品。信息写清改了什么；提交只在本地，推送按 §7 手动进行。
- 模板新增、隐藏、恢复、改代码或介绍后必须构建；隐藏不免除自测，算法修复还要跑受影响的同章测试。
- 模板按人手抄需要排版，不用 clang-format 批量重排（`format.sh` 已删）。自测选项看 `./check.sh -h`；只列不跑用 `-l`。
- `./build.sh --watch` 只在启动前生成一次章节，再监视 Typst 依赖；**增删文件、修改 manifest 或隐藏标记后须重启**。只生成章节用 `python3 gen.py`。
- 自测用 `g++ -std=c++17 -O2`，在 check 所在目录运行（含小节子目录）；每个运行默认限时 60 秒，超时按 FAIL。`XCPC_CHECK_TIMEOUT` 调整限时，`CXX` 覆盖编译器。
- 产物、日志、调试脚本放已忽略的 `tmp/`；提交/推送前检查 `git diff --check`，只暂存本次改动，不夹带用户改动、临时文件或参考库。
- 附加 worktree 统一放在 `~/XCPC/.worktrees/`，该目录须被 Git 忽略；不在 `~` 下另建并列容器。

## 3. 文件与生成规则

| 文件/目录 | 职责 |
|---|---|
| `templates/<章>/` | 源码、同名介绍 `.typ`、同名自测 `.check.cpp` |
| `templates/_check_base.hpp` | 自测用 house header 兼容层，不进 PDF |
| `.manifest.json` | 章节顺序、条目标题与历史状态 |
| `gen.py` → `sections.typ` | 扫描、对账、生成正文，跳过 check 与隐藏项 |
| `xcpc.typ` → `xcpc.pdf` | 主排版与成品，正文 include 生成章节 |
| `templates/_table.typ` | 紧凑资料表的公共样式 |
| `tests/`、`build.sh`、`check.sh` | 工具回归、构建与自测入口 |
| `old_versions/xcpc.typ.bak` | 迁移前内联版，用于恢复顺序与标题 |

- 章节是 `templates/` 的一级目录；新章节/模板自动追加，已有标题**以 manifest 为准**，修改头注释不会覆盖它。
- `missing: true` 记忆曾存在的条目，同路径文件回来会原位恢复；`hidden: true` 来自头部 `// hide` 或 `// 隐藏`，去掉标记即可恢复上书。manifest 是状态，不是待清理垃圾。
- 子目录表示二级小节，内部文件是三级条目；与子目录同名的 `.typ` 是小节正文，不单独占子条（如 `数据结构/四边形不等式/`）。
- 有同基底名代码的 `.typ` 是介绍，放在代码前，不进 manifest；代码隐藏时介绍一起跳过，介绍自身的标记只隐藏介绍。
- 没有同名代码的 `.typ` 是独立资料页：大质数表、常数速查表、四边形不等式、积性函数、反射容斥；不要漏掉 `数学/反射容斥.typ`。
- **同名项按路径区分**：隐藏 Dinic 只指 `图论/dinic.cpp`，不含费用流 `mcmf`；`图论/支配树.typ` 单独隐藏，代码保留。
- **两份 Pollard-Rho 都保留**：`数学/pollard-rho.cpp` 是隐藏旧版（自带 MR，因子升序去重），上书的是 `数论/Pollard-Rho.cpp`（依赖 `chkp`，批量 gcd，支持 u64，因子升序带重数）。旧版 MR 已并入 `数论/Miller-Rabin.cpp`；两份分解不是改名残留，不合并、不删其中一份。

## 4. 模板写法与 house 风格

```cpp
// 接口名(): 章节 · 短标题 (备注)

// 模板共用的前置条件与接口契约。
```

- 第一行注释用于新条目的标题，用短名词；尾部备注按配平括号剥离，无标题注释则回退文件名。**标题后留空行**；渲染只丢第一个非空行，后面的契约照常打印。
- 隐藏标记放在开头连续行注释区（允许空行）；放在正文中不会隐藏。
- 片段不写 `main`、`#include`、`using namespace std;`，不重复定义读者 base header 的宏和基础类型；算法自身的结构体正常保留。
- 使用 `For/rFor`、`vect<T>`、`ll/db/poly/ksm`；`pii` 是 `array<int,2>`，不是 pair。几何沿用 `p2/seg/line` 与 `eps/cmp/sgn/cross/det`。
- `vect` 追加用 `+=`，不假设用户容器有 `push_back`；需要标准容器就明确用 `std::vector`。自测兼容层继承 vector，不能反推用户 header 的接口。
- 普通多项式用 `poly = vector<int>`、`M = 998244353`，系数规范到 `[0, M)`；乘法升 `ll` 并逐步取模，非零逆元用 `ksm(x, M - 2, M)`。
- `ksm` 内部是 signed `ll` 乘法；指数非负、底数规范到 `[0, p)`，正模数须满足 `(p-1)^2 <= LLONG_MAX`（`p <= 3037000500`）。大模数参照 Miller-Rabin、二次剩余模板使用安全模乘。
- **完整 u64×u64 范围用 `unsigned __int128`**；signed 溢出是未定义行为。有界的带符号计算仍可用 `i128`，不要求所有 `__int128` 都改成 unsigned。
- 补容器形状用 `resize`，不用会清空调用方数据的 `assign`；可重复调用的入口须重置自己维护的状态，不能默默沿用上次结果。单次构建等例外必须明确注释。

移植优先顺序：`~/0/Code/`（用户板子）→ `~/0/Lib/atcoder/`（官方参考）→ `skip2004-ICPC-Templates/contents/`（第三方参考）；只移植需要的部分，不改来源。

### 4.1 函数前注释：功能与调用约束

- **每个独立函数前写简短中文 `//` 注释**，包括成员函数和独立辅助函数；局部 lambda 不必机械添加。已有准确注释直接保留或精简，不重复添加。
- 格式为 **“功能；要求……；返回/修改……”**，默认 1 行，必要时 2～3 行；只写实际需要的信息，不机械列参数或复述函数名。
- **输入与数学前提**：规模、值域、素数/互质、有序/单调、非零等；下标从 0 还是 1 开始，区间闭开也要在有歧义时说明。
- **状态与依赖**：是否须初始化、先调用其他接口、能否重复调用；递归栈等资源要求按需说明。
- **结果与副作用**：无解返回值、最小解还是任意解、是否有序/去重、追加还是覆盖、是否修改输入或内部状态。
- **数值限制**：中间结果的溢出上界、模数范围、精度要求；能给具体条件就不只写“注意溢出”。特殊复杂度按需写，不逐函数机械重复。
- 共用约束集中在模板头部，各函数只补自己的特殊约束；不能为了凑一行省略必要条件。算法原理放 `.typ`，代码内部只解释不直观的步骤。
- **注释必须与实现一致**：先读实现与自测，不凭算法名称猜契约，也不修改代码来迎合注释；发现不一致按 §1 先报告。

```cpp
// 分解质因数；要求 x >= 1；升序带重数追加到 res。
void fact(u64 x, vect<u64> &res);
```

## 5. 自测标准与关键契约

check 首行写 `// X 自测:<覆盖内容>`，随后 include `../_check_base.hpp` 和 `X.cpp`；小节子目录改为 `../../_check_base.hpp`，成功用 `PASSED("X")`。兼容层接口直接读头文件，不在本文复制实现。

- **独立参考实现 + 大规模随机对拍 + 退化/极端用例 + 性质断言**；失败打印具体输入，不能只测“不崩”。
- 新写/强化 check 后做变异测试：临时故意改错模板，确认测试变红再恢复；不留下故障注入，不覆盖用户已有改动。
- 可能卡死的测试须有超时保护；`check.sh` 的外层 timeout 是统一兜底，内部可用子进程限时或 alarm。**fork 隔离不是看门狗**，不能把阻塞 waitpid 当超时机制。
- 内存/整数疑点用 ASAN+UBSAN；编译产物放 `tmp/`，运行时仍切到 check 所在目录：

```bash
mkdir -p tmp
g++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -o "$PWD/tmp/asan_x" templates/图论/支配树.check.cpp
(cd templates/图论 && ../../tmp/asan_x)
```

以下只保留容易被误改的契约；完整边界和覆盖以源码、自测为准：

| 模板 | 必须保留的契约 |
|---|---|
| `字符串/sa.cpp` | 只构造 `sa/rk/height`，任意后缀 LCP 查询需另配 RMQ。 |
| `数据结构/全局平衡二叉树.cpp` | 静态计数不重置，单进程只建一次，多组另开进程；非负点权，中间和在 int 内，递归须有足够栈。按头注释初始化。 |
| `计算几何/geo.cpp` | eps 是绝对容差，极小尺度须先缩放；`proj/reflect` 要求非零长直线，没有通用于所有运算的精确尺度分界。 |
| `数论/Pollard-Rho.cpp` | `fact(x, res)` 升序带重数追加；`rho` 依赖同时 include Miller-Rabin 的 `chkp`，x 为素数时原样返回。初筛只到 20 与完全平方，其余交给 rho。 |
| `数论/Miller-Rabin.cpp` | `chkp` 对 u64 全范围确定性判定且为 `constexpr`，保留 7 底数确定集；`a % p == 0` 必须跳过，否则 73 等素数会被误判。安全模乘与快速幂沿用 lambda，保持编译期求值兼容，不换成普通 `ksm`。 |

## 6. Typst 介绍与排版

- 介绍首行用 `// 介绍: <标题>`，不再造与模板同级的标题；小标题单独一行 `*加粗*`，不用会进入目录和编号的 `===`。粗体用单 `*`/`_`，不用 Markdown 双星。
- 数学不用 LaTeX 反斜杠命令或花括号分组：写 `sqrt(n)`、`log(n)`、`10^(-5)`，不写 `\sqrt`、`10^{-5}`；相邻变量留空格，多字母名称用文本或 `op("gcd")(a,b)`。
- 正负号用 `plus.minus`，共轭用 `macron(pi)`，`bar(pi)` 是模长竖线；避免不成对的反引号，内容块中的 `_ ^ # $ [ ]` 等特殊字符须转义。
- 两栏表格复用 `_table.typ` 的字体、间距与网格；整表连标题/介绍放进 `#block(breakable: false)[...]`，长表可拆成左右两半并排。
- 主字体固定为 Libertinus Serif / DejaVu Sans Mono + Noto Sans CJK JP，不依赖系统 fallback。
- 调试小样放 `tmp/`（需要时放项目目录），可用 `#panic(...)` 输出值；构建错误按 `templates/xxx.typ:行号` 修复。**版面看 PNG，不用 pdftotext 判断换行/上标**：

```bash
pdftoppm -png -r 150 -f N -l N xcpc.pdf tmp/page
pdfinfo xcpc.pdf | grep '^Pages:'
```

## 7. 环境与 GitHub

- Typst 0.15.1 位于 `/usr/bin/typst`；优先用 `./build.sh`，支持 `TYPST` 覆盖，优先普通 PATH 二进制，并设置可写的 `tmp/xdg`。
- 远端为 private 仓库 `Z3O1/xcpc-templates`，origin 使用 SSH。22 端口不可靠时，用已登录 gh 的 HTTPS 凭据助手；**token 不写进 URL、日志或仓库文件**：

```bash
gh auth status
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  fetch https://github.com/Z3O1/xcpc-templates.git main:refs/remotes/origin/main
# 检查远端是否领先并确认提交内容，再 push；不自动强推
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  push https://github.com/Z3O1/xcpc-templates.git main:main
```

以上配置只用于单次命令，不永久修改 origin/helper。SSH 可走 `ssh.github.com:443`，但首次信任须按 GitHub 官方指纹核验，不未经核验追加 host key。拉取失败时 `origin/main` 可能陈旧；推送后核实远端 SHA，不凭旧缓存宣称已同步。
