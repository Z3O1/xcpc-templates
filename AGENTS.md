# AGENTS.md — XCPC 模板集

算法竞赛纸质板子书：`templates/` 放代码与介绍，`gen.py` 生成章节，Typst 编译 `xcpc.pdf`（A4 横向，正文双栏、目录三栏）。代码片段各配同目录的 `.check.cpp`。

本文件只保留执行规则与必要环境信息；模板清单/顺序看 `.manifest.json`，待办看 `TODO.md`，测试覆盖看各 `.check.cpp` 首行。通常只更新变动处，不复制命令 help、代码注释或提交历史；用户明确要求重整时可重写。

## 1. 优先遵守

- **“删/去掉板子”默认隐藏，不物理删除**：在模板头部加 `// hide`，保留源码、自测及介绍；只有明确要求物理删除才删文件。同名 `.typ` 介绍也可单独加标记，代码继续显示。
- **`sections.typ` 不手改**；标题与顺序改 `.manifest.json`。历史 `missing` 条目不要清理；保留 `old_versions/xcpc.typ.bak`，用于状态丢失时恢复。
- **模板本体发现 bug，先报告再修**。改完跑相关 check 和 `./build.sh`；不能把编译成功当作正确性验证。
- **来源目录只读**：`~/0/**` 与 `skip2004-ICPC-Templates/**` 只作参考；移植写入本仓库，第三方参考目录不入库。无关残留 `a.cpp` 不动、不删。
- **不重新引入 ModInt/mint**。普通模运算直接用整数，任意模数卷积保留 `Mul::mul(a, b, P)`。
- 介绍不写“坑:”说明段。`数论/高斯整数.typ` 由用户自行编辑，改前必须读盘上真实内容，不能照记忆覆盖。

## 2. 工作流程与验证

在仓库根目录运行：

```bash
./build.sh                         # gen.py + Typst → sections.typ / xcpc.pdf
./check.sh -x 原根                  # 改一个模板，选跑它的自测
./check.sh -j 4                     # 收尾/推送前全量自测
python3 -m unittest discover -s tests  # 改生成、构建、选测、格式工具时必跑
./format.sh --check                 # 格式验证；选项看 ./format.sh -h
```

- 模板新增、隐藏、恢复、改代码或介绍后，必须构建；隐藏不免除自测。算法修复还要跑受影响的同章测试（例如 `./check.sh 图论`）。
- `./format.sh` 默认不改 check；避免与任务无关的批量格式改动。自测选项看 `./check.sh -h`，只列不跑用 `-l`。
- `./build.sh --watch` 仅启动前跑一次生成器，再监视 Typst 依赖；**加删文件、修改 manifest 或隐藏标记后要重启**，不能指望它重新扫描目录。只更新章节用 `python3 gen.py`。
- 自测用 `g++ -std=c++17 -O2`，运行目录为该 check 所在目录（含小节子目录）。默认每个运行限时 60 秒，超时按 FAIL；用 `XCPC_CHECK_TIMEOUT` 调整，`CXX` 可覆盖编译器。
- 产物、日志、调试脚本统一放已忽略的 `tmp/`。推送前检查 `git diff --check`，不要顺手提交临时文件或参考库。

最近验证快照（版面/模板变动后重算，不是永久保证）：

| 项目 | 当前值 |
|---|---|
| 章节 | 8 章：字符串 / 数据结构 / 图论 / 数学 / 数论 / 多项式 / 计算几何 / 通用 |
| 源码与自测 | 69 份源码、69 份 check；源码隐藏 6 份，隐藏模板仍测试 |
| PDF | 38 页；68 个生成条目 = 63 份可见源码 + 5 份独立资料页；另渲染 15 份介绍 |
| 验证 | 69 C++ passed、29 Python passed；格式检查与构建通过 |

隐藏源码为 `sam`、`pam`、`线性筛`、`李超树`、`全局平衡二叉树`、`dinic`。**“网络流”只隐藏 Dinic，费用流 `mcmf` 保留**；`支配树.typ` 单独隐藏，支配树代码保留。

## 3. 文件与生成规则

| 文件/目录 | 职责 |
|---|---|
| `templates/<章>/` | 代码片段、同名介绍 `.typ`、同名自测 `.check.cpp` |
| `templates/_check_base.hpp` | 自测用 house header 兼容层；不进 PDF |
| `.manifest.json` | 章节顺序、条目标题与历史状态 |
| `gen.py` → `sections.typ` | 扫描、对账、生成正文；跳过 check、隐藏文件/目录与隐藏条目 |
| `xcpc.typ` → `xcpc.pdf` | 主排版与成品；正文 include 生成章节 |
| `templates/_table.typ` | 紧凑资料表的公共样式 |
| `tests/`、`build.sh/check.sh/format.sh` | 工具回归与执行入口 |
| `old_versions/xcpc.typ.bak` | 迁移前内联版，恢复章节顺序/标题的依据 |

- 章节是 `templates/` 的一级目录；新章节/模板自动追加，已有标题**以 manifest 为准**，修改头注释不会覆盖它。调顺序/标题直接改 manifest。
- `missing: true` 记忆曾存在的条目；同路径文件回来会原位恢复。`hidden: true` 来自头部 `// hide` 或 `// 隐藏`，去掉标记即可恢复上书。manifest 是状态，不是待清理垃圾。
- 子目录表示二级小节，里面的文件是三级条目；与子目录同名的 `.typ` 是小节正文，不单独占子条（例：`数据结构/四边形不等式/`）。
- 有同基底名代码的 `.typ` 是介绍，放在代码前，不进 manifest、不占条目；代码隐藏时介绍一起跳过。介绍自身的头部隐藏标记只隐藏介绍。
- 没有同名代码的 `.typ` 是独立资料页。当前为大质数表、常数速查表、四边形不等式、积性函数、反射容斥；不要漏掉 `数学/反射容斥.typ`。
- `数学/pollard-rho.cpp` 与 `数论/Pollard-Rho.cpp` 是两份仍存在的独立模板，不能当成改名残留合并。

## 4. 模板写法与 house 风格

```cpp
// 接口名(): 章节 · 短标题 (备注)

// 前置条件、接口契约；这些行会印进 PDF。
```

- 第一行注释用于新条目的标题，标题用短名词，避免复杂括号/标记；会原样进入 Typst 标题。尾部备注按配平括号剥离，无标题注释则回退文件名。
- **标题后留空行**，供代码匹配/历史恢复使用。渲染只丢第一个非空行，第 2 行起的契约照常打印。隐藏标记须在开头连续行注释区，允许空行；放到正文中不会隐藏。
- 片段不写 `main`、`#include`、`using namespace std;`，不重复定义读者 base header 的宏和基础类型；算法自己的结构体正常保留。
- 使用 `For/rFor/ForD`、`vect<T>`、`ll/db/poly/ksm`。`pii` 是 `array<int,2>`，不是 pair；几何沿用 `p2/seg/line` 与 `eps/cmp/sgn/cross/det`。
- 片段中 `vect` 追加用 `+=`，不要假设用户 house 容器有 `push_back`；需要标准容器就明确用 `std::vector`。自测兼容层继承 vector，不能反推用户 header 的接口。
- 普通多项式用 `poly = vector<int>`、`MOD = 998244353`，系数规范到 `[0, MOD)`；乘法升 `ll` 并逐步取模，非零逆元用 `ksm(x, MOD - 2, MOD)`。
- `ksm` 内部是 signed `ll` 乘法；非负指数、底数按模数约束，统一安全条件为 `(p-1)^2 <= LLONG_MAX`（正模数 `p <= 3037000500`）。大模数自行用安全模乘，参照 `MR::mul`、`QR::mul`。
- **完整 u64×u64 范围用 `unsigned __int128`**，signed 溢出是未定义行为；有界的带符号计算仍可用 `i128`，不是所有 `__int128` 都必须 unsigned。
- 补容器形状用 `resize`，不要用会清空调用方数据的 `assign`。多次调用入口重置自己维护的状态，不能默默沿用上次结果。

移植优先顺序：`~/0/Code/`（用户板子，最接近 house 风格）→ `~/0/Lib/atcoder/`（官方算法参考）→ `skip2004-ICPC-Templates/contents/`（第三方参考）。只移植需要的部分，不改来源。

## 5. 自测标准与重要契约

```cpp
// X 自测:<覆盖内容>
#include "../_check_base.hpp"
#include "X.cpp"
int main() { /* 独立参考 + 对拍 + 断言 */ PASSED("X"); }
```

小节子目录中 base 头改为 `../../_check_base.hpp`。兼容层提供常用类型/宏、`ksm`、`MOD`、`rng/rnd`、`CHECK/ok/PASSED`；具体接口读头文件，不在本文复制实现。

- 标准是**独立参考实现 + 大规模随机对拍 + 退化/极端用例 + 性质断言**，失败打印具体输入，不能只测“不崩”。
- 新写/强化 check 后做变异测试：临时把模板故意改错，确认测试变红，再恢复；不要把故障注入留在工作区。
- 可能卡死的测试必须有超时保护。`check.sh` 的外层 timeout 是统一兜底；内部可用子进程限时或 alarm。**fork 隔离本身不是看门狗**，不能把阻塞 waitpid 当超时机制。
- 内存/整数疑点用 ASAN+UBSAN，编译产物写 `tmp/`：

```bash
mkdir -p tmp
# 替换成目标 check 的真实路径；执行时 cd 到该 check 所在目录
# 输出路径用绝对路径，避免运行目录不同找不到二进制
g++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -o "$PWD/tmp/asan_x" templates/图论/支配树.check.cpp
(cd templates/图论 && ../../tmp/asan_x)
```

| 模板 | 必须保留的契约 |
|---|---|
| `字符串/sa.cpp` | 只构造 `sa/rk/height`，没有任意后缀 LCP 查询；需要查询时由使用者另配 RMQ。输入边界等看头注释。 |
| `数据结构/全局平衡二叉树.cpp` | 静态计数不重置，单进程只能建一次；非负点权，中间和须在 int 内；递归需要足够栈。按头注释初始化，多组另开进程。 |
| `计算几何/geo.cpp` | eps 为绝对容差，极小尺度需先缩放；`proj/reflect` 要求非零长直线。没有与所有运算通用的精确尺度分界。 |

具体测试覆盖与已修缺陷以 `.check.cpp` 和提交记录为准，不在这里维护重复历史。

## 6. Typst 介绍与排版

- 介绍首行习惯 `// 介绍: <标题>`；不要再造与模板同级的标题。小标题用单独一行 `*加粗*`，不要用 `===`（会进目录、参与编号）。粗体用单 `*`/`_`，Typst 不认 Markdown 的双星粗体。
- 用户偏好说人话：“4k+1 型素数”“模 4 余 1”，少堆同余符号；说明什么算同一个，不用含混的“有序/无序解数”。公式组织参考 `数学/lagrange.typ`：前提 → 公式 → 怎么做。
- 数学不用 LaTeX 反斜杠命令或花括号分组：写 `sqrt(n)`、`log(n)`、`10^(-5)`，不要 `\sqrt`、`10^{-5}`。相邻变量留空格，如 `a + b i`；多字母名称用文本或 `op("gcd")(a,b)`。
- 正负号写 `plus.minus`，共轭用 `macron(pi)`；`bar(pi)` 是模长竖线，不是共轭。避免不成对的反引号；内容块中 `_ ^ # $ [ ]` 等特殊字符需要转义。
- 两栏表格沿用 `_table.typ` 的字体、间距和网格，不复制另一套参数。整表连标题/介绍一起放进 `#block(breakable: false)[...]`；长表可拆成左右两半并排。
- 主字体在 `xcpc.typ` 固定为 Libertinus Serif / DejaVu Sans Mono + Noto Sans CJK JP，不能依赖系统 fallback。
- Typst 0.15 空数组 `.join()` 返回 none，需要显式给空串。调试小样放项目目录或 `tmp/`，可用 `#panic(...)` 输出值。
- 构建出错按 `templates/xxx.typ:行号` 修复；版面用 PNG 判断，不用 `pdftotext` 判断换行/上标：

```bash
pdftoppm -png -r 150 -f N -l N xcpc.pdf tmp/page
pdfinfo xcpc.pdf | grep '^Pages:'
```

## 7. 环境与 GitHub

- 本机当前 Typst 0.15.1 位于 `/usr/bin/typst`。优先用 `./build.sh`：支持 `TYPST` 覆盖，优先普通 PATH 二进制，snap 载荷只作后备，并设置可写 `tmp/xdg`。
- 如仅有 snap，不直接调用 `/snap/bin/typst` 启动器；PID 隔离下可能无法创建 host scope。载荷为 `/snap/typst/current/bin/typst`，其 `/tmp` 可见性受限，所以统一使用项目 `tmp/`；这不是所有 Typst 的限制。
- 渲染依赖 `@preview/zebraw:0.6.3`；缺包先查 Typst 缓存/下载，不擅自删渲染层。
- 远端是 private 仓库 `Z3O1/xcpc-templates`，origin 使用 SSH；22 端口不可靠时，用已登录 gh 的 HTTPS 凭据助手。**不把 token 写进 URL、日志或仓库文件**：

```bash
gh auth status
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  fetch https://github.com/Z3O1/xcpc-templates.git main:refs/remotes/origin/main
# 检查远端是否领先，确认提交内容后再 push；不自动强推
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  push https://github.com/Z3O1/xcpc-templates.git main:main
```

以上是单次命令配置，不是永久改了 origin/helper。SSH 也可走 `ssh.github.com:443`，但首次信任必须按 GitHub 官方指纹核验，不能未经核验追加 host key。拉取失败时本地 origin/main 可能陈旧，推送后核实远端 SHA；不要凭旧缓存宣称已同步。
