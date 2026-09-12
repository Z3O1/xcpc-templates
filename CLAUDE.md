# CLAUDE.md

算法竞赛模板集:XCPC 板子库。`templates/` 下按章目录放代码模板,`gen.py` 扫描生成 `sections.typ`,
再由 Typst 编译成 `xcpc.pdf`(横向 A4、双栏正文、右上+右下大页码)。

**现状**:8 章 66 个模板(图论 4 / 多项式 13 / 字符串 6 / 数学 5 / 数据结构 7 / 数论 13 / 计算几何 17 / 通用 1),
其中 65 个配了 `X.check.cpp` 自测且全绿(`./check.sh`);PDF 35 页。待补模板见 `TODO.md`,
历史 bug 与踩坑教训见 `BUGS.md`。

## 常用命令

```bash
./build.sh                   # gen.py + typst compile → xcpc.pdf(自动选 typst 二进制、设 XDG_RUNTIME_DIR)
./build.sh --watch           # 改 templates/ 后自动重编译
python3 gen.py               # 只看生成(sections.typ)、不编译
./check.sh                   # 编译并运行 templates/**/*.check.cpp 全部自测
./check.sh -v 数论           # 只跑某章并打印每条断言的输出
```

- 改模板后跑 `./build.sh`,Typst 报错会直接给出 `templates/xxx.typ:行号`,照报错修
- 单个 check 有 60s 墙钟上限(`XCPC_CHECK_TIMEOUT` 可调);卡死会报成 FAIL 而不是挂住整轮
- 视觉抽查:`pdftoppm -png -r 150 -f N -l N xcpc.pdf tmp/x`(注意别输出到 /tmp:snap 版 typst/工具读不到)
- 环境坑(snap typst、XDG_RUNTIME_DIR、沙箱网络)见 `BUGS.md` 末尾「环境」

## 架构:生成链路

```
templates/<章目录>/<模板>.cpp|sh|py  ─┐
templates/<章目录>/<模板>.typ(介绍)  ─┼─→ gen.py ─→ sections.typ ─→ #include 进 xcpc.typ ─→ xcpc.pdf
templates/<章目录>/<模板>.check.cpp  ─┘(被 gen.py 跳过,只给 check.sh 用)
```

- Typst 的 `read()` 不能遍历目录,所以"扫描"全部由 `gen.py` 承担:扫磁盘 → 与 `.manifest.json`
  对账(记忆章节顺序/标题,增删改都是增量) → 输出 `sections.typ` → 确保 `xcpc.typ` 正文是 `#include "sections.typ"`
- `sections.typ` 是生成物,**不要手改**;`.manifest.json` 也是状态文件,只手工改它来调整顺序/标题
- `old_versions/xcpc.typ.bak` 是迁移前原始内联版(39499 B,勿删):首次无 manifest 时用来还原顺序与标题
- 外部依赖:`sections.typ` 首行 `#import "@preview/zebraw:0.6.3"`(代码高亮),报包找不到先看 `~/.cache/typst`

## 模板约定(写模板必须遵守)

- **章节**:`templates/` 下每个目录是一章,顺序由 manifest 决定(新目录追加到最后)
- **标题**:取首行注释 `// 接口名(): 章节 · 标题 (备注)` 里「· 标题」部分;冒号前可为空(允许 `接口()` 这种写法),
  尾部括号备注会被**配平剥离**(能处理嵌套括号);无注释则回退文件名。标题存 manifest,与注释解耦
- **改标题**:gen.py 只给新条目从注释取标题,已有条目以 manifest 为准 → 想换标题先从 `.manifest.json` 摘掉该条再 `python3 gen.py`
- **隐藏**:文件头注释区含 `// hide` 或 `// 隐藏` → 整块不进 PDF(`字符串/sam.cpp` 现处于此状态),去掉标记即恢复
- **介绍(推荐)**:同名 `.typ` 渲染在该代码前面,不进 manifest、不占条目;**判断规则是"存在同基底名非 `.typ` 文件"**——
  没有同名代码的 `.typ` 是独立可渲染模板(`通用/常数速查表.typ` 就是这种)
- **头注释区**:第一个空行之前的内容不渲染。
  - **头注释后必须留一个空行**,否则 `codelines()` 取不到正文、PDF 里渲染成**空代码块**(实测踩过 11 个文件)
  - 头部注释之后、正文里的 `//` 行与空行都会被 `readcode` 跳过,所以正文中间可以随意写注释
- **不写 `main`**:模板是给读者抄的片段。自带 main 的模板会把输入输出调度整块渲染进 PDF(已有 3 个
  这样的模板被抽象掉了),新模板请直接不写
- **不写 `#include` / `using namespace std;`**:片段默认读者有 base header(`多项式/多项式复合.cpp` 历史上带着这两个,已清)
- 介绍 `.typ` 里不要再写与 `== 标题` 同级的标题;可用 `*文本*`/`_文本_` 加粗、`` `代码` `` 等 Typst 语法
  (注意 Typst **不认 `**双星粗体**`**,会渲染成"星号吞掉")

## 代码风格(house 宏,写模板必看)

模板是**竞赛代码片段**,不是能单独编译的单元:文件内不定义宏/类型,沿用读者 base header 的那套命名。
拿不准就看 `skip2004-ICPC-Templates/` 与 `~/0/Code/` 里既有写法:

- `For(i, l, r)` / `rFor(i, r, l)` / `ForD` 代替 `for`;`vect<T>` 代替 `vector<T>`
- `ll` = `long long`、`db` = 浮点(常用 `long double`)、`mint` = 模数类(带 `.inv()`)、`poly` = 多项式/`vector<mint>`、`ksm` = 快速幂
- `vect<T>` 是带 `+=` 的自定义容器(**没有 `push_back`**,追加用 `p += x`);要标准语义就用 `std::vector`
- 计算几何统一 `struct p2`(点)、`seg`(线段)、`line`;`eps`/`cmp`/`sgn`/`cross`/`det` 常规
- **`ksm` 的模数限制**:它内部按 `ll` 相乘,模数 > 2^32 时 `a * a` 溢出算错 → 大模数场合**自带模乘**
  (见 `数论/Miller-Rabin.cpp` 的 `MR::mul`、`数论/二次剩余.cpp` 的 `QR::mul`)
- **补形状用 `resize` 不要用 `assign`**:`assign` 会清空调用方传进来的数据(踩过,见 `BUGS.md`)

## 模板自测(check,仿 skip2004)

每个模板配一份同目录 `X.check.cpp`,写法固定:

```cpp
// X 自测:<测什么>
#include "../_check_base.hpp"   // house 宏 + CHECK/PASSED/ok/rnd/rng
#include "X.cpp"                // 直接 include 被测模板本体
int main() { ...; PASSED("X"); }
```

- `_check_base.hpp` 提供 `For/rFor/ForD/vect/ll/db/mint/ksm` 与断言宏;`mint` 的成员是 `int v`(还有 union 别名 `x`)。
  注意它的 `ksm` 同样有上面的模数限制
- **标准是"内置独立参考实现 + 大规模随机对拍 + 退化/极端用例 + 性质断言"**,失败要打印具体输入。
  别写成"跑一遍不崩"——用变异测试(故意改错模板,看 check 是否变红)验证 check 有牙
- 遇到会死循环的用例:用 `alarm`+`sigsetjmp` 或 `fork`+`waitpid` 看门狗,让它报 FAIL 而不是挂住整轮
- **模板本体出 bug 时先报告再改**:改完必须跑该章 check 与 `./build.sh`

## 移植来源(补模板时优先搬,别从零写)

- **`~/0/Code/`** —— 用户自己的板子(接近 house 风格):`一般图最大匹配.cpp`/`上下界.cpp`/`Hopcroft–Karp.cpp`/`FHQ.cpp`/
  `FFT.cpp`/`poly 1..3.cpp`/`Bostan-Mori.cpp`/`杜教筛.cpp`/`Gause.cpp`/`det.cpp`/`Matrix.cpp`/`BigInt.cpp`/`Hash.cpp`/`Geo.cpp` 等
- **`~/0/Lib/atcoder/`** —— 官方 atcoder 库(`maxflow`/`mincostflow`/`scc`/`twosat`/`segtree`/`lazysegtree`/`convolution`/`math`/`string`),
  当"该算法该怎么做"的权威参考
- **`skip2004-ICPC-Templates/contents/<分类>/`** —— 第三方完整模板集(173 个文件;`numbertheory`/`geometry`/`math`/`string` 等,
  还带 `tables/Addon.tex`、`tables/dn2.md` 这类常数表)。**只当参考来源,目录本身不要改、也不入库**

## Typst 陷阱(本环境实测)

- `$...$` 数学模式里**不要用带反斜杠的多字母符号**:`\log`、`\sqrt`、`\alpha` 报 `unknown variable`;
  多字母标识符(`len`、`mcf`)会被拆成变量序列 → 复杂度等写成反引号文本或 Unicode 文本,数学里只留单字母与基本运算
- **LaTeX 源不能直接抄**:`\times`、`\cdot`、`^{0}`、`\multicolumn` 这些要手工转成 Typst 写法
- 强调用**单** `*` 或 `_`;`**` 无效(编译 warning "no text within stars",排版乱)
- 反引号必须成对,不成对报 `unclosed raw text`(报错跨度可能跨行,很迷惑);内容块 `[...]` 里的 `_`、`^`、`#`、`$` 要转义
- **表格排版**(两栏很窄):`#table(inset: 1.5pt, row-gutter: 0pt, column-gutter: 3pt)` + `#set text(size: 6.5pt)` +
  `#set par(leading: 0.35em)` 才够紧;想让整张表留在同一页用 `#block(breakable: false)[...]`(会整体推到下一页),
  长表可以改成**左右两半并排**(6 列布局)把行数减半
- 标题要 Typst 安全:标题会原样进 `== 标题`,带 `(`/`)` 会被当函数调用、带 `*`/反引号会被当标记 → 标题写短名词
- Typst 0.15:`().join("\n")` 返回 `none`(空数组 join),需要显式 `if ... else { "" }`
- 调试单点小样:写 `_t.typ` 到**项目目录**(snap 版读不到 /tmp),用 `#panic("...")` 把值打到 stderr
  —— 别用 `pdftotext` 调字符串(页面折行会伪造换行;上标还会被压平,如 `10^6` 显示成 `106`)

## 其它

- 已是 git 仓库(远端 `git@github.com:Z3O1/xcpc-templates.git`,private)。推送若报 `Connection closed ... port 22`,
  是网络层断 SSH,改用 gh token 走 HTTPS:`git push https://x-access-token:$TOK@github.com/Z3O1/xcpc-templates.git main:main`
  (`$TOK` 取 `~/.config/gh/hosts.yml` 里的 `oauth_token`)
- `tmp/` 是编译夹具与调试产物(已 gitignore,里面 `tmp/数论/` 有早期测试驱动、`tmp/xdg/` 是 typst 运行目录)
- `a.cpp`、`in2.txt` 等是无关残留文件,不动也不删
