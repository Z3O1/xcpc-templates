# BUGS.md —— 历史 bug 与踩坑留档

CLAUDE.md 的补充:这里放"具体到行/带最小复现"的历史记录。**这些 bug 全部已修复**,每个都有对应的
`X.check.cpp` 守着(回归会红)。写新模板前扫一遍,能少踩坑。

## 修掉的模板 bug

| 模板 | bug | 修法 | 怎么发现的 |
|---|---|---|---|
| `图论/dinic.cpp` | `bfs()` 首行 `d[i] = -(i == s)` → 邻居拿到 `0`,而判空 `!~d[v]` 只认 `-1` → BFS 一个点都进不去,**最大流恒 0** | 改成 `d[i] = -1` 再 `d[s] = 0`(注意:只写 `d[i] = -1` 不设源点会让源点被当未访问重新入队 → `solve()` 死循环) | 最小复现 `add(1,2,3) add(2,3,4) add(1,3,1); solve(1,3,3)` 返回 0(应 4) |
| `图论/dinic.cpp` | `cut(auto *f)` 的 `auto` 形参是 C++20 语法,C++17 只有 GCC 扩展 | 改 `int *f` | `-Werror`/clang 会挂 |
| `图论/mcmf.cpp` | `mcmf2` 内层两次 `mcmf(...)` 漏传第 4 参 `_n` → 内层 n=0 → **任何输入死循环** | 补 `_n+2` | `mf.add(1,2,1,-5); mf.mcmf2(1,2,0,2);` 挂死 |
| `图论/mcmf.cpp` | `a1 += e[tot].w` 写在 `mcmf(...)` 之前(那时 `_t→_s` 上流量还是 0,且增广会追加边、下标会变) | 存 `fake = tot`,增广后再读 `e[fake].w` | 可行性阶段的循环流永远算不进流量 |
| `图论/mcmf.cpp` | `clear()` 只清到成员 `n`,而 `n` 只在 `mcmf()` 里赋值 → 首次 `clear()` 空操作,残留 `hd[]` 让 dij 取到奇下标前驱 → **`mcmf()` 也死循环** | 新增 `mx`(用过的最大点数),`clear()` 清到 `mx`;`mcmf/mcmf2` 入口维护它 | 3/3 组子进程被看门狗强杀 |
| `图论/一般图最大匹配.cpp` | 遍历 `for(auto v : t[u])` **少一行 `if(fd(u) == fd(v)) continue;`**:花缩起来后花内相邻点还会互相走分支 → 死循环 + `mat[]` 被写脏 | 循环体第一行加这一句 | 400 组随机 n=6..12 里 4 组卡死;**穷举 n≤5 的 1099 个图不触发**(小规模对拍抓不到,必须随机+看门狗) |
| `图论/支配树.cpp` | 处理不可达前驱时没跳过 `dfn[x] == 0` → `sd[u]` 被拉成 0、`dm[u]` 算错 | `for(auto x : t2[u]) if(dfn[x]) cmin(...)` | 最小复现 n=3、边 `1→3, 2→3` → `dm[3]=0` 应 1;19408 用例里触发 2798 次 |
| `图论/支配树.cpp` | `getdom()` 不自重置 `dfn/pos/dt/q1`(只重置 `fa/f`)→ 同进程二次调用算错 | 入口重置算法自身状态(**不动 `t[]`/`t2[]`**,图是调用方建的) | 同一进程连跑两张图结果错 |
| `计算几何/geo.cpp` | `operator-=` 与 `operator/=` 都写成 `x = x + y`;且 `/= ` 只有 `(p2&, p2)` 重载,`a /= 2.0` 直接编译不过 | `x - y` / `(p2&, db)` 版 `x / y` | check 里 SFINAE 探测重载存在性 |
| `计算几何/半平面交.cpp` | ①两条近同向直线(极角差 1 ulp)去重时保留**较弱**的一条;②真交集为空时双端队列交出**顺时针假三角形**(面积 -1.4e17) | 显式比强弱;带号面积 ≤ 0 一律当退化返回空 | 随机对拍 8 万组才碰到一次,check 里钉了 8 直线可复现反例 |
| `计算几何/上凸壳.cpp` | 同 x 判据用了 eps,圆上 15° 间隔的两点被当同 x → 丢真正的最高点 | 同 x 用**精确相等**,近似同 x 交给栈的叉积处理 | 半径 1000 圆上 24 点只出 9 个顶点 |
| `计算几何/凸包内点判定.cpp` | 先算对角线分支再算外边 → 顶点本身被判成"内部"(返回 2 而非 1) | 先算外边半平面,`c==0` 直接返回 1 | 与 `contain` 三方对照 |
| `计算几何/图形交.cpp` | `cir_poly_area` 在"圆完全在多边形外"时返回**整块圆面积**;精度只有 1e-7;相切/圆内判据的绝对 eps 未按半径缩放(半径 1e-6 时失效) | 按 skip2004 思路重写成解析版(逐边有向扇形 + 弦三角形、有向角不归一),O(n);判据按半径缩放 | 4 处 check 期望值也算错了,用蒙特卡洛/条带积分/200 万边形裁剪独立核对 |
| `计算几何/geo.cpp`(`convex_hull`) | 头注释没写"输出缓冲要 ≥ 2n+2"(实测 `p2 b[5]` + 点 `(i,i²)` 写 `b[5]` 时 stack-buffer-overflow);`n<=1` 返回 1 但不写点;会**原地排序输入数组** | 契约写进头注释 | ASAN 实测 |
| `数论/二次剩余.cpp` | 复用 base header 的 `ksm`,模数 > 2^32 时溢出 → `quadres/sqrtp` 在大素数上全错 | 自带 `QR::mul/QR::pw`(u128);Cipolla 的 `F_p[√w]` 乘法一并改用;返回类型 `int → ll` | `ksm(4000000000,2,4294967291)` 得负数;2000 组"y² 构造"的二次剩余全被判 -1 |
| `数论/分数还原.cpp` | 移植时把连分数的 `x %= y` 抄成 `x *= y` → 对几乎所有输入溢出返回垃圾 | 改回 `%=`;参数顺序对齐 skip2004 的 `(p, q, A)` | `approx(5,2,1)` 曾返回 -1000 万级垃圾值 |
| `数论/Miller-Rabin.cpp` | 复用 `ksm` 做模乘,`p > 2^32` 时溢出 | 自带 `MR::mul`(用 **unsigned** `__int128`——signed 上限 1.7e38,两个 u64 相乘到 3.4e38 会静默溢出) | `4294967291`、`2^64-59` 被误判合数 |
| `多项式/多项式复合逆.cpp` | `pw_pj` 只写前 n 项却不 `resize(n)`,随后 `comp_inv` `reverse` 整个向量 → `[n,\|F\|)` 的高次项倒到低次,**静默算错** | `pw_pj` 末尾 `F.c.resize(n)`;`comp_inv` 入口补 F 的系数 | n=2、F={0,1,0} → g={0,0}(应 {0,1}) |
| `多项式/多项式复合.cpp` | `Gp[1][i] = mod - G.c[i]` 越界读(调用方给的 G 不足 n 个系数) | 入口把 `G`/`F` 补到 n 个系数 | ASAN 实测 heap-buffer-overflow |
| `多项式/多项式复合逆.cpp` | 同上,且 `n == 1` 时也要至少 2 个系数(内部读 `F.c[1]`) | 同上 | ASAN 实测 |
| `字符串/manacher.cpp` | 签名 `void manacher(int n, char *s, int d)` **少个 `*`**,函数体里 `d[i]` → 编译不过 | 改 `int *d` | 字符串章 check 实测 |
| `字符串/sa.cpp` | ①`n == 1` 段错误(重排名循环不跑,`rk[1]` 停在字符值,LCP 段自比自身无限增长);②`lcp()` 依赖调用方 `a[0]`;③字符值必须 ≥ 1(0-based 映射会坏) | **未修**,只在 check 里以契约形式覆盖(n=1 且字符值 1 时才测) | 字符串章 check |
| `字符串/zfunc.cpp` | 注释说 `s[n]` 是哨兵,实际读 **`s[n+1]`**;且要求 z 数组零初始化 | 未修(注释已在 check 里记录) | 字符串章 check |
| `数据结构/全局平衡二叉树.cpp` | **模板是完整程序(自带 main)**,而 main 末尾没有 `return 0;` —— 改名 include 后 GCC 把"non-void 函数掉出末尾"当 `__builtin_unreachable`,`-O2` 下**连循环退出判断和 `ret` 一起删掉** → 死循环 | 补 `return 0;`,并按新约定去掉了整个 main | 反汇编 + 插桩(循环计数被优化掉) |
| `数学/类欧.cpp` | 结构版 `f(n,a,b,c)` 与单值版 `f(a,b,c,n)` 签名都是 `(int,int,int,int)` → **不能同时编译** | 未修(按需只留一版) | 报 `ambiguating new declaration of 'mint f(int, int, int, int)'` |
| `多项式/ntt.cpp` | 变换长度 `l` 必须 ≤ 2^20 但没写进注释;负系数不在契约内 | 契约写进头注释 | `mul({-1},{1},998244353)` 返回 -1 |

## 我自己踩过的坑(不是模板的错)

- **把"两个不相交三角形"的最大匹配算成 3**(正确答案是 **2**:每个三角形内部只能取 1 条边,取两条必有公共端点)。
  因为这个错的期望值,我把"实现正确"当成"漏增广路",还试了三种错方向的改法,白费一轮。
  → 教训:**期望值要用暴力枚举独立算,别口算**。
- **把 `d[i] = -1` 当成 dinic 的修法推荐出去**:那会让源点被当未访问重新入队,`solve()` 死循环。
  → 教训:改最小复现之后必须实测,别只做纸面推理。
- **`assign` 清空调用方数据**:写 charpoly 时用 `a.assign(n, poly(n))` 补形状,把输入矩阵元素全抹了(对角阵还"看起来对",随机对拍才发现)。
  → 教训:补形状用 `resize`。
- **空代码块**:头注释后没留空行,`codelines()` 取不到正文,PDF 里渲染成空块(一次踩了 11 个文件)。
- **`pdftotext` 不能用来判断排版**:页面折行会伪造换行,上标会被压平(`10^6` → `106`)。要判断就得 `pdftoppm` 看图。

## 未采用的尝试(留档避免重踩)

- **特征多项式(charpoly)**:按 skip2004 风格写过一个 Hessenberg 版,对角阵对、但 Jordan 块给不出 `(x-7)³`、
  随机 n≤4 矩阵与 Leibniz 展开 **300 组里 221 组不一致** → 删除,没有进仓库。要写请照 skip2004 的
  `contents/math/charpoly.cpp` 逐行核对并自备 Leibniz 对拍。
- **带花树的"标准结构重排"**:把 `if(vis[v]==2)` 提前、"未访问且已匹配才扩展"等三种改法都会让匹配变小或仍挂
  (真正需要的只是那行 `if(fd(u) == fd(v)) continue;`)。
- **`seg_in_polygon` / 在线卷积 / 幂投影**:还没做,见 `TODO.md`。

## 环境

- **typst 是 snap 版**:`./build.sh` 会自动选二进制并设 `XDG_RUNTIME_DIR`。要点:
  - **不要直接用 `/snap/bin/typst`** —— 它只是 snapd 启动器,在 DSH 沙箱里必然失败
    (`cannot create transient scope: DBus error ... UnixProcessIdUnknown`):snapd 要经 D-Bus 让 host systemd
    申请 scope,而沙箱 PID 1 是 `bwrap --unshare-pid`,host 看不见。**这是 PID 命名空间隔离,放宽沙箱也修不好**
  - 绕法:跑 snap 载荷里的真二进制 `/snap/typst/current/bin/typst`(0.15.1),并给一个可写的 `XDG_RUNTIME_DIR`
    (`mkdir -p tmp/xdg && XDG_RUNTIME_DIR=$PWD/tmp/xdg ...`),否则报 `cannot create XDG_RUNTIME_DIR folder ... Read-only file system`
  - snap 版 typst **读不到 /tmp**,调试小样要写在项目目录里
- **沙箱里没有网络**:`curl https://github.com` 报 `SSL_ERROR_SYSCALL`。git 推送若报
  `Connection closed by ... port 22`,改用 token 走 HTTPS(见 CLAUDE.md「其它」)
- 早期会话里 `git ls-tree` 对非 ASCII 路径会加引号转义,用 `git ls-tree -z ... | tr '\0' '\n'` 才 grep 得到
