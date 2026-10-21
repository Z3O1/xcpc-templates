# 开发操作参考

按需查阅的调试与连接命令；项目约束见 [AGENTS.md](AGENTS.md)。以下命令均从仓库根目录运行。

## 内存与整数检查

内存或整数计算有疑点时，用 ASAN+UBSAN 编译相应自测。产物放 `tmp/`，运行时切到 check 所在目录；以下以支配树为例：

```bash
mkdir -p tmp
g++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -o "$PWD/tmp/asan_x" templates/图论/支配树.check.cpp
(cd templates/图论 && timeout 60s ../../tmp/asan_x)
```

## 样式快速预览

共享样式在 `style.typ`，正式入口 `xcpc.typ` 与小样 `tests/style-preview.typ` 都使用它。调整样式可先运行：

```bash
./preview.sh          # 持续监听,Ctrl+C 停止
./preview.sh --once   # 只编译一次
./preview.sh --open   # 持续监听并打开 PDF 阅读器
```

保持命令运行，修改并保存 `style.typ` 后会自动更新 `tmp/style-preview.pdf`；第 1 页是目录，后面是正文示例。脚本自动创建输出目录、设置项目根，不生成章节，也不覆盖正式 PDF；支持 `TYPST` 指定编译器，其他参数原样传给 Typst。

`style.typ` 只定义样式，不应单独编译成 `style.pdf`；预览统一使用脚本，避免 PDF 散落在根目录或 `tests/` 中。

常调的配置都在 `style.typ` 中：

- `set text(..., size: ...)`：正文字号。
- `show raw: set text(..., size: ...)`：代码字号。
- `set par(leading: ..., spacing: ...)`：文字行间空隙、段落间距。
- `show heading: set block(above: ..., below: ...)`：标题前后间距。
- `zebraw-init.with(inset: (top: ..., bottom: ...))`：每行代码上下留白，不是纸张页边距。

确认后用 `./build.sh` 更新整本 `xcpc.pdf`。

## 版面检查

版面检查用 PNG，不用 `pdftotext` 判断换行或上标。将 `N` 改为要查看的页码：

```bash
mkdir -p tmp
N=1
pdftoppm -png -r 150 -f "$N" -l "$N" xcpc.pdf tmp/page
pdfinfo xcpc.pdf | grep '^Pages:'
```

一次性 Typst 调试小样放 `tmp/`（需要时放项目目录），可用 `#panic(...)` 输出值；长期共用的版式示例见 `tests/style-preview.typ`。
