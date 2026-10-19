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

## 版面检查

版面检查用 PNG，不用 `pdftotext` 判断换行或上标。将 `N` 改为要查看的页码：

```bash
mkdir -p tmp
N=1
pdftoppm -png -r 150 -f "$N" -l "$N" xcpc.pdf tmp/page
pdfinfo xcpc.pdf | grep '^Pages:'
```

Typst 调试小样放 `tmp/`（需要时放项目目录），可用 `#panic(...)` 输出值。
