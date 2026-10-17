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

## GitHub 连接排障

远端为 private 仓库 `Z3O1/xcpc-templates`，origin 使用 SSH。22 端口不可靠时，可用已登录 gh 的 HTTPS 凭据助手；token 不写进 URL、日志或仓库文件，以下配置只对单次命令生效，不永久修改 origin/helper。

先检查登录并获取远端状态：

```bash
gh auth status
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  fetch https://github.com/Z3O1/xcpc-templates.git main:refs/remotes/origin/main
```

拉取失败时，`origin/main` 可能陈旧；不能据此认定已同步。检查远端是否领先并确认提交内容后，才手动推送，不自动强推：

```bash
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  push https://github.com/Z3O1/xcpc-templates.git main:main
```

推送后查询实际远端 SHA，与本地 `git rev-parse main` 比较，不凭旧缓存宣称成功：

```bash
git -c credential.helper= -c credential.helper='!gh auth git-credential' \
  ls-remote https://github.com/Z3O1/xcpc-templates.git refs/heads/main
```

SSH 也可走 `ssh.github.com:443`，但首次信任须按 GitHub 官方指纹核验，不未经核验追加 host key。
