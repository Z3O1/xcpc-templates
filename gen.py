#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 templates/ 自动生成 sections.typ, 并让 xcpc.typ 引用它。

Typst 的 read() 只能按显式路径读单文件(无 glob/目录遍历), 所以"遍历"
由本脚本承担: 每次运行时扫描 templates/, 子目录=一级标题, 目录内文件=代码块。

自动识别 Typst 模板: 后缀为 .typ 的模板以 #include 引入(其中 markup/code
会被 Typst 真正执行, 适合嵌可渲染的 Typst 片段); 其余后缀(.cpp/.sh/...)
作为代码块 raw 原样显示。.typ 模板文件头的 // 注释不渲染, 且其内容不要再
写与生成文件里 "== 标题" 同级的标题。

增删处理(通过 .manifest.json 状态文件记忆):
- 新增模板: 自动追加到所属章节末尾, 标题取文件首行注释中的 "· 标题" 部分
     // xxx(): 通用 · 快读 (备注写在括号里会被去掉)
  无注释则用文件名(去扩展名)。
- 删除模板: 从 sections.typ 移除(即 PDF 不再包含), 但记忆保留; 以后加回来
  时按原位置原位出现, 标题不变。
- 删除章节目录: 整章移除; 目录重建后整章恢复。

介绍文件: 每份代码可附同名 .typ 介绍(如 ntt.cpp 配 ntt.typ)。判断规则:
目录内 x.typ 与某个非 .typ 文件(x.cpp 等)同名时, 该 .typ 视为那
份代码的介绍——不单独成块、不进 manifest, 渲染在标题/简介之后、代
码之前; 介绍了同名代码后 .typ 自动退回"独立可渲染模板"的语义。

标题层级(靠目录结构): 章内文件是二级条目(==);章内的子目录是一个二级小节,
子目录里的文件是它的三级子条(===, 编号如 2.5.1)。子目录里与目录同名的 .typ
(如 四边形不等式/四边形不等式.typ)是该小节的正文, 不算子条。

顺序/标题第一次生成时取自 xcpc.typ.bak(原始内联版), 之后以 manifest 为准。

用法: python3 gen.py && typst compile xcpc.typ   (或直接 ./build.sh)
"""
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).parent
TPL = ROOT / 'templates'
SEC = ROOT / 'sections.typ'
MANIFEST = ROOT / '.manifest.json'
MAIN = ROOT / 'xcpc.typ'
BAK = ROOT / 'xcpc.typ.bak'

LANG = {'cpp': 'cpp', 'sh': 'bash', 'py': 'python', 'rs': 'rust', 'typ': 'typst'}
FENCE = re.compile(r'^```(\w+)?\s*$')


def lang_of(p: pathlib.Path) -> str:
    return LANG.get(p.suffix.lstrip('.'), p.suffix.lstrip('.'))


def norm(s: str) -> str:
    return '\n'.join(l for l in s.split('\n') if l.strip())


def codelines(p: pathlib.Path) -> list[str]:
    """readcode 的等价物: 去掉文件头注释(第一个空行之前的内容)。"""
    ls = p.read_text().splitlines()
    i = 0
    while i < len(ls) and ls[i] != '':
        i += 1
    return ls[i + 1:]


def strip_trailing_parens(t: str) -> str:
    """去掉标题尾部的括号备注,支持嵌套括号(如 "... (单调栈 O(n log n))")"""
    t = t.rstrip()
    while t.endswith(')') or t.endswith('）'):
        close, open_ = (')', '(') if t.endswith(')') else ('）', '（')
        depth = 0
        k = len(t) - 1
        while k >= 0:
            if t[k] == close:
                depth += 1
            elif t[k] == open_:
                depth -= 1
                if depth == 0:
                    break
            k -= 1
        if k <= 0 or depth != 0:
            break          # 括号不配平,别硬剥
        t = t[:k].rstrip()
    return t


def title_from_comment(p: pathlib.Path) -> str:
    """从首行注释取标题: // 符号: 章节 · 标题 (括号备注)"""
    for line in p.read_text().splitlines()[:1]:
        # 注意:接口名后面常常直接跟一对空括号(如 `// poly_inv(): 多项式 · 求逆 (牛顿迭代)`),
        # 所以冒号前允许为空;备注可能带嵌套括号(如 "(有向直线左侧求交 O(n log n))"),
        # 所以剥离尾部括号要用配平而不是 [^)]*。
        m = re.match(r'^[#/]+\s*[^:]*:\s*[^·]*·\s*(.+)$', line)
        if m:
            t = strip_trailing_parens(m.group(1)).strip()
            if t:
                return t
        break
    return p.stem


def header_lines(p: pathlib.Path) -> list[str]:
    """文件头注释区(第一个空行之前的所有行)"""
    out = []
    for l in p.read_text().splitlines():
        if l.strip() == '':
            break
        out.append(l)
    return out


def hidden_from_header(p: pathlib.Path) -> bool:
    """打印隐藏标记: 头注释区含 // hide(或 // 隐藏) 则该模板不进入 PDF"""
    return any(re.match(r'^[/#]+\s*(hide|隐藏)\b', l.strip()) for l in header_lines(p))


# ---------- 1) 扫描磁盘 ----------
# 章 = templates/ 下的一级目录; 章内文件 = 二级条目; 章内的子目录 = 一个二级小节,
# 子目录里的文件 = 它的三级子条。子目录里与目录同名的 .typ(如 四边形不等式/四边形不等式.typ)
# 是这个小节的正文, 不算子条。
files = {}   # 章名 -> {相对路径(可含子目录): Path}, 不含介绍 .typ
intros = {}  # 章名 -> {相对路径去扩展名: Path}: 同名 .typ 介绍(渲染在代码前)
for d in sorted(p for p in TPL.iterdir() if p.is_dir()):
    fs = {p.relative_to(d).as_posix(): p for p in sorted(d.rglob('*'))
          if p.is_file() and not p.name.startswith('.')
          and not p.name.endswith('.check.cpp')}   # check 自测不入 PDF
    if not fs:
        continue
    code = {n[:-4] for n in fs if not n.endswith('.typ')}   # 有同名代码的基底路径
    files[d.name], intros[d.name] = {}, {}
    for n, p in fs.items():
        if n.endswith('.typ') and n[:-4] in code:
            intros[d.name][n[:-4]] = n   # 同目录同名 .typ = 代码介绍
        else:
            files[d.name][n] = p


# ---------- 2) 状态: manifest ----------
def parse_backup() -> dict:
    """{'通用/快读.cpp': ['通用', '快读'], ...} —— 备份里可匹配到的文件"""
    lines = BAK.read_text().splitlines()
    by_norm = {}
    for d, fs in files.items():
        for name, p in fs.items():
            by_norm.setdefault(norm('\n'.join(codelines(p))), []).append((d, name))
    out, cur_sec, cur_title = {}, None, None
    i = 0
    while i < len(lines):
        m = re.match(r'^= (.+)$', lines[i])
        if m:
            cur_sec, cur_title = m.group(1).strip(), None
            i += 1
            continue
        m = re.match(r'^== (.+)$', lines[i])
        if m:
            cur_title = m.group(1).strip()
            i += 1
            continue
        if FENCE.match(lines[i]) and cur_sec:
            j = i + 1
            while j < len(lines) and not FENCE.match(lines[j]):
                j += 1
            hits = by_norm.get(norm('\n'.join(lines[i + 1:j])), [])
            if len(hits) == 1:
                d, name = hits[0]
                out[f'{d}/{name}'] = (d, cur_title)
            elif not hits:
                print(f'! 原文档代码块无法匹配模板: {cur_title}')
            else:
                print(f'! 代码块匹配到多个模板, 跳过: {cur_title} {hits}')
            i = j + 1
            continue
        i += 1
    return out


def parse_state() -> dict:
    """{'通用/快读.cpp': ['通用', '快读'], ...} —— 现有 sections.typ 的状态(标题)"""
    secs, cur_sec, cur_title = {}, None, None
    for line in SEC.read_text().splitlines():
        m = re.match(r'^= (.+)$', line)
        if m:
            cur_sec, cur_title = m.group(1).strip(), None
            continue
        m = re.match(r'^== (.+)$', line)
        if m:
            cur_title = m.group(1).strip()
            continue
        m = re.match(r'^#raw\(readcode\("templates/[^/]+/(.+)"\), lang: "(\w+)"\)$', line)
        if m and cur_sec:
            secs[f'{cur_sec}/{m.group(1)}'] = (cur_sec, cur_title)
            continue
        m = re.match(r'^#include "templates/[^/]+/(.+)"$', line)
        if m and cur_sec:
            name = m.group(1)
            # 同名 .typ 是代码的介绍, 不单独成条目(标题随代码块)
            if name.endswith('.typ') and intros.get(cur_sec, {}).get(name[:-4]):
                continue
            secs[f'{cur_sec}/{name}'] = (cur_sec, cur_title)
    return secs


if MANIFEST.exists():
    data = json.loads(MANIFEST.read_text())
    print('使用 .manifest.json 作为状态来源')
else:
    if BAK.exists():
        state = parse_backup()
        print('首次运行: 从 xcpc.typ.bak 还原章节顺序与标题')
    elif SEC.exists():
        state = parse_state()
        print('无 xcpc.typ.bak: 从现有 sections.typ 还原标题')
    else:
        state = {}
        print('无任何状态来源: 将按文件名排序')
    data = {'sections': [], 'entries': []}
    # 先按备份/状态顺序 + 其章节名
    for key, (d, title) in state.items():
        if d not in data['sections']:
            data['sections'].append(d)
        data['entries'].append({'key': key, 'title': title})
    # 磁盘上有但状态里没有的(如 fgcd/rmq), 设为 missing=False 并按目录追加
    state_keys = set(state)
    for d in sorted(files):
        if d not in data['sections']:
            data['sections'].append(d)
        for name in files[d]:
            key = f'{d}/{name}'
            if key not in state_keys:
                data['entries'].append({'key': key, 'title': title_from_comment(files[d][name])})
    # 状态里有但磁盘上已删除的 -> 打上 missing(首次生成时提醒)
    for e in data['entries']:
        d, name = e['key'].split('/', 1)
        e['missing'] = d not in files or name not in files[d]
    for e in data['entries']:
        if e['missing']:
            print(f'- 已移除: {e["key"]}')

# ---------- 3) 与磁盘对账 ----------
# 曾作为独立模板的 .typ 如今有了同名代码 -> 身份转为介绍, 从条目里移除
intro_keys = {f'{d}/{stem}.typ' for d, m in intros.items() for stem in m}
old = [e['key'] for e in data['entries'] if e['key'] in intro_keys]
if old:
    data['entries'] = [e for e in data['entries'] if e['key'] not in intro_keys]
    for k in old:
        print(f'- 已转为介绍: {k} (不再独立显示)')

entry_keys = {e['key'] for e in data['entries']}
# 磁盘新增 -> 追加到其章节末尾
for d in sorted(files):
    if d not in data['sections']:
        data['sections'].append(d)
    for name in files[d]:
        key = f'{d}/{name}'
        if key not in entry_keys:
            data['entries'].append({'key': key, 'title': title_from_comment(files[d][name]),
                                    'missing': False})
            print(f'+ 新增: {key} (标题: {data["entries"][-1]["title"]})')

# 打标记 & 报告状态迁移
new_entries = []
for e in data['entries']:
    d, name = e['key'].split('/', 1)
    on_disk = d in files and name in files[d]
    was_missing = e.get('missing', False)
    if on_disk and was_missing:
        print(f'^ 已恢复: {e["key"]} -> 原位 {d} / {e["title"] or name}')
    e['missing'] = not on_disk
    was_hidden = e.get('hidden', False)
    hidden_now = False
    if on_disk:
        hidden_now = hidden_from_header(files[d][name])
        if hidden_now and not was_hidden:
            print(f'- 已隐藏: {e["key"]} (// hide 生效)')
        elif not hidden_now and was_hidden:
            print(f'^ 已恢复显示: {e["key"]} (去掉 // hide 标记)')
    e['hidden'] = hidden_now
    new_entries.append(e)
data['entries'] = new_entries
MANIFEST.write_text(json.dumps(data, ensure_ascii=False, indent=1) + '\n')

# ---------- 4) 输出 sections.typ ----------
out = [
    '// 由 gen.py 自动生成 —— 请勿手动修改; 改动 templates/ 后运行: ./build.sh',
    '// 读取模板文件: 去掉文件头的说明注释(第一个空行之前)与末尾空行',
    '#import "@preview/zebraw:0.6.3": zebraw',
    '#let readcode(path) = {',
    '  let lines = read(path).split("\\n")',
    '  let i = 0',
    '  // 跳过头部注释行(// 开头)与空行; 无注释的文件从第一行代码开始',
    '  while i < lines.len() and (lines.at(i) == "" or lines.at(i).starts-with("//")) {',
    '    i += 1',
    '  }',
    '  while i < lines.len() and lines.at(i) == "" {',
    '    i += 1',
    '  }',
    '  let body0 = lines.slice(i)',
    '  let body = if body0.len() > 0 and body0.last() == "" {',
    '    body0.slice(0, body0.len() - 1)',
    '  } else {',
    '    body0',
    '  }',
    '  // 空数组 join() 返回 none(Typst 0.15), 这里显式给空串',
    '  if body.len() > 0 { body.join("\\n") } else { "" }',
    '}',
    '',
]
total = 0
for d in data['sections']:
    items = [e for e in data['entries']
             if not e['missing'] and not e.get('hidden')
             and e['key'].startswith(d + '/')]
    if not items:
        continue
    # <小节>/<小节>.typ 是小节正文, 其余 <小节>/xxx 是它的子条
    def group_body(rel):
        parts = rel.split('/')
        return len(parts) == 2 and parts[1] == parts[0] + '.typ'

    groups = {}
    for e in items:
        rel = e['key'].split('/', 1)[1]
        if '/' in rel and not group_body(rel):
            groups.setdefault(rel.split('/', 1)[0], []).append(e)

    def emit(e, lvl):
        if e['title'] is not None:
            out.append('')
            out.append(f'{"=" * lvl} {e["title"]}')
        name = e['key'].split('/', 1)[1]
        path = files[d][name]
        # 同目录同名 .typ 介绍: 排在正文代码前面
        intro = intros.get(d, {}).get(name.rsplit('.', 1)[0])
        if intro:
            out.append('')
            out.append(f'#include "templates/{d}/{intro}"')
        out.append('')
        if path.suffix == '.typ':
            out.append(f'#include "templates/{e["key"]}"')
        else:
            out.append(f'#zebraw(lang: false)[#raw(readcode("templates/{e["key"]}"), lang: "{lang_of(path)}", block: true)]')

    out.append(f'= {d}')
    done = set()
    for e in items:
        rel = e['key'].split('/', 1)[1]
        if '/' in rel and not group_body(rel):
            continue                       # 子条, 跟着小节正文一起出
        emit(e, 2)
        g = rel.split('/', 1)[0] if '/' in rel else rel.rsplit('.', 1)[0]
        if g in groups and g not in done:  # 小节正文后紧跟它的三级子条
            for ce in groups[g]:
                emit(ce, 3)
            done.add(g)
    for g, ces in groups.items():          # 没有正文的小节: 兜底只出标题
        if g not in done:
            out += ['', f'== {g}']
            for ce in ces:
                emit(ce, 3)
    out.append('')
    total += len(items)
SEC.write_text('\n'.join(out).rstrip('\n') + '\n')
print(f'\n已生成 sections.typ: {sum(1 for d in data["sections"] if any(not e["missing"] and e["key"].startswith(d + "/") for e in data["entries"]))} 章 / {total} 个代码块')

# ---------- 5) xcpc.typ 切换为 include 模式 ----------
main = MAIN.read_text()
if '#include "sections.typ"' not in main:
    lines = main.splitlines()
    idx = next(i for i, l in enumerate(lines) if l.startswith('= '))
    assert '#include' not in '\n'.join(lines[:idx])
    main = '\n'.join(lines[:idx]).rstrip('\n') + '\n\n#include "sections.typ"\n'
    MAIN.write_text(main)
    print('xcpc.typ: 正文已替换为 #include "sections.typ"')
else:
    print('xcpc.typ: 已是 include 模式, 未改动')
