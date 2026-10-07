#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""扫描 templates/ → .manifest.json → sections.typ → xcpc.typ。

章节顺序与标题以 manifest 为准;删除只标 missing,同名文件回来时原位恢复。
同名 .typ 是代码介绍,不占条目;子目录是小节,同名 .typ 是小节正文。
首行注释取新条目的标题,头部 // hide 或 // 隐藏使代码及介绍不进 PDF。
manifest 丢失时优先从 old_versions/xcpc.typ.bak,再从 sections.typ 恢复。

用法: python3 gen.py (或 ./build.sh)。导入本模块不扫描、不写文件。
"""
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent
LANG = {'cpp': 'cpp', 'sh': 'bash', 'py': 'python', 'rs': 'rust', 'typ': 'typst'}
FENCE = re.compile(r'^```(\w+)?\s*$')
CODE_TOKEN = re.compile(
    r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|'
    r'>>=|<<=|\.\.\.|->\*|::|\+\+|--|->|==|!=|<=|>=|&&|\|\||<<|>>|'
    r'\+=|-=|\*=|/=|%=|&=|\|=|\^=|##|[^\s]', re.S,
)


def lang_of(path: pathlib.Path) -> str:
    suffix = path.suffix.lstrip('.')
    return LANG.get(suffix, suffix)


def norm(text: str) -> str:
    """备份代码按 token 匹配:忽略格式/注释,但保留字符串与运算符的区别。"""
    return '\0'.join(token for token in CODE_TOKEN.findall(text)
                     if not token.startswith(('//', '/*')))


def codelines(path: pathlib.Path) -> list[str]:
    """迁移前的代码匹配规则:去掉第一个空行之前的头部。"""
    lines = path.read_text(encoding='utf-8').splitlines()
    index = next((i for i, line in enumerate(lines) if not line.strip()), len(lines))
    return lines[index + 1:]


def strip_trailing_parens(title: str) -> str:
    """配平剥离标题尾部的中英文括号备注,支持嵌套。"""
    title = title.rstrip()
    while title.endswith((')', '）')):
        close, open_ = (')', '(') if title.endswith(')') else ('）', '（')
        depth, index = 0, len(title) - 1
        while index >= 0:
            if title[index] == close:
                depth += 1
            elif title[index] == open_:
                depth -= 1
                if depth == 0:
                    break
            index -= 1
        if index <= 0 or depth != 0:
            break
        title = title[:index].rstrip()
    return title


def title_from_comment(path: pathlib.Path) -> str:
    """首行格式:// 接口(): 章节 · 标题 (备注);不匹配则用文件名。"""
    with path.open(encoding='utf-8') as source:
        line = source.readline()
    match = re.match(r'^[#/]+\s*[^:]*:\s*[^·]*·\s*(.+)$', line)
    if match:
        title = strip_trailing_parens(match.group(1)).strip()
        if title:
            return title
    return path.stem


def header_lines(path: pathlib.Path) -> list[str]:
    """读取开头连续的行注释,允许标题与说明之间留空行。"""
    lines = []
    with path.open(encoding='utf-8') as source:
        for line in source:
            line = line.strip()
            if not line:
                continue
            if not line.startswith(('//', '#')):
                break
            lines.append(line)
    return lines


def hidden_from_header(path: pathlib.Path) -> bool:
    return any(re.match(r'^[/#]+\s*(hide|隐藏)\b', line) for line in header_lines(path))


def scan_templates(directory: pathlib.Path) -> tuple[dict, dict]:
    """返回章→代码路径与章→介绍路径;忽略自测和隐藏目录。"""
    files, intros = {}, {}
    for chapter in sorted(path for path in directory.iterdir() if path.is_dir() and not path.name.startswith('.')):
        entries = {}
        for path in sorted(chapter.rglob('*')):
            relative = path.relative_to(chapter)
            if (not path.is_file() or any(part.startswith('.') for part in relative.parts)
                    or path.name.endswith('.check.cpp')):
                continue
            entries[relative.as_posix()] = path
        if not entries:
            continue
        code_stems = {name.rsplit('.', 1)[0] for name, path in entries.items() if path.suffix != '.typ'}
        files[chapter.name], intros[chapter.name] = {}, {}
        for name, path in entries.items():
            if path.suffix == '.typ' and name[:-4] in code_stems:
                intros[chapter.name][name[:-4]] = name
            else:
                files[chapter.name][name] = path
    return files, intros


def parse_backup(path: pathlib.Path, files: dict) -> dict:
    """按旧内联版代码内容恢复条目顺序与标题。"""
    lines = path.read_text(encoding='utf-8').splitlines()
    by_code = {}
    for chapter, entries in files.items():
        for name, source in entries.items():
            by_code.setdefault(norm('\n'.join(codelines(source))), []).append((chapter, name))
    state, chapter, title = {}, None, None
    index = 0
    while index < len(lines):
        match = re.match(r'^= (.+)$', lines[index])
        if match:
            chapter, title = match.group(1).strip(), None
        else:
            match = re.match(r'^== (.+)$', lines[index])
            if match:
                title = match.group(1).strip()
            elif FENCE.match(lines[index]) and chapter:
                end = index + 1
                while end < len(lines) and not FENCE.match(lines[end]):
                    end += 1
                hits = by_code.get(norm('\n'.join(lines[index + 1:end])), [])
                if len(hits) == 1:
                    actual_chapter, name = hits[0]
                    state[f'{actual_chapter}/{name}'] = (actual_chapter, title)
                else:
                    print(f'! 原文档代码块无法唯一匹配模板: {title} {hits}')
                index = end
        index += 1
    return state


def parse_state(path: pathlib.Path, intros: dict) -> dict:
    """兼容旧 raw 与现有 zebraw,恢复二级/三级标题,跳过代码介绍。"""
    state, title = {}, None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('= '):
            title = None
        match = re.match(r'^={2,3} (.+)$', line)
        if match:
            title = match.group(1).strip()
        match = re.search(r'readcode\("templates/([^/]+)/([^"\n]+)"\)', line)
        if not match:
            match = re.match(r'^#include "templates/([^/]+)/([^"\n]+)"$', line)
        if match:
            chapter, name = match.groups()
            if name.endswith('.typ') and name[:-4] in intros.get(chapter, {}):
                continue
            state[f'{chapter}/{name}'] = (chapter, title)
    return state


def validate_manifest(data: dict) -> None:
    """拒绝损坏的状态文件,避免对账时覆盖掉原有顺序。"""
    if not isinstance(data, dict) or not isinstance(data.get('sections'), list) or not isinstance(data.get('entries'), list):
        raise ValueError('manifest 需要 sections 与 entries 数组')
    chapters = data['sections']
    if any(not isinstance(chapter, str) or not chapter for chapter in chapters) or len(set(chapters)) != len(chapters):
        raise ValueError('manifest 章节必须是非空且不重复的字符串')
    keys = set()
    for entry in data['entries']:
        if not isinstance(entry, dict) or not isinstance(entry.get('key'), str):
            raise ValueError('manifest 条目需要字符串 key')
        key = entry['key']
        parts = key.split('/')
        if len(parts) < 2 or any(part in ('', '.', '..') for part in parts) or key in keys:
            raise ValueError(f'manifest 路径无效或重复: {key}')
        if 'title' not in entry or (entry['title'] is not None and not isinstance(entry['title'], str)):
            raise ValueError(f'manifest 标题必须是字符串或 null: {key}')
        keys.add(key)


def load_manifest(root: pathlib.Path, files: dict, intros: dict) -> dict:
    manifest = root / '.manifest.json'
    if manifest.exists():
        data = json.loads(manifest.read_text(encoding='utf-8'))
        validate_manifest(data)
        print('使用 .manifest.json 作为状态来源')
        return data
    backups = (root / 'old_versions/xcpc.typ.bak', root / 'xcpc.typ.bak')
    backup = next((path for path in backups if path.exists()), None)
    if backup:
        state = parse_backup(backup, files)
        print(f'首次运行: 从 {backup.relative_to(root)} 还原章节顺序与标题')
    elif (root / 'sections.typ').exists():
        state = parse_state(root / 'sections.typ', intros)
        print('无备份: 从现有 sections.typ 还原标题')
    else:
        state = {}
        print('无任何状态来源: 将按文件名排序')
    data = {'sections': [], 'entries': []}
    for key, (chapter, title) in state.items():
        if chapter not in data['sections']:
            data['sections'].append(chapter)
        data['entries'].append({'key': key, 'title': title})
    return data


def reconcile(data: dict, files: dict, intros: dict) -> dict:
    intro_keys = {f'{chapter}/{stem}.typ' for chapter, entries in intros.items() for stem in entries}
    for entry in data['entries']:
        if entry['key'] in intro_keys:
            print(f'- 已转为介绍: {entry["key"]} (不再独立显示)')
    data['entries'] = [entry for entry in data['entries'] if entry['key'] not in intro_keys]
    keys = {entry['key'] for entry in data['entries']}
    for chapter, entries in files.items():
        if chapter not in data['sections']:
            data['sections'].append(chapter)
        for name, path in entries.items():
            key = f'{chapter}/{name}'
            if key not in keys:
                title = title_from_comment(path)
                data['entries'].append({'key': key, 'title': title})
                keys.add(key)
                print(f'+ 新增: {key} (标题: {title})')
    for entry in data['entries']:
        chapter, name = entry['key'].split('/', 1)
        path = files.get(chapter, {}).get(name)
        missing = path is None
        hidden = not missing and hidden_from_header(path)
        if missing and not entry.get('missing'):
            print(f'- 已移除: {entry["key"]}')
        elif not missing and entry.get('missing'):
            print(f'^ 已恢复: {entry["key"]} -> 原位 {chapter} / {entry["title"] or name}')
        if hidden != entry.get('hidden', False):
            print(f'{"- 已隐藏" if hidden else "^ 已恢复显示"}: {entry["key"]}')
        entry['missing'], entry['hidden'] = missing, hidden
    return data


PREAMBLE = r'''// 由 gen.py 自动生成 —— 请勿手动修改; 改动 templates/ 后运行: ./build.sh
// 读取模板文件: 只丢掉第 1 个非空行(标题行), 其余内容含说明注释一律渲染
#import "@preview/zebraw:0.6.3": zebraw
#let readcode(path) = {
  let lines = read(path).split("\n")
  let i = 0
  // 第 1 个非空行是元数据行(标题), 不渲染; 其后原样输出
  while i < lines.len() and lines.at(i) == "" {
    i += 1
  }
  if i < lines.len() and (lines.at(i).starts-with("//") or lines.at(i).starts-with("#")) {
    i += 1
  }
  while i < lines.len() and lines.at(i) == "" {
    i += 1
  }
  let body0 = lines.slice(i)
  let body = if body0.len() > 0 and body0.last() == "" {
    body0.slice(0, body0.len() - 1)
  } else {
    body0
  }
  // 空数组 join() 返回 none(Typst 0.15), 这里显式给空串
  if body.len() > 0 { body.join("\n") } else { "" }
}
'''


def group_body(relative: str) -> bool:
    parts = relative.split('/')
    return len(parts) == 2 and parts[1] == parts[0] + '.typ'


def render_sections(data: dict, files: dict, intros: dict) -> tuple[str, int, int]:
    out = PREAMBLE.splitlines() + ['']
    total, chapters = 0, 0
    for chapter in data['sections']:
        items = [entry for entry in data['entries'] if not entry['missing'] and not entry['hidden']
                 and entry['key'].startswith(chapter + '/')]
        if not items:
            continue
        groups = {}
        for entry in items:
            relative = entry['key'].split('/', 1)[1]
            if '/' in relative and not group_body(relative):
                groups.setdefault(relative.split('/', 1)[0], []).append(entry)

        def emit(entry: dict, level: int) -> None:
            if entry['title'] is not None:
                out.extend(['', f'{"=" * level} {entry["title"]}'])
            name = entry['key'].split('/', 1)[1]
            path = files[chapter][name]
            intro = intros.get(chapter, {}).get(name.rsplit('.', 1)[0])
            if intro:
                out.extend(['', f'#include "templates/{chapter}/{intro}"'])
            out.append('')
            if path.suffix == '.typ':
                out.append(f'#include "templates/{entry["key"]}"')
            else:
                out.append(f'#zebraw(lang: false)[#raw(readcode("templates/{entry["key"]}"), '
                           f'lang: "{lang_of(path)}", block: true)]')

        out.append(f'= {chapter}')
        done = set()
        for entry in items:
            relative = entry['key'].split('/', 1)[1]
            if '/' in relative and not group_body(relative):
                continue
            emit(entry, 2)
            group = relative.split('/', 1)[0] if '/' in relative else relative.rsplit('.', 1)[0]
            if group in groups and group not in done:
                for child in groups[group]:
                    emit(child, 3)
                done.add(group)
        for group, children in groups.items():
            if group not in done:
                out.extend(['', f'== {group}'])
                for child in children:
                    emit(child, 3)
        out.append('')
        total += len(items)
        chapters += 1
    return '\n'.join(out).rstrip('\n') + '\n', chapters, total


def include_main(text: str) -> str:
    if '#include "sections.typ"' in text:
        return text
    lines = text.splitlines()
    index = next((i for i, line in enumerate(lines) if line.startswith('= ')), None)
    if index is None:
        raise ValueError('xcpc.typ 没有正文一级标题,无法安全切换 include 模式')
    return '\n'.join(lines[:index]).rstrip('\n') + '\n\n#include "sections.typ"\n'


def write_if_changed(path: pathlib.Path, text: str) -> None:
    """内容不变时不触碰时间戳,避免 watch 无谓重编译。"""
    if not path.exists() or path.read_text(encoding='utf-8') != text:
        path.write_text(text, encoding='utf-8')


def main(root: pathlib.Path = ROOT) -> None:
    files, intros = scan_templates(root / 'templates')
    data = reconcile(load_manifest(root, files, intros), files, intros)
    sections, chapters, total = render_sections(data, files, intros)
    main_path = root / 'xcpc.typ'
    old_main = main_path.read_text(encoding='utf-8')
    new_main = include_main(old_main)  # 所有内容先算好再落盘,主文件出错时不改状态
    write_if_changed(root / '.manifest.json', json.dumps(data, ensure_ascii=False, indent=1) + '\n')
    write_if_changed(root / 'sections.typ', sections)
    write_if_changed(main_path, new_main)
    print(f'\n已生成 sections.typ: {chapters} 章 / {total} 个代码块')
    print('xcpc.typ: 已是 include 模式, 未改动' if new_main == old_main else 'xcpc.typ: 正文已替换为 include 模式')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as error:
        print(f'gen.py: {error}', file=sys.stderr)
        sys.exit(1)
