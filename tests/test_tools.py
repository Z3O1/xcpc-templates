"""生成器、选测、构建与格式工具的回归;运行 python3 -m unittest discover -s tests -v。"""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import gen

ROOT = Path(__file__).resolve().parents[1]


class Workspace(unittest.TestCase):
    def setUp(self):
        (ROOT / 'tmp').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix='tools test ', dir=ROOT / 'tmp')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'templates').mkdir()
        self.put('xcpc.typ', '#include "sections.typ"\n')

    def put(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding='utf-8')
        return path

    def script(self, name):
        shutil.copy2(ROOT / name, self.root / name)

    def run_script(self, script, *args, env=None):
        return subprocess.run(['bash', script, *args], cwd=self.root,
                              env={**os.environ, 'XCPC_CHECK_TIMEOUT': '60', **(env or {})},
                              capture_output=True, text=True, timeout=15)


class GeneratorTests(Workspace):
    def generate(self):
        with contextlib.redirect_stdout(io.StringIO()):
            gen.main(self.root)
        return json.loads((self.root / '.manifest.json').read_text(encoding='utf-8'))

    def source(self, relative, title='示例', body='int x;'):
        return self.put('templates/' + relative, f'// f(): 章 · {title}\n\n{body}\n')

    def test_titles_nested_parens_and_unbalanced_suffix(self):
        self.assertEqual(gen.strip_trailing_parens('标题 (O(n log n))（附注）'), '标题')
        self.assertEqual(gen.strip_trailing_parens('标题 (未配平))'), '标题 (未配平))')
        path = self.source('章/a.cpp', '标题 (备注 (嵌套))')
        self.assertEqual(gen.title_from_comment(path), '标题')
        path.write_text('int x;\n', encoding='utf-8')
        self.assertEqual(gen.title_from_comment(path), 'a')

    def test_backup_matching_ignores_format_but_preserves_literals_and_operators(self):
        self.assertEqual(gen.norm('int x = f(1, 2); // note'), gen.norm('int x=f(1,2);'))
        self.assertNotEqual(gen.norm('"a b"'), gen.norm('"ab"'))
        self.assertNotEqual(gen.norm('a + +b'), gen.norm('a++b'))

    def test_import_has_no_generation_side_effects(self):
        self.script('gen.py')
        self.source('章/a.cpp')
        spec = importlib.util.spec_from_file_location('isolated_gen', self.root / 'gen.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertFalse((self.root / '.manifest.json').exists())
        self.assertFalse((self.root / 'sections.typ').exists())

    def test_scan_ignores_checks_hidden_files_and_hidden_directories(self):
        self.source('章/a.cpp')
        self.source('章/a.check.cpp')
        self.source('章/.private.cpp')
        self.source('章/.cache/private.cpp')
        self.source('.cache/private.cpp')
        self.put('templates/章/a.typ', '// 介绍\n正文\n')
        files, intros = gen.scan_templates(self.root / 'templates')
        self.assertEqual(list(files), ['章'])
        self.assertEqual(list(files['章']), ['a.cpp'])
        self.assertEqual(intros, {'章': {'a': 'a.typ'}})

    def test_missing_restore_hidden_and_manifest_titles(self):
        path = self.source('章/a.cpp', '原题')
        data = self.generate()
        data['entries'][0]['title'] = '自订标题'
        self.put('.manifest.json', json.dumps(data, ensure_ascii=False))
        path.unlink()
        self.assertTrue(self.generate()['entries'][0]['missing'])
        self.source('章/a.cpp', '新题', '// hide\n\nint x;')
        data = self.generate()
        self.assertEqual(data['entries'][0]['title'], '自订标题')
        self.assertFalse(data['entries'][0]['missing'])
        self.assertTrue(data['entries'][0]['hidden'])
        files, intros = gen.scan_templates(self.root / 'templates')
        _, chapters, total = gen.render_sections(data, files, intros)
        self.assertEqual((chapters, total), (0, 0))
        self.source('章/a.cpp', '新题')
        self.assertFalse(self.generate()['entries'][0]['hidden'])

    def test_hide_marker_in_body_does_not_hide_template(self):
        path = self.source('章/a.cpp', body='int x;\n// hide')
        self.assertFalse(gen.hidden_from_header(path))
        path = self.source('章/a.cpp', body='// 隐藏\nint x;')
        self.assertTrue(gen.hidden_from_header(path))

    def test_intro_transition_and_hidden_intro(self):
        self.put('templates/章/a.typ', '独立资料\n')
        self.generate()
        self.source('章/a.cpp')
        data = self.generate()
        self.assertEqual([e['key'] for e in data['entries']], ['章/a.cpp'])
        sections = (self.root / 'sections.typ').read_text(encoding='utf-8')
        self.assertLess(sections.index('#include "templates/章/a.typ"'), sections.index('readcode("templates/章/a.cpp")'))
        self.source('章/a.cpp', body='// hide\nint x;')
        self.generate()
        self.assertNotIn('templates/章/a.typ', (self.root / 'sections.typ').read_text(encoding='utf-8'))
        (self.root / 'templates/章/a.cpp').unlink()
        data = self.generate()
        self.assertEqual([e['key'] for e in data['entries'] if not e['missing']], ['章/a.typ'])

    def test_intro_can_be_hidden_without_hiding_code(self):
        for name in ('a', '小节/a'):
            with self.subTest(name=name):
                self.source(f'章/{name}.cpp')
                self.put(f'templates/章/{name}.typ', '// 介绍\n// hide\n原理正文\n')
                data = self.generate()
                entry = next(e for e in data['entries'] if e['key'] == f'章/{name}.cpp')
                self.assertFalse(entry['hidden'])
                sections = (self.root / 'sections.typ').read_text(encoding='utf-8')
                self.assertIn(f'readcode("templates/章/{name}.cpp")', sections)
                self.assertNotIn(f'#include "templates/章/{name}.typ"', sections)
                self.put(f'templates/章/{name}.typ', '// 介绍\n原理正文\n')
                self.generate()
                sections = (self.root / 'sections.typ').read_text(encoding='utf-8')
                self.assertIn(f'#include "templates/章/{name}.typ"', sections)

    def test_hide_marker_in_intro_body_does_not_hide_it(self):
        self.source('章/a.cpp')
        self.put('templates/章/a.typ', '// 介绍\n原理正文\n// hide\n')
        self.generate()
        sections = (self.root / 'sections.typ').read_text(encoding='utf-8')
        self.assertIn('#include "templates/章/a.typ"', sections)

    def test_groups_with_and_without_body(self):
        self.source('章/小节/a.cpp', '子条')
        self.put('templates/章/小节/小节.typ', '小节正文\n')
        self.source('章/无正文/b.cpp', '另一子条')
        self.generate()
        sections = (self.root / 'sections.typ').read_text(encoding='utf-8')
        self.assertIn('=== 子条', sections)
        self.assertIn('== 无正文', sections)
        self.assertLess(sections.index('#include "templates/章/小节/小节.typ"'), sections.index('=== 子条'))

    def test_recover_from_archived_backup(self):
        self.source('A/a.cpp', body='int a = 1; // formatted')
        self.source('Z/z.cpp', body='int z;')
        self.put('old_versions/xcpc.typ.bak', '= Z\n== 旧 Z\n```cpp\nint z;\n```\n= A\n== 旧 A\n```cpp\nint a=1;\n```\n')
        data = self.generate()
        self.assertEqual(data['sections'], ['Z', 'A'])
        self.assertEqual([e['title'] for e in data['entries']], ['旧 Z', '旧 A'])

    def test_recover_zebraw_and_third_level_titles(self):
        self.source('章/小节/a.cpp', '代码标题')
        self.put('templates/章/小节/a.typ', '介绍\n')
        self.put('sections.typ', '= 章\n== 小节\n=== 保存标题\n#include "templates/章/小节/a.typ"\n'
                 '#zebraw(lang: false)[#raw(readcode("templates/章/小节/a.cpp"), lang: "cpp", block: true)]\n')
        data = self.generate()
        self.assertEqual([e['title'] for e in data['entries']], ['保存标题'])
        self.assertEqual([e['key'] for e in data['entries']], ['章/小节/a.cpp'])

    def test_unchanged_generation_preserves_mtimes(self):
        self.source('章/a.cpp')
        self.generate()
        paths = [self.root / name for name in ('.manifest.json', 'sections.typ', 'xcpc.typ')]
        for path in paths:
            os.utime(path, ns=(1000000000, 1000000000))
        before = [path.stat().st_mtime_ns for path in paths]
        self.generate()
        self.assertEqual(before, [path.stat().st_mtime_ns for path in paths])

    def test_invalid_manifest_is_not_overwritten(self):
        for data in [[], {'sections': [], 'entries': [{}]},
                     {'sections': ['章', '章'], 'entries': []},
                     {'sections': ['章'], 'entries': [{'key': '章/../a.cpp', 'title': 'a'}]},
                     {'sections': ['章'], 'entries': [{'key': '章/a.cpp', 'title': 1}]}]:
            with self.subTest(data=data):
                text = json.dumps(data)
                self.put('.manifest.json', text)
                with self.assertRaises(ValueError):
                    self.generate()
                self.assertEqual((self.root / '.manifest.json').read_text(), text)
                self.assertFalse((self.root / 'sections.typ').exists())

    def test_main_migration_and_failure_before_writes(self):
        text = '#include "header.typ"\n#set text(size: 9pt)\n= 旧正文\n内容\n'
        self.assertEqual(gen.include_main(text), '#include "header.typ"\n#set text(size: 9pt)\n\n#include "sections.typ"\n')
        self.source('章/a.cpp')
        self.put('xcpc.typ', '#set text(size: 9pt)\n')
        with self.assertRaises(ValueError):
            self.generate()
        self.assertFalse((self.root / '.manifest.json').exists())
        self.assertFalse((self.root / 'sections.typ').exists())


class CheckRunnerTests(Workspace):
    def setUp(self):
        super().setUp()
        self.script('check.sh')
        self.put('templates/A/item.cpp', '// 模板,不应被当成 check\n')
        self.check('A/item', 'puts("A item"); return 0;')
        self.check('B/item', 'puts("B item"); return 0;')

    def check(self, name, body):
        return self.put(f'templates/{name}.check.cpp', '#include <cstdio>\nint main(){' + body + '}\n')

    def test_exact_paths_resolve_to_checks_not_sources(self):
        for name in ('A/item', 'templates/A/item.cpp', './templates/A/item.check.cpp', str(self.root / 'templates/A/item.cpp')):
            with self.subTest(name=name):
                result = self.run_script('check.sh', '-l', '-x', name)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('templates/A/item.check.cpp', result.stdout)
                self.assertIn('共 1 个', result.stdout)
                self.assertNotIn('templates/B/', result.stdout)

    def test_ambiguous_and_unknown_paths_fail(self):
        ambiguous = self.run_script('check.sh', '-l', '-x', 'item')
        self.assertNotEqual(ambiguous.returncode, 0)
        self.assertIn('歧义', ambiguous.stderr)
        self.assertNotEqual(self.run_script('check.sh', '-l', '-x', 'wrong/item').returncode, 0)

    def test_normalized_files_are_deduplicated(self):
        result = self.run_script('check.sh', '-l', '-f', './templates/A/item.check.cpp', '-f', 'templates/A/item.check.cpp')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('共 1 个', result.stdout)

    def test_invalid_options_and_timeouts(self):
        for args in (('-j', '0'), ('-j',), ('-f',), ('--bad',)):
            self.assertEqual(self.run_script('check.sh', *args).returncode, 2)
        for timeout in ('0', '-1', 'bad', '1s'):
            self.assertEqual(self.run_script('check.sh', '-l', env={'XCPC_CHECK_TIMEOUT': timeout}).returncode, 2)

    @unittest.skipUnless(shutil.which('g++'), '需要 g++')
    def test_parallel_order_verbose_and_chapter_working_directory(self):
        self.put('templates/A/relative.txt', 'input')
        self.check('A/cwd', 'FILE *f=fopen("relative.txt","r"); if(!f) return 1; fclose(f); puts("chapter cwd ok"); return 0;')
        result = self.run_script('check.sh', '-j', '2', '-v', '-x', 'A/cwd', 'B/item')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('2 passed, 0 failed', result.stdout)
        self.assertIn('chapter cwd ok', result.stdout)
        self.assertLess(result.stdout.index('ok    templates/A/cwd'), result.stdout.index('ok    templates/B/item'))
        self.assertFalse(list((self.root / 'tmp').glob('check.*')))

    @unittest.skipUnless(shutil.which('g++'), '需要 g++')
    def test_compile_runtime_and_sigterm_resistant_timeout_failures(self):
        self.put('templates/A/broken.check.cpp', 'not C++\n')
        self.check('A/failed', 'puts("runtime failure"); return 3;')
        self.put('templates/A/hung.check.cpp', '#include <csignal>\nint main(){signal(SIGTERM,SIG_IGN); for(;;){}}\n')
        result = self.run_script('check.sh', '-j', '3', '-x', 'broken', 'failed', 'hung', env={'XCPC_CHECK_TIMEOUT': '0.1'})
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn('COMPILE FAIL', result.stdout)
        self.assertIn('runtime failure', result.stdout)
        self.assertIn('超时', result.stdout)
        self.assertIn('0 passed, 3 failed', result.stdout)


class BuildTests(Workspace):
    def setUp(self):
        super().setUp()
        self.script('build.sh')
        self.script('gen.py')
        self.put('templates/章/a.cpp', '// f(): 章 · 示例\n\nint x;\n')
        self.compiler = self.put('fake typst', '#!/usr/bin/env bash\nprintf "%s\\n" "$@" >> "$CALL_LOG"\n'
                                 'if [[ "$1" == --version ]]; then echo "typst test"; fi\n')
        self.compiler.chmod(0o755)
        self.log = self.root / 'calls.log'
        self.env = {'TYPST': str(self.compiler), 'CALL_LOG': str(self.log)}

    def test_build_without_arguments(self):
        result = self.run_script('build.sh', env=self.env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.log.read_text().splitlines(), ['--version', 'compile', 'xcpc.typ', 'xcpc.pdf'])
        self.assertTrue((self.root / 'tmp/xdg').is_dir())

    def test_watch_consumes_option_and_forwards_remaining_arguments(self):
        result = self.run_script('build.sh', '--watch', '--font-path', 'font dir', env=self.env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.log.read_text().splitlines(), ['--version', 'watch', 'xcpc.typ', 'xcpc.pdf', '--font-path', 'font dir'])

    def test_missing_compiler_does_not_modify_generated_state(self):
        result = self.run_script('build.sh', env={'TYPST': str(self.root / 'missing')})
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / '.manifest.json').exists())


@unittest.skipUnless(shutil.which('typst') and shutil.which('pdffonts'), '需要 Typst 与 pdffonts')
class TypographyTests(Workspace):
    def test_document_keeps_original_body_and_code_fonts(self):
        self.put('xcpc.typ', (ROOT / 'xcpc.typ').read_text(encoding='utf-8'))
        self.put('sections.typ', '= 字符串\n中文正文\n```cpp\n// 中文注释\nint main() {}\n```\n')
        runtime = self.root / 'tmp/xdg'
        runtime.mkdir(parents=True)
        result = subprocess.run(['typst', 'compile', str(self.root / 'xcpc.typ'), str(self.root / 'fonts.pdf')],
                                env={**os.environ, 'XDG_RUNTIME_DIR': str(runtime)},
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        fonts = subprocess.check_output(['pdffonts', str(self.root / 'fonts.pdf')], text=True)
        self.assertIn('NotoSansCJKjp-Regular', fonts)
        self.assertIn('NotoSansCJKjp-Bold', fonts)
        self.assertIn('DejaVuSansMono', fonts)
        self.assertNotRegex(fonts, r'\+(KaiTi|NSimSun|SimSun)\b')


@unittest.skipUnless(shutil.which('clang-format'), '需要 clang-format 14+')
class FormatTests(Workspace):
    def test_default_excludes_checks_and_check_mode_is_read_only(self):
        self.script('format.sh')
        self.script('.clang-format')
        path = self.put('templates/章/a.cpp', '// 标题\n\nint f(){return 1;}\n')
        check = self.put('templates/章/a.check.cpp', 'int main(){return 0;}\n')
        before = path.read_text()
        self.assertNotEqual(self.run_script('format.sh', '--check').returncode, 0)
        self.assertEqual(path.read_text(), before)
        self.assertEqual(self.run_script('format.sh').returncode, 0)
        self.assertEqual(self.run_script('format.sh', '--check').returncode, 0)
        self.assertEqual(check.read_text(), 'int main(){return 0;}\n')


if __name__ == '__main__':
    unittest.main()
