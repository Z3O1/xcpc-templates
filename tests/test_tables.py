"""资料表回归:锁住整理前的数据,并检查双栏排版与字体。"""
from decimal import Decimal, InvalidOperation
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def table_blocks(text):
    blocks = []
    for match in re.finditer(r'#(?:ref-)?table\(', text):
        start = match.end()
        depth, quoted, escape = 1, False, False
        for index in range(start, len(text)):
            char = text[index]
            if quoted:
                if escape:
                    escape = False
                elif char == '\\':
                    escape = True
                elif char == '"':
                    quoted = False
            elif char == '"':
                quoted = True
            elif char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
                if not depth:
                    blocks.append(text[start:index])
                    break
        else:
            raise ValueError('表格括号未闭合')
    return blocks


def cells(block):
    return re.findall(r'\[([^\[\]\n]*)\]', block)


def growth_rows(block):
    rows = {}
    for line in block.splitlines():
        if not line.lstrip().startswith('['):
            continue
        values = cells(line)
        for chunk in (values[:4], values[4:8]):
            if len(chunk) != 4:
                continue
            label = chunk[0].strip('$ ')
            try:
                bound = 10 ** int(label[3:]) if label.startswith('10^') else int(Decimal(label))
            except (ValueError, InvalidOperation):
                continue
            if bound in rows:
                raise ValueError(f'增长上界重复行: {bound}')
            rows[bound] = [chunk[1], chunk[2], None if chunk[3] in ('', '—') else chunk[3]]
    return sorted(rows.items())


def digest(value):
    encoded = json.dumps(value, ensure_ascii=False, separators=(',', ':')).encode()
    return hashlib.sha256(encoded).hexdigest()


class TableTests(unittest.TestCase):
    def setUp(self):
        directory = ROOT / 'templates/通用'
        self.primes = table_blocks((directory / '大质数表.typ').read_text(encoding='utf-8'))
        self.constants = table_blocks((directory / '常数速查表.typ').read_text(encoding='utf-8'))

    def test_prime_sequence_catalan_and_factorization_cells_unchanged(self):
        self.assertEqual(len(self.primes), 1)
        self.assertEqual(len(self.constants), 4)
        expected = [
            (self.primes[0], 198, 'b993fc5e68dd68e638bec426669de668740f06d30ae486db477a26ad1bb5e025'),
            (self.constants[0], 198, 'c600380438cf4c1e03d5fa002905a887b7dbb1a4cc7363081e7f440f67ceb1f9'),
            (self.constants[1], 88, '028af6e2af73b0fcbbe97654db40e121d41a3e2eadf1cb70aef7502694f8ea0c'),
            (self.constants[3], 216, '5a3a09f4ae216b1ca482071e0934b172419cd14a3e2a16d87de0fdcaf5e4bad8'),
        ]
        for block, count, checksum in expected:
            with self.subTest(checksum=checksum):
                values = cells(block)
                self.assertEqual(len(values), count)
                self.assertEqual(digest(values), checksum)

    def test_growth_bounds_preserve_values_and_missing_prime_counts(self):
        rows = growth_rows(self.constants[2])
        self.assertEqual([bound for bound, _ in rows], [10 ** k for k in range(1, 19)])
        self.assertEqual(digest(rows), 'abe19dff8562442f92aa97478a2bb4ccc272b9371f6acfc77dc9c022cfc2fcdc')
        self.assertEqual(sum(values[2] is None for _, values in rows), 6)
        self.assertNotIn('[], []', self.constants[2])

    @unittest.skipUnless(all(shutil.which(tool) for tool in ('typst', 'pdfinfo', 'pdffonts')), '需要 Typst 与 Poppler')
    def test_reference_tables_fit_one_page_with_original_fonts(self):
        (ROOT / 'tmp').mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as directory:
            directory = Path(directory)
            (directory / 'xcpc.typ').write_text((ROOT / 'xcpc.typ').read_text(encoding='utf-8'), encoding='utf-8')
            (directory / 'sections.typ').write_text(
                '= 通用\n== 大质数表\n#include "../../templates/通用/大质数表.typ"\n'
                '== 常数速查表\n#include "../../templates/通用/常数速查表.typ"\n', encoding='utf-8')
            pdf = directory / 'tables.pdf'
            runtime = directory / 'xdg'
            runtime.mkdir()
            result = subprocess.run(['typst', 'compile', '--root', str(ROOT), str(directory / 'xcpc.typ'), str(pdf)],
                                    env={**os.environ, 'XDG_RUNTIME_DIR': str(runtime)},
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            info = subprocess.check_output(['pdfinfo', str(pdf)], text=True, stderr=subprocess.DEVNULL)
            # 主文件还会单独生成目录页,资料表正文应仅占一页。
            self.assertRegex(info, r'Pages:\s+2\b')
            fonts = subprocess.check_output(['pdffonts', str(pdf)], text=True)
            self.assertIn('NotoSansCJKjp-Regular', fonts)
            self.assertIn('NotoSansCJKjp-Bold', fonts)
            self.assertNotRegex(fonts, r'\+(KaiTi|NSimSun|SimSun)\b')


if __name__ == '__main__':
    unittest.main()
