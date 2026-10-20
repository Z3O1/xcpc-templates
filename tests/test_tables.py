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


# 从受控 Typst 文本提取 table/ref-table 参数块；识别双引号与括号，未闭合抛 ValueError。
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


# 按原顺序提取单行、不含嵌套方括号的单元格内容。
def cells(block):
    return re.findall(r'\[([^\[\]\n]*)\]', block)


# 解析四列/八列增长表并按上界排序；缺失质数计数记 None，重复上界抛 ValueError。
# 无法解析的行标签跳过；非缺失的数据列保留原文。
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


# 返回紧凑 UTF-8 JSON 编码的 SHA-256 指纹；value 须可 JSON 序列化。
def digest(value):
    encoded = json.dumps(value, ensure_ascii=False, separators=(',', ':')).encode()
    return hashlib.sha256(encoded).hexdigest()


class TableTests(unittest.TestCase):
    # 读取真实资料表并缓存表块，不修改源文件。
    def setUp(self):
        directory = ROOT / 'templates/通用'
        self.primes = table_blocks((directory / '大质数表.typ').read_text(encoding='utf-8'))
        self.constants = table_blocks((directory / '常数速查表.typ').read_text(encoding='utf-8'))

    # 锁定四张固定数据表的单元格数量与内容摘要。
    def test_prime_sequence_catalan_and_factorization_cells_unchanged(self):
        self.assertEqual(len(self.primes), 1)
        self.assertEqual(len(self.constants), 4)
        expected = [
            (self.primes[0], 198, 'b993fc5e68dd68e638bec426669de668740f06d30ae486db477a26ad1bb5e025'),
            (self.constants[0], 198, 'c600380438cf4c1e03d5fa002905a887b7dbb1a4cc7363081e7f440f67ceb1f9'),
            (self.constants[1], 88, '028af6e2af73b0fcbbe97654db40e121d41a3e2eadf1cb70aef7502694f8ea0c'),
            (self.constants[3], 216, 'd2b8894c63498debe6df6645d2b5d21adff41fd8c80d2ddc3285d0c531b248b8'),
        ]
        for block, count, checksum in expected:
            with self.subTest(checksum=checksum):
                values = cells(block)
                self.assertEqual(len(values), count)
                self.assertEqual(digest(values), checksum)
        self.assertNotIn('[], []', self.constants[3])

    # 锁定 10^1..10^18 的增长上界数据及六处质数计数缺项。
    def test_growth_bounds_preserve_values_and_missing_prime_counts(self):
        rows = growth_rows(self.constants[2])
        self.assertEqual([bound for bound, _ in rows], [10 ** k for k in range(1, 19)])
        self.assertEqual(digest(rows), 'abe19dff8562442f92aa97478a2bb4ccc272b9371f6acfc77dc9c022cfc2fcdc')
        self.assertEqual(sum(values[2] is None for _, values in rows), 6)
        self.assertNotIn('[], []', self.constants[2])

    # 独立编译资料页，核对正文一页及指定中文字体；工具缺失时跳过。
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
