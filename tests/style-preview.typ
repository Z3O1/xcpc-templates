// 在仓库根目录运行: mkdir -p tmp && typst watch --root . tests/style-preview.typ tmp/style-preview.pdf
// 调整 ../style.typ,保存后此小样和正式板子书使用同一套样式。
#import "../style.typ": book
#import "@preview/zebraw:0.6.3": zebraw
#show: book

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

= 样式预览
== 正文与段落
这是一段中文说明,用于观察字号、文字行距与自动折行。English text、数字 0123456789 和行内代码 `vector<int>` 用来检查中英文混排。修改 `style.typ` 中的 `leading` 调整同一段内相邻文字行之间的空隙。

这是第二段。修改 `spacing` 可调整它与上一段的距离;修改标题的 `above`、`below` 可观察标题前后的留白。代码每行的上下留白则由 `inset` 控制,不是这里的段落间距。

=== 公式与列表
前缀和满足 $s_i = sum_(j=1)^i a_j$,区间和为 $s_r - s_(l-1)$。
$ sum_(i=1)^n i = frac(n(n+1), 2) $

- 第一项:检查列表的缩进与行间距。
- 第二项:包含较长的中文说明,用于观察列表项自动折行后的行距以及与下一段的距离。

== 代码:短行与空行
#zebraw(lang: false)[```cpp
// 这是中文注释,用于观察代码字体与每行上下的留白。
int sum(int n) {
    int ans = 0;
    for(int i = 1; i <= n; ++i) {
        ans += i;
    }
    return ans;
}

vector<int> a = {1, 2, 3, 4};
int result = sum((int)a.size());
```]

// 固定换栏,方便并排比较正文、短行和长行代码。
#colbreak()
== 代码:长行自动折行
#zebraw(lang: false)[```cpp
// 这是一条故意写长的中文注释,用于检查代码自动折行后是否清楚,并比较 top、bottom 和 leading 对代码原始行与折行片段之间留白的影响。
vector<int> collect(const vector<int> &values, int lower_bound, int upper_bound) {
    vector<int> result;
    for(int value : values) {
        if(lower_bound <= value && value <= upper_bound) result.push_back(value);
    }
    sort(result.begin(), result.end(), [](int left, int right) { return left < right; });
    return result;
}

vector<pair<int, int>> intervals = {{1, 3}, {5, 8}, {13, 21}, {34, 55}, {89, 144}, {233, 377}, {610, 987}, {1597, 2584}};
```]

== 标题前后的间距
这是二级标题后的第一段,用来检查标题与正文的距离。

=== 三级标题
这是三级标题后的正文。标题应与后续内容保持在同一栏,避免单独留在栏尾。

=== 相邻代码块
#zebraw(lang: false)[```cpp
int minimum(int a, int b) {
    return a < b ? a : b;
}
```]

#zebraw(lang: false)[```cpp
int maximum(int a, int b) {
    return a > b ? a : b;
}
```]

=== 长代码

#zebraw(lang: false)[#raw(readcode("../templates/计算几何/平面剖分/Delaunay.cpp"), lang: "cpp", block: true)]