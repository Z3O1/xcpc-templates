// 板子书与 tests/style-preview.typ 共用; 修改后用小样预览。
#import "@preview/zebraw:0.6.3": zebraw-init

#let book(body) = {
  set page(
    paper: "a4",
    flipped: true, // 横向
    margin: (top: 0.8cm, bottom: 0.8cm, left: 0.6cm, right: 0.6cm),
    header: context {
      let p = here().page()
      let chs = query(heading.where(level: 1)).filter(it => it.numbering != none and it.location().page() <= p)
      let secs = query(heading.where(level: 2)).filter(it => it.numbering != none and it.location().page() <= p)
      let label = if chs.len() > 0 {
        let ch = chs.last()
        let cnum = counter(heading).at(ch.location()).first()
        let sec = if secs.len() > 0 { secs.last() } else { none }
        if sec != none and counter(heading).at(sec.location()).first() == cnum {
          let snum = counter(heading).at(sec.location()).slice(0, 2).map(str).join(".")
          [#cnum #ch.body · #snum #sec.body #h(0.6em) *#p*]
        } else {
          [#cnum #ch.body #h(0.6em) *#p*]
        }
      } else {
        [*#p*]
      }
      place(right + top, dx: 0cm, dy: 0.25cm, text(size: 8.5pt, fill: luma(70), label))
    },
  )

  // 正文字号; 固定字体,避免系统新增字体改变中文 fallback。
  set text(font: ("Libertinus Serif", "Noto Sans CJK JP"), size: 8pt)
  show raw: set text(font: ("DejaVu Sans Mono", "Noto Sans CJK JP"), size: 6pt)

  set heading(numbering: "1.1")
  show heading: set block(sticky: true)
  show outline.entry.where(level: 1): it => v(1.4em, weak: false) + it

  // 目录三栏,独占页面。
  set page(columns: 3)
  outline()
  pagebreak()

  // 正文双栏; leading 是文字行间空隙,spacing 是段落间距。
  set page(columns: 2)
  set par(leading: 0.5em, spacing: 0.6em)
  show heading: set block(above: 0.7em, below: 0.35em)

  // 每行代码上、下的内边距; 数值越小越紧凑,不是纸张页边距。
  show: zebraw-init.with(inset: (top: 0.25em, bottom: 0.25em))

  body
}
