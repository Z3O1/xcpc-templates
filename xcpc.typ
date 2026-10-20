#set page(
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

// 固定原版字体,避免系统新增字体改变中文 fallback。
#set text(font: ("Libertinus Serif", "Noto Sans CJK JP"), size: 9pt)
#show raw: set text(font: ("DejaVu Sans Mono", "Noto Sans CJK JP"))

// 标题编号: 大点 (Level 1) = 1, 小点 (Level 2) = 1.1; sticky 防止小节标题孤立在栏尾
#set heading(numbering: "1.1")
#show heading: set block(sticky: true)

// 目录: 每个一级标题(大块)后面加一个空行
#show outline.entry.where(level: 1): it => v(1.4em, weak: false) + it

// ===== 目录: Typst 原生 outline() =====
#set page(columns: 3)
#outline()
#pagebreak()

// ===== 正文: 两列 =====
#set page(columns: 2)

#include "sections.typ"
