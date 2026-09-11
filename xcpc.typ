#set page(
  paper: "a4",
  flipped: true, // 横向
  margin: (top: 0.8cm, bottom: 0.8cm, left: 0.6cm, right: 0.6cm),
  header: [
    #place(
      right + top,
      dx: 0.5cm, // 距右边界的偏移
      dy: 0.75cm, // 距顶边的偏移
      [#text(size: 12pt)[#context here().page()]],
    )
  ],
  footer: [
    #place(
      right + bottom,
      dx: 0.5cm, // 距离右边界的偏移
      dy: -0.75cm, // 距离底边的偏移
      [#text(size: 12pt)[#context here().page()]],
    )
  ],
)

#set text(size: 9pt)

// 标题编号: 大点 (Level 1) = 1, 小点 (Level 2) = 1.1
#set heading(numbering: "1.1")

// 目录: 每个一级标题(大块)后面加一个空行
#show outline.entry.where(level: 1): it => v(1.4em, weak: false) + it

// ===== 目录: Typst 原生 outline() =====
#set page(columns: 3)
#outline()
#pagebreak()

// ===== 正文: 两列 =====
#set page(columns: 2)

#include "sections.typ"
