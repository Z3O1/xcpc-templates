// 通用资料表:满栏、连续细网格与灰色表头;不改变文档字体。
#let ref-table(columns: (), header-rows: (0,), align: right, ..cells) = {
  block(width: 100%, breakable: false)[
    #set text(size: 6.5pt, number-type: "lining", number-width: "tabular")
    #set par(leading: 0.35em)
    #table(
      columns: columns,
      inset: (x: 2pt, y: 1pt),
      column-gutter: 0pt,
      row-gutter: 0pt,
      stroke: 0.3pt + luma(40%),
      fill: (x, y) => if header-rows.contains(y) { luma(92%) } else { none },
      align: (x, y) => {
        let side = if type(align) == array { align.at(x) } else { align }
        if header-rows.contains(y) { center + horizon } else { side + horizon }
      },
      ..cells.pos(),
    )
  ]
}
