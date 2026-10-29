

#let conf(doc) = {
	import "@preview/zebraw:0.6.1": *
	show: zebraw
  ////////////////////
  // Document Settings
  ////////////////////
  let fntsz = 11pt
  let after = fntsz * 1
  set page(margin: .5in, width: 8.5in, height: 11in, numbering: "1")
  set columns(gutter: fntsz)
  set text(font: "Georgia", size: fntsz)
  set par(leading:  0.45em, spacing: 1em, justify: true)

  ////////////////
  // Set Defaults
  ////////////////

  show math.equation: set text(size: 12pt)

  show title: it => {
    set align(center + horizon)
    set text(size: 20pt, weight: "bold")
    v(-5em)
    it
  }

  show outline: it => {
    show heading: set text(size: 15pt)
    show outline.entry: it => pad(left: 1em, right: 0em, it)
    it
  }

  show link: it => {
    set text(blue)
    it
  }

  show ref: it => {
    if it.element != none {
      set text(baseline: 0em, size: 1em, weight: "bold")
      it
    } else {
      set text(baseline: -.11em, size: .8em)
      it
    }
  }

  show heading: it => {
    block(below: after)[#it]
  }
  show figure.where(kind: "code"): it => {
    // let col = rgb("#2746e397")
    // block(fill: none, stroke: col, inset: 8pt,below:5pt, radius: 20pt, width: 100%)[
    //     #it.body
    // ]
    block[
      #it.body
    ]
    v(-1em)
    it.caption
  }

  let fig_width = 200pt
  show figure.caption: it => {
    let sz = fntsz - 3pt
    set text(sz)
    set par(leading: sz / 2)
    set align(left)
    h(10pt)
    box(width: fig_width)[
      *#it.supplement #it.counter.display(it.numbering)#it.separator* #it.body
    ]
  }

  set heading(numbering: none)

  let theme_col = rgb("#8c2fc185")

  show "XXX": text.with(fill: red)

  show "GH": box[
    #box(
      link("https://github.com/canderson318/BIOS6011_midterm", image("github.svg", height: .8em)),
    )
  ]

  show heading.where(level: 1): it => {
    set text(fntsz + 1pt, weight: "bold")
    it
  }

  show heading.where(level: 2): it => {
    set text(fntsz, weight: "bold")
    it
  }

  doc
}

////////////////////////
// Function Definitions
////////////////////////

//************************************************************************\\
// Creates a filled rect with text
// Make a Todo note
//
// Arguments:
//   body: Text to show
//   color: bg fill, (default: yellow)
//   title: rect label (default: TODO)
//
// Returns: a colored box with body text labeled TODO
#let TODO(body, color: yellow, title: "TODO") = {
  rect(
    width: 100%,
    radius: 3pt,
    stroke: 0.5pt,
    fill: color,
  )[
    #text(weight: 700)[#title]: #body
  ]
}

//************************************************************************\\

#let today = datetime.today()

#let fmt-num(n, decimals: 2, format: "f") = {
  if format == "e" {
    let exp = if n == 0 { 0 } else { calc.floor(calc.log(calc.abs(n), base: 10)) }
    let mantissa = n / calc.pow(10, exp)
    let m = str(calc.round(mantissa, digits: decimals))
    let parts = m.split(".")
    let frac = if parts.len() > 1 { parts.at(1) } else { "" }
    while frac.len() < decimals { frac = frac + "0" }
    parts.at(0) + "." + frac + "e+" + str(exp)
  } else {
    let s = str(calc.round(n, digits: decimals))
    let parts = s.split(".")
    let whole = parts.at(0)
    let frac = if parts.len() > 1 { parts.at(1) } else { "" }
    while frac.len() < decimals { frac = frac + "0" }
    let chars = whole.clusters().rev()
    let grouped = ()
    for (i, c) in chars.enumerate() {
      if i > 0 and calc.rem(i, 3) == 0 { grouped.push(",") }
      grouped.push(c)
    }
    grouped.rev().join("") + if decimals > 0 { "." + frac } else { "" }
  }
}
