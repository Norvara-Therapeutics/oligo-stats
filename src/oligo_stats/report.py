"""Self-contained HTML reporting.

The document embeds its own style rules, uses no JavaScript and loads nothing
from the network, so a report stays readable from a local file or an email
attachment years later. Output is byte-deterministic for a given input and set
of parameters: nothing in it depends on the clock or the filesystem beyond the
source name that is passed in.
"""

from __future__ import annotations

from dataclasses import fields
from html import escape
from typing import Sequence

from oligo_stats.oligo import CLAMP_WINDOW, Analysis, Criteria, Flag, OligoStats
from oligo_stats.structure import Dimer, Hairpin

_STYLE = """
html { -webkit-text-size-adjust: 100%; }
body { font: 16px/1.55 Georgia, 'Times New Roman', serif; color: #1a1a1a;
       background: #fff; margin: 0 auto; padding: 3rem 1.25rem 4rem;
       max-width: 46rem; }
h1 { font-size: 1.5rem; font-weight: normal; margin: 0 0 .25rem; }
h2, h3 { font-size: .8rem; font-variant: small-caps; letter-spacing: .08em;
     font-weight: normal; margin: 2.5rem 0 .6rem; border-bottom: 1px solid #1a1a1a;
     padding-bottom: .2rem; }
h3 { margin-top: 1.8rem; }
h2.oligo { font-size: 1.1rem; font-variant: normal; letter-spacing: normal;
     font-style: italic; border-bottom: none; padding-bottom: 0; margin: 2.8rem 0 .2rem; }
p { margin: .5rem 0; }
.lede { color: #555; font-size: .9rem; margin-bottom: .5rem; }
.note { color: #555; font-size: .85rem; }
code, .seq, .struct, td.num, th.num { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
table { border-collapse: collapse; width: 100%; margin: .6rem 0 .2rem;
        font-size: .9rem; }
th { text-align: left; font-weight: normal; font-variant: small-caps;
     letter-spacing: .04em; border-bottom: 1px solid #1a1a1a; padding: .25rem .5rem .25rem 0; }
td { padding: .22rem .5rem .22rem 0; border-bottom: 1px solid #e6e6e6;
     vertical-align: baseline; }
td.num, th.num { text-align: right; padding-right: 1rem; font-size: .85rem; }
tr:last-child td { border-bottom: none; }
dl.metrics { display: grid; grid-template-columns: max-content 1fr; gap: .15rem 1.2rem;
             margin: .6rem 0; font-size: .9rem; }
dl.metrics dt { color: #555; }
dl.metrics dd { margin: 0; font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
.seq { font-size: .8rem; line-height: 1.45; white-space: pre-wrap; word-break: break-all;
       margin: .2rem 0 1rem; }
.struct { font-size: .8rem; line-height: 1.35; white-space: pre; overflow-x: auto;
          margin: .2rem 0 .9rem; }
.bar { display: block; height: 2px; background: #1a1a1a; margin: .3rem 0 .1rem; }
.tag { font-variant: small-caps; letter-spacing: .04em; font-size: .8rem; color: #0f6f6a; }
.tag.warn { color: #8a5a00; }
ul.flags { margin: .4rem 0 .2rem; padding-left: 1.1rem; font-size: .9rem; }
ul.flags li { margin: .15rem 0; }
footer { margin-top: 3rem; padding-top: .5rem; border-top: 1px solid #e6e6e6;
         color: #555; font-size: .8rem; }
footer p { margin: .35rem 0; }
@media print { body { padding: 0; max-width: none; } h2, h3 { page-break-after: avoid; } }
"""

_SEVERITY_LABEL = {"ok": "pass", "note": "note", "warn": "review"}


def _tag(severity: str) -> str:
    label = _SEVERITY_LABEL[severity]
    css = "tag warn" if severity == "warn" else "tag"
    return f'<span class="{css}">{label}</span>'


def _number(value: float | None, digits: int = 2, suffix: str = "") -> str:
    if value is None:
        return "n/a"
    return f"{value:,.{digits}f}{suffix}"


def _wrap(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def _conditions(analysis: Analysis) -> list[str]:
    ions = analysis.ions
    rows = [
        ("Oligonucleotide", f"{ions.strand * 1e6:,.3f} uM total strand"),
        ("Sodium", f"{ions.sodium * 1e3:,.1f} mM"),
    ]
    if ions.potassium:
        rows.append(("Potassium", f"{ions.potassium * 1e3:,.1f} mM"))
    if ions.tris:
        rows.append(("Tris", f"{ions.tris * 1e3:,.1f} mM, counted as half monovalent"))
    rows.append(("Monovalent total", f"{ions.monovalent * 1e3:,.1f} mM"))
    if ions.magnesium or ions.dntp:
        rows.append(("Magnesium", f"{ions.magnesium * 1e3:,.2f} mM"))
        rows.append(("dNTP, total", f"{ions.dntp * 1e3:,.2f} mM"))
        rows.append(("Free magnesium", f"{ions.free_magnesium * 1e3:,.2f} mM"))
    rows.append(("Ion correction", f"{ions.regime} regime"))
    rows.append(
        ("5' terminus", "hydroxyl" if analysis.terminus == "hydroxyl" else "phosphate")
    )
    parts = ['<dl class="metrics">']
    for name, value in rows:
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    return parts


def _exact(value: float) -> str:
    """Render a threshold as configured, without rounding it away.

    The report is the record of which rules produced its flags, so a limit of
    42.5% must not appear as 42%.
    """
    if isinstance(value, int) or value == int(value):
        return str(int(value))
    return f"{value:g}"


def _criteria_rows(criteria: Criteria) -> list[tuple[str, str]]:
    """Every design threshold, paired with the field it came from.

    Built field by field from :class:`Criteria` and checked against
    ``dataclasses.fields`` so that a threshold added to the rules cannot be
    quietly left out of the document that reports them.
    """
    rows = [
        (
            "Length",
            f"{_exact(criteria.length_min)}-{_exact(criteria.length_max)} nt",
            ("length_min", "length_max"),
        ),
        (
            "GC content",
            f"{_exact(criteria.gc_min)}-{_exact(criteria.gc_max)} %",
            ("gc_min", "gc_max"),
        ),
        (
            "Melting temperature",
            f"{_exact(criteria.melting_min)}-{_exact(criteria.melting_max)} C",
            ("melting_min", "melting_max"),
        ),
        (
            f"3' GC clamp, last {CLAMP_WINDOW} bases",
            f"{_exact(criteria.clamp_min)}-{_exact(criteria.clamp_max)} of G or C",
            ("clamp_min", "clamp_max"),
        ),
        (
            "3' end stacking",
            f"at or below {_exact(criteria.end_stability)} kcal/mol",
            ("end_stability",),
        ),
        (
            "Single-base run",
            f"{_exact(criteria.max_run)} bases or more",
            ("max_run",),
        ),
        (
            "Hairpin",
            f"at or below {_exact(criteria.hairpin_free_energy)} kcal/mol",
            ("hairpin_free_energy",),
        ),
        (
            "Dimer",
            f"at or below {_exact(criteria.dimer_free_energy)} kcal/mol",
            ("dimer_free_energy",),
        ),
        (
            "Dimer reaching a 3' end",
            f"at or below {_exact(criteria.dimer_three_prime_free_energy)} kcal/mol",
            ("dimer_three_prime_free_energy",),
        ),
    ]
    covered = {name for _, _, names in rows for name in names}
    missing = {field.name for field in fields(criteria)} - covered
    if missing:  # pragma: no cover - guarded by a test over Criteria
        raise AssertionError(f"criteria absent from the report: {sorted(missing)}")
    return [(label, value) for label, value, _ in rows]


def _criteria_block(analysis: Analysis) -> list[str]:
    parts = ['<dl class="metrics">']
    for label, value in _criteria_rows(analysis.criteria):
        parts.append(f"<dt>{escape(label)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    parts.append(
        '<p class="note">These are design conventions for a standard '
        "amplification primer, not physical constants. An oligonucleotide "
        "designed for another purpose, such as a probe or a cloning adapter, "
        "will fail some of them by intent.</p>"
    )
    return parts


def _summary_table(oligos: Sequence[OligoStats]) -> list[str]:
    parts = [
        "<table>",
        "<tr><th>Oligo</th>"
        '<th class="num">nt</th><th class="num">GC %</th>'
        '<th class="num">Tm &deg;C</th><th class="num">dG37</th>'
        '<th class="num">MW</th><th class="num">e260</th>'
        '<th class="num">Hairpin</th><th class="num">Self-dimer</th>'
        "<th></th></tr>",
    ]
    for item in oligos:
        hairpin = item.best_hairpin
        dimer = item.best_self_dimer
        if not item.structures_screened:
            # "none" here would read as "nothing found", which is the opposite
            # of what an unscreened oligonucleotide is entitled to.
            hairpin_cell = dimer_cell = "not screened"
        else:
            hairpin_cell = _number(hairpin.free_energy, 1) if hairpin else "none"
            dimer_cell = _number(dimer.free_energy, 1) if dimer else "none"
        status = (
            _tag(item.worst_severity)
            if item.structures_screened or item.worst_severity == "warn"
            else '<span class="tag">part screened</span>'
        )
        parts.append(
            f"<tr><td>{escape(item.name)}</td>"
            f'<td class="num">{item.length}</td>'
            f'<td class="num">{_number(item.composition.gc_percent, 1)}</td>'
            f'<td class="num">{_number(item.melting, 1)}</td>'
            f'<td class="num">{_number(item.free_energy, 1)}</td>'
            f'<td class="num">{item.optics.molecular_weight:,.1f}</td>'
            f'<td class="num">{item.optics.extinction:,.0f}</td>'
            f'<td class="num">{hairpin_cell}</td>'
            f'<td class="num">{dimer_cell}</td>'
            f"<td>{status}</td></tr>"
        )
    parts.append("</table>")
    parts.append(
        '<p class="note">Free energies are kcal/mol at 37 &deg;C, molecular '
        "weight g/mol for the single strand, and e260 the molar extinction "
        "coefficient at 260 nm in 1/(M cm).</p>"
    )
    return parts


def _composition_table(item: OligoStats) -> list[str]:
    counts = item.composition.counts
    widest = max(counts.values()) or 1
    parts = ["<table>", '<tr><th>Base</th><th class="num">Count</th>'
             '<th class="num">%</th><th></th></tr>']
    for base in ("A", "C", "G", "T"):
        count = counts[base]
        share = 100.0 * count / item.length if item.length else 0.0
        width = 0.0 if not count else max(1.0, 100.0 * count / widest)
        bar = f'<span class="bar" style="width: {width:.1f}%"></span>' if count else ""
        parts.append(
            f"<tr><td>{base}</td>"
            f'<td class="num">{count}</td>'
            f'<td class="num">{share:.1f}</td><td>{bar}</td></tr>'
        )
    parts.append("</table>")
    return parts


def _metrics_block(item: OligoStats) -> list[str]:
    rows: list[tuple[str, str]] = [
        ("Length", f"{item.length} nt"),
        ("GC content", _number(item.composition.gc_percent, 2, " %")),
        ("Tm, intended duplex", _number(item.melting, 2, " C")),
        ("dH", _number(item.thermo.enthalpy, 1, " kcal/mol")),
        ("dS", _number(item.thermo.entropy, 1, " cal/(mol K)")),
        ("dG37", _number(item.free_energy, 2, " kcal/mol")),
        ("Molecular weight", _number(item.optics.molecular_weight, 2, " g/mol")),
        ("Extinction, 260 nm", f"{item.optics.extinction:,.0f} 1/(M cm)"),
        ("At A260 = 1.000", f"{item.optics.micrograms_per_od:,.2f} ug/mL"),
        ("", f"{item.optics.nanomoles_per_od:,.2f} nmol/mL"),
        ("3' end stacking", _number(item.end_stability, 2, " kcal/mol")),
        ("3' GC clamp", f"{item.gc_clamp} of the last {CLAMP_WINDOW} bases"),
    ]
    parts = ['<dl class="metrics">']
    for name, value in rows:
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    return parts


def _flag_list(flags: Sequence[Flag]) -> list[str]:
    if not flags:
        return ['<p>No design criterion was breached.</p>']
    parts = ['<ul class="flags">']
    for flag in flags:
        parts.append(
            f"<li>{_tag(flag.severity)} {escape(flag.message)} "
            f"<code>{escape(flag.code)}</code></li>"
        )
    parts.append("</ul>")
    return parts


def _notable(items: Sequence, threshold) -> list:
    """The most stable structure, plus any other that reaches its threshold.

    A report that lists every register found buries the one that matters under
    near-neutral alternatives, so the runners-up are shown only where they too
    would be flagged.
    """
    if not items:
        return []
    kept = [items[0]]
    kept.extend(
        item for item in items[1:] if item.free_energy <= threshold(item)
    )
    return kept


def _hairpin_block(hairpins: Sequence[Hairpin]) -> list[str]:
    if not hairpins:
        return [
            "<p>No stem-loop of three or more base pairs can form, so no "
            "hairpin was scored.</p>"
        ]
    parts: list[str] = []
    for hairpin in hairpins:
        start, end = hairpin.loop
        parts.append(
            f'<p class="note">dG37 {hairpin.free_energy:,.2f} kcal/mol '
            f"&middot; Tm {_number(hairpin.melting, 1)} &deg;C "
            f"&middot; {hairpin.stem_length} bp stem "
            f"&middot; {hairpin.loop_length} nt loop at {start + 1}&ndash;{end}</p>"
        )
        parts.append(f'<div class="struct">{escape(hairpin.diagram)}</div>')
    return parts


def _dimer_block(dimers: Sequence[Dimer], label: str) -> list[str]:
    if not dimers:
        return [f"<p>No favourable {label} register was found.</p>"]
    parts: list[str] = []
    for dimer in dimers:
        if dimer.terminal_paired:
            extendable = ' <span class="tag warn">3\u2032 terminus paired</span>'
        elif dimer.three_prime_paired:
            extendable = ' <span class="tag warn">3\u2032 overhang</span>'
        else:
            extendable = ""
        helices = "" if dimer.helices == 1 else f" in {dimer.helices} helices"
        parts.append(
            f'<p class="note">dG37 {dimer.free_energy:,.2f} kcal/mol '
            f"&middot; Tm {_number(dimer.melting, 1)} &deg;C "
            f"&middot; {dimer.pairs} bp{helices}{extendable}</p>"
        )
        parts.append(f'<div class="struct">{escape(dimer.diagram)}</div>')
    return parts


def _dimer_threshold(criteria: Criteria):
    """Pick the threshold a dimer is judged against, by its 3'-end involvement."""

    def threshold(dimer: Dimer) -> float:
        return (
            criteria.dimer_three_prime_free_energy
            if dimer.three_prime_paired
            else criteria.dimer_free_energy
        )

    return threshold


def _oligo_section(item: OligoStats, criteria: Criteria) -> list[str]:
    parts = [f'<h2 class="oligo">{escape(item.name)}</h2>']
    if item.description:
        parts.append(f'<p class="lede">{escape(item.description)}</p>')
    parts.append(f'<div class="seq">5\u2032-{escape(_wrap(item.sequence))}-3\u2032</div>')

    parts.append("<h3>Measurements</h3>")
    parts.extend(_metrics_block(item))

    parts.append("<h3>Base composition</h3>")
    parts.extend(_composition_table(item))

    if not item.structures_screened:
        parts.append("<h3>Structures</h3>")
        parts.append(
            "<p>The hairpin and dimer screens were not run for this report.</p>"
        )
    else:
        parts.append("<h3>Hairpins</h3>")
        parts.extend(
            _hairpin_block(
                _notable(item.hairpins, lambda _: criteria.hairpin_free_energy)
            )
        )
        parts.append("<h3>Self-dimers</h3>")
        parts.extend(
            _dimer_block(
                _notable(item.self_dimers, _dimer_threshold(criteria)), "self-dimer"
            )
        )

    parts.append("<h3>Design flags</h3>")
    parts.extend(_flag_list(item.flags))
    return parts


def _cross_dimer_section(analysis: Analysis) -> list[str]:
    if not analysis.cross_dimers_screened:
        return []
    parts = ["<h2>Cross-dimers</h2>"]
    pairs = len(analysis.oligos) * (len(analysis.oligos) - 1) // 2
    if not analysis.cross_dimers:
        parts.append(
            f"<p>{pairs:,} pair{'' if pairs == 1 else 's'} screened; none reaches "
            "the dimer threshold.</p>"
        )
        return parts
    parts.append(
        f"<p>{pairs:,} pair{'' if pairs == 1 else 's'} screened; "
        f"{len(analysis.cross_dimers):,} at or below the threshold, most stable "
        "first.</p>"
    )
    for entry in analysis.cross_dimers:
        parts.append(
            f'<h3>{escape(entry.first)} &middot; {escape(entry.second)}</h3>'
        )
        parts.extend(_dimer_block([entry.dimer], "cross-dimer"))
    return parts


def render_html(analysis: Analysis, source: str, version: str) -> str:
    """Render a complete, self-contained HTML document.

    Args:
        analysis: the measured oligonucleotides, their cross-dimers and the
            conditions and criteria the run applied.
        source: name of the input shown in the header.
        version: tool version recorded in the footer.
    """
    oligos = analysis.oligos
    plural = "" if len(oligos) == 1 else "s"
    flagged = sum(1 for item in oligos if item.worst_severity == "warn")
    span = analysis.melting_span
    melting = (
        f" &middot; Tm {span[0]:,.1f}&ndash;{span[1]:,.1f} &deg;C"
        if span is not None and len(oligos) > 1
        else ""
    )

    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>Oligonucleotide analysis of {escape(source)}</title>",
        f"<style>{_STYLE}</style>",
        "</head>",
        "<body>",
        "<h1>Primer and oligonucleotide analysis</h1>",
        f'<p class="lede">{escape(source)} &middot; {len(oligos)} '
        f"oligonucleotide{plural}{melting} &middot; {flagged} flagged for review</p>",
        "<h2>Reaction conditions</h2>",
    ]
    parts.extend(_conditions(analysis))

    parts.append("<h2>Summary</h2>")
    parts.extend(_summary_table(oligos))
    parts.append("<h2>Design criteria</h2>")
    parts.extend(_criteria_block(analysis))

    parts.extend(_cross_dimer_section(analysis))

    for item in oligos:
        parts.extend(_oligo_section(item, analysis.criteria))

    parts.extend(
        [
            "<footer>",
            "<p>Melting temperature is the two-state nearest-neighbour value for "
            "the oligonucleotide against its perfect complement, using the unified "
            "DNA parameters of SantaLucia (1998), PNAS 95:1460, corrected to the "
            "stated ion composition by Owczarzy et al. (2004), Biochemistry "
            "43:3537 and, where magnesium is present, Owczarzy et al. (2008), "
            "Biochemistry 47:5336.</p>",
            "<p>The extinction coefficient is the single-strand nearest-neighbour "
            "value of Warshaw &amp; Cantor (1970), Biopolymers 9:1079; an annealed "
            "duplex absorbs less than the sum of its strands, so these factors "
            "must not be used to quantify annealed material. Molecular weight is "
            "anhydrous and excludes counterions.</p>",
            "<p>Hairpins are single contiguous stems closing one loop of at least "
            "three nucleotides; dimers are single antiparallel registers scored "
            "helix by helix. Bulges, internal loops and multibranch structures "
            "are not modelled, so these free energies rank candidates rather than "
            "predict a fold. In a structure diagram, <code>|</code> marks a scored "
            "base pair and <code>.</code> a complementary position whose helix was "
            "too short or too weak to count.</p>",
            "<p>A dimer marked 3&prime; terminus paired has a 3&prime;-terminal "
            "base within a scored helix and can be extended as it stands. One "
            "marked 3&prime; overhang has a scored pair near a 3&prime; end but "
            "not at it, so extension needs the terminal bases to pair as well; "
            "both are held to the stricter of the two dimer thresholds, because "
            "both are screened as extendable.</p>",
            f"<p>oligo-stats {escape(version)}</p>",
            "</footer>",
            "</body>",
            "</html>",
            "",
        ]
    )
    return "\n".join(parts)
