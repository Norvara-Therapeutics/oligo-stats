"""The self-contained HTML report."""

from oligo_stats import __version__
from oligo_stats.oligo import Analysis, Criteria, analyse_oligo, screen_cross_dimers
from oligo_stats.report import render_html
from oligo_stats.thermo import IonConditions

IONS = IonConditions(strand=0.25e-6, sodium=0.05, magnesium=3e-3, dntp=0.8e-3)


def _analysis(sequences, criteria=None, **kwargs):
    rules = Criteria() if criteria is None else criteria
    oligos = [
        analyse_oligo(sequence, name=name, ions=IONS, criteria=rules, **kwargs)
        for name, sequence in sequences
    ]
    return Analysis(
        oligos=tuple(oligos),
        ions=IONS,
        criteria=rules,
        terminus="hydroxyl",
        cross_dimers=tuple(screen_cross_dimers(oligos, IONS, rules)),
        cross_dimers_screened=len(oligos) > 1,
    )


PANEL = [
    ("phiX_F", "TCATGAAGGATGGTGTTAATGC"),
    ("phiX_hairpin", "CATGCCTCCAAATCTTGGAGGC"),
]


def test_report_is_a_complete_html_document():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")
    assert html.count("<html") == 1
    assert html.count("</body>") == 1


def test_report_loads_nothing_from_the_network_and_runs_no_script():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    assert "<script" not in html
    assert "http://" not in html
    assert "https://" not in html
    assert "<link" not in html
    assert "<style>" in html


def test_report_is_byte_deterministic():
    first = render_html(_analysis(PANEL), source="primers.fasta", version="1.0.0")
    second = render_html(_analysis(PANEL), source="primers.fasta", version="1.0.0")
    assert first == second


def test_report_records_the_source_the_version_and_the_conditions():
    html = render_html(_analysis(PANEL), source="primers.fasta", version="9.9.9")
    assert "primers.fasta" in html
    assert "oligo-stats 9.9.9" in html
    assert "mixed regime" in html
    assert "Free magnesium" in html


def test_report_carries_every_oligonucleotide_and_its_measurements():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    for name, sequence in PANEL:
        assert name in html
        assert sequence in html
    assert "Extinction, 260 nm" in html
    assert "Base composition" in html


def test_report_shows_structures_and_their_legend():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    assert 'class="struct"' in html
    assert "marks a scored " in html
    assert "Hairpins" in html
    assert "Self-dimers" in html


def test_report_marks_a_flagged_oligonucleotide_for_review():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    assert "review" in html
    assert "pass" in html


def test_report_states_how_many_pairs_the_cross_dimer_screen_covered():
    html = render_html(_analysis(PANEL), source="primers.fasta", version=__version__)
    assert "1 pair screened" in html


def test_report_lists_a_cross_dimer_when_one_reaches_the_threshold():
    analysis = _analysis(PANEL, criteria=Criteria(dimer_free_energy=-1.0))
    html = render_html(analysis, source="primers.fasta", version=__version__)
    assert "at or below the threshold" in html
    assert "phiX_F &middot; phiX_hairpin" in html


def test_report_records_every_criterion_exactly_as_configured():
    """The report is the record of which rules produced its flags."""
    criteria = Criteria(
        gc_min=42.5,
        gc_max=60.5,
        melting_min=55.5,
        melting_max=65.5,
        end_stability=-7.5,
        max_run=6,
    )
    html = render_html(
        _analysis(PANEL, criteria=criteria), source="p.fasta", version=__version__
    )
    for text in ("42.5-60.5 %", "55.5-65.5 C", "-7.5 kcal/mol", "6 bases or more"):
        assert text in html


def test_every_criterion_appears_in_the_report():
    """A threshold added to the rules cannot be left out of the document."""
    import dataclasses

    from oligo_stats.report import _criteria_rows

    rows = _criteria_rows(Criteria())
    assert len(rows) == 9
    # _criteria_rows raises when a Criteria field is unaccounted for.
    assert {field.name for field in dataclasses.fields(Criteria)}


def test_a_skipped_screen_is_not_summarised_as_a_clean_one():
    html = render_html(
        _analysis(PANEL[:1], screen_structures=False),
        source="p.fasta",
        version=__version__,
    )
    assert html.count("not screened") == 2
    assert "part screened" in html
    assert '<span class="tag">pass</span>' not in html


def test_a_screened_clean_oligonucleotide_still_reads_as_passing():
    html = render_html(
        _analysis([PANEL[0]]), source="p.fasta", version=__version__
    )
    assert '<span class="tag">pass</span>' in html
    assert "not screened" not in html


def test_report_escapes_names_and_descriptions():
    oligos = (
        analyse_oligo(
            "TCATGAAGGATGGTGTTAATGC",
            name="<script>bad</script>",
            description="a & b",
            ions=IONS,
        ),
    )
    analysis = Analysis(oligos=oligos, ions=IONS, criteria=Criteria(), terminus="hydroxyl")
    html = render_html(analysis, source="<b>x</b>", version=__version__)
    assert "<script>bad" not in html
    assert "&lt;script&gt;bad" in html
    assert "a &amp; b" in html
    assert "&lt;b&gt;x&lt;/b&gt;" in html


def test_report_handles_a_single_oligonucleotide_without_a_cross_dimer_section():
    analysis = _analysis(PANEL[:1])
    html = render_html(analysis, source="one.fasta", version=__version__)
    assert "1 oligonucleotide " in html
    assert "<h2>Cross-dimers</h2>" not in html


def test_report_reports_an_absent_structure_screen():
    analysis = _analysis(PANEL[:1], screen_structures=False)
    html = render_html(analysis, source="one.fasta", version=__version__)
    assert "screens were not run" in html
    assert "<h3>Hairpins</h3>" not in html
