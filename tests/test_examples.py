"""The bundled example panel."""

from pathlib import Path

import pytest

from oligo_stats.inputs import read
from oligo_stats.oligo import analyse_oligo
from oligo_stats.thermo import AMPLIFICATION

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
#: The example panel is documented against the default amplification buffer.
PCR = AMPLIFICATION


@pytest.fixture(scope="module")
def panel():
    return read(EXAMPLES / "primers.fasta")


def test_the_example_panel_is_six_named_twenty_two_mers(panel):
    assert len(panel) == 6
    assert [item.name for item in panel] == [
        "phiX_F",
        "phiX_R",
        "phiX_hairpin",
        "phiX_dimer",
        "phiX_gcrich",
        "phiX_atrich",
    ]
    assert all(len(item.sequence) == 22 for item in panel)
    assert all("NC_001422.1" in item.description for item in panel)


def test_the_clean_pair_is_matched_and_unflagged(panel):
    forward, reverse = panel[0], panel[1]
    first = analyse_oligo(forward.sequence, ions=PCR)
    second = analyse_oligo(reverse.sequence, ions=PCR)
    assert first.flags == ()
    assert second.flags == ()
    assert abs(first.melting - second.melting) < 0.1


def test_each_demonstration_oligonucleotide_shows_what_it_is_named_for(panel):
    measured = {item.name: analyse_oligo(item.sequence, ions=PCR) for item in panel}
    assert measured["phiX_hairpin"].best_hairpin.stem_length == 8
    assert "hairpin" in [flag.code for flag in measured["phiX_hairpin"].flags]
    assert measured["phiX_dimer"].best_self_dimer.three_prime_paired is True
    assert "self_dimer" in [flag.code for flag in measured["phiX_dimer"].flags]
    assert measured["phiX_gcrich"].composition.gc_percent > 70
    assert measured["phiX_atrich"].gc_clamp == 0


def test_the_order_sheet_example_parses_as_a_table():
    oligos = read(EXAMPLES / "order_sheet.tsv")
    assert [item.name for item in oligos] == ["phiX_F", "phiX_R"]
    assert all(len(item.sequence) == 22 for item in oligos)
