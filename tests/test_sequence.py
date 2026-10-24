"""Normalisation, composition and degeneracy."""

import pytest

from oligo_stats.sequence import (
    SequenceError,
    composition,
    degenerate_positions,
    expand,
    expansion_count,
    has_uracil,
    is_self_complementary,
    longest_runs,
    normalise,
    reverse_complement,
    validate,
)


def test_reverse_complement_round_trips():
    assert reverse_complement("ACGTACGT") == "ACGTACGT"
    assert reverse_complement("AAACCC") == "GGGTTT"
    assert reverse_complement(reverse_complement("GATTACA")) == "GATTACA"


def test_reverse_complement_handles_ambiguity_and_case():
    assert reverse_complement("acgtRY") == "RYacgt"
    assert reverse_complement("NNN") == "NNN"


def test_normalise_strips_whitespace_and_transcribes_uracil():
    assert normalise("acg t\nacg") == "ACGTACG"
    assert normalise("AUGC") == "ATGC"
    assert has_uracil("augc") is True
    assert has_uracil("ATGC") is False


def test_validate_rejects_empty_and_unknown_characters():
    with pytest.raises(SequenceError, match="empty"):
        validate("   ")
    with pytest.raises(SequenceError, match="not a nucleotide base"):
        validate("ACGTX")


def test_validate_rejects_ambiguity_by_default_and_points_at_expand():
    with pytest.raises(SequenceError, match="expand"):
        validate("ACGTN")
    assert validate("ACGTN", allow_degenerate=True) == "ACGTN"


def test_degenerate_positions_are_zero_based():
    assert degenerate_positions("ACGTN") == [4]
    assert degenerate_positions("RACGT") == [0]
    assert degenerate_positions("ACGT") == []


def test_expansion_count_multiplies_degeneracies():
    assert expansion_count("ACGT") == 1
    assert expansion_count("N") == 4
    assert expansion_count("RY") == 4
    assert expansion_count("NNN") == 64


def test_expand_enumerates_in_input_order():
    assert expand("AR") == ["AA", "AG"]
    assert expand("NN")[:4] == ["AA", "AC", "AG", "AT"]
    assert len(expand("NNN")) == 64


def test_expand_refuses_to_exceed_its_limit():
    with pytest.raises(SequenceError, match="above the limit"):
        expand("NNNN", limit=10)


def test_composition_counts_and_gc():
    result = composition("AACCGGTT")
    assert result.length == 8
    assert result.counts == {"A": 2, "C": 2, "G": 2, "T": 2}
    assert result.gc_percent == 50.0
    assert result.gc_count == 4
    assert result.degenerate == 0


def test_gc_percent_is_stored_unrounded():
    """It is compared against a caller's threshold, so it must not be rounded."""
    result = composition("GGGCCCAAAAAAAA")
    assert result.gc_percent == pytest.approx(600.0 / 14.0, abs=1e-12)
    assert result.gc_percent != round(result.gc_percent, 2)


def test_composition_excludes_ambiguity_from_gc():
    result = composition("GCNN")
    assert result.degenerate == 2
    assert result.gc_percent == 100.0


def test_composition_of_fully_ambiguous_sequence_has_no_gc():
    assert composition("NNNN").gc_percent is None


def test_self_complementarity():
    assert is_self_complementary("GAATTC") is True
    assert is_self_complementary("ACGT") is True
    assert is_self_complementary("ACGTA") is False
    assert is_self_complementary("AAAA") is False


def test_longest_runs_reports_every_base_separately():
    assert longest_runs("AAGGGT") == {"A": 2, "G": 3, "T": 1}
    assert longest_runs("ACGT") == {"A": 1, "C": 1, "G": 1, "T": 1}
    assert longest_runs("") == {}
    assert longest_runs("TTTTAAAAA") == {"T": 4, "A": 5}


def test_longest_runs_keeps_the_longest_run_of_a_base_that_recurs():
    """A base that runs twice is credited with its longer run, not its last."""
    assert longest_runs("GGGGATTGG")["G"] == 4
    assert longest_runs("GGATTGGGG")["G"] == 4
