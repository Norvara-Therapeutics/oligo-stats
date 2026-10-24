"""Assembly of the measurement set and the design flags."""

import pytest

from oligo_stats.oligo import (
    Criteria,
    analyse_oligo,
    gc_clamp,
    screen_cross_dimers,
    worst_dimer,
)
from oligo_stats.sequence import SequenceError
from oligo_stats.thermo import IonConditions

IONS = IonConditions(strand=0.25e-6, sodium=0.05, magnesium=3e-3, dntp=0.8e-3)
CLEAN = "TCATGAAGGATGGTGTTAATGC"


def _codes(sequence, **kwargs):
    return [flag.code for flag in analyse_oligo(sequence, ions=IONS, **kwargs).flags]


def test_a_clean_primer_raises_no_flag():
    stats = analyse_oligo(CLEAN, name="phiX_F", ions=IONS)
    assert stats.flags == ()
    assert stats.worst_severity == "ok"
    assert stats.name == "phiX_F"
    assert stats.length == 22


def test_measurements_are_carried_through_from_each_module():
    stats = analyse_oligo(CLEAN, ions=IONS)
    assert stats.composition.gc_percent == pytest.approx(40.91, abs=0.01)
    assert stats.optics.molecular_weight > 6000
    assert stats.melting is not None
    assert stats.free_energy < 0
    assert stats.end_stability < 0
    assert stats.gc_clamp == 2


def test_input_is_normalised_before_measurement():
    spaced = analyse_oligo("tcatg aagga tggtg ttaatgc", ions=IONS)
    assert spaced.sequence == CLEAN
    assert spaced.rna_input is False


def test_uracil_is_accepted_and_noted():
    stats = analyse_oligo("ACGUACGUACGUACGUACGU", ions=IONS)
    assert stats.rna_input is True
    assert "rna_input" in [flag.code for flag in stats.flags]
    assert stats.sequence == "ACGTACGTACGTACGTACGT"


def test_ambiguity_codes_are_refused():
    with pytest.raises(SequenceError, match="expand"):
        analyse_oligo("ACGTNACGTACGTACGTACG", ions=IONS)


def test_gc_content_outside_the_window_is_flagged():
    assert "gc_content" in _codes("AATATCTATAACAACTATTTTA")
    assert "gc_content" in _codes("ACCACCGCCCCGAAGGGGACGA")


def test_gc_content_is_judged_on_the_exact_value_not_the_displayed_one():
    """42.857% must not clear a 42.859% minimum because it prints as 42.86."""
    sequence = "GGGCCCAAAAAAAA"
    just_above = Criteria(gc_min=42.859, gc_max=60.0)
    just_below = Criteria(gc_min=42.855, gc_max=60.0)
    assert "gc_content" in [
        flag.code
        for flag in analyse_oligo(sequence, ions=IONS, criteria=just_above).flags
    ]
    assert "gc_content" not in [
        flag.code
        for flag in analyse_oligo(sequence, ions=IONS, criteria=just_below).flags
    ]


def test_an_absent_melting_temperature_is_flagged_rather_than_ignored():
    from oligo_stats.thermo import IonConditions

    # 1e-31 M: far below any real reaction, which is the point - the two-state
    # relation has no root above freezing there.
    dilute = IonConditions(strand=1e-31)
    stats = analyse_oligo("TCATGAAGGATGGTGTTAATGC", ions=dilute)
    assert stats.melting is None
    assert stats.worst_severity == "warn"
    flag = [item for item in stats.flags if item.code == "melting_temperature"][0]
    assert "could not be applied" in flag.message


def test_melting_temperature_outside_the_window_is_flagged():
    assert "melting_temperature" in _codes("AATATCTATAACAACTATTTTA")


def test_a_missing_gc_clamp_is_flagged_and_a_heavy_one_is_noted():
    stats = analyse_oligo("AATATCTATAACAACTATTTTA", ions=IONS)
    clamp = [flag for flag in stats.flags if flag.code == "gc_clamp"]
    assert clamp and clamp[0].severity == "warn"
    heavy = analyse_oligo("CATGCCTCCAAATCTTGGAGGC", ions=IONS)
    note = [flag for flag in heavy.flags if flag.code == "gc_clamp"]
    assert note and note[0].severity == "note"


def test_a_poly_g_run_warns_while_other_runs_only_note():
    runs = [
        flag for flag in analyse_oligo("ACCACCGCCCCGAAGGGGACGA", ions=IONS).flags
        if flag.code == "base_run"
    ]
    # Both a CCCC and a GGGG run are present; each is reported on its own.
    assert len(runs) == 2
    assert [flag.severity for flag in runs] == ["note", "warn"]
    assert "G bases" in runs[1].message
    adenine = [
        flag for flag in analyse_oligo("GAGGCGGTCAAAAAGCCGCCTC", ions=IONS).flags
        if flag.code == "base_run"
    ]
    assert adenine and adenine[0].severity == "note"


def test_a_strong_hairpin_and_self_dimer_are_flagged():
    codes = _codes("CATGCCTCCAAATCTTGGAGGC")
    assert "hairpin" in codes
    assert "self_dimer" in codes


def test_a_qualifying_three_prime_register_is_flagged_though_it_is_not_the_strongest():
    """The strongest register passes its threshold; a weaker 3' one does not."""
    sequence = "GGGCCCATGCTGGCATCTGCGA"
    stats = analyse_oligo(sequence, ions=IONS)
    assert stats.best_self_dimer.three_prime_paired is False
    assert stats.best_self_dimer.free_energy > Criteria().dimer_free_energy
    flag = [item for item in stats.flags if item.code == "self_dimer"]
    assert flag, "a register breaching the 3' threshold was not flagged"
    assert flag[0].message.startswith("self-dimer of -5.73")


def test_a_flag_for_a_proximal_register_does_not_claim_direct_extension():
    stats = analyse_oligo("GGGCCCATGCTGGCATCTGCGA", ions=IONS)
    message = [item for item in stats.flags if item.code == "self_dimer"][0].message
    assert "unpaired 3' overhang" in message
    assert "polymerase can extend" not in message


def test_a_flag_for_a_terminal_register_does_claim_direct_extension():
    stats = analyse_oligo("CATGCCTCCAAATCTTGGAGGC", ions=IONS)
    message = [item for item in stats.flags if item.code == "self_dimer"][0].message
    assert "3'-terminal base" in message


def test_worst_dimer_picks_the_breaching_register_in_either_class():
    from oligo_stats.structure import find_self_dimers

    criteria = Criteria()
    registers = find_self_dimers("GGGCCCATGCTGGCATCTGCGA", IONS, limit=3)
    offender = worst_dimer(registers, criteria)
    assert offender is not None
    assert offender.three_prime_paired is True
    assert worst_dimer([], criteria) is None
    assert worst_dimer(registers, Criteria(dimer_free_energy=-100.0,
                                           dimer_three_prime_free_energy=-100.0)) is None


def test_cross_dimer_screen_sees_a_qualifying_three_prime_register_too():
    oligos = [
        analyse_oligo("GGGCCCATGCTGGCATCTGCGA", name="a", ions=IONS),
        analyse_oligo("TCGCAGATGCCAGCATGGGCCC", name="b", ions=IONS),
    ]
    found = screen_cross_dimers(oligos, IONS, Criteria())
    assert found
    assert found[0].dimer.free_energy <= Criteria().dimer_three_prime_free_energy


def test_a_self_complementary_oligo_is_flagged():
    codes = _codes("GGGGGGCCCCCC")
    assert "self_complementary" in codes


def test_a_length_outside_the_window_is_noted_not_warned():
    stats = analyse_oligo("ACGTCCTGATCA", ions=IONS)
    length = [flag for flag in stats.flags if flag.code == "length"]
    assert length and length[0].severity == "note"


def test_criteria_are_honoured():
    permissive = Criteria(gc_min=0.0, gc_max=100.0, melting_min=0.0, melting_max=100.0)
    codes = [
        flag.code
        for flag in analyse_oligo(
            "ACCACCGCCCCGAAGGGGACGA", ions=IONS, criteria=permissive
        ).flags
    ]
    assert "gc_content" not in codes
    assert "melting_temperature" not in codes


def test_turning_off_the_structure_screen_removes_its_flags():
    stats = analyse_oligo("CATGCCTCCAAATCTTGGAGGC", ions=IONS, screen_structures=False)
    assert stats.hairpins == ()
    assert stats.self_dimers == ()
    assert stats.best_hairpin is None
    assert stats.best_self_dimer is None
    assert "hairpin" not in [flag.code for flag in stats.flags]


def test_the_structure_limit_bounds_what_is_retained():
    stats = analyse_oligo("CATGCCTCCAAATCTTGGAGGC", ions=IONS, structures=1)
    assert len(stats.hairpins) == 1
    assert len(stats.self_dimers) == 1


def test_a_phosphate_terminus_changes_the_mass_only():
    plain = analyse_oligo(CLEAN, ions=IONS)
    phosphorylated = analyse_oligo(CLEAN, ions=IONS, terminus="phosphate")
    assert phosphorylated.optics.molecular_weight > plain.optics.molecular_weight
    assert phosphorylated.optics.extinction == plain.optics.extinction
    assert phosphorylated.melting == plain.melting


def test_gc_clamp_counts_only_the_last_bases():
    assert gc_clamp("AAAAAGGGGG") == 5
    assert gc_clamp("GGGGGAAAAA") == 0
    assert gc_clamp("AAAAAGGGGG", window=2) == 2


def test_cross_dimer_screen_reports_pairs_below_threshold_most_stable_first():
    oligos = [
        analyse_oligo("CATGCCTCCAAATCTTGGAGGC", name="a", ions=IONS),
        analyse_oligo("GAGGCGGTCAAAAAGCCGCCTC", name="b", ions=IONS),
        analyse_oligo(CLEAN, name="c", ions=IONS),
    ]
    found = screen_cross_dimers(oligos, IONS, Criteria(dimer_free_energy=-2.0))
    assert found
    energies = [entry.dimer.free_energy for entry in found]
    assert energies == sorted(energies)
    assert all(entry.first != entry.second for entry in found)


def test_every_criterion_is_applied_somewhere():
    """A declared threshold that nothing reads is a promise the tool breaks."""
    import dataclasses
    import pathlib

    source = (
        pathlib.Path(__file__).resolve().parent.parent / "src" / "oligo_stats"
    )
    applied = "\n".join(path.read_text() for path in source.glob("*.py"))
    for field in dataclasses.fields(Criteria):
        assert applied.count(f"criteria.{field.name}") + applied.count(
            f"rules.{field.name}"
        ), f"Criteria.{field.name} is never read"


def test_a_duplicated_sequence_in_a_panel_scores_as_identical_strands():
    oligos = [
        analyse_oligo("CATGCCTCCAAATCTTGGAGGC", name="a", ions=IONS),
        analyse_oligo("CATGCCTCCAAATCTTGGAGGC", name="b", ions=IONS),
    ]
    found = screen_cross_dimers(oligos, IONS, Criteria())
    assert found
    # The same register scored as one molecule melts above the two-strand form.
    assert found[0].dimer.melting == oligos[0].best_self_dimer.melting


def test_cross_dimer_screen_is_empty_when_nothing_reaches_the_threshold():
    oligos = [
        analyse_oligo(CLEAN, name="a", ions=IONS),
        analyse_oligo("AATATCTATAACAACTATTTTA", name="b", ions=IONS),
    ]
    assert screen_cross_dimers(oligos, IONS, Criteria()) == []
