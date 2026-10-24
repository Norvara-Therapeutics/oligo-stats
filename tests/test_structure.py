"""Hairpin and dimer searches."""

import pytest

from oligo_stats.structure import (
    find_dimers,
    find_hairpins,
    find_self_dimers,
    perfect_complement,
)
from oligo_stats.thermo import IonConditions

IONS = IonConditions(strand=0.25e-6, sodium=0.05)


def test_a_designed_stem_loop_is_found_with_the_expected_stem_and_loop():
    # Six complementary bases either side of a four-base loop.
    hairpins = find_hairpins("GCGCGCTTTTGCGCGC", IONS)
    best = hairpins[0]
    assert best.stem_length == 6
    assert best.loop_length == 4
    assert best.loop == (6, 10)
    assert best.free_energy < -4.0


def test_a_sequence_with_no_possible_stem_yields_nothing():
    assert find_hairpins("AAAAAAAAAAAAAAA", IONS) == []


def test_hairpin_stems_shorter_than_the_minimum_are_not_reported():
    assert find_hairpins("GCGCGCTTTTGCGCGC", IONS, min_stem=7) == []
    with pytest.raises(ValueError, match="min_stem"):
        find_hairpins("GCGCGCTTTTGCGCGC", IONS, min_stem=0)


def test_hairpins_are_returned_most_stable_first_and_one_per_loop():
    hairpins = find_hairpins("CATGCCTCCAAATCTTGGAGGC", IONS, limit=5)
    energies = [item.free_energy for item in hairpins]
    assert energies == sorted(energies)
    loops = [item.loop for item in hairpins]
    assert len(loops) == len(set(loops))


def test_hairpin_loop_never_closes_below_three_nucleotides():
    for hairpin in find_hairpins("CATGCCTCCAAATCTTGGAGGC", IONS, limit=20):
        assert hairpin.loop_length >= 3


def test_hairpin_diagram_rows_describe_the_structure():
    hairpin = find_hairpins("GCGCGCTTTTGCGCGC", IONS)[0]
    rows = hairpin.diagram.splitlines()
    assert len(rows) == 3
    assert rows[0].startswith("5'-") and rows[0].endswith("-\\")
    assert rows[2].startswith("3'-") and rows[2].endswith("-/")
    assert rows[1].count("|") == hairpin.stem_length
    assert rows[1].endswith("TTTT")


def test_a_longer_stem_is_more_stable_than_a_shorter_one():
    short = find_hairpins("GCGTTTTGCGC", IONS)[0]
    long = find_hairpins("GCGCGCTTTTGCGCGC", IONS)[0]
    assert long.free_energy < short.free_energy


def test_a_palindrome_self_dimerises_over_its_whole_length():
    dimer = find_self_dimers("GAATTC", IONS)[0]
    assert dimer.pairs == 6
    assert dimer.shift == 0
    assert dimer.three_prime_paired is True


def test_self_dimer_of_a_sequence_with_no_complementarity_is_absent():
    assert find_self_dimers("AAAAAAAA", IONS) == []


def test_a_complementary_three_prime_end_is_detected_as_extendable():
    # The last six bases are a palindrome, so two copies pair 3' end to 3' end.
    dimers = find_self_dimers("TTTTTTTTTTGAATTC", IONS)
    assert dimers[0].three_prime_paired is True
    assert dimers[0].pairs == 6


def test_a_dimer_away_from_the_three_prime_end_is_not_extendable():
    dimers = find_self_dimers("GAATTCTTTTTTTTTT", IONS)
    assert dimers[0].pairs == 6
    assert dimers[0].three_prime_paired is False


def test_dimers_are_returned_most_stable_first():
    dimers = find_dimers("CATGCCTCCAAATCTTGGAGGC", "GAGGCGGTCAAAAAGCCGCCTC", IONS, limit=5)
    energies = [item.free_energy for item in dimers]
    assert energies == sorted(energies)
    assert all(item.free_energy < 0 for item in dimers)


def test_every_reported_register_pairs_at_least_two_bases():
    for dimer in find_dimers("ACGTCCTGATCAAGCTTAGG", "CCTAAGCTTGATCAGGACGT", IONS, limit=20):
        assert dimer.pairs >= 2
        assert dimer.helices >= 1


def test_the_perfect_complement_is_the_most_stable_register():
    sequence = "TCATGAAGGATGGTGTTAATGC"
    best = find_dimers(sequence, perfect_complement(sequence), IONS)[0]
    assert best.pairs == len(sequence)
    assert best.helices == 1
    assert best.shift == 0


def test_dimer_diagram_aligns_the_two_strands():
    sequence = "TCATGAAGGATGGTGTTAATGC"
    dimer = find_dimers(sequence, perfect_complement(sequence), IONS)[0]
    rows = dimer.diagram.splitlines()
    assert len(rows) == 3
    assert rows[0] == f"5'-{sequence}-3'"
    # The lower row is drawn 3'->5', so it reads as the base-by-base complement.
    complement = perfect_complement(sequence)[::-1]
    assert rows[2] == f"3'-{complement}-5'"
    assert rows[1].strip() == "|" * len(sequence)


def test_an_unscored_contact_is_marked_separately_from_a_scored_pair():
    """A lone complementary base cannot stack, so it is shown but not counted."""
    dimers = find_self_dimers("CATGCCTCCAAATCTTGGAGGC", IONS, limit=10)
    marks = "".join(item.diagram.splitlines()[1] for item in dimers)
    assert "." in marks
    assert "|" in marks


#: A sequence whose most stable self-dimer leaves both 3' ends free while a
#: weaker register pairs near one of them. The weaker register is the one that
#: breaches, because it is judged against the stricter threshold.
THREE_PRIME_TRAP = "GGGCCCATGCTGGCATCTGCGA"


def test_the_most_stable_three_prime_register_survives_the_limit():
    """It is judged against its own threshold, so it must not be truncated away."""
    registers = find_self_dimers(THREE_PRIME_TRAP, IONS, limit=1)
    assert registers[0].three_prime_paired is False
    assert any(item.three_prime_paired for item in registers)
    strongest_three_prime = min(
        (item for item in find_self_dimers(THREE_PRIME_TRAP, IONS, limit=200)
         if item.three_prime_paired),
        key=lambda item: item.free_energy,
    )
    assert strongest_three_prime.shift in [
        item.shift for item in registers if item.three_prime_paired
    ]


def test_a_register_without_a_three_prime_contact_adds_nothing():
    sequence = "GAATTCTTTTTTTTTT"
    registers = find_self_dimers(sequence, IONS, limit=1)
    assert len(registers) == 1
    assert registers[0].three_prime_paired is False


def test_terminal_pairing_is_distinguished_from_mere_proximity():
    """The screen is wider than the claim that a polymerase can extend."""
    proximal = next(
        item for item in find_self_dimers(THREE_PRIME_TRAP, IONS, limit=200)
        if item.three_prime_paired
    )
    assert proximal.three_prime_paired is True
    assert proximal.terminal_paired is False
    sequence = "TCATGAAGGATGGTGTTAATGC"
    blunt = find_dimers(sequence, perfect_complement(sequence), IONS)[0]
    assert blunt.terminal_paired is True


def test_a_palindrome_pairs_its_terminal_bases():
    dimer = find_self_dimers("GAATTC", IONS)[0]
    assert dimer.terminal_paired is True


def test_identical_strands_raise_the_dimer_melting_temperature():
    sequence = "GCGCGCGCGCGCGCGC"
    together = find_dimers(sequence, sequence, IONS, identical_strands=True)[0]
    apart = find_dimers(sequence, sequence, IONS, identical_strands=False)[0]
    assert together.melting > apart.melting


def test_perfect_complement_is_the_reverse_complement():
    assert perfect_complement("ACGTT") == "AACGT"
