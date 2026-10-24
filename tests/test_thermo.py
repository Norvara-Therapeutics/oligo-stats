"""Nearest-neighbour parameters, melting temperature and ion corrections."""

import math

import pytest

from oligo_stats.sequence import reverse_complement
from oligo_stats.thermo import (
    AMPLIFICATION,
    INIT_ENTHALPY,
    INIT_ENTROPY,
    NN_ENTHALPY,
    NN_ENTROPY,
    SUPPLIER,
    SYMMETRY_ENTROPY,
    T37,
    IonConditions,
    ThermoError,
    duplex,
    loop_free_energy,
    melting_temperature,
    stack,
    three_prime_stability,
)

#: Free energies at 37 C published alongside the enthalpies and entropies in
#: SantaLucia (1998), Table 1. They are an independent column of that table, so
#: deriving them from the tabulated enthalpy and entropy checks both for
#: transcription errors. The tolerance is the paper's own rounding.
PUBLISHED_FREE_ENERGY = {
    "AA": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "GT": -1.44,
    "CT": -1.28,
    "GA": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
}


def _monovalent(**ions) -> IonConditions:
    """Ion conditions holding only what the caller names.

    The package defaults are an amplification buffer, so a test about one ion
    has to clear the others to vary it alone.
    """
    settings = {"potassium": 0.0, "tris": 0.0, "magnesium": 0.0, "dntp": 0.0}
    settings.update(ions)
    return IonConditions(**settings)


def _free_energy(enthalpy: float, entropy: float) -> float:
    return enthalpy - T37 * entropy / 1000.0


@pytest.mark.parametrize("step,published", sorted(PUBLISHED_FREE_ENERGY.items()))
def test_stacking_parameters_reproduce_published_free_energies(step, published):
    derived = _free_energy(NN_ENTHALPY[step], NN_ENTROPY[step])
    assert derived == pytest.approx(published, abs=0.02)


def test_initiation_and_symmetry_reproduce_published_free_energies():
    assert _free_energy(INIT_ENTHALPY["GC"], INIT_ENTROPY["GC"]) == pytest.approx(
        0.98, abs=0.02
    )
    assert _free_energy(INIT_ENTHALPY["AT"], INIT_ENTROPY["AT"]) == pytest.approx(
        1.03, abs=0.02
    )
    assert _free_energy(0.0, SYMMETRY_ENTROPY) == pytest.approx(0.43, abs=0.02)


@pytest.mark.parametrize("step", sorted(NN_ENTHALPY))
def test_stacking_tables_are_strand_symmetric(step):
    """A duplex read from either strand must be scored identically."""
    mirror = reverse_complement(step)
    assert NN_ENTHALPY[step] == NN_ENTHALPY[mirror]
    assert NN_ENTROPY[step] == NN_ENTROPY[mirror]


def test_stacking_tables_cover_every_dinucleotide():
    expected = {first + second for first in "ACGT" for second in "ACGT"}
    assert set(NN_ENTHALPY) == expected
    assert set(NN_ENTROPY) == expected


def test_stack_sums_steps_and_both_initiations():
    helix = stack("GC")
    expected_enthalpy = (
        NN_ENTHALPY["GC"] + 2 * INIT_ENTHALPY["GC"]
    )
    assert helix.enthalpy == pytest.approx(expected_enthalpy)
    assert helix.base_pairs == 2
    assert helix.gc_fraction == 1.0


def test_stack_of_one_base_has_no_interaction():
    helix = stack("A")
    assert helix.enthalpy == 0.0
    assert helix.entropy == 0.0
    assert helix.base_pairs == 1


def test_stack_rejects_ambiguity():
    with pytest.raises(ThermoError, match="no nearest-neighbour parameter"):
        stack("ACGN")


def test_duplex_applies_the_symmetry_correction_only_when_symmetric():
    palindrome = duplex("GAATTC")
    assert palindrome.symmetric is True
    assert palindrome.entropy == pytest.approx(stack("GAATTC").entropy + SYMMETRY_ENTROPY)
    ordinary = duplex("GAATTCA")
    assert ordinary.symmetric is False
    assert ordinary.entropy == pytest.approx(stack("GAATTCA").entropy)


def test_duplex_free_energy_is_negative_for_a_real_primer():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    assert thermo.free_energy() < 0
    assert thermo.stable is True


def test_melting_temperature_matches_the_two_state_relation_at_one_molar_sodium():
    """At 1 M Na+ the ion correction is the identity, so the raw relation shows."""
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    ions = _monovalent(strand=0.25e-6, sodium=1.0)
    gas_constant = 1.9872
    expected = thermo.enthalpy / (
        thermo.entropy / 1000.0 + gas_constant * math.log(ions.strand / 4.0) / 1000.0
    ) - 273.15
    assert melting_temperature(thermo, ions) == pytest.approx(expected, abs=1e-9)


def test_identical_strands_drop_the_factor_of_four():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    ions = _monovalent(strand=0.25e-6, sodium=1.0)
    apart = melting_temperature(thermo, ions)
    together = melting_temperature(thermo, ions, identical_strands=True)
    assert together > apart


def test_hairpin_melting_is_independent_of_concentration():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    dilute = _monovalent(strand=1e-9, sodium=0.05)
    concentrated = _monovalent(strand=1e-5, sodium=0.05)
    assert melting_temperature(thermo, dilute, unimolecular=True) == pytest.approx(
        melting_temperature(thermo, concentrated, unimolecular=True)
    )


def test_melting_temperature_rises_with_salt():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    low = melting_temperature(thermo, _monovalent(sodium=0.01))
    high = melting_temperature(thermo, _monovalent(sodium=0.5))
    assert high > low


def test_melting_temperature_rises_with_gc_content():
    ions = AMPLIFICATION
    at_rich = melting_temperature(duplex("ATATATATATATATATATAT"), ions)
    gc_rich = melting_temperature(duplex("GCGCGCGCGCGCGCGCGCGC"), ions)
    assert gc_rich > at_rich


def test_melting_temperature_rises_with_length():
    ions = AMPLIFICATION
    short = melting_temperature(duplex("ACGTCCTGATCA"), ions)
    long = melting_temperature(duplex("ACGTCCTGATCAAGCTTAGG"), ions)
    assert long > short


def test_magnesium_raises_melting_temperature_and_dntp_offsets_it():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    plain = melting_temperature(thermo, _monovalent(sodium=0.05))
    magnesium = melting_temperature(
        thermo, _monovalent(sodium=0.05, magnesium=3e-3)
    )
    chelated = melting_temperature(
        thermo, _monovalent(sodium=0.05, magnesium=3e-3, dntp=2.5e-3)
    )
    assert magnesium > plain
    assert plain < chelated < magnesium


def test_structure_that_cannot_melt_above_freezing_has_no_melting_point():
    ions = _monovalent(strand=1e-12, sodium=0.05)
    assert melting_temperature(duplex("ATAT"), ions) is None


def test_ion_regimes_follow_the_published_boundaries():
    assert _monovalent(sodium=0.05).regime == "monovalent"
    assert _monovalent(sodium=0.05, magnesium=1e-5).regime == "monovalent"
    assert _monovalent(sodium=0.05, magnesium=3e-3).regime == "mixed"
    assert _monovalent(sodium=1e-4, magnesium=3e-3).regime == "divalent"
    assert _monovalent(sodium=0.0, magnesium=3e-3).regime == "divalent"


def test_free_magnesium_never_goes_negative():
    ions = _monovalent(sodium=0.05, magnesium=1e-3, dntp=5e-3)
    assert ions.free_magnesium == 0.0
    assert ions.regime == "monovalent"


def test_tris_counts_as_half_monovalent():
    assert _monovalent(sodium=0.0, tris=0.02).monovalent == pytest.approx(0.01)
    assert _monovalent(sodium=0.01, potassium=0.04, tris=0.02).monovalent == (
        pytest.approx(0.06)
    )


def test_ion_conditions_reject_impossible_input():
    with pytest.raises(ThermoError, match="negative"):
        _monovalent(sodium=-1.0)
    with pytest.raises(ThermoError, match="positive"):
        _monovalent(strand=0.0, sodium=0.05)
    with pytest.raises(ThermoError, match="no cation"):
        _monovalent(sodium=0.0)


def test_the_named_presets_describe_their_buffers():
    assert AMPLIFICATION == IonConditions()
    assert AMPLIFICATION.regime == "mixed"
    assert AMPLIFICATION.monovalent == pytest.approx(0.055)
    assert AMPLIFICATION.free_magnesium == pytest.approx(0.7e-3)
    assert SUPPLIER.monovalent == pytest.approx(0.05)
    assert SUPPLIER.free_magnesium == 0.0
    assert SUPPLIER.regime == "monovalent"


def test_an_amplification_buffer_melts_higher_than_the_supplier_convention():
    thermo = duplex("TCATGAAGGATGGTGTTAATGC")
    assert melting_temperature(thermo, AMPLIFICATION) > melting_temperature(
        thermo, SUPPLIER
    )


def test_loop_free_energy_is_tabulated_and_extrapolated():
    assert loop_free_energy(3) == 3.5
    assert loop_free_energy(10) == 4.6
    assert loop_free_energy(30) == 6.3
    assert loop_free_energy(11) > loop_free_energy(10)
    assert loop_free_energy(40) > loop_free_energy(30)
    with pytest.raises(ThermoError, match="cannot close"):
        loop_free_energy(2)


def test_three_prime_stability_uses_only_the_last_bases():
    whole = three_prime_stability("AAAAAAAAAAGGGGG")
    tail = three_prime_stability("GGGGG")
    assert whole == pytest.approx(tail)
    assert three_prime_stability("A") is None


def test_three_prime_stability_is_more_negative_for_a_gc_end():
    assert three_prime_stability("AAAAAGGCGC") < three_prime_stability("GGCGCAAAAA")
