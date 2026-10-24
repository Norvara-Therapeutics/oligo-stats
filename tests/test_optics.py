"""Mass, extinction coefficient and the conversions that follow."""

import pytest

from oligo_stats.optics import (
    DINUCLEOTIDE_EXTINCTION,
    HYDROXYL_ADJUSTMENT,
    MONONUCLEOTIDE_EXTINCTION,
    PHOSPHATE_MASS,
    RESIDUE_MASS,
    OpticsError,
    extinction_coefficient,
    molecular_weight,
    optics,
)


def test_molecular_weight_matches_the_residue_sum():
    sequence = "ACGT"
    expected = sum(RESIDUE_MASS[base] for base in sequence) + HYDROXYL_ADJUSTMENT
    assert molecular_weight(sequence) == pytest.approx(round(expected, 2))


def test_molecular_weight_of_a_known_primer():
    """The arithmetic a supplier's specification sheet reports, spelled out."""
    sequence = "TCATGAAGGATGGTGTTAATGC"
    counts = {base: sequence.count(base) for base in "ACGT"}
    expected = round(
        counts["A"] * 313.21
        + counts["C"] * 289.18
        + counts["G"] * 329.21
        + counts["T"] * 304.20
        - 61.96,
        2,
    )
    assert molecular_weight(sequence) == expected


def test_a_five_prime_phosphate_adds_its_mass():
    assert molecular_weight("ACGT", "phosphate") == pytest.approx(
        molecular_weight("ACGT") + PHOSPHATE_MASS
    )


def test_molecular_weight_rejects_unknown_input():
    with pytest.raises(OpticsError, match="residue mass"):
        molecular_weight("ACGN")
    with pytest.raises(OpticsError, match="terminus"):
        molecular_weight("ACGT", "amine")


def test_extinction_of_a_single_base_is_the_mononucleotide_value():
    for base, value in MONONUCLEOTIDE_EXTINCTION.items():
        assert extinction_coefficient(base) == value


@pytest.mark.parametrize("step", sorted(DINUCLEOTIDE_EXTINCTION))
def test_extinction_of_a_dinucleotide_is_the_tabulated_value(step):
    """With no internal base to subtract, the model reduces to the table."""
    assert extinction_coefficient(step) == DINUCLEOTIDE_EXTINCTION[step]


def test_extinction_follows_the_nearest_neighbour_model():
    sequence = "ACG"
    expected = (
        DINUCLEOTIDE_EXTINCTION["AC"]
        + DINUCLEOTIDE_EXTINCTION["CG"]
        - MONONUCLEOTIDE_EXTINCTION["C"]
    )
    assert extinction_coefficient(sequence) == pytest.approx(expected)


def test_extinction_grows_with_length_and_is_lower_for_gc():
    assert extinction_coefficient("AAAAAAAAAA") > extinction_coefficient("AAAAA")
    assert extinction_coefficient("GCGCGCGCGC") < extinction_coefficient("AAAAAAAAAA")


def test_extinction_rejects_empty_and_unknown_input():
    with pytest.raises(OpticsError, match="empty"):
        extinction_coefficient("")
    with pytest.raises(OpticsError, match="extinction coefficient"):
        extinction_coefficient("ACGN")


def test_optical_density_conversions_invert_the_extinction_coefficient():
    result = optics("TCATGAAGGATGGTGTTAATGC")
    assert result.nanomoles_per_od == pytest.approx(1e6 / result.extinction)
    assert result.micrograms_per_od == pytest.approx(
        result.nanomoles_per_od * result.molecular_weight / 1000.0
    )
    assert result.terminus == "hydroxyl"
