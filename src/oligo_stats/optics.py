"""Molecular weight, extinction coefficient and the conversions they support.

Residue masses are those of the deoxynucleoside monophosphates as free acids;
the terminal adjustment removes one phosphate group, giving an oligonucleotide
with a 5'-hydroxyl, which is what standard synthesis delivers. A 5'-phosphate
can be requested instead.

The extinction coefficient at 260 nm uses the nearest-neighbour model

    e = sum(e of each dinucleotide) - sum(e of each internal nucleotide)

with the dinucleotide and mononucleotide values of Warshaw & Cantor (1970),
Biopolymers 9:1079, as tabulated by Cavaluzzi & Borer (2004), Nucleic Acids Res
32:e13. Hypochromicity makes a duplex absorb less than the sum of its strands,
so these single-strand values must not be used to quantify annealed material.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Residue mass of each deoxynucleoside monophosphate in a chain, g/mol.
RESIDUE_MASS: dict[str, float] = {
    "A": 313.21,
    "C": 289.18,
    "G": 329.21,
    "T": 304.20,
}

#: Removal of the 5'-terminal phosphate, leaving a 5'-hydroxyl, g/mol.
HYDROXYL_ADJUSTMENT = -61.96

#: Mass of a 5'-monophosphate relative to a 5'-hydroxyl terminus, g/mol.
PHOSPHATE_MASS = 79.98

#: Molar extinction coefficient at 260 nm of each mononucleotide, 1/(M cm).
MONONUCLEOTIDE_EXTINCTION: dict[str, float] = {
    "A": 15400.0,
    "C": 7400.0,
    "G": 11500.0,
    "T": 8700.0,
}

#: Molar extinction coefficient at 260 nm of each dinucleotide, 1/(M cm).
DINUCLEOTIDE_EXTINCTION: dict[str, float] = {
    "AA": 27400.0, "AC": 21200.0, "AG": 25000.0, "AT": 22800.0,
    "CA": 21200.0, "CC": 14600.0, "CG": 18000.0, "CT": 15200.0,
    "GA": 25200.0, "GC": 17600.0, "GG": 21600.0, "GT": 20000.0,
    "TA": 23400.0, "TC": 16200.0, "TG": 19000.0, "TT": 16800.0,
}

TERMINI = ("hydroxyl", "phosphate")


class OpticsError(ValueError):
    """Raised when mass or absorbance cannot be computed for a sequence."""


@dataclass(frozen=True)
class Optics:
    """Mass and absorbance properties of one single-stranded oligonucleotide.

    ``extinction`` is the molar extinction coefficient at 260 nm in 1/(M cm).
    ``micrograms_per_od`` is the mass in one millilitre of a solution with an
    absorbance of 1.000 at 260 nm in a 1 cm path, and ``nanomoles_per_od`` the
    corresponding amount of substance.
    """

    molecular_weight: float
    extinction: float
    terminus: str

    @property
    def micrograms_per_od(self) -> float:
        """Micrograms per millilitre at an absorbance of 1.000 at 260 nm."""
        return 1000.0 * self.molecular_weight / self.extinction

    @property
    def nanomoles_per_od(self) -> float:
        """Nanomoles per millilitre at an absorbance of 1.000 at 260 nm."""
        return 1e6 / self.extinction


def molecular_weight(sequence: str, terminus: str = "hydroxyl") -> float:
    """Anhydrous molecular weight of a single strand, g/mol.

    Args:
        sequence: a normalised, unambiguous sequence.
        terminus: ``hydroxyl`` for the 5'-OH of standard synthesis, or
            ``phosphate`` for a 5'-phosphorylated oligonucleotide.

    Raises:
        OpticsError: the sequence holds a base with no tabulated mass, or the
            terminus is not recognised.
    """
    if terminus not in TERMINI:
        raise OpticsError(f"unknown 5' terminus: {terminus!r}")
    total = 0.0
    for base in sequence:
        if base not in RESIDUE_MASS:
            raise OpticsError(f"no residue mass for base {base!r}")
        total += RESIDUE_MASS[base]
    total += HYDROXYL_ADJUSTMENT
    if terminus == "phosphate":
        total += PHOSPHATE_MASS
    return round(total, 2)


def extinction_coefficient(sequence: str) -> float:
    """Molar extinction coefficient of a single strand at 260 nm, 1/(M cm).

    Raises:
        OpticsError: the sequence is empty or holds an untabulated base.
    """
    if not sequence:
        raise OpticsError("sequence is empty")
    for base in sequence:
        if base not in MONONUCLEOTIDE_EXTINCTION:
            raise OpticsError(f"no extinction coefficient for base {base!r}")
    if len(sequence) == 1:
        return MONONUCLEOTIDE_EXTINCTION[sequence]
    nearest_neighbour = sum(
        DINUCLEOTIDE_EXTINCTION[sequence[index : index + 2]]
        for index in range(len(sequence) - 1)
    )
    internal = sum(MONONUCLEOTIDE_EXTINCTION[base] for base in sequence[1:-1])
    return round(nearest_neighbour - internal, 1)


def optics(sequence: str, terminus: str = "hydroxyl") -> Optics:
    """Mass and absorbance properties of one oligonucleotide."""
    return Optics(
        molecular_weight=molecular_weight(sequence, terminus),
        extinction=extinction_coefficient(sequence),
        terminus=terminus,
    )
