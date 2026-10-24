"""Oligonucleotide normalisation, composition and degeneracy handling."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

#: Complement map covering the standard bases and the IUPAC ambiguity codes.
_COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)

#: IUPAC ambiguity codes expanded to the DNA bases they stand for.
IUPAC_CODES: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

BASES = ("A", "C", "G", "T")
_WHITESPACE = frozenset(" \t\r\n")


class SequenceError(ValueError):
    """Raised for input that is not usable as an oligonucleotide."""


def reverse_complement(sequence: str) -> str:
    """Reverse-complement a sequence, preserving case and ambiguity codes."""
    return sequence.translate(_COMPLEMENT)[::-1]


def normalise(sequence: str) -> str:
    """Upper-case a sequence, drop whitespace, and read U as T.

    Internal whitespace is common in sequences pasted from a vendor
    specification or a publication, and carries no meaning here. Uracil is
    accepted so an RNA oligonucleotide can be measured, but every
    thermodynamic parameter in this package is for a DNA duplex.
    """
    return "".join(
        character for character in sequence.upper() if character not in _WHITESPACE
    ).replace("U", "T")


def has_uracil(sequence: str) -> bool:
    """True when the sequence as given contains U, before normalisation."""
    return "U" in sequence.upper()


def degenerate_positions(sequence: str) -> list[int]:
    """Zero-based positions holding an ambiguity code, in a normalised sequence."""
    return [index for index, base in enumerate(sequence) if base not in BASES]


def validate(sequence: str, allow_degenerate: bool = False) -> str:
    """Normalise a sequence and reject anything that is not nucleic acid.

    Args:
        sequence: the sequence as supplied.
        allow_degenerate: accept IUPAC ambiguity codes as well as A, C, G, T.

    Returns:
        The normalised sequence.

    Raises:
        SequenceError: the sequence is empty, or holds a character that is not
            a recognised base, or holds an ambiguity code that is not allowed.
    """
    normalised = normalise(sequence)
    if not normalised:
        raise SequenceError("sequence is empty")
    unknown = sorted(set(normalised) - set(IUPAC_CODES))
    if unknown:
        listed = ", ".join(repr(character) for character in unknown)
        raise SequenceError(f"not a nucleotide base: {listed}")
    degenerate = degenerate_positions(normalised)
    if degenerate and not allow_degenerate:
        codes = sorted({normalised[index] for index in degenerate})
        raise SequenceError(
            f"sequence contains the ambiguity code(s) {', '.join(codes)}; "
            "nearest-neighbour parameters are defined for A, C, G and T only. "
            "Resolve them first with 'oligo-stats expand'"
        )
    return normalised


def expansion_count(sequence: str) -> int:
    """Number of unambiguous sequences a degenerate sequence stands for."""
    total = 1
    for base in sequence:
        total *= len(IUPAC_CODES[base])
    return total


def expand(sequence: str, limit: int = 256) -> list[str]:
    """Enumerate every unambiguous sequence a degenerate sequence stands for.

    Args:
        sequence: a normalised sequence, which may hold ambiguity codes.
        limit: the largest number of sequences to produce.

    Returns:
        The expansion in lexicographic order of the substituted bases, which is
        the order of the input positions read left to right.

    Raises:
        SequenceError: the expansion would exceed ``limit``.
    """
    total = expansion_count(sequence)
    if total > limit:
        raise SequenceError(
            f"sequence expands to {total:,} sequences, above the limit of {limit:,}"
        )
    choices = [IUPAC_CODES[base] for base in sequence]
    return ["".join(combination) for combination in product(*choices)]


@dataclass(frozen=True)
class Composition:
    """Base composition of one oligonucleotide.

    ``gc_percent`` is computed over unambiguous bases only, unrounded, and is
    ``None`` when the sequence holds none. It is unrounded because it is
    compared against a threshold the caller sets; the output formatters round
    it for display. ``degenerate`` counts positions holding an
    ambiguity code, which are excluded from the four base counts.
    """

    length: int
    counts: dict[str, int] = field(default_factory=dict)
    degenerate: int = 0
    gc_percent: float | None = None

    @property
    def gc_count(self) -> int:
        """Number of unambiguous G and C bases."""
        return self.counts["G"] + self.counts["C"]


def composition(sequence: str) -> Composition:
    """Count the bases of a normalised sequence."""
    counts = {base: sequence.count(base) for base in BASES}
    unambiguous = sum(counts.values())
    # Stored unrounded: this value is compared against a configurable threshold,
    # and a value rounded for display passes limits the sequence does not meet.
    gc = 100.0 * (counts["G"] + counts["C"]) / unambiguous if unambiguous else None
    return Composition(
        length=len(sequence),
        counts=counts,
        degenerate=len(sequence) - unambiguous,
        gc_percent=gc,
    )


def is_self_complementary(sequence: str) -> bool:
    """True when a sequence is its own reverse complement.

    Such an oligonucleotide forms a symmetric duplex with itself, which changes
    the concentration term and the entropy of duplex formation.
    """
    return len(sequence) % 2 == 0 and sequence == reverse_complement(sequence)


def longest_runs(sequence: str) -> dict[str, int]:
    """The longest run of each base present, keyed by base.

    A run of guanine matters wherever it falls, not only where it happens to be
    the longest run in the sequence, so every base is reported separately.
    """
    longest: dict[str, int] = {}
    current_base, current_length = "", 0
    for base in sequence:
        current_length = current_length + 1 if base == current_base else 1
        current_base = base
        if current_length > longest.get(base, 0):
            longest[base] = current_length
    return longest
