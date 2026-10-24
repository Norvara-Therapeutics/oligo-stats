"""Nearest-neighbour duplex thermodynamics and melting temperature.

Stacking parameters are the unified DNA/DNA set of SantaLucia (1998), PNAS
95:1460, Table 1, reproduced in SantaLucia & Hicks (2004), Annu Rev Biophys
Biomol Struct 33:415. Only enthalpy and entropy are tabulated; free energy is
derived from them at the requested temperature, so a single source of truth
sets every quantity the module reports.

Melting temperature at 1 M Na+ follows the two-state relation

    Tm = dH / (dS + R ln(Ct / x))

with ``x`` 4 for a duplex of two different strands at equal concentration and 1
where the two strands are chemically identical, and is then corrected to the
reaction's own ion composition with the empirical relations of Owczarzy et al.
(2004), Biochemistry 43:3537 for monovalent cations and Owczarzy et al. (2008),
Biochemistry 47:5336 where magnesium is present.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from oligo_stats.sequence import is_self_complementary

#: Gas constant in cal/(mol K), the unit of the tabulated entropies.
R = 1.9872

#: Reference temperature of the tabulated free energies, in kelvin.
T37 = 310.15

ZERO_CELSIUS = 273.15

#: Stacking enthalpy, kcal/mol, keyed by the 5'->3' dinucleotide of either
#: strand of a Watson-Crick duplex.
NN_ENTHALPY: dict[str, float] = {
    "AA": -7.9, "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5, "TG": -8.5,
    "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8,
    "GA": -8.2, "TC": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0, "CC": -8.0,
}

#: Stacking entropy, cal/(mol K), keyed as :data:`NN_ENTHALPY`.
NN_ENTROPY: dict[str, float] = {
    "AA": -22.2, "TT": -22.2,
    "AT": -20.4,
    "TA": -21.3,
    "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0,
    "GA": -22.2, "TC": -22.2,
    "CG": -27.2,
    "GC": -24.4,
    "GG": -19.9, "CC": -19.9,
}

#: Helix initiation, applied once per helix end, by the base pair at that end.
INIT_ENTHALPY: dict[str, float] = {"GC": 0.1, "AT": 2.3}
INIT_ENTROPY: dict[str, float] = {"GC": -2.8, "AT": 4.1}

#: Entropy correction for a duplex with two-fold rotational symmetry.
SYMMETRY_ENTROPY = -1.4

#: Hairpin loop initiation free energy at 37 C, kcal/mol, by loop length.
#: SantaLucia & Hicks (2004), Table 4. Lengths between tabulated entries and
#: above 30 are interpolated by the Jacobson-Stockmayer relation.
_LOOP_FREE_ENERGY: dict[int, float] = {
    3: 3.5,
    4: 3.5,
    5: 3.3,
    6: 4.0,
    7: 4.2,
    8: 4.3,
    9: 4.5,
    10: 4.6,
    12: 5.0,
    14: 5.1,
    16: 5.3,
    18: 5.5,
    20: 5.7,
    25: 6.1,
    30: 6.3,
}

#: Penalty added to a triloop closed by an A-T pair (SantaLucia & Hicks 2004).
TRILOOP_AT_CLOSURE = 0.5

MIN_LOOP = 3


class ThermoError(ValueError):
    """Raised when a duplex cannot be scored with the tabulated parameters."""


@dataclass(frozen=True)
class IonConditions:
    """Reaction ion and strand concentrations.

    Concentrations are molar. ``tris`` is the total concentration of Tris
    buffer, of which half is counted as monovalent cation: at the pH of a
    standard amplification buffer about half the buffer is protonated.
    ``dntp`` is the total concentration of all four deoxynucleotides, which
    chelate magnesium one for one and so reduce the free magnesium available to
    the duplex. ``strand`` is the total concentration of oligonucleotide
    strands.

    The defaults are :data:`AMPLIFICATION`, a standard amplification reaction,
    because that is the buffer a primer anneals in and the one an annealing
    temperature has to be derived from. Setting one ion leaves the rest of that
    buffer in place; :data:`SUPPLIER`, the convention a specification sheet is
    quoted against, is spelled out separately.
    """

    strand: float = 0.25e-6
    sodium: float = 0.0
    potassium: float = 0.05
    tris: float = 0.01
    magnesium: float = 1.5e-3
    dntp: float = 0.8e-3

    def __post_init__(self) -> None:
        for name in ("strand", "sodium", "potassium", "tris", "magnesium", "dntp"):
            if getattr(self, name) < 0:
                raise ThermoError(f"{name} concentration must not be negative")
        if self.strand <= 0:
            raise ThermoError("strand concentration must be positive")
        if self.monovalent <= 0 and self.free_magnesium <= 0:
            raise ThermoError(
                "no cation is present: give a sodium, potassium, Tris or "
                "magnesium concentration"
            )

    @property
    def monovalent(self) -> float:
        """Total monovalent cation concentration, molar."""
        return self.sodium + self.potassium + self.tris / 2.0

    @property
    def free_magnesium(self) -> float:
        """Magnesium left after deoxynucleotide chelation, molar."""
        return max(0.0, self.magnesium - self.dntp)

    @property
    def regime(self) -> str:
        """Which empirical ion correction applies: monovalent, mixed or divalent.

        The ratio of the square root of free magnesium to monovalent cation
        concentration decides which term dominates, following Owczarzy et al.
        (2008): below 0.22 the monovalent relation is used alone, above 6.0 the
        divalent relation is used alone, and between them both contribute.
        """
        if self.free_magnesium <= 0:
            return "monovalent"
        if self.monovalent <= 0:
            return "divalent"
        ratio = math.sqrt(self.free_magnesium) / self.monovalent
        if ratio < 0.22:
            return "monovalent"
        if ratio < 6.0:
            return "mixed"
        return "divalent"


@dataclass(frozen=True)
class Duplex:
    """Thermodynamics of one duplex or helix.

    ``enthalpy`` is kcal/mol and ``entropy`` cal/(mol K), both for duplex
    formation, so both are negative for a stabilising interaction.
    ``base_pairs`` is the number of Watson-Crick pairs and ``gc_fraction`` the
    fraction of them that are G-C, which the ion corrections require.
    """

    enthalpy: float
    entropy: float
    base_pairs: int
    gc_fraction: float
    symmetric: bool = False

    def free_energy(self, temperature: float = T37) -> float:
        """Free energy of formation at ``temperature`` kelvin, kcal/mol."""
        return self.enthalpy - temperature * self.entropy / 1000.0

    @property
    def stable(self) -> bool:
        """True when the duplex is favourable at 37 C."""
        return self.free_energy() < 0


def _pairs(top: str) -> list[str]:
    return [top[index : index + 2] for index in range(len(top) - 1)]


def _initiation(pair: str) -> str:
    return "GC" if pair in ("G", "C") else "AT"


def stack(top: str) -> Duplex:
    """Score a contiguous Watson-Crick helix from its 5'->3' top strand.

    The complementary strand is implied, so only the top strand is needed. Both
    helix ends are charged an initiation term. A helix of fewer than two base
    pairs has no stacking interaction and is returned with zero enthalpy and
    entropy.

    Raises:
        ThermoError: the sequence holds a base outside A, C, G and T.
    """
    unknown = sorted(set(top) - set("ACGT"))
    if unknown:
        listed = ", ".join(repr(base) for base in unknown)
        raise ThermoError(f"no nearest-neighbour parameter for: {listed}")
    if len(top) < 2:
        return Duplex(0.0, 0.0, len(top), 0.0)

    enthalpy = sum(NN_ENTHALPY[step] for step in _pairs(top))
    entropy = sum(NN_ENTROPY[step] for step in _pairs(top))
    for end in (top[0], top[-1]):
        enthalpy += INIT_ENTHALPY[_initiation(end)]
        entropy += INIT_ENTROPY[_initiation(end)]
    gc = sum(1 for base in top if base in "GC")
    return Duplex(
        enthalpy=enthalpy,
        entropy=entropy,
        base_pairs=len(top),
        gc_fraction=gc / len(top),
    )


def duplex(sequence: str, symmetric: bool | None = None) -> Duplex:
    """Score the perfect duplex of an oligonucleotide with its complement.

    Args:
        sequence: a normalised, unambiguous sequence.
        symmetric: whether the duplex has two-fold symmetry. Determined from
            the sequence when omitted.
    """
    helix = stack(sequence)
    is_symmetric = is_self_complementary(sequence) if symmetric is None else symmetric
    if not is_symmetric:
        return Duplex(
            helix.enthalpy, helix.entropy, helix.base_pairs, helix.gc_fraction, False
        )
    return Duplex(
        enthalpy=helix.enthalpy,
        entropy=helix.entropy + SYMMETRY_ENTROPY,
        base_pairs=helix.base_pairs,
        gc_fraction=helix.gc_fraction,
        symmetric=True,
    )


def loop_free_energy(length: int) -> float:
    """Hairpin loop initiation free energy at 37 C, kcal/mol.

    Lengths absent from the table are extrapolated from the nearest shorter
    tabulated loop by the Jacobson-Stockmayer relation
    ``dG(n) = dG(x) + 2.44 R T ln(n / x)``.

    Raises:
        ThermoError: a loop shorter than three nucleotides cannot close.
    """
    if length < MIN_LOOP:
        raise ThermoError(f"a hairpin loop of {length} nt cannot close")
    if length in _LOOP_FREE_ENERGY:
        return _LOOP_FREE_ENERGY[length]
    anchor = max(size for size in _LOOP_FREE_ENERGY if size < length)
    return _LOOP_FREE_ENERGY[anchor] + 2.44 * R * T37 * math.log(length / anchor) / 1000.0


def melting_temperature(
    thermo: Duplex,
    ions: IonConditions,
    unimolecular: bool = False,
    identical_strands: bool = False,
) -> float | None:
    """Melting temperature in degrees Celsius, corrected for ion composition.

    Args:
        thermo: enthalpy and entropy of formation of the structure.
        ions: reaction ion and strand concentrations.
        unimolecular: a hairpin, whose melting point does not depend on
            concentration.
        identical_strands: the two strands of a bimolecular structure are
            chemically identical, as in a self-dimer or a self-complementary
            duplex, which removes the factor of four from the concentration
            term.

    Returns:
        The melting temperature, or ``None`` when the structure does not melt
        above 0 C. A structure that unfolds below the freezing point of the
        buffer does not exist under any usable condition, and both the
        two-state model and the ion corrections are extrapolations there.
    """
    if thermo.entropy == 0.0 or thermo.enthalpy == 0.0:
        return None
    denominator = thermo.entropy / 1000.0
    if not unimolecular:
        divisor = 1.0 if (identical_strands or thermo.symmetric) else 4.0
        denominator += R * math.log(ions.strand / divisor) / 1000.0
    if denominator == 0.0:
        return None
    kelvin = thermo.enthalpy / denominator
    if kelvin <= 0.0:
        return None
    corrected = _correct_for_ions(kelvin, thermo.gc_fraction, thermo.base_pairs, ions)
    if corrected <= ZERO_CELSIUS:
        return None
    return corrected - ZERO_CELSIUS


def _correct_for_ions(
    kelvin: float, gc_fraction: float, base_pairs: int, ions: IonConditions
) -> float:
    """Shift a 1 M Na+ melting temperature to the reaction's ion composition."""
    regime = ions.regime
    if regime == "monovalent":
        return _monovalent_correction(kelvin, gc_fraction, ions.monovalent)
    return _divalent_correction(kelvin, gc_fraction, base_pairs, ions, regime)


def _monovalent_correction(kelvin: float, gc_fraction: float, monovalent: float) -> float:
    """Owczarzy et al. (2004) equation 22, valid from 69 mM to 1.02 M Na+.

    The relation is applied outside that range as well, which is what every
    implementation does; it is the best available estimate for the dilute
    buffers used in amplification.
    """
    if monovalent <= 0:
        return kelvin
    logarithm = math.log(monovalent)
    inverse = (
        1.0 / kelvin
        + (4.29 * gc_fraction - 3.95) * 1e-5 * logarithm
        + 9.40e-6 * logarithm * logarithm
    )
    return 1.0 / inverse if inverse > 0 else 0.0


def _divalent_correction(
    kelvin: float,
    gc_fraction: float,
    base_pairs: int,
    ions: IonConditions,
    regime: str,
) -> float:
    """Owczarzy et al. (2008) equation 16, for magnesium-containing buffers.

    In the mixed regime the three sequence-independent coefficients are
    themselves functions of the monovalent concentration, following the paper's
    equations 18 to 20.
    """
    magnesium = ions.free_magnesium
    if magnesium <= 0:
        return _monovalent_correction(kelvin, gc_fraction, ions.monovalent)
    log_magnesium = math.log(magnesium)
    a, d, g = 3.92e-5, 1.42e-5, 8.31e-5
    if regime == "mixed":
        monovalent = ions.monovalent
        log_monovalent = math.log(monovalent)
        a *= 0.843 - 0.352 * math.sqrt(monovalent) * log_monovalent
        d *= 1.279 - 4.03e-3 * log_monovalent - 8.03e-3 * log_monovalent**2
        g *= 0.486 - 0.258 * log_monovalent + 5.25e-3 * log_monovalent**3
    length_term = 0.0
    if base_pairs > 1:
        length_term = (
            -4.82e-4 + 5.25e-4 * log_magnesium + g * log_magnesium * log_magnesium
        ) / (2.0 * (base_pairs - 1))
    inverse = (
        1.0 / kelvin
        + a
        - 9.11e-6 * log_magnesium
        + gc_fraction * (6.26e-5 + d * log_magnesium)
        + length_term
    )
    return 1.0 / inverse if inverse > 0 else 0.0


#: The default: 0.25 uM primer in 50 mM KCl, 10 mM Tris-HCl, 1.5 mM MgCl2 and
#: 0.2 mM of each deoxynucleotide, which is a standard amplification reaction.
AMPLIFICATION = IonConditions()

#: The convention a supplier's specification sheet is quoted against: 0.25 uM
#: oligonucleotide in 50 mM Na+ with no magnesium. Melting temperatures
#: computed here are several degrees below those of the same oligonucleotide in
#: an amplification buffer, and the two must not be compared.
SUPPLIER = IonConditions(
    strand=0.25e-6, sodium=0.05, potassium=0.0, tris=0.0, magnesium=0.0, dntp=0.0
)


def three_prime_stability(sequence: str, window: int = 5) -> float | None:
    """Stacking free energy of the last ``window`` bases at 37 C, kcal/mol.

    Initiation is excluded, so the value measures only how firmly the 3' end of
    a primer stacks on a template. A strongly stabilised 3' end tolerates
    mismatches further in and so primes off-target sites. Returns ``None`` for a
    sequence with fewer than two bases.
    """
    tail = sequence[-window:]
    if len(tail) < 2:
        return None
    return sum(
        NN_ENTHALPY[step] - T37 * NN_ENTROPY[step] / 1000.0 for step in _pairs(tail)
    )
