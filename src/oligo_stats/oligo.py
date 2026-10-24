"""Assembly of the per-oligonucleotide measurement set and its quality flags.

The thresholds in :class:`Criteria` are design conventions for a standard
amplification primer, not physical constants. They are the defaults a reviewer
would apply by eye, collected in one place so that a run states which rules it
applied. An oligonucleotide designed for another purpose, such as a sequencing
primer, a probe or a cloning adapter, will fail some of them by intent.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Sequence

from oligo_stats.optics import Optics, optics
from oligo_stats.sequence import (
    Composition,
    composition,
    has_uracil,
    is_self_complementary,
    longest_runs,
    validate,
)
from oligo_stats.structure import (
    THREE_PRIME_WINDOW,
    Dimer,
    Hairpin,
    find_dimers,
    find_hairpins,
    find_self_dimers,
)
from oligo_stats.thermo import (
    Duplex,
    IonConditions,
    duplex,
    melting_temperature,
    three_prime_stability,
)

#: How many 3'-terminal bases the GC clamp is counted over.
CLAMP_WINDOW = 5


#: Stem lengths whose spoken form begins with a vowel. A stem cannot reach the
#: other such numbers, since it is bounded by half the length of the oligo.
_VOWEL_INITIAL = frozenset({8, 11, 18})


def _article(count: int) -> str:
    """The indefinite article that precedes a number read aloud."""
    return "an" if count in _VOWEL_INITIAL else "a"


@dataclass(frozen=True)
class Criteria:
    """Design thresholds applied when flagging an oligonucleotide.

    Free-energy thresholds are upper bounds: a structure is flagged when its
    free energy at 37 C is at or below the threshold. ``dimer_free_energy``
    applies to a dimer that leaves both 3' ends unpaired, and
    ``dimer_three_prime_free_energy``, which is deliberately less permissive, to
    one that pairs a 3' end and can therefore be extended.
    """

    length_min: int = 18
    length_max: int = 30
    gc_min: float = 40.0
    gc_max: float = 60.0
    melting_min: float = 55.0
    melting_max: float = 65.0
    clamp_min: int = 1
    clamp_max: int = 3
    end_stability: float = -9.0
    hairpin_free_energy: float = -2.0
    dimer_free_energy: float = -9.0
    dimer_three_prime_free_energy: float = -5.0
    max_run: int = 4


@dataclass(frozen=True)
class Flag:
    """One design observation about an oligonucleotide.

    ``severity`` is ``warn`` where the observation would normally stop a primer
    going to synthesis and ``note`` where it is worth seeing but not
    disqualifying.
    """

    code: str
    severity: str
    message: str


@dataclass(frozen=True)
class OligoStats:
    """Every measurement this package makes for one oligonucleotide."""

    name: str
    description: str
    sequence: str
    composition: Composition
    thermo: Duplex
    melting: float | None
    optics: Optics
    end_stability: float | None
    gc_clamp: int
    self_complementary: bool
    rna_input: bool
    structures_screened: bool = True
    hairpins: Sequence[Hairpin] = ()
    self_dimers: Sequence[Dimer] = ()
    flags: Sequence[Flag] = ()

    @property
    def length(self) -> int:
        """Length in nucleotides."""
        return len(self.sequence)

    @property
    def free_energy(self) -> float:
        """Free energy of the intended duplex at 37 C, kcal/mol."""
        return self.thermo.free_energy()

    @property
    def best_hairpin(self) -> Hairpin | None:
        """The most stable hairpin found, if any stem-loop can form."""
        return self.hairpins[0] if self.hairpins else None

    @property
    def best_self_dimer(self) -> Dimer | None:
        """The most stable self-dimer register found, if any helix forms."""
        return self.self_dimers[0] if self.self_dimers else None

    @property
    def worst_severity(self) -> str:
        """``warn``, ``note`` or ``ok``, whichever the flags justify."""
        severities = {flag.severity for flag in self.flags}
        if "warn" in severities:
            return "warn"
        return "note" if "note" in severities else "ok"


@dataclass(frozen=True)
class CrossDimer:
    """The most stable register between two different oligonucleotides."""

    first: str
    second: str
    dimer: Dimer


@dataclass(frozen=True)
class Analysis:
    """A whole run: the oligonucleotides, their cross-dimers and the settings."""

    oligos: Sequence[OligoStats]
    ions: IonConditions
    criteria: Criteria
    terminus: str
    cross_dimers: Sequence[CrossDimer] = field(default_factory=tuple)
    cross_dimers_screened: bool = False

    @property
    def melting_span(self) -> tuple[float, float] | None:
        """Lowest and highest duplex melting temperature across the run."""
        values = [item.melting for item in self.oligos if item.melting is not None]
        if not values:
            return None
        return min(values), max(values)


def gc_clamp(sequence: str, window: int = CLAMP_WINDOW) -> int:
    """Count G and C among the last ``window`` bases of a sequence."""
    return sum(1 for base in sequence[-window:] if base in "GC")


def _flag_composition(stats: OligoStats, criteria: Criteria) -> list[Flag]:
    flags: list[Flag] = []
    if stats.length < criteria.length_min or stats.length > criteria.length_max:
        flags.append(
            Flag(
                "length",
                "note",
                f"{stats.length} nt is outside the usual "
                f"{criteria.length_min}-{criteria.length_max} nt range",
            )
        )
    gc = stats.composition.gc_percent
    if gc is not None and (gc < criteria.gc_min or gc > criteria.gc_max):
        flags.append(
            Flag(
                "gc_content",
                "warn",
                f"GC content {gc:.1f}% is outside "
                f"{criteria.gc_min:.0f}-{criteria.gc_max:.0f}%",
            )
        )
    runs = longest_runs(stats.sequence)
    for base in ("A", "C", "G", "T"):
        run = runs.get(base, 0)
        if run < criteria.max_run:
            continue
        flags.append(
            Flag(
                "base_run",
                "warn" if base == "G" else "note",
                f"a run of {run} consecutive {base} bases"
                + (
                    "; four or more consecutive G can form a quadruplex and "
                    "couples poorly during synthesis"
                    if base == "G"
                    else ""
                ),
            )
        )
    if stats.rna_input:
        flags.append(
            Flag(
                "rna_input",
                "note",
                "input contained uracil; it was read as thymine and DNA "
                "nearest-neighbour parameters were applied",
            )
        )
    if stats.self_complementary:
        flags.append(
            Flag(
                "self_complementary",
                "warn",
                "the sequence is its own reverse complement, so it anneals to "
                "itself in preference to any template",
            )
        )
    return flags


def _flag_thermodynamics(stats: OligoStats, criteria: Criteria) -> list[Flag]:
    flags: list[Flag] = []
    if stats.melting is None:
        # Silence here would let an oligonucleotide with no computable melting
        # point pass every melting-range check and leave a gate reporting that
        # it had nothing to say.
        flags.append(
            Flag(
                "melting_temperature",
                "warn",
                "no melting temperature above 0 C under these conditions, so "
                "the melting range could not be applied; check the strand and "
                "ion concentrations",
            )
        )
    elif (
        stats.melting < criteria.melting_min or stats.melting > criteria.melting_max
    ):
        flags.append(
            Flag(
                "melting_temperature",
                "warn",
                f"Tm {stats.melting:.1f} C is outside "
                f"{criteria.melting_min:.0f}-{criteria.melting_max:.0f} C",
            )
        )
    if stats.gc_clamp < criteria.clamp_min:
        flags.append(
            Flag(
                "gc_clamp",
                "warn",
                f"no G or C in the last {CLAMP_WINDOW} bases, so the 3' end "
                "does not clamp to the template",
            )
        )
    elif stats.gc_clamp > criteria.clamp_max:
        flags.append(
            Flag(
                "gc_clamp",
                "note",
                f"{stats.gc_clamp} of the last {CLAMP_WINDOW} bases are G or C, "
                "which stabilises mismatched 3' ends as well as matched ones",
            )
        )
    if stats.end_stability is not None and stats.end_stability <= criteria.end_stability:
        flags.append(
            Flag(
                "end_stability",
                "warn",
                f"the last {CLAMP_WINDOW} bases stack at "
                f"{stats.end_stability:.1f} kcal/mol, which can prime a "
                "mismatched site",
            )
        )
    return flags


def dimer_threshold(dimer: Dimer, criteria: Criteria) -> float:
    """The free energy at which this register is flagged.

    A register that pairs a 3' end can be extended, so it is held to the less
    permissive of the two thresholds.
    """
    return (
        criteria.dimer_three_prime_free_energy
        if dimer.three_prime_paired
        else criteria.dimer_free_energy
    )


def worst_dimer(dimers: Sequence[Dimer], criteria: Criteria) -> Dimer | None:
    """The most stable register that breaches the threshold for its own class.

    Judging only the most stable register misses a qualifying one: a weaker
    register that pairs a 3' end is held to a less permissive threshold, so it
    can breach while the strongest register of all passes.
    """
    offenders = [
        dimer for dimer in dimers if dimer.free_energy <= dimer_threshold(dimer, criteria)
    ]
    if not offenders:
        return None
    return min(offenders, key=lambda dimer: dimer.free_energy)


def _describe_three_prime(dimer: Dimer) -> str:
    """How the register sits against the 3' ends, in the terms it earns."""
    if dimer.terminal_paired:
        return (
            "pairing a 3'-terminal base, so a polymerase can extend it into a "
            "primer dimer"
        )
    if dimer.three_prime_paired:
        return (
            f"pairing within {THREE_PRIME_WINDOW} bases of a 3' end, which "
            "leaves a short unpaired 3' overhang; extension needs the terminal "
            "bases to pair as well"
        )
    return "leaving both 3' ends unpaired"


def _flag_structures(stats: OligoStats, criteria: Criteria) -> list[Flag]:
    flags: list[Flag] = []
    hairpin = stats.best_hairpin
    if hairpin is not None and hairpin.free_energy <= criteria.hairpin_free_energy:
        melting = (
            f", melting at {hairpin.melting:.1f} C" if hairpin.melting is not None else ""
        )
        flags.append(
            Flag(
                "hairpin",
                "warn",
                f"hairpin of {hairpin.free_energy:.2f} kcal/mol with "
                f"{_article(hairpin.stem_length)} {hairpin.stem_length} bp "
                f"stem{melting}",
            )
        )
    dimer = worst_dimer(stats.self_dimers, criteria)
    if dimer is not None:
        flags.append(
            Flag(
                "self_dimer",
                "warn",
                f"self-dimer of {dimer.free_energy:.2f} kcal/mol over "
                f"{dimer.pairs} bp, {_describe_three_prime(dimer)}",
            )
        )
    return flags


def analyse_oligo(
    sequence: str,
    name: str = "oligo",
    description: str = "",
    ions: IonConditions | None = None,
    criteria: Criteria | None = None,
    terminus: str = "hydroxyl",
    structures: int = 3,
    screen_structures: bool = True,
) -> OligoStats:
    """Measure one oligonucleotide and flag it against the design criteria.

    Args:
        sequence: the sequence as supplied; whitespace and case are tolerated
            and uracil is read as thymine.
        name: an identifier for the oligonucleotide.
        description: free text carried through to the outputs, such as the
            header of the FASTA record the sequence came from.
        ions: reaction conditions. Defaults to
            :data:`~oligo_stats.thermo.AMPLIFICATION`, a standard amplification
            reaction: 0.25 uM strand in 50 mM KCl, 10 mM Tris-HCl and 1.5 mM
            MgCl2 with 0.2 mM of each deoxynucleotide.
        criteria: design thresholds. Defaults to :class:`Criteria`.
        terminus: ``hydroxyl`` or ``phosphate``, the 5' end as synthesised.
        structures: how many hairpins and self-dimer registers to retain.
        screen_structures: run the hairpin and self-dimer searches. Turning the
            screen off removes the only part of the analysis whose cost grows
            faster than the length of the sequence.

    Raises:
        SequenceError: the sequence is empty, or holds a character that is not
            an unambiguous nucleotide base.
    """
    conditions = IonConditions() if ions is None else ions
    rules = Criteria() if criteria is None else criteria
    rna = has_uracil(sequence)
    normalised = validate(sequence)

    thermo = duplex(normalised)
    stats = OligoStats(
        name=name,
        description=description,
        sequence=normalised,
        composition=composition(normalised),
        thermo=thermo,
        melting=melting_temperature(thermo, conditions),
        optics=optics(normalised, terminus),
        end_stability=three_prime_stability(normalised, CLAMP_WINDOW),
        gc_clamp=gc_clamp(normalised),
        self_complementary=is_self_complementary(normalised),
        rna_input=rna,
        structures_screened=screen_structures,
        hairpins=tuple(find_hairpins(normalised, conditions, limit=structures))
        if screen_structures
        else (),
        self_dimers=tuple(find_self_dimers(normalised, conditions, limit=structures))
        if screen_structures
        else (),
    )
    flags = (
        _flag_composition(stats, rules)
        + _flag_thermodynamics(stats, rules)
        + (_flag_structures(stats, rules) if screen_structures else [])
    )
    return replace(stats, flags=tuple(flags))


def screen_cross_dimers(
    oligos: Sequence[OligoStats],
    ions: IonConditions,
    criteria: Criteria | None = None,
) -> list[CrossDimer]:
    """Score every pair of oligonucleotides and keep the registers that matter.

    A pair is reported when any of its registers is at or below the dimer
    threshold for that register's 3'-end involvement, which is the same rule
    applied to self-dimers, and the register reported is the most stable of
    those that breach. Pairs are returned most stable first.
    """
    rules = Criteria() if criteria is None else criteria
    found: list[CrossDimer] = []
    for first_index in range(len(oligos)):
        for second_index in range(first_index + 1, len(oligos)):
            first = oligos[first_index]
            second = oligos[second_index]
            registers = find_dimers(
                first.sequence,
                second.sequence,
                ions,
                limit=1,
                # Two entries of a panel may carry the same sequence under
                # different names, and then the two strands really are the
                # same molecule.
                identical_strands=first.sequence == second.sequence,
            )
            if not registers:
                continue
            dimer = worst_dimer(registers, rules)
            if dimer is not None:
                found.append(CrossDimer(first.name, second.name, dimer))
    found.sort(key=lambda item: item.dimer.free_energy)
    return found
