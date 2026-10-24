"""Hairpin and dimer screening.

Both searches are exhaustive over the structures they model and both score with
the nearest-neighbour parameters of :mod:`oligo_stats.thermo`, so a reported
free energy is comparable with the free energy of the intended duplex. The
models are deliberately narrow:

* A hairpin is one contiguous Watson-Crick stem closing one loop of at least
  three nucleotides. Bulges, internal loops and multibranch structures are not
  considered, and no terminal-mismatch stabilisation is credited to the loop.
* A dimer is one antiparallel register of two strands. Every contiguous
  Watson-Crick run of at least two base pairs in that register is scored as its
  own helix, paying initiation at both of its ends, and a run is counted only
  where it is favourable at 37 C. Bulges are not modelled, so a register whose
  complementarity is interrupted is scored as separate short helices rather
  than one long one.

What follows from this is that the free energies are screening estimates for
ranking candidates, not a substitute for a folding algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from oligo_stats.sequence import reverse_complement
from oligo_stats.thermo import (
    INIT_ENTHALPY,
    INIT_ENTROPY,
    MIN_LOOP,
    NN_ENTHALPY,
    NN_ENTROPY,
    T37,
    TRILOOP_AT_CLOSURE,
    Duplex,
    IonConditions,
    loop_free_energy,
    melting_temperature,
    stack,
)

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

#: Smallest stem and smallest helix the searches will report.
MIN_STEM = 3
MIN_RUN = 2

#: How many 3'-terminal bases count as "the 3' end" when judging a dimer.
THREE_PRIME_WINDOW = 3


def _complementary(first: str, second: str) -> bool:
    return _COMPLEMENT.get(first) == second


@dataclass(frozen=True)
class Hairpin:
    """One stem-loop structure of a single oligonucleotide.

    Positions are zero-based into the sequence. The stem occupies
    ``[stem_start, stem_start + stem_length)`` on the 5' side and pairs with the
    equally long block ending at ``stem_end`` on the 3' side.
    """

    stem_start: int
    stem_length: int
    loop_start: int
    loop_length: int
    stem_end: int
    thermo: Duplex
    melting: float | None
    diagram: str

    @property
    def free_energy(self) -> float:
        """Free energy of folding at 37 C, kcal/mol."""
        return self.thermo.free_energy()

    @property
    def loop(self) -> tuple[int, int]:
        """Loop span as a zero-based half-open interval."""
        return self.loop_start, self.loop_start + self.loop_length


@dataclass(frozen=True)
class Dimer:
    """One antiparallel register of two strands.

    ``shift`` positions the second strand against the first: the base at index
    ``index`` of the first strand faces index ``len(second) - 1 - index + shift``
    of the second. ``pairs`` is the number of base pairs scored.

    Two flags describe the 3' ends, and they are not the same thing.
    ``three_prime_paired`` is true when a scored pair falls within
    :data:`THREE_PRIME_WINDOW` bases of either 3' terminus; that is the
    screening criterion, and it is deliberately the wider of the two.
    ``terminal_paired`` is true only when a 3'-terminal base is itself paired,
    which is the register a polymerase can extend without the duplex shifting.
    A register that is ``three_prime_paired`` but not ``terminal_paired`` leaves
    a short unpaired 3' overhang and is the weaker risk of the two.
    """

    shift: int
    pairs: int
    helices: int
    thermo: Duplex
    melting: float | None
    diagram: str
    three_prime_paired: bool
    terminal_paired: bool = False

    @property
    def free_energy(self) -> float:
        """Free energy of association at 37 C, kcal/mol."""
        return self.thermo.free_energy()


def _loop_entropy(loop_length: int, closing_base: str) -> float:
    """Entropy equivalent of the loop initiation penalty, cal/(mol K).

    Loop initiation is tabulated as a free energy at 37 C and treated as purely
    entropic, which is the usual convention: a loop costs conformational
    freedom rather than base-pairing enthalpy.
    """
    penalty = loop_free_energy(loop_length)
    if loop_length == MIN_LOOP and closing_base in ("A", "T"):
        penalty += TRILOOP_AT_CLOSURE
    return -1000.0 * penalty / T37


def _hairpin_thermo(sequence: str, stem_start: int, stem_length: int, loop_length: int) -> Duplex:
    """Enthalpy and entropy of one stem-loop.

    The stem contributes its stacks. The exterior end of the stem is a free
    helix end and is charged initiation; the loop end is closed by the loop
    penalty instead.
    """
    stem = sequence[stem_start : stem_start + stem_length]
    enthalpy = sum(NN_ENTHALPY[stem[i : i + 2]] for i in range(len(stem) - 1))
    entropy = sum(NN_ENTROPY[stem[i : i + 2]] for i in range(len(stem) - 1))
    exterior = "GC" if stem[0] in "GC" else "AT"
    enthalpy += INIT_ENTHALPY[exterior]
    entropy += INIT_ENTROPY[exterior]
    entropy += _loop_entropy(loop_length, stem[-1])
    gc = sum(1 for base in stem if base in "GC")
    return Duplex(
        enthalpy=enthalpy,
        entropy=entropy,
        base_pairs=stem_length,
        gc_fraction=gc / stem_length,
    )


def _hairpin_diagram(sequence: str, stem_start: int, stem_length: int, loop_length: int) -> str:
    """Render a stem-loop as three monospaced rows."""
    stem_five = sequence[stem_start : stem_start + stem_length]
    loop_start = stem_start + stem_length
    stem_three_start = loop_start + loop_length
    stem_three = sequence[stem_three_start : stem_three_start + stem_length]
    loop = sequence[loop_start:stem_three_start]
    tail_five = sequence[:stem_start]
    tail_three = sequence[stem_three_start + stem_length :]
    width = max(len(tail_five), len(tail_three))
    bars = "|" * stem_length
    return "\n".join(
        [
            f"5'-{tail_five:>{width}} {stem_five}-\\",
            f"   {'':>{width}} {bars}  {loop}",
            f"3'-{tail_three[::-1]:>{width}} {stem_three[::-1]}-/",
        ]
    )


def find_hairpins(
    sequence: str,
    ions: IonConditions,
    min_stem: int = MIN_STEM,
    limit: int = 3,
) -> list[Hairpin]:
    """Find the most stable stem-loop structures of one oligonucleotide.

    Every stem-loop with a stem of at least ``min_stem`` contiguous base pairs
    and a loop of at least three nucleotides is enumerated. Structures sharing a
    loop are the same hairpin with the stem trimmed, so only the most stable of
    each loop is kept.

    Args:
        sequence: a normalised, unambiguous sequence.
        ions: reaction conditions, used for the melting temperature.
        min_stem: shortest stem to consider.
        limit: how many structures to return, most stable first.

    Returns:
        Up to ``limit`` hairpins sorted by free energy, most stable first. The
        list is empty when no stem-loop can form at all; it is not filtered by
        stability, so a caller can report the best structure even where that
        structure is unfavourable.
    """
    if min_stem < 1:
        raise ValueError("min_stem must be at least 1")
    best: dict[tuple[int, int], Hairpin] = {}
    length = len(sequence)
    for start in range(length):
        for end in range(start + 2 * min_stem + MIN_LOOP - 1, length):
            if not _complementary(sequence[start], sequence[end]):
                continue
            stem_length = 1
            while True:
                loop_length = end - start - 2 * stem_length + 1
                if loop_length < MIN_LOOP:
                    break
                if stem_length >= min_stem:
                    thermo = _hairpin_thermo(sequence, start, stem_length, loop_length)
                    loop_start = start + stem_length
                    key = (loop_start, loop_length)
                    candidate = Hairpin(
                        stem_start=start,
                        stem_length=stem_length,
                        loop_start=loop_start,
                        loop_length=loop_length,
                        stem_end=end,
                        thermo=thermo,
                        melting=melting_temperature(thermo, ions, unimolecular=True),
                        diagram=_hairpin_diagram(sequence, start, stem_length, loop_length),
                    )
                    previous = best.get(key)
                    if previous is None or candidate.free_energy < previous.free_energy:
                        best[key] = candidate
                if not _complementary(
                    sequence[start + stem_length], sequence[end - stem_length]
                ):
                    break
                stem_length += 1
    ordered = sorted(best.values(), key=lambda item: (item.free_energy, item.loop_start))
    return ordered[:limit]


def _runs(paired: Sequence[bool]) -> list[tuple[int, int]]:
    """Half-open spans of consecutive paired columns."""
    spans: list[tuple[int, int]] = []
    start: int | None = None
    for index, flag in enumerate(paired):
        if flag and start is None:
            start = index
        elif not flag and start is not None:
            spans.append((start, index))
            start = None
    if start is not None:
        spans.append((start, len(paired)))
    return spans


def _three_prime_state(
    first: str, second: str, shift: int, scored: Sequence[bool], window: int
) -> tuple[bool, bool]:
    """Whether a scored pair is near, and whether one is at, either 3' terminus.

    The second strand is written 3'->5' beneath the first, so its 3'-terminal
    base faces index ``shift`` of the first strand.
    """
    near = any(scored[-window:])
    for index in range(len(first)):
        partner = len(second) - 1 - index + shift
        if scored[index] and partner >= len(second) - window:
            near = True
            break
    terminal = bool(scored) and scored[-1]
    if 0 <= shift < len(first) and scored[shift]:
        terminal = True
    return near or terminal, terminal


def _dimer_diagram(first: str, second: str, shift: int, scored: Sequence[bool],
                   contacts: Sequence[bool]) -> str:
    """Render one register as three monospaced rows.

    A scored base pair is marked ``|`` and a complementary position that was not
    scored, because its helix is too short or unfavourable, is marked ``.``.
    """
    pad_first = max(0, -shift)
    pad_second = max(0, shift)
    marks = []
    for index in range(len(first)):
        if scored[index]:
            marks.append("|")
        elif contacts[index]:
            marks.append(".")
        else:
            marks.append(" ")
    return "\n".join(
        [
            f"5'-{' ' * pad_first}{first}-3'",
            f"   {' ' * pad_first}{''.join(marks).rstrip()}",
            f"3'-{' ' * pad_second}{second[::-1]}-5'",
        ]
    )


def _score_register(first: str, second: str, shift: int) -> tuple[Duplex, int, int, list[bool], list[bool]]:
    """Score one antiparallel register of two strands."""
    contacts: list[bool] = []
    for index in range(len(first)):
        partner = len(second) - 1 - index + shift
        contacts.append(
            0 <= partner < len(second)
            and _complementary(first[index], second[partner])
        )
    enthalpy = 0.0
    entropy = 0.0
    pairs = 0
    helices = 0
    gc = 0
    scored = [False] * len(first)
    for start, end in _runs(contacts):
        if end - start < MIN_RUN:
            continue
        helix = stack(first[start:end])
        if helix.free_energy() >= 0.0:
            continue
        enthalpy += helix.enthalpy
        entropy += helix.entropy
        pairs += end - start
        helices += 1
        gc += sum(1 for base in first[start:end] if base in "GC")
        for index in range(start, end):
            scored[index] = True
    thermo = Duplex(
        enthalpy=enthalpy,
        entropy=entropy,
        base_pairs=pairs,
        gc_fraction=gc / pairs if pairs else 0.0,
    )
    return thermo, pairs, helices, scored, contacts


def find_dimers(
    first: str,
    second: str,
    ions: IonConditions,
    limit: int = 3,
    identical_strands: bool = False,
) -> list[Dimer]:
    """Find the most stable registers in which two strands associate.

    Args:
        first: a normalised, unambiguous sequence, drawn 5'->3'.
        second: the partner strand, drawn 3'->5' beneath it.
        ions: reaction conditions, used for the melting temperature.
        limit: how many registers to return, most stable first.
        identical_strands: the two strands are the same molecule, as in a
            self-dimer, which removes the factor of four from the
            concentration term of the melting temperature.

    Returns:
        Up to ``limit`` registers sorted by free energy, most stable first,
        followed by the most stable register that pairs a 3' end where that is
        not already among them. Registers in which no helix is favourable are
        omitted entirely.
    """
    found: list[Dimer] = []
    for shift in range(-(len(second) - 1), len(first)):
        thermo, pairs, helices, scored, contacts = _score_register(first, second, shift)
        if not pairs:
            continue
        near, terminal = _three_prime_state(
            first, second, shift, scored, THREE_PRIME_WINDOW
        )
        found.append(
            Dimer(
                shift=shift,
                pairs=pairs,
                helices=helices,
                thermo=thermo,
                melting=melting_temperature(
                    thermo, ions, identical_strands=identical_strands
                ),
                diagram=_dimer_diagram(first, second, shift, scored, contacts),
                three_prime_paired=near,
                terminal_paired=terminal,
            )
        )
    found.sort(key=lambda item: (item.free_energy, item.shift))
    kept = found[:limit]
    # A register that pairs a 3' end is judged against its own, less permissive
    # threshold, so the most stable register overall is not a sufficient
    # summary: a weaker 3' register can breach while the strongest one passes.
    if not any(item.three_prime_paired for item in kept):
        extendable = next((item for item in found if item.three_prime_paired), None)
        if extendable is not None:
            kept.append(extendable)
    return kept


def find_self_dimers(sequence: str, ions: IonConditions, limit: int = 3) -> list[Dimer]:
    """Find the registers in which an oligonucleotide pairs with a copy of itself."""
    return find_dimers(sequence, sequence, ions, limit=limit, identical_strands=True)


def perfect_complement(sequence: str) -> str:
    """The exact reverse complement, the partner of the intended duplex."""
    return reverse_complement(sequence)
