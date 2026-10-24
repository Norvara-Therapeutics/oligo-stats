"""Readers for the two ways an oligonucleotide list arrives.

A list of primers is kept either as FASTA or as the name-and-sequence table that
a synthesis order sheet uses. Both are accepted from a path or an open handle
and the format is decided by the first non-blank line, so neither needs to be
declared. Record order is preserved in both cases, because a primer list is
usually written in the order the primers are used.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Iterator

#: Characters permitted in sequence text: the IUPAC nucleotide codes.
VALID_RESIDUES: frozenset[str] = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn")

#: Column separators recognised in a table of oligonucleotides.
_SEPARATORS = ("\t", ",", ";")

#: Column names that identify a header row. Recognising a header by name rather
#: than by "the second field does not parse as sequence" matters: the loose test
#: silently swallows a mistyped base in the first data row, which is exactly the
#: row a reader would trust.
_HEADER_FIELDS = frozenset(
    {"sequence", "seq", "oligo", "oligonucleotide", "primer", "bases", "dna"}
)


class InputError(ValueError):
    """Raised when input cannot be parsed as a list of oligonucleotides."""


@dataclass(frozen=True)
class Oligo:
    """One named oligonucleotide as it was supplied.

    Attributes:
        name: identifier, taken from the FASTA header or the table's first
            column, or generated when the input gives none.
        description: any text after the identifier on a FASTA header line, with
            internal whitespace collapsed to single spaces.
        sequence: sequence text with whitespace removed and case preserved.
    """

    name: str
    description: str
    sequence: str


def _check_residues(text: str, lineno: int) -> None:
    """Reject text that is not sequence.

    Without this, a file that is not a sequence at all, such as prose or a
    saved error page, is absorbed as bases and reported with a melting
    temperature, which is worse than an error.
    """
    invalid = sorted({character for character in text if character not in VALID_RESIDUES})
    if invalid:
        shown = ", ".join(repr(character) for character in invalid[:5])
        if len(invalid) > 5:
            shown += ", ..."
        raise InputError(
            f"line {lineno}: not nucleotide sequence; unexpected character(s): {shown}"
        )


def parse_fasta(lines: Iterable[str]) -> Iterator[Oligo]:
    """Yield oligonucleotides from FASTA text."""
    name: str | None = None
    description = ""
    chunks: list[str] = []
    seen_header = False

    for lineno, raw in enumerate(lines, start=1):
        if lineno == 1:
            raw = raw.lstrip("\ufeff")
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield Oligo(name, description, "".join(chunks))
            seen_header = True
            head = line[1:].strip()
            if not head:
                raise InputError(f"line {lineno}: header has no identifier")
            fields = head.split(None, 1)
            name = fields[0]
            # Collapse internal whitespace: a description may carry a tab, which
            # would add a column to a tab-separated row written downstream.
            description = " ".join(fields[1].split()) if len(fields) > 1 else ""
            chunks = []
        else:
            if not seen_header:
                raise InputError(f"line {lineno}: sequence data before any '>' header")
            residues = "".join(line.split())
            _check_residues(residues, lineno)
            chunks.append(residues)

    if name is not None:
        yield Oligo(name, description, "".join(chunks))
    elif not seen_header:
        raise InputError("input contains no FASTA records")


def _split_row(line: str) -> list[str]:
    for separator in _SEPARATORS:
        if separator in line:
            return [field.strip() for field in line.split(separator)]
    return line.split()


def parse_table(lines: Iterable[str]) -> Iterator[Oligo]:
    """Yield oligonucleotides from a name-and-sequence table.

    Each non-blank line holds either a bare sequence or a name and a sequence
    separated by a tab, a comma, a semicolon or whitespace. A line beginning
    with ``#`` is a comment. A first row whose second field is one of the
    recognised column names, such as ``sequence`` or ``seq``, is a header and is
    skipped; any other first row is data, so a mistyped base in it is an error
    rather than a silently dropped record. An unnamed row is called ``oligo<N>``
    after its position in the input. Cells are read by position: a trailing
    separator is ignored, but an empty sequence column is an error rather than
    an invitation to read the next column instead.
    """
    index = 0
    for lineno, raw in enumerate(lines, start=1):
        if lineno == 1:
            raw = raw.lstrip("\ufeff")
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = _split_row(line)
        # A trailing separator is harmless, but an empty cell between two filled
        # ones is a hole in the row: dropping it would slide the sequence in
        # from whatever column happened to come next.
        while fields and not fields[-1]:
            fields.pop()
        if not fields:
            continue
        index += 1
        if len(fields) == 1:
            sequence = "".join(fields[0].split())
            _check_residues(sequence, lineno)
            yield Oligo(f"oligo{index}", "", sequence)
            continue
        name, sequence = fields[0], "".join(fields[1].split())
        if not sequence:
            raise InputError(
                f"line {lineno}: the sequence column is empty for {name!r}"
            )
        if index == 1 and sequence.lower() in _HEADER_FIELDS:
            # A header row such as "name,sequence": drop it and renumber.
            index = 0
            continue
        _check_residues(sequence, lineno)
        yield Oligo(name, " ".join(fields[2:]), sequence)


def parse(lines: Iterable[str]) -> list[Oligo]:
    """Parse FASTA or a table, deciding from the first non-blank line."""
    materialised = list(lines)
    for line in materialised:
        stripped = line.strip().lstrip("\ufeff")
        if not stripped:
            continue
        if stripped.startswith(">"):
            return list(parse_fasta(materialised))
        break
    oligos = list(parse_table(materialised))
    if not oligos:
        raise InputError("input contains no oligonucleotides")
    return oligos


def read(source: str | Path | IO[str]) -> list[Oligo]:
    """Read every oligonucleotide from a path or an open text handle."""
    if hasattr(source, "read"):
        return parse(source)  # type: ignore[arg-type]
    path = Path(source)
    if path.is_dir():
        raise IsADirectoryError(f"not a file: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"no such file: {path}")
    with path.open("r", encoding="utf-8-sig") as handle:
        return parse(handle)


def unique_names(oligos: Iterable[Oligo]) -> list[str]:
    """Name each oligonucleotide uniquely, suffixing only names that repeat.

    The second and later occurrences of a name become ``<name>.2``, ``<name>.3``
    and so on; a name that appears once is unchanged. A generated name is
    advanced past any collision with a name already in the input.
    """
    items = list(oligos)
    totals: dict[str, int] = {}
    for item in items:
        totals[item.name] = totals.get(item.name, 0) + 1
    taken = set(totals)
    seen: dict[str, int] = {}
    names: list[str] = []
    for item in items:
        seen[item.name] = seen.get(item.name, 0) + 1
        if totals[item.name] == 1 or seen[item.name] == 1:
            names.append(item.name)
            continue
        suffix = seen[item.name]
        candidate = f"{item.name}.{suffix}"
        while candidate in taken:
            suffix += 1
            candidate = f"{item.name}.{suffix}"
        taken.add(candidate)
        names.append(candidate)
    return names
