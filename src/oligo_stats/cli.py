"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import IO, Sequence

from oligo_stats import __version__
from oligo_stats.inputs import InputError, Oligo, read, unique_names
from oligo_stats.oligo import (
    CLAMP_WINDOW,
    Analysis,
    Criteria,
    OligoStats,
    analyse_oligo,
    screen_cross_dimers,
)
from oligo_stats.optics import OpticsError, TERMINI
from oligo_stats.report import render_html
from oligo_stats.sequence import SequenceError, expand, expansion_count, validate
from oligo_stats.structure import Dimer, Hairpin
from oligo_stats.thermo import IonConditions, ThermoError

STDOUT_FORMATS = ("table", "tsv", "json", "none")

_TSV_COLUMNS = (
    "name",
    "description",
    "length",
    "gc_percent",
    "tm_c",
    "dh_kcal_mol",
    "ds_cal_mol_k",
    "dg37_kcal_mol",
    "molecular_weight",
    "extinction_260",
    "ug_per_od260",
    "nmol_per_od260",
    "three_prime_stacking",
    "gc_clamp",
    "hairpin_dg",
    "hairpin_tm",
    "hairpin_stem_bp",
    "self_dimer_dg",
    "self_dimer_tm",
    "self_dimer_bp",
    "self_dimer_three_prime",
    "severity",
    "flags",
    "sequence",
)


def _finite(text: str) -> float:
    """Parse a real number, refusing NaN and the infinities.

    A non-finite concentration does not announce itself downstream: every
    comparison against NaN is false, so the ion-regime test falls through to
    its last branch and the run reports a confident melting temperature for a
    buffer that was never specified. JSON has no spelling for these values
    either.
    """
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if value != value or value in (float("inf"), float("-inf")):
        raise argparse.ArgumentTypeError(f"not a finite number: {text!r}")
    return value


def _concentration(text: str) -> float:
    """Parse a non-negative, finite concentration."""
    value = _finite(text)
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative: {value}")
    return value


def _positive(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or greater: {value}")
    return value


def _range(text: str) -> tuple[float, float]:
    """Parse a MIN:MAX or MIN-MAX window."""
    separator = ":" if ":" in text else "-"
    fields = text.split(separator)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError(f"not a MIN:MAX range: {text!r}")
    low, high = _finite(fields[0]), _finite(fields[1])
    if low > high:
        raise argparse.ArgumentTypeError(f"range is inverted: {text!r}")
    return low, high


def _integer_range(text: str) -> tuple[int, int]:
    """Parse a MIN:MAX window over whole numbers.

    Length and the GC clamp are counts of bases. Accepting ``18.9:30.1`` and
    silently applying ``18:30`` reports a rule the caller did not ask for.
    """
    low, high = _range(text)
    for value in (low, high):
        if value != int(value):
            raise argparse.ArgumentTypeError(
                f"must be a whole number of bases: {value:g}"
            )
    return int(low), int(high)


def _free_energy(text: str) -> float:
    """Parse a free-energy threshold, which must be finite and not positive."""
    value = _finite(text)
    if value > 0:
        raise argparse.ArgumentTypeError(
            f"a free-energy threshold must be zero or negative: {value}"
        )
    return value


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="oligo-stats",
        description=(
            "Measure PCR primers and other short oligonucleotides: melting "
            "temperature from nearest-neighbour thermodynamics, base "
            "composition, molecular weight, extinction coefficient, and a "
            "screen for hairpins and dimers. Ion concentrations default to a "
            "standard amplification reaction; set them explicitly to compare "
            "against a supplier's specification sheet, which is quoted at "
            "50 mM Na+ with no magnesium."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"oligo-stats {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyse = subparsers.add_parser(
        "analyse",
        aliases=["analyze"],
        help="measure and screen one or more oligonucleotides",
        description=(
            "Measure every oligonucleotide given on the command line or read "
            "from a file. Sequences must be unambiguous; resolve ambiguity "
            "codes first with 'oligo-stats expand'."
        ),
    )
    analyse.add_argument(
        "sequences",
        nargs="*",
        metavar="SEQUENCE",
        help="oligonucleotide sequence, repeatable",
    )
    analyse.add_argument(
        "-i",
        "--input",
        metavar="PATH",
        help=(
            "FASTA file, or a table of name and sequence columns, or '-' to "
            "read standard input"
        ),
    )
    analyse.add_argument(
        "-c",
        "--conc",
        type=_concentration,
        default=0.25,
        metavar="UM",
        help="total oligonucleotide strand concentration, uM (default: 0.25)",
    )
    analyse.add_argument(
        "--na",
        type=_concentration,
        default=0.0,
        metavar="MM",
        help="sodium concentration, mM (default: 0)",
    )
    analyse.add_argument(
        "--k",
        type=_concentration,
        default=50.0,
        metavar="MM",
        help="potassium concentration, mM (default: 50)",
    )
    analyse.add_argument(
        "--tris",
        type=_concentration,
        default=10.0,
        metavar="MM",
        help=(
            "Tris buffer concentration, mM, of which half is counted as "
            "monovalent cation (default: 10)"
        ),
    )
    analyse.add_argument(
        "--mg",
        type=_concentration,
        default=1.5,
        metavar="MM",
        help="magnesium concentration, mM (default: 1.5)",
    )
    analyse.add_argument(
        "--dntp",
        type=_concentration,
        default=0.8,
        metavar="MM",
        help=(
            "total concentration of all four dNTPs, mM, which chelate "
            "magnesium one for one (default: 0.8, i.e. 0.2 mM of each)"
        ),
    )
    analyse.add_argument(
        "--five-prime",
        choices=TERMINI,
        default="hydroxyl",
        help="5' terminus as synthesised, which sets the mass (default: hydroxyl)",
    )
    analyse.add_argument(
        "--tm-range",
        type=_range,
        default=(55.0, 65.0),
        metavar="MIN:MAX",
        help="acceptable melting temperature window, Celsius (default: 55:65)",
    )
    analyse.add_argument(
        "--gc-range",
        type=_range,
        default=(40.0, 60.0),
        metavar="MIN:MAX",
        help="acceptable GC content window, per cent (default: 40:60)",
    )
    analyse.add_argument(
        "--length-range",
        type=_integer_range,
        default=(18, 30),
        metavar="MIN:MAX",
        help="acceptable length window, nucleotides (default: 18:30)",
    )
    analyse.add_argument(
        "--clamp-range",
        type=_integer_range,
        default=(1, 3),
        metavar="MIN:MAX",
        help=(
            "acceptable count of G and C among the last five bases, the 3' GC "
            "clamp (default: 1:3)"
        ),
    )
    analyse.add_argument(
        "--end-stability",
        type=_free_energy,
        default=-9.0,
        metavar="DG",
        help=(
            "flag a 3' end whose last five bases stack at or below this free "
            "energy, kcal/mol (default: -9.0)"
        ),
    )
    analyse.add_argument(
        "--max-run",
        type=_positive,
        default=4,
        metavar="N",
        help="flag a run of this many identical bases or more (default: 4)",
    )
    analyse.add_argument(
        "--hairpin-threshold",
        type=_free_energy,
        default=-2.0,
        metavar="DG",
        help=(
            "flag a hairpin at or below this free energy at 37 C, kcal/mol "
            "(default: -2.0)"
        ),
    )
    analyse.add_argument(
        "--dimer-threshold",
        type=_free_energy,
        default=-9.0,
        metavar="DG",
        help=(
            "flag a dimer at or below this free energy at 37 C, kcal/mol "
            "(default: -9.0)"
        ),
    )
    analyse.add_argument(
        "--dimer-3p-threshold",
        type=_free_energy,
        default=-5.0,
        metavar="DG",
        help=(
            "flag a dimer that pairs a 3' end at or below this free energy, "
            "kcal/mol (default: -5.0)"
        ),
    )
    analyse.add_argument(
        "--structures",
        type=_positive,
        default=3,
        metavar="N",
        help="hairpins and self-dimer registers to report each (default: 3)",
    )
    analyse.add_argument(
        "--no-structures",
        action="store_true",
        help="skip the hairpin and dimer screens entirely",
    )
    analyse.add_argument(
        "--no-cross-dimers",
        action="store_true",
        help="skip the pairwise cross-dimer screen",
    )
    analyse.add_argument(
        "-f",
        "--format",
        dest="stdout_format",
        choices=STDOUT_FORMATS,
        default="table",
        help="what to write to standard output (default: table)",
    )
    analyse.add_argument(
        "-o", "--output", metavar="PATH", help="write the measurements as TSV"
    )
    analyse.add_argument(
        "-r", "--report", metavar="PATH", help="write the self-contained HTML report"
    )
    analyse.add_argument(
        "--strict",
        action="store_true",
        help=(
            "exit non-zero when any oligonucleotide is flagged for review or "
            "any pair reaches the cross-dimer threshold"
        ),
    )

    expander = subparsers.add_parser(
        "expand",
        help="resolve IUPAC ambiguity codes into unambiguous sequences",
        description=(
            "Write out every unambiguous sequence a degenerate oligonucleotide "
            "stands for, as FASTA, so that each can be measured in its own "
            "right."
        ),
    )
    expander.add_argument(
        "sequences", nargs="*", metavar="SEQUENCE", help="sequence, repeatable"
    )
    expander.add_argument(
        "-i", "--input", metavar="PATH", help="FASTA or table input, or '-' for stdin"
    )
    expander.add_argument(
        "--limit",
        type=_positive,
        default=256,
        metavar="N",
        help="largest expansion to produce per sequence (default: 256)",
    )
    expander.add_argument(
        "-o", "--output", metavar="PATH", help="write FASTA here instead of stdout"
    )
    return parser


def _collect(args: argparse.Namespace) -> list[Oligo]:
    """Gather oligonucleotides from positional arguments and from the input."""
    oligos: list[Oligo] = []
    for index, sequence in enumerate(args.sequences, start=1):
        oligos.append(Oligo(f"oligo{index}", "", sequence))
    if args.input:
        source = sys.stdin if args.input == "-" else args.input
        oligos.extend(read(source))
    if not oligos:
        raise InputError(
            "no oligonucleotide given: pass a sequence, or --input to read a file"
        )
    return oligos


def _ions(args: argparse.Namespace) -> IonConditions:
    return IonConditions(
        strand=args.conc * 1e-6,
        sodium=args.na * 1e-3,
        potassium=args.k * 1e-3,
        tris=args.tris * 1e-3,
        magnesium=args.mg * 1e-3,
        dntp=args.dntp * 1e-3,
    )


def _criteria(args: argparse.Namespace) -> Criteria:
    return Criteria(
        length_min=args.length_range[0],
        length_max=args.length_range[1],
        gc_min=args.gc_range[0],
        gc_max=args.gc_range[1],
        melting_min=args.tm_range[0],
        melting_max=args.tm_range[1],
        clamp_min=args.clamp_range[0],
        clamp_max=args.clamp_range[1],
        end_stability=args.end_stability,
        max_run=args.max_run,
        hairpin_free_energy=args.hairpin_threshold,
        dimer_free_energy=args.dimer_threshold,
        dimer_three_prime_free_energy=args.dimer_3p_threshold,
    )


def _tsv_field(text: str) -> str:
    """Make text safe to place in a tab-separated column.

    A FASTA description may legitimately contain a tab, which would otherwise
    add a column to that row alone and silently misalign every field after it.
    """
    return " ".join(text.split())


def _optional(value: float | None, digits: int = 2) -> str:
    """Format a value for a machine-read column, leaving it empty when absent."""
    return "" if value is None else f"{value:.{digits}f}"


def _shown(value: float | None, digits: int = 2) -> str:
    """Format a value for a person to read, naming an absent one."""
    return "n/a" if value is None else f"{value:.{digits}f}"


def format_tsv(analysis: Analysis) -> str:
    """Render the measurements as a tab-separated table."""
    lines = ["\t".join(_TSV_COLUMNS)]
    for item in analysis.oligos:
        hairpin = item.best_hairpin
        dimer = item.best_self_dimer
        row = [
            _tsv_field(item.name),
            _tsv_field(item.description),
            str(item.length),
            _optional(item.composition.gc_percent),
            _optional(item.melting),
            f"{item.thermo.enthalpy:.2f}",
            f"{item.thermo.entropy:.2f}",
            f"{item.free_energy:.2f}",
            f"{item.optics.molecular_weight:.2f}",
            f"{item.optics.extinction:.1f}",
            f"{item.optics.micrograms_per_od:.2f}",
            f"{item.optics.nanomoles_per_od:.2f}",
            _optional(item.end_stability),
            str(item.gc_clamp),
            _optional(hairpin.free_energy) if hairpin else "",
            _optional(hairpin.melting, 1) if hairpin else "",
            str(hairpin.stem_length) if hairpin else "",
            _optional(dimer.free_energy) if dimer else "",
            _optional(dimer.melting, 1) if dimer else "",
            str(dimer.pairs) if dimer else "",
            str(dimer.three_prime_paired).lower() if dimer else "",
            item.worst_severity,
            ";".join(flag.code for flag in item.flags),
            item.sequence,
        ]
        assert len(row) == len(_TSV_COLUMNS)
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def _hairpin_json(hairpin: Hairpin) -> dict:
    return {
        "free_energy": round(hairpin.free_energy, 3),
        "melting_temperature": None
        if hairpin.melting is None
        else round(hairpin.melting, 2),
        "enthalpy": round(hairpin.thermo.enthalpy, 2),
        "entropy": round(hairpin.thermo.entropy, 2),
        "stem_start": hairpin.stem_start + 1,
        "stem_length": hairpin.stem_length,
        "loop_start": hairpin.loop_start + 1,
        "loop_length": hairpin.loop_length,
        "diagram": hairpin.diagram,
    }


def _dimer_json(dimer: Dimer) -> dict:
    return {
        "free_energy": round(dimer.free_energy, 3),
        "melting_temperature": None
        if dimer.melting is None
        else round(dimer.melting, 2),
        "enthalpy": round(dimer.thermo.enthalpy, 2),
        "entropy": round(dimer.thermo.entropy, 2),
        "shift": dimer.shift,
        "base_pairs": dimer.pairs,
        "helices": dimer.helices,
        "three_prime_paired": dimer.three_prime_paired,
        "terminal_paired": dimer.terminal_paired,
        "diagram": dimer.diagram,
    }


def _oligo_json(item: OligoStats) -> dict:
    return {
        "name": item.name,
        "description": item.description,
        "sequence": item.sequence,
        "length": item.length,
        "composition": {
            "counts": item.composition.counts,
            "gc_percent": item.composition.gc_percent,
        },
        "thermodynamics": {
            "enthalpy": round(item.thermo.enthalpy, 2),
            "entropy": round(item.thermo.entropy, 2),
            "free_energy": round(item.free_energy, 3),
            "melting_temperature": None
            if item.melting is None
            else round(item.melting, 2),
            "self_complementary": item.self_complementary,
        },
        "optics": {
            "molecular_weight": item.optics.molecular_weight,
            "extinction_260": item.optics.extinction,
            "micrograms_per_od260": round(item.optics.micrograms_per_od, 3),
            "nanomoles_per_od260": round(item.optics.nanomoles_per_od, 3),
            "terminus": item.optics.terminus,
        },
        "three_prime": {
            "stacking_free_energy": None
            if item.end_stability is None
            else round(item.end_stability, 3),
            "gc_clamp": item.gc_clamp,
            "window": CLAMP_WINDOW,
        },
        "structures_screened": item.structures_screened,
        "hairpins": [_hairpin_json(hairpin) for hairpin in item.hairpins],
        "self_dimers": [_dimer_json(dimer) for dimer in item.self_dimers],
        "severity": item.worst_severity,
        "flags": [
            {"code": flag.code, "severity": flag.severity, "message": flag.message}
            for flag in item.flags
        ],
    }


def format_json(analysis: Analysis) -> str:
    """Render the whole analysis as JSON."""
    ions = analysis.ions
    document = {
        "version": __version__,
        "conditions": {
            "strand_molar": ions.strand,
            "sodium_molar": ions.sodium,
            "potassium_molar": ions.potassium,
            "tris_molar": ions.tris,
            "magnesium_molar": ions.magnesium,
            "dntp_molar": ions.dntp,
            "monovalent_molar": round(ions.monovalent, 9),
            "free_magnesium_molar": round(ions.free_magnesium, 9),
            "regime": ions.regime,
            "terminus": analysis.terminus,
        },
        "criteria": {
            "length": [analysis.criteria.length_min, analysis.criteria.length_max],
            "gc_percent": [analysis.criteria.gc_min, analysis.criteria.gc_max],
            "melting_temperature": [
                analysis.criteria.melting_min,
                analysis.criteria.melting_max,
            ],
            "gc_clamp": [analysis.criteria.clamp_min, analysis.criteria.clamp_max],
            "end_stability": analysis.criteria.end_stability,
            "max_run": analysis.criteria.max_run,
            "hairpin_free_energy": analysis.criteria.hairpin_free_energy,
            "dimer_free_energy": analysis.criteria.dimer_free_energy,
            "dimer_three_prime_free_energy": (
                analysis.criteria.dimer_three_prime_free_energy
            ),
        },
        "oligos": [_oligo_json(item) for item in analysis.oligos],
        "cross_dimers_screened": analysis.cross_dimers_screened,
        "cross_dimers": [
            {
                "first": entry.first,
                "second": entry.second,
                **_dimer_json(entry.dimer),
            }
            for entry in analysis.cross_dimers
        ],
    }
    return json.dumps(document, indent=2, sort_keys=False) + "\n"


def format_table(analysis: Analysis) -> str:
    """Render a readable per-oligonucleotide summary for the terminal."""
    lines: list[str] = []
    for item in analysis.oligos:
        header = item.name
        if item.description:
            header = f"{header}  {item.description}"
        lines.append(f"> {header}")
        lines.append(
            f"  {item.length} nt  GC {_shown(item.composition.gc_percent)}%  "
            f"Tm {_shown(item.melting)} C  dG37 {item.free_energy:.2f} kcal/mol"
        )
        lines.append(
            f"  MW {item.optics.molecular_weight:,.2f} g/mol  "
            f"e260 {item.optics.extinction:,.0f} /(M cm)  "
            f"{item.optics.micrograms_per_od:,.2f} ug/OD  "
            f"{item.optics.nanomoles_per_od:,.2f} nmol/OD"
        )
        lines.append(
            f"  3' clamp {item.gc_clamp}/{CLAMP_WINDOW}  "
            f"3' stacking {_shown(item.end_stability)} kcal/mol"
        )
        hairpin = item.best_hairpin
        dimer = item.best_self_dimer
        if item.structures_screened:
            lines.append(
                "  hairpin "
                + (
                    f"{hairpin.free_energy:+.2f} kcal/mol "
                    f"({hairpin.stem_length} bp stem)"
                    if hairpin
                    else "none possible"
                )
                + "  self-dimer "
                + (
                    f"{dimer.free_energy:+.2f} kcal/mol ({dimer.pairs} bp"
                    + (", 3' paired)" if dimer.three_prime_paired else ")")
                    if dimer
                    else "none favourable"
                )
            )
        if item.flags:
            for flag in item.flags:
                lines.append(f"  {flag.severity:<5} {flag.message}")
        else:
            lines.append("  pass")
        lines.append("")
    span = analysis.melting_span
    if span is not None and len(analysis.oligos) > 1:
        lines.append(
            f"Tm across the panel: {span[0]:.2f} to {span[1]:.2f} C "
            f"(span {span[1] - span[0]:.2f} C)"
        )
    if analysis.cross_dimers_screened:
        pairs = len(analysis.oligos) * (len(analysis.oligos) - 1) // 2
        lines.append(
            f"{pairs} pair{'' if pairs == 1 else 's'} screened for cross-dimers; "
            f"{len(analysis.cross_dimers)} at or below threshold"
        )
        for entry in analysis.cross_dimers:
            lines.append(
                f"  {entry.first} + {entry.second}  "
                f"{entry.dimer.free_energy:+.2f} kcal/mol "
                f"({entry.dimer.pairs} bp"
                + (", 3' paired)" if entry.dimer.three_prime_paired else ")")
            )
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n" if lines else ""


def _write(path: str, text: str) -> None:
    Path(path).write_text(text, encoding="utf-8")


def _same_file(first: Path, second: Path) -> bool:
    """Whether two paths name one file on disk.

    Comparing resolved names is not enough: a hard link gives one file two
    names that resolve differently, and writing through either truncates it.
    """
    if first == second:
        return True
    if not (first.exists() and second.exists()):
        return False
    return first.samefile(second)


def _check_output_paths(args: argparse.Namespace) -> None:
    """Refuse to overwrite the input or write two outputs to one file.

    Each output is written independently, so sharing a path would silently
    leave only whichever was written last, and a path that reaches the input by
    any route destroys the source once it has been read.
    """
    requested = {"--output": args.output, "--report": getattr(args, "report", None)}
    claimed: dict[Path, str] = {}
    input_path = (
        None
        if not args.input or args.input == "-"
        else Path(args.input).expanduser().resolve()
    )
    for flag, value in requested.items():
        if value is None:
            continue
        resolved = Path(value).expanduser().resolve()
        if input_path is not None and _same_file(resolved, input_path):
            raise ValueError(f"{flag} would overwrite the input file: {value}")
        for other_path, other_flag in claimed.items():
            if _same_file(resolved, other_path):
                raise ValueError(f"{other_flag} and {flag} would both write to {value}")
        claimed[resolved] = flag


def run_analyse(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the analyse subcommand. Returns a process exit code."""
    _check_output_paths(args)
    oligos = _collect(args)
    names = unique_names(oligos)
    ions = _ions(args)
    criteria = _criteria(args)
    measured = [
        analyse_oligo(
            item.sequence,
            name=name,
            description=item.description,
            ions=ions,
            criteria=criteria,
            terminus=args.five_prime,
            structures=args.structures,
            screen_structures=not args.no_structures,
        )
        for item, name in zip(oligos, names)
    ]
    screen = (
        not args.no_cross_dimers and not args.no_structures and len(measured) > 1
    )
    analysis = Analysis(
        oligos=tuple(measured),
        ions=ions,
        criteria=criteria,
        terminus=args.five_prime,
        cross_dimers=tuple(screen_cross_dimers(measured, ions, criteria))
        if screen
        else (),
        cross_dimers_screened=screen,
    )

    if args.report:
        source = (
            "standard input"
            if args.input == "-"
            else Path(args.input).name
            if args.input
            else "command line"
        )
        _write(args.report, render_html(analysis, source=source, version=__version__))
    if args.output:
        _write(args.output, format_tsv(analysis))

    if args.stdout_format == "table":
        stdout.write(format_table(analysis))
    elif args.stdout_format == "tsv":
        stdout.write(format_tsv(analysis))
    elif args.stdout_format == "json":
        stdout.write(format_json(analysis))

    if args.strict and (
        any(item.worst_severity == "warn" for item in analysis.oligos)
        or analysis.cross_dimers
    ):
        return 1
    return 0


def run_expand(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the expand subcommand. Returns a process exit code."""
    _check_output_paths(args)
    oligos = _collect(args)
    names = unique_names(oligos)
    lines: list[str] = []
    for item, name in zip(oligos, names):
        sequence = validate(item.sequence, allow_degenerate=True)
        total = expansion_count(sequence)
        if total == 1:
            lines.append(f">{name}")
            lines.append(sequence)
            continue
        width = len(str(total))
        for index, variant in enumerate(expand(sequence, limit=args.limit), start=1):
            lines.append(f">{name}.{index:0{width}d} of {total}")
            lines.append(variant)
    text = "\n".join(lines) + "\n"
    if args.output:
        _write(args.output, text)
    else:
        stdout.write(text)
    return 0


def main(argv: Sequence[str] | None = None, stdout: IO[str] | None = None) -> int:
    """Entry point. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    out = sys.stdout if stdout is None else stdout
    try:
        if args.command in ("analyse", "analyze"):
            return run_analyse(args, out)
        if args.command == "expand":
            return run_expand(args, out)
    except (
        InputError,
        SequenceError,
        ThermoError,
        OpticsError,
        FileNotFoundError,
        IsADirectoryError,
        OSError,
        ValueError,
    ) as error:
        parser.exit(2, f"oligo-stats: error: {error}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
