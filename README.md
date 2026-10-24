# oligo-stats

PCR primer and oligonucleotide analysis from the command line. Melting
temperature from nearest-neighbour thermodynamics with explicit ion
corrections, base composition, molecular weight, extinction coefficient, and a
screen for hairpins and dimers. Pure Python with no third-party dependencies,
and a self-contained HTML report that opens from a local disk without a network
connection.

## Install

Install the prerelease from PyPI:

```bash
python -m pip install --pre norvara-oligo-stats
```

The Python package and command keep the shorter project name:

```python
import oligo_stats
```

```bash
oligo-stats --version
```

## Quickstart

```bash
oligo-stats analyse TCATGAAGGATGGTGTTAATGC
```

```
> oligo1
  22 nt  GC 40.91%  Tm 57.47 C  dG37 -25.58 kcal/mol
  MW 6,829.53 g/mol  e260 220,500 /(M cm)  30.97 ug/OD  4.54 nmol/OD
  3' clamp 2/5  3' stacking -5.58 kcal/mol
  hairpin +1.94 kcal/mol (3 bp stem)  self-dimer -4.36 kcal/mol (6 bp)
  pass
```

A whole panel, with the HTML report:

```bash
oligo-stats analyse -i examples/primers.fasta --report primers.html
```

```
> phiX_hairpin  NC_001422.1 3953-3974 forward, folds into an 8 bp stem
  22 nt  GC 54.55%  Tm 63.32 C  dG37 -28.27 kcal/mol
  MW 6,695.40 g/mol  e260 203,200 /(M cm)  32.95 ug/OD  4.92 nmol/OD
  3' clamp 4/5  3' stacking -6.66 kcal/mol
  hairpin -6.00 kcal/mol (8 bp stem)  self-dimer -17.93 kcal/mol (16 bp, 3' paired)
  note  4 of the last 5 bases are G or C, which stabilises mismatched 3' ends as well as matched ones
  warn  hairpin of -6.00 kcal/mol with an 8 bp stem, melting at 60.2 C
  warn  self-dimer of -17.93 kcal/mol over 16 bp, pairing a 3'-terminal base, so a polymerase can extend it into a primer dimer
```

The report adds the full composition breakdown, the thermodynamic terms, the
design criteria exactly as they were configured, the pairwise cross-dimer
screen, and drawn structures: the most stable hairpin and dimer, plus any other
retained structure that reaches its threshold. `--structures N` sets how many
are retained per oligonucleotide in the first place.

## Input

Sequences are given on the command line, or in a file with `--input`, or on
standard input with `--input -`. Positional sequences and a file may be
combined. A file is read as FASTA when its first non-blank line begins with
`>`, and otherwise as a table of name and sequence columns separated by a tab,
a comma, a semicolon or whitespace, with `#` comments. A first row is treated
as a header only when its second field is a recognised column name such as
`sequence` or `seq`; any other first row is data, so a mistyped base in it is an
error rather than a silently dropped primer. Case, internal whitespace and
uracil are all tolerated.

```bash
oligo-stats analyse -i examples/order_sheet.tsv
printf 'F1\tTCATGAAGGATGGTGTTAATGC\n' | oligo-stats analyse -i -
```

Nearest-neighbour parameters are defined for A, C, G and T only, so a sequence
carrying an IUPAC ambiguity code is refused. `oligo-stats expand` resolves one
into the sequences it stands for, which pipes straight back in:

```bash
oligo-stats expand -i degenerate.fasta | oligo-stats analyse -i -
```

## Options

```
-i, --input PATH          FASTA or name-and-sequence table, or '-' for stdin
-c, --conc UM             total oligonucleotide strand concentration
                          (default: 0.25)
    --na MM               sodium (default: 0)
    --k MM                potassium (default: 50)
    --tris MM             Tris buffer, half counted as cation (default: 10)
    --mg MM               magnesium (default: 1.5)
    --dntp MM             total dNTP, which chelates magnesium (default: 0.8)
    --five-prime {hydroxyl,phosphate}
                          5' terminus as synthesised, which sets the mass
    --tm-range MIN:MAX    acceptable Tm window, Celsius (default: 55:65)
    --gc-range MIN:MAX    acceptable GC window, per cent (default: 40:60)
    --length-range MIN:MAX
                          acceptable length window, nt, whole numbers
                          (default: 18:30)
    --clamp-range MIN:MAX acceptable G and C count among the last five bases,
                          whole numbers (default: 1:3)
    --end-stability DG    flag a 3' end stacking at or below, kcal/mol
                          (default: -9.0)
    --max-run N           flag a run of this many identical bases
                          (default: 4)
    --hairpin-threshold DG      flag at or below, kcal/mol (default: -2.0)
    --dimer-threshold DG        flag at or below, kcal/mol (default: -9.0)
    --dimer-3p-threshold DG     the same where a 3' end is paired
                                (default: -5.0)
    --structures N        hairpins and dimer registers to keep each (default: 3)
    --no-structures       skip the hairpin and dimer screens
    --no-cross-dimers     skip the pairwise cross-dimer screen
    --strict              exit non-zero when anything is flagged for review,
                          including a pair reaching the cross-dimer threshold
```

The ion defaults are a standard amplification reaction, because that is the
buffer a primer anneals in and the one an annealing temperature has to be
derived from. A supplier's specification sheet is quoted at 50 mM Na+ with no
magnesium, which melts the same oligonucleotide several degrees lower; to
reproduce it, set the ions explicitly:

```bash
oligo-stats analyse TCATGAAGGATGGTGTTAATGC --na 50 --k 0 --tris 0 --mg 0 --dntp 0
```

## Outputs

| Flag | Output |
| --- | --- |
| `--format table` (default) | the summary above, on stdout |
| `--format tsv` \| `json` \| `none` | the same on stdout, or nothing |
| `--output PATH` | one row per oligonucleotide as TSV |
| `--report PATH` | the self-contained HTML report |

Any combination may be requested in a single pass. JSON carries everything: the
thermodynamic terms, the drawn structures, every threshold that was applied,
and whether the structure screens were run at all, so a flag in the output can
always be traced to the rule that raised it. TSV leaves an unmeasurable value
empty, where the terminal summary and the report both name it `n/a`, and no
output is ever non-finite, so the JSON is always strict JSON. `--strict` makes
the exit code usable as a gate in an ordering pipeline: it fails on an
individual review flag and on any pair that reaches the cross-dimer threshold,
which is the finding a per-oligo check cannot see.

## What is computed

- **Melting temperature** is the two-state nearest-neighbour value for the
  oligonucleotide against its perfect complement, using the unified DNA
  parameters of SantaLucia (1998), PNAS 95:1460. Only enthalpy and entropy are
  tabulated; free energy is derived from them, so one source sets every
  reported quantity. The 1 M Na+ result is corrected to the stated ion
  composition by Owczarzy et al. (2004), Biochemistry 43:3537, and where
  magnesium is present by Owczarzy et al. (2008), Biochemistry 47:5336, whose
  ratio of divalent to monovalent cation selects the monovalent, mixed or
  divalent relation. Deoxynucleotides chelate magnesium one for one and only
  the free magnesium counts. Tris is counted as half monovalent cation.
- **Molecular weight** is anhydrous and excludes counterions. The terminal
  adjustment gives a 5'-hydroxyl, which is what standard synthesis delivers;
  `--five-prime phosphate` adds the phosphate.
- **Extinction coefficient** at 260 nm is the single-strand nearest-neighbour
  model over the dinucleotide and mononucleotide values of Warshaw & Cantor
  (1970), Biopolymers 9:1079. The micrograms and nanomoles per unit absorbance
  follow from it.
- **3' end stacking** is the free energy of the last five bases with initiation
  excluded. A strongly stacked 3' end tolerates mismatches further in and so
  primes off-target sites; Primer3 gates on the same quantity.
- **Hairpins** are enumerated exhaustively over every contiguous stem of three
  or more base pairs closing a loop of three or more nucleotides, scored as the
  stem's stacks plus the loop initiation of SantaLucia & Hicks (2004), Annu Rev
  Biophys Biomol Struct 33:415, treated as purely entropic and extrapolated
  between tabulated loop sizes by the Jacobson-Stockmayer relation.
- **Dimers** are enumerated over every antiparallel register of the two
  strands. Each contiguous Watson-Crick run of two or more pairs in a register
  is scored as its own helix and pays initiation at both ends, and a run counts
  only where it is favourable at 37 C, so a fragmented register scores as the
  short helices it really is. A register that pairs either 3' end is marked,
  because that is the one that extends into a primer dimer, and it is judged
  against a less permissive threshold.

## Reading the output

- Free energies are kcal/mol at 37 C and enthalpies kcal/mol, both for
  formation, so both are negative for a stabilising interaction. Entropies are
  cal/(mol K).
- A melting temperature is reported only above 0 C. Below the freezing point of
  the buffer the two-state model and the ion corrections are extrapolations and
  the structure does not exist under any usable condition. An oligonucleotide
  with no melting temperature in range is flagged for review rather than passed:
  the melting window could not be applied to it, which is not the same as
  meeting it.
- GC content is compared against the threshold unrounded and displayed to two
  decimal places, so a sequence at 42.857% does not pass a minimum of 42.859%
  by virtue of how it is printed.
- Columns in a table input are read by position. A trailing separator is
  ignored, but an empty sequence column is an error rather than a licence to
  read the next column instead.
- A self-dimer and a self-complementary duplex are made of chemically identical
  strands, which removes the factor of four from the concentration term; a
  cross-dimer between two different oligonucleotides keeps it.
- In a structure diagram `|` marks a scored base pair and `.` a complementary
  position whose helix was too short or too weak to count.
- Three windows are called "the 3' end" and they differ. The GC clamp and the
  stacking free energy are counted over the last five bases, which is about how
  firmly a primer sits down. A dimer is screened as `3' overhang` when a scored
  pair falls within three bases of either terminus, and as `3' terminus paired`
  only when a terminal base is itself in a helix. Both are held to the stricter
  dimer threshold, because both are worth seeing, but only the second can be
  extended as it stands: an overhang needs the terminal bases to pair first.
- A dimer is judged register by register, not by its strongest register alone.
  A weaker register that reaches a 3' end is held to the stricter threshold, so
  it can breach while the most stable register of all passes; the flag names
  whichever register actually breached, and the summary column still reports
  the most stable one.
- `pass`, `note` and `review` summarise the flags. The thresholds are design
  conventions for a standard amplification primer, not physical constants: a
  probe, a sequencing primer or a cloning adapter will fail some of them by
  intent. Every threshold is a flag on the command line.

## Limitations

- Hairpins are single contiguous stems and dimers are single registers. Bulges,
  internal loops, multibranch structures and coaxial stacking are not modelled,
  and no terminal-mismatch stabilisation is credited to a loop. These free
  energies rank candidates; they are not a folding prediction, and a hairpin
  this screen reports as marginal may be real.
- Only Watson-Crick pairs are scored. Wobble and mismatch parameters, dangling
  ends, modified bases, locked nucleic acids and non-natural backbones are all
  outside the model, so a probe carrying a dye, a quencher or a minor-groove
  binder is measured as though it were bare DNA.
- RNA input is accepted and read as DNA. RNA/RNA and RNA/DNA duplexes have
  their own stacking parameters, which are not bundled.
- The extinction coefficient is for the single strand. An annealed duplex
  absorbs less than the sum of its strands, so these factors must not be used
  to quantify annealed material.
- The ion corrections were fitted over roughly 69 mM to 1.02 M monovalent
  cation and 0.5 to 125 mM magnesium. They are applied outside those ranges,
  which is what every implementation does, but a very dilute or very
  concentrated buffer is an extrapolation.
- A dimer register is judged on Watson-Crick contacts alone, so "within three
  bases of a 3' end" is a screening proxy for extendability rather than a
  measurement of it. Whether a polymerase actually extends a given register
  depends on the mismatch identity and the enzyme, neither of which is modelled.
- No template is considered. Specificity, amplicon secondary structure and
  mispriming against a genome need a search tool, not this one.

## Development

```bash
pip install -e ".[dev]"
pytest -q
python -m pyflakes src/oligo_stats tests
```

## License

Licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0).
See [LICENSE](LICENSE) for the complete terms.
