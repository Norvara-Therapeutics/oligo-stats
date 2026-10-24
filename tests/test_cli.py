"""The command-line interface."""

import io
import json
import os

import pytest

from oligo_stats.cli import main

CLEAN = "TCATGAAGGATGGTGTTAATGC"
PROBLEM = "CATGCCTCCAAATCTTGGAGGC"
#: The defaults are already an amplification buffer; this raises the magnesium
#: so a test can show the conditions reaching the output.
PCR = ["--mg", "3"]
#: What a supplier's specification sheet is quoted against.
SUPPLIER = ["--na", "50", "--k", "0", "--tris", "0", "--mg", "0", "--dntp", "0"]


def _run(argv):
    out = io.StringIO()
    code = main(argv, stdout=out)
    return code, out.getvalue()


def test_version_is_reported(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert "oligo-stats" in capsys.readouterr().out


def test_a_sequence_on_the_command_line_is_measured():
    code, out = _run(["analyse", CLEAN])
    assert code == 0
    assert "22 nt" in out
    assert "pass" in out


def test_the_american_spelling_is_accepted_as_an_alias():
    assert _run(["analyze", CLEAN])[1] == _run(["analyse", CLEAN])[1]


def test_several_sequences_are_numbered_and_screened_as_a_set():
    code, out = _run(["analyse", CLEAN, PROBLEM])
    assert code == 0
    assert "> oligo1" in out and "> oligo2" in out
    assert "1 pair screened for cross-dimers" in out


def test_reaction_conditions_change_the_melting_temperature():
    """An amplification buffer melts an oligonucleotide well above 50 mM Na+."""
    amplification = _run(["analyse", CLEAN, "-f", "json"])[1]
    supplier = _run(["analyse", CLEAN, "-f", "json"] + SUPPLIER)[1]
    hotter = json.loads(amplification)["oligos"][0]["thermodynamics"]
    cooler = json.loads(supplier)["oligos"][0]["thermodynamics"]
    assert hotter["melting_temperature"] > cooler["melting_temperature"] + 3
    assert hotter["free_energy"] == cooler["free_energy"]


def test_tsv_output_has_a_header_and_one_row_per_oligonucleotide():
    code, out = _run(["analyse", CLEAN, PROBLEM, "-f", "tsv"])
    lines = out.strip().splitlines()
    assert code == 0
    assert len(lines) == 3
    assert lines[0].split("\t")[0] == "name"
    assert len(lines[1].split("\t")) == len(lines[0].split("\t"))


def test_json_output_is_valid_and_carries_the_conditions_and_structures():
    code, out = _run(["analyse", PROBLEM, "-f", "json"] + PCR)
    document = json.loads(out)
    assert code == 0
    assert document["conditions"]["regime"] == "mixed"
    assert document["conditions"]["free_magnesium_molar"] == pytest.approx(2.2e-3)
    oligo = document["oligos"][0]
    assert oligo["sequence"] == PROBLEM
    assert oligo["hairpins"][0]["stem_length"] == 8
    assert "diagram" in oligo["self_dimers"][0]
    assert oligo["severity"] == "warn"


def test_the_none_format_writes_nothing_to_standard_output():
    assert _run(["analyse", CLEAN, "-f", "none"]) == (0, "")


def test_an_html_report_and_a_tsv_are_written(tmp_path):
    report = tmp_path / "report.html"
    table = tmp_path / "table.tsv"
    code, _ = _run(
        ["analyse", CLEAN, "-f", "none", "-r", str(report), "-o", str(table)]
    )
    assert code == 0
    assert report.read_text().startswith("<!DOCTYPE html>")
    assert table.read_text().startswith("name\t")


def test_input_is_read_from_a_file(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(f">F1 forward\n{CLEAN}\n", encoding="utf-8")
    code, out = _run(["analyse", "-i", str(path)])
    assert code == 0
    assert "> F1  forward" in out


def test_input_is_read_from_standard_input(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO(f">F1\n{CLEAN}\n"))
    code, out = _run(["analyse", "-i", "-"])
    assert code == 0
    assert "> F1" in out


def test_positional_sequences_and_a_file_are_combined(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(f">F1\n{CLEAN}\n", encoding="utf-8")
    code, out = _run(["analyse", PROBLEM, "-i", str(path), "-f", "tsv"])
    assert code == 0
    assert len(out.strip().splitlines()) == 3


def test_no_input_at_all_is_an_error(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse"])
    assert exit_info.value.code == 2
    assert "no oligonucleotide given" in capsys.readouterr().err


def test_an_ambiguity_code_is_refused_with_a_pointer_at_expand(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", "ACGTNACGTACGTACGTACG"])
    assert exit_info.value.code == 2
    assert "expand" in capsys.readouterr().err


def test_a_missing_input_file_is_reported(capsys, tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", "-i", str(tmp_path / "absent.fasta")])
    assert exit_info.value.code == 2
    assert "no such file" in capsys.readouterr().err


def test_writing_over_the_input_is_refused(capsys, tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(f">F1\n{CLEAN}\n", encoding="utf-8")
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", "-i", str(path), "-o", str(path)])
    assert exit_info.value.code == 2
    assert "overwrite the input" in capsys.readouterr().err


def test_an_output_hard_linked_to_the_input_is_refused(capsys, tmp_path):
    """Resolved names differ for a hard link, but writing truncates the file."""
    source = tmp_path / "primers.fasta"
    source.write_text(f">F1\n{CLEAN}\n", encoding="utf-8")
    alias = tmp_path / "out.tsv"
    os.link(source, alias)
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", "-i", str(source), "-o", str(alias)])
    assert exit_info.value.code == 2
    assert "overwrite the input" in capsys.readouterr().err
    assert source.read_text() == f">F1\n{CLEAN}\n"


def test_two_outputs_hard_linked_together_are_refused(capsys, tmp_path):
    first = tmp_path / "a.tsv"
    first.write_text("", encoding="utf-8")
    second = tmp_path / "b.html"
    os.link(first, second)
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "-o", str(first), "-r", str(second)])
    assert exit_info.value.code == 2
    assert "would both write to" in capsys.readouterr().err


def test_an_oligonucleotide_with_no_melting_temperature_does_not_pass_strict():
    """A window that could not be applied is not a window that was met."""
    code, out = _run(["analyse", CLEAN, "--conc", "1e-25", "--strict"])
    assert "Tm n/a C" in out
    assert "no melting temperature above 0 C" in out
    assert code == 1


def test_an_empty_sequence_column_is_an_error_not_a_column_shift(
    capsys, monkeypatch
):
    monkeypatch.setattr("sys.stdin", io.StringIO(f"F1,,{CLEAN}\n"))
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", "-i", "-"])
    assert exit_info.value.code == 2
    assert "sequence column is empty" in capsys.readouterr().err


def test_two_outputs_may_not_share_one_path(capsys, tmp_path):
    shared = tmp_path / "out"
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "-o", str(shared), "-r", str(shared)])
    assert exit_info.value.code == 2
    assert "would both write to" in capsys.readouterr().err


def test_strict_mode_exits_non_zero_only_when_something_is_flagged():
    assert _run(["analyse", CLEAN, "--strict", "-f", "none"])[0] == 0
    assert _run(["analyse", PROBLEM, "--strict", "-f", "none"])[0] == 1


#: A primer and its own reverse complement: each is individually clean, and
#: together they are the worst cross-dimer there is.
REVERSE_COMPLEMENT = "GCATTAACACCATCCTTCATGA"


def test_strict_mode_fails_on_a_cross_dimer_between_individually_clean_oligos():
    """An ordering gate that passes a primer and its complement is no gate."""
    code, out = _run(["analyse", CLEAN, REVERSE_COMPLEMENT, "--strict", "-f", "table"])
    assert "pass" in out
    assert "1 at or below threshold" in out
    assert code == 1


def test_strict_mode_still_passes_a_panel_with_no_finding():
    assert _run(["analyse", CLEAN, "--strict", "-f", "none"])[0] == 0
    assert _run(
        ["analyse", CLEAN, REVERSE_COMPLEMENT, "--strict", "--no-cross-dimers", "-f", "none"]
    )[0] == 0


@pytest.mark.parametrize(
    "option",
    [
        ["--na", "nan"],
        ["--na", "inf"],
        # A leading minus must be attached, or argparse claims the token as an
        # option before any type conversion runs.
        ["--mg=-inf"],
        ["--conc", "nan"],
        ["--dntp", "inf"],
        ["--tm-range", "nan:inf"],
        ["--gc-range", "40:inf"],
        ["--hairpin-threshold=-inf"],
        ["--end-stability=nan"],
        ["--dimer-threshold=-inf"],
    ],
)
def test_a_non_finite_number_is_refused(capsys, option):
    """NaN does not propagate visibly: it routes the ion regime to its last
    branch and a confident melting temperature comes out of a buffer that was
    never given."""
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN] + option)
    assert exit_info.value.code == 2
    assert "not a finite number" in capsys.readouterr().err


def test_a_negative_concentration_is_refused(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "--mg=-5"])
    assert exit_info.value.code == 2
    assert "must not be negative" in capsys.readouterr().err


def test_a_json_document_is_always_strict_json():
    code, out = _run(["analyse", CLEAN, "-f", "json"])
    assert code == 0

    def reject(constant):
        raise AssertionError(f"non-finite constant in JSON: {constant}")

    json.loads(out, parse_constant=reject)


def test_a_tab_in_a_description_does_not_add_a_tsv_column(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(f">F1 desc\twith tab\n{CLEAN}\n", encoding="utf-8")
    code, out = _run(["analyse", "-i", str(path), "-f", "tsv"])
    rows = [line.split("\t") for line in out.strip().splitlines()]
    assert code == 0
    assert len(rows[1]) == len(rows[0])
    assert dict(zip(rows[0], rows[1]))["description"] == "desc with tab"


def test_every_tsv_row_has_exactly_the_header_columns(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(
        f">F1 a\tb\tc\n{CLEAN}\n>F2  spaced   out \n{PROBLEM}\n", encoding="utf-8"
    )
    out = _run(["analyse", "-i", str(path), "-f", "tsv"])[1]
    rows = [line.split("\t") for line in out.strip().splitlines()]
    assert len({len(row) for row in rows}) == 1


@pytest.mark.parametrize(
    "option", [["--length-range", "18.9:30.1"], ["--clamp-range", "0.5:3.5"]]
)
def test_a_fractional_count_range_is_refused_rather_than_truncated(capsys, option):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN] + option)
    assert exit_info.value.code == 2
    assert "whole number of bases" in capsys.readouterr().err


def test_a_whole_number_count_range_is_still_accepted():
    document = json.loads(
        _run(
            ["analyse", CLEAN, "--length-range", "20:24", "--clamp-range", "0:5",
             "-f", "json"]
        )[1]
    )
    assert document["criteria"]["length"] == [20, 24]
    assert document["criteria"]["gc_clamp"] == [0, 5]


def test_json_separates_a_paired_terminus_from_a_three_prime_overhang():
    document = json.loads(
        _run(["analyse", "GGGCCCATGCTGGCATCTGCGA", "-f", "json"])[1]
    )
    registers = document["oligos"][0]["self_dimers"]
    proximal = [item for item in registers if item["three_prime_paired"]]
    assert proximal
    assert proximal[0]["terminal_paired"] is False


def test_thresholds_are_configurable():
    relaxed = _run(["analyse", PROBLEM, "--hairpin-threshold", "-20", "-f", "json"])[1]
    assert "hairpin" not in [
        flag["code"] for flag in json.loads(relaxed)["oligos"][0]["flags"]
    ]


def test_a_positive_free_energy_threshold_is_refused(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "--hairpin-threshold", "2"])
    assert exit_info.value.code == 2
    assert "zero or negative" in capsys.readouterr().err


def test_an_inverted_range_is_refused(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "--tm-range", "65:55"])
    assert exit_info.value.code == 2
    assert "inverted" in capsys.readouterr().err


def test_a_buffer_with_no_cation_is_refused(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["analyse", CLEAN, "--na", "0", "--k", "0", "--tris", "0", "--mg", "0"])
    assert exit_info.value.code == 2
    assert "no cation" in capsys.readouterr().err


def test_the_structure_screen_can_be_turned_off():
    code, out = _run(["analyse", PROBLEM, "--no-structures"])
    assert code == 0
    assert "hairpin" not in out


def test_the_cross_dimer_screen_can_be_turned_off():
    code, out = _run(["analyse", CLEAN, PROBLEM, "--no-cross-dimers"])
    assert code == 0
    assert "cross-dimers" not in out


def test_a_five_prime_phosphate_changes_the_reported_mass():
    hydroxyl = _run(["analyse", CLEAN, "-f", "json"])[1]
    phosphate = _run(["analyse", CLEAN, "--five-prime", "phosphate", "-f", "json"])[1]
    first = json.loads(hydroxyl)["oligos"][0]["optics"]["molecular_weight"]
    second = json.loads(phosphate)["oligos"][0]["optics"]["molecular_weight"]
    assert second == pytest.approx(first + 79.98)


def test_an_unmeasurable_value_is_named_rather_than_left_blank():
    """A one-base oligonucleotide has no stack, so it has no Tm to report."""
    code, out = _run(["analyse", "A"])
    assert code == 0
    assert "Tm n/a C" in out
    assert "3' stacking n/a kcal/mol" in out


def test_tsv_leaves_an_unmeasurable_value_empty_for_a_reader_of_columns():
    code, out = _run(["analyse", "A", "-f", "tsv"])
    header, row = out.strip().splitlines()
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert code == 0
    assert fields["tm_c"] == ""
    assert fields["three_prime_stacking"] == ""


def test_tsv_carries_the_description_like_the_other_formats(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(f">F1 forward primer\n{CLEAN}\n", encoding="utf-8")
    code, out = _run(["analyse", "-i", str(path), "-f", "tsv"])
    header, row = out.strip().splitlines()
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert code == 0
    assert fields["description"] == "forward primer"


def test_json_records_every_threshold_that_was_applied():
    """The JSON has to say which rules produced its flags, not just some."""
    import dataclasses

    from oligo_stats.oligo import Criteria

    document = json.loads(_run(["analyse", CLEAN, "-f", "json"])[1])
    recorded = document["criteria"]
    assert set(recorded) == {
        "length",
        "gc_percent",
        "melting_temperature",
        "gc_clamp",
        "end_stability",
        "max_run",
        "hairpin_free_energy",
        "dimer_free_energy",
        "dimer_three_prime_free_energy",
    }
    # Every Criteria field is represented, counting the paired ranges once.
    paired = sum(len(value) for value in recorded.values() if isinstance(value, list))
    singles = sum(1 for value in recorded.values() if not isinstance(value, list))
    assert paired + singles == len(dataclasses.fields(Criteria))


def test_json_distinguishes_an_unscreened_oligonucleotide_from_a_clean_one():
    screened = json.loads(_run(["analyse", CLEAN, "-f", "json"])[1])["oligos"][0]
    skipped = json.loads(
        _run(["analyse", CLEAN, "--no-structures", "-f", "json"])[1]
    )["oligos"][0]
    assert screened["structures_screened"] is True
    assert skipped["structures_screened"] is False
    assert skipped["hairpins"] == []


def test_the_clamp_end_stability_and_run_thresholds_are_configurable():
    strict = json.loads(
        _run(["analyse", CLEAN, "--clamp-range", "4:5", "-f", "json"])[1]
    )["oligos"][0]
    assert "gc_clamp" in [flag["code"] for flag in strict["flags"]]
    lax = json.loads(
        _run(["analyse", "GAGGCGGTCAAAAAGCCGCCTC", "--max-run", "6", "-f", "json"])[1]
    )["oligos"][0]
    assert "base_run" not in [flag["code"] for flag in lax["flags"]]
    ends = json.loads(
        _run(["analyse", CLEAN, "--end-stability", "-3", "-f", "json"])[1]
    )["oligos"][0]
    assert "end_stability" in [flag["code"] for flag in ends["flags"]]


def test_the_panel_melting_span_is_reported_for_more_than_one_oligonucleotide():
    one = _run(["analyse", CLEAN])[1]
    many = _run(["analyse", CLEAN, "AATATCTATAACAACTATTTTA"])[1]
    assert "Tm across the panel" not in one
    assert "Tm across the panel" in many


def test_expand_writes_every_variant_as_fasta():
    code, out = _run(["expand", "ACGTRY"])
    assert code == 0
    assert out.count(">") == 4
    assert ">oligo1.1 of 4" in out
    assert "ACGTAC" in out


def test_expand_leaves_an_unambiguous_sequence_alone():
    code, out = _run(["expand", "ACGT"])
    assert code == 0
    assert out == ">oligo1\nACGT\n"


def test_expand_refuses_an_expansion_above_its_limit(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["expand", "NNNNN", "--limit", "10"])
    assert exit_info.value.code == 2
    assert "above the limit" in capsys.readouterr().err


def test_expand_refuses_to_overwrite_its_own_input(capsys, tmp_path):
    """expand reads the file and then writes it, so sharing a path loses it."""
    path = tmp_path / "degenerate.fasta"
    path.write_text(">F1\nACGTRCCTGATCAAGCTTAGG\n", encoding="utf-8")
    with pytest.raises(SystemExit) as exit_info:
        main(["expand", "-i", str(path), "-o", str(path)])
    assert exit_info.value.code == 2
    assert "overwrite the input" in capsys.readouterr().err
    assert path.read_text() == ">F1\nACGTRCCTGATCAAGCTTAGG\n"


def test_expand_writes_to_a_file(tmp_path):
    path = tmp_path / "expanded.fasta"
    code, out = _run(["expand", "ACGTR", "-o", str(path)])
    assert code == 0
    assert out == ""
    assert path.read_text().count(">") == 2


def test_expand_output_feeds_straight_back_into_analyse(monkeypatch):
    expanded = _run(["expand", "ACGTRCCTGATCAAGCTTAGG"])[1]
    monkeypatch.setattr("sys.stdin", io.StringIO(expanded))
    code, out = _run(["analyse", "-i", "-", "-f", "tsv"])
    assert code == 0
    assert len(out.strip().splitlines()) == 3
