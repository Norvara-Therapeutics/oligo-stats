"""FASTA and table readers."""

import pytest

from oligo_stats.inputs import (
    InputError,
    Oligo,
    parse,
    parse_fasta,
    read,
    unique_names,
)


def test_fasta_is_parsed_with_identifier_and_description():
    oligos = parse([">F1 forward primer", "ACGTACGT"])
    assert oligos == [Oligo("F1", "forward primer", "ACGTACGT")]


def test_fasta_tolerates_blank_lines_wrapping_and_a_byte_order_mark():
    oligos = parse(["\ufeff>F1", "ACGT", "", "ACGT", ">F2", "TTTT"])
    assert [item.sequence for item in oligos] == ["ACGTACGT", "TTTT"]


def test_fasta_rejects_an_empty_header():
    with pytest.raises(InputError, match="no identifier"):
        parse([">", "ACGT"])
    with pytest.raises(InputError, match="no identifier"):
        parse([">F1", "ACGT", ">", "ACGT"])


def test_fasta_rejects_sequence_before_any_header():
    with pytest.raises(InputError, match="before any"):
        list(parse_fasta(["ACGT", ">F1", "ACGT"]))


def test_fasta_rejects_text_that_is_not_sequence():
    with pytest.raises(InputError, match="not nucleotide sequence"):
        parse([">F1", "the quick brown fox"])


def test_a_tab_separated_table_is_parsed():
    oligos = parse(["F1\tACGTACGT", "R1\tTTTTAAAA"])
    assert [(item.name, item.sequence) for item in oligos] == [
        ("F1", "ACGTACGT"),
        ("R1", "TTTTAAAA"),
    ]


@pytest.mark.parametrize("separator", ["\t", ",", ";", " "])
def test_every_recognised_separator_works(separator):
    oligos = parse([f"F1{separator}ACGTACGT"])
    assert oligos[0].name == "F1"
    assert oligos[0].sequence == "ACGTACGT"


@pytest.mark.parametrize("header", ["sequence", "Seq", "OLIGO", "primer", "bases"])
def test_a_header_row_is_dropped_and_rows_are_renumbered(header):
    oligos = parse([f"name,{header}", "F1,ACGTACGT"])
    assert len(oligos) == 1
    assert oligos[0].name == "F1"


def test_a_mistyped_base_in_the_first_row_is_an_error_not_a_dropped_record():
    """Treating any unparseable first row as a header hides the worst typo."""
    with pytest.raises(InputError, match="line 1"):
        parse(["F1\tACGTXACGT", "R1\tACGTACGT"])
    with pytest.raises(InputError, match="line 1"):
        parse(["F1\tACGTXACGT"])


def test_bare_sequences_are_named_by_position():
    oligos = parse(["ACGTACGT", "TTTTAAAA"])
    assert [item.name for item in oligos] == ["oligo1", "oligo2"]


def test_comments_and_blank_lines_are_skipped():
    oligos = parse(["# a comment", "", "ACGTACGT"])
    assert [item.name for item in oligos] == ["oligo1"]


def test_an_empty_input_is_an_error():
    with pytest.raises(InputError, match="no oligonucleotides"):
        parse(["", "   "])


def test_an_empty_sequence_cell_is_an_error():
    """Dropping the hole would slide the sequence in from the next column."""
    with pytest.raises(InputError, match="sequence column is empty"):
        parse(["F1,,ACGTACGT"])
    with pytest.raises(InputError, match="sequence column is empty"):
        parse(["F1\t\tACGTACGT"])


def test_a_trailing_separator_is_ignored():
    assert parse(["F1,ACGTACGT,"])[0].sequence == "ACGTACGT"
    assert parse(["F1,ACGTACGT,"])[0].description == ""
    assert parse(["F1,ACGTACGT"])[0].sequence == "ACGTACGT"


def test_extra_columns_become_the_description():
    oligos = parse(["F1\tACGT\tforward\tprimer"])
    assert oligos[0].description == "forward primer"


def test_read_reports_a_missing_path_and_a_directory(tmp_path):
    with pytest.raises(FileNotFoundError):
        read(tmp_path / "absent.fasta")
    with pytest.raises(IsADirectoryError):
        read(tmp_path)


def test_read_accepts_a_path(tmp_path):
    path = tmp_path / "primers.fasta"
    path.write_text(">F1\nACGTACGT\n", encoding="utf-8")
    assert read(path)[0].name == "F1"


def test_unique_names_suffixes_only_repeats():
    oligos = [Oligo("F1", "", "A"), Oligo("F1", "", "C"), Oligo("R1", "", "G")]
    assert unique_names(oligos) == ["F1", "F1.2", "R1"]


def test_unique_names_steps_past_a_collision_with_an_existing_name():
    oligos = [Oligo("F1", "", "A"), Oligo("F1", "", "C"), Oligo("F1.2", "", "G")]
    assert unique_names(oligos) == ["F1", "F1.3", "F1.2"]
