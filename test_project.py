import pytest
from project import (fetch_gene_summary, fetch_sequence_data, reverse_complement, calculate_cg_percentage,
                     calculate_melting_temps, count_nucleotides)


def test_fetch_gene_summary_true():
    accession = 'NM_003979'
    gene_summary_data, gene_id = fetch_gene_summary(accession)
    assert gene_summary_data is not None
    assert gene_id is not None


def test_fetch_gene_summary_false():
    accession = 'invalid_id'
    with pytest.raises(ValueError):
        fetch_gene_summary(accession)


def test_count_nucleotides():
    flanks = ['ACCGTTATATGGACGAATCGATCGATCGATCGA']
    expected_output = [{'A': 10, 'C': 7, 'G': 8, 'T': 8}]
    assert count_nucleotides(flanks) == expected_output


def test_fetch_sequence_data():
    gene_id = "1722698043"
    sequence = fetch_sequence_data(gene_id)
    assert sequence is not None


def test_reverse_complement():
    # Test with a valid sequence
    assert reverse_complement("ATGC") == "GCAT"

    # Test with an invalid sequence
    with pytest.raises(KeyError):
        reverse_complement("invalid_sequence")


def test_calculate_cg_percentage():
    # Test with a valid primer list
    assert calculate_cg_percentage(["ATGC", "CGTA"], []) == [50.0, 50.0]

    # Test with an empty primer list
    assert calculate_cg_percentage([], []) == []


def test_calculate_melting_temps():
    # Test with a valid primer list
    assert calculate_melting_temps(["ATGC", "CGTA"], []) == [12, 12]

    # Test with an empty primer list
    assert calculate_melting_temps([], []) == []


