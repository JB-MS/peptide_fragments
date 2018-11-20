import pytest
import numpy as np
from pandas import DataFrame

from peptide_fragmentor import PeptideFragment0r


def test_fragment_one_aa_peptide_y_series():
    """Test y1 fragmentation"""
    fragger = PeptideFragment0r('K')
    fragments = fragger.fragment_peptide(ion_series=['y'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 1

    row = fragments.iloc[0]
    assert row['name'] == 'y1'
    assert row['cc'] == 'C(6)H(14)N(2)O(2)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 147.11280646


def test_fragment_two_aa_peptide_y_series():
    """Test y2 fragmentation"""
    fragger = PeptideFragment0r('KK')
    fragments = fragger.fragment_peptide(ion_series=['y'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 2

    row = fragments.iloc[1]
    print(row)
    assert row['name'] == 'y2'
    assert row['cc'] == 'C(12)H(26)N(4)O(3)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        274.2004907132
    )


def test_fragment_two_aa_peptide_b_series():
    """Test b2 fragmentation"""
    fragger = PeptideFragment0r('KK')
    fragments = fragger.fragment_peptide(ion_series=['b'])
    assert isinstance(fragments, DataFrame)

    row = fragments.iloc[1]
    assert row['name'] == 'b2'
    assert row['cc'] == 'C(12)H(24)N(4)O(2)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 257.19720249556997


def test_fragment_three_aa_peptide_b_series():
    """Test b3 fragmentation"""
    fragger = PeptideFragment0r('KKK')
    fragments = fragger.fragment_peptide(ion_series=['b'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 3
    row = fragments.iloc[2]
    assert row['name'] == 'b3'
    assert row['cc'] == 'C(18)H(36)N(6)O(3)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 385.29216550997


def test_fragment_two_aa_peptide_a_series():
    """Test a1 fragmentation"""
    fragger = PeptideFragment0r('KK')
    fragments = fragger.fragment_peptide(ion_series=['a'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 2
    row = fragments.iloc[0]
    assert row['name'] == 'a1'
    assert row['cc'] == 'C(5)H(12)N(2)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 101.10732486116999


def test_fragment_three_aa_peptide_a_series():
    """Test a2 fragmentation"""
    fragger = PeptideFragment0r('KKK')
    fragments = fragger.fragment_peptide(ion_series=['a'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 3
    row = fragments.iloc[1]
    assert row['name'] == 'a2'
    assert row['cc'] == 'C(11)H(24)N(4)O(1)'
    assert row['charge'] == 1
    assert pytest.approx(
        row['mz'],
        5e-6
    ) == 229.20228787557


def test_fragment_two_aa_peptide_y_series_with_mod():
    """Test y2 fragmentation"""
    fragger = PeptideFragment0r('MK#Oxidation:1')
    fragments = fragger.fragment_peptide(ion_series=['y'])
    assert isinstance(fragments, DataFrame)
    assert len(fragments) == 2

    row = fragments.iloc[1]
    print(row)
    assert row['name'] == 'y2'
    assert row['cc'] == 'C(11)H(23)N(3)O(4)S(1)'
    # assert row['charge'] == 1
    # assert pytest.approx(
    #     row['mz'],
    #     5e-6
    # ) == 275.207767179


# def test_fragment_one_aa_peptide_internal():
#     fragger = PeptideFragment0r('RKKR')
#     fragments = fragger.fragment_peptide(ion_series=['internal'])
#     assert isinstance(fragments, DataFrame)
#     assert len(fragments) == 1
