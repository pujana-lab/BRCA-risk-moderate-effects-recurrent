import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from common import allele_key, is_indel_alleles, is_palindromic


def test_allele_key_order_independent():
    assert allele_key("A", "G") == allele_key("G", "A")


def test_palindromic():
    assert is_palindromic("A", "T")
    assert is_palindromic("C", "G")
    assert not is_palindromic("A", "G")


def test_indel():
    assert is_indel_alleles("A", "AT")
    assert not is_indel_alleles("A", "G")
