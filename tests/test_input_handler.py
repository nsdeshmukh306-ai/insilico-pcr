"""Unit tests for input_handler module."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from modules.input_handler import (
    validate_primer,
    parse_primers_from_strings,
    PCRParams,
    default_params,
)


class TestValidatePrimer:
    def test_valid_primer(self):
        assert validate_primer("ATCGATCGATCGATCGATCG") == "ATCGATCGATCGATCGATCG"

    def test_lowercase_converted(self):
        assert validate_primer("atcgatcg atcgatcg") == "ATCGATCGATCGATCG"

    def test_iupac_allowed(self):
        seq = validate_primer("ATCGRYSWKMBDHVN", allow_iupac=True)
        assert "R" in seq

    def test_iupac_rejected_strict(self):
        with pytest.raises(ValueError):
            validate_primer("ATCGR", allow_iupac=False)

    def test_invalid_base(self):
        with pytest.raises(ValueError):
            validate_primer("ATCGX")

    def test_too_short(self):
        with pytest.raises(ValueError):
            validate_primer("ATCG")

    def test_too_long_warning(self):
        with pytest.raises(ValueError):
            validate_primer("A" * 61)

    def test_rna_converted(self):
        # U should be converted to T
        result = validate_primer("AUCGAUCGAUCGAUCGAUCG")
        assert "U" not in result
        assert "T" in result


class TestParsePrimersFromStrings:
    def test_basic_parse(self):
        pair = parse_primers_from_strings("GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC")
        assert pair.forward == "GCTAGCTAGCTAGCTAGCTA"
        assert pair.reverse == "TAGCTAGCTAGCTAGCTAGC"
        assert pair.name == "primer_pair_1"

    def test_custom_name(self):
        pair = parse_primers_from_strings(
            "GCTAGCTAGCTAGCTAGCTA", "TAGCTAGCTAGCTAGCTAGC", name="myPair"
        )
        assert pair.name == "myPair"


class TestPCRParams:
    def test_defaults(self):
        p = default_params()
        assert p.max_mismatches == 3
        assert p.min_amplicon_size == 50
        assert p.max_amplicon_size == 3000
        assert p.three_prime_strict is True

    def test_custom(self):
        p = PCRParams(max_mismatches=2, min_amplicon_size=100)
        assert p.max_mismatches == 2
        assert p.min_amplicon_size == 100
