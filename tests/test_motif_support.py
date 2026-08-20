import numpy as np

from edge_causality.motif_support import pwm_from_pfm, scan_sequence


def test_pwm_has_expected_shape():
    pfm = {"A": [9, 1], "C": [1, 1], "G": [1, 1], "T": [1, 9]}
    assert pwm_from_pfm(pfm).shape == (4, 2)


def test_scan_finds_forward_and_reverse_hits():
    pfm = {"A": [100, 0], "C": [0, 0], "G": [0, 0], "T": [0, 100]}
    pwm = pwm_from_pfm(pfm)
    forward = scan_sequence("CCATCC", pwm, 0.90)
    reverse = scan_sequence("CCATCC".translate(str.maketrans("ACGT", "TGCA"))[::-1], pwm, 0.90)
    assert forward["motif_support"]
    assert reverse["motif_support"]


def test_scan_rejects_low_scoring_sequence():
    pfm = {"A": [100, 0], "C": [0, 0], "G": [0, 0], "T": [0, 100]}
    result = scan_sequence("CCCCCC", pwm_from_pfm(pfm), 0.90)
    assert not result["motif_support"]
