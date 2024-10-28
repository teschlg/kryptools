import pytest
from fractions import Fraction
from kryptools import Matrix, gram_schmidt


def test_gram_schmidt():
	V = Matrix([[5, 8], [0, 1]], ring=Fraction)
	Vs, M = gram_schmidt(V)
	assert V == Vs * M
