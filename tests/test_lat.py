# pragma pylint: disable=C0114,C0116
import pytest  # pylint: disable=W0611
from random import seed  # pylint: disable=C0411
from fractions import Fraction  # pylint: disable=C0411
from kryptools import Matrix, gram_schmidt, random_unimodular_matrix
from kryptools import Zmod, sis_search, sis_lll, isis_search, isis_lll


def test_gram_schmidt():
    V = Matrix([[5, 8], [0, 1]], ring=Fraction)
    Vs, M = gram_schmidt(V)
    assert V == Vs * M


def test_random_unimodular_matrix():
    for i in range(2, 5):
        U = random_unimodular_matrix(i)
        U.map(Fraction)
        assert U.det() in (-1, 1)

def test_sis():
    seed(0)
    for ring in Zmod(3), Zmod(4), Zmod(5), Zmod(6):
        for _ in range(25):
            A = Matrix([[ring.random() for _ in range(6)] for _ in range(3)])
            xs = sis_search(A)
            assert not A * xs
            x = sis_lll(A)
            assert not A * x
            assert x.norm() == xs.norm()
            b = Matrix([[ring.random()] for _ in range(3)])
            xs = isis_search(A, b)
            if xs is not None:
                assert A * xs == b
            x = isis_lll(A, b)
            if x is not None:
                assert A * x == b
            if xs is not None:
                assert xs.norm() <= x.norm() <= 1.5 * xs.norm()
