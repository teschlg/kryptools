# pragma pylint: disable=C0114,C0116
import pytest  # pylint: disable=W0611
from fractions import Fraction  # pylint: disable=C0411
from kryptools import Matrix, Zmod, GF2


def test_Matrix():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 12]])
    Z = M.zeros()
    assert M
    assert M == M + Z
    assert not Z
    Z[0] = 1
    assert Z
    assert Z != Z.zeros()
    assert len(M) == 3 * 3
    assert M[0] == 1
    assert M[:] == [1, 2, 3, 4, 5, 6, 7, 8, 12]
    assert 2 * M - M == M
    M.map(Fraction)
    Mi = M.inv()
    assert M * Mi == M.eye()

    Z_11 = Zmod(11)
    assert 2 * M - M == M
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 12]], ring=Z_11)
    Mi = M.inv()
    assert M * Mi == M.eye()

def test_kernel():
    for gf in Zmod(2), Zmod(7), GF2(1), GF2(8):
        for m,n in ((3, 5), (5, 3), (3, 3)):
            for _ in range(25):
                A = Matrix([[gf.random() for _ in range(m)] for _ in range(n)], ring = gf)
                K = A.kernel()
                assert not A * K
                assert K.rows == A.cols
                if K:
                    assert K.cols == A.nullity()
                else:
                    assert K.cols == 1
