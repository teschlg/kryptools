# pragma pylint: disable=C0114,C0116
import pytest  # pylint: disable=W0611
from fractions import Fraction  # pylint: disable=C0411
from random import randint, seed  # pylint: disable=C0411
from kryptools import Matrix, Zmod, GF2, eye, circulant, BinaryMatrix


def test_Matrix():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 12]])
    Z = M.zeros()
    assert M
    assert M == M + Z
    assert not Z
    Z[0] = 1
    assert Z
    assert Z != Z.zeros()
    assert circulant([ 1, 0, 0]) == eye(3)
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

def test_BinaryMatrix():
    num_tests = 100
    Z_2 = Zmod(2)
    seed(0)

    def matrix2bitmatrix(M: Matrix) -> list:
        return M.applyfunc(int).matrix

    def Matrix2(M: list) -> Matrix:
        return Matrix(M.bitmatrix(), ring = Z_2)

    n = 5

    for _ in range(num_tests):
        m = randint(n-1,n+1)
        mat1 = BinaryMatrix([randint(0,2**n-1) for _ in range(m)], cols = n)
        mat2 = BinaryMatrix([randint(0,2**n-1) for _ in range(m)], cols = n)
        mat3 = BinaryMatrix([randint(0,2**m-1) for _ in range(n)], cols = n)
        Mat1 = Matrix2(mat1)
        Mat2 = Matrix2(mat2)
        Mat3 = Matrix2(mat3)
        assert matrix2bitmatrix(Mat1 + Mat2) == (mat1 + mat2).bitmatrix()
        assert matrix2bitmatrix(Mat1 - Mat2) == (mat1 + mat2).bitmatrix()
        assert matrix2bitmatrix(Mat1 * Mat3) == (mat1 * mat3).bitmatrix()

    for _ in range(num_tests):
        m = randint(n-1,n+1)
        mat = BinaryMatrix([randint(0,2**n-1) for _ in range(m)], cols = n)
        b = [randint(0,1) for _ in range(m)]
        x = [randint(0,1) for _ in range(n)]
        Mat = Matrix2(mat)
        X = Matrix(x, ring = Z_2)
        assert matrix2bitmatrix(Mat.rref()) == mat.rref().bitmatrix()
        assert matrix2bitmatrix(Mat.kernel()) == mat.kernel().bitmatrix()
        assert matrix2bitmatrix(Mat.transpose()) == mat.transpose().bitmatrix()
        if m == n:
            d = mat.det()
            assert int(Mat.det()) == d
            if d:
                assert matrix2bitmatrix(Mat.inv()) == mat.inv().bitmatrix()
        assert Mat.rank() == mat.rank()
        sol = mat.solve(b)
        Sol = Mat.solve(b, ring = Z_2)
        assert (sol is None and Sol is None) or Sol.applyfunc(int).matrix == sol.bitmatrix()
        if mat.rows == 1:
            assert int(Mat * X) == mat.dot(x)[0]
        else:
            assert list((Mat * X).applyfunc(int)) == mat.dot(x)
