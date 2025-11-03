# pragma pylint: disable=C0114,C0116
import pytest  # pylint: disable=W0611
from random import seed, randint  # pylint: disable=C0411
from fractions import Fraction  # pylint: disable=C0411
from kryptools import Matrix, Lattice, gram_schmidt, random_unimodular_matrix
from kryptools import Zmod, sis_search, sis, isis_search, isis


def test_gram_schmidt():
    seed(0)
    for n in range(2, 5):
        for _ in range(10):
            found = False
            while not found:
                V = Matrix([[ randint(-10,10) for _ in range(n)] for _ in range(n) ], ring=Fraction)
                found = bool(V.rank() == n)
            Vs, M = gram_schmidt(V)
            assert V == Vs * M
            D = V.eye()
            for i in range(Vs.cols):
                assert Vs[:,i].dot(Vs[:,i]) == Vs[:,i].norm2()
                D[i,i] = Vs[:,i].dot(Vs[:,i])
            assert Vs.transpose() * Vs == D

def test_random_unimodular_matrix():
    for i in range(2, 5):
        U = random_unimodular_matrix(i)
        U.map(Fraction)
        assert U.det() in (-1, 1)

def check(self, j: int = 0) -> bool:
    "Check Gram-Schmit decomposition and weigts."
    if self.Nu is None:
        self.Nu = list(map(self.norm2, self.U))
    if not all( n > 0 for n in self.Nu):
        raise ValueError(f"Zero vector in basis {j}")
    if self.Mu is not None:
        oldMu = self.Mu
        oldNu = self.Nu
        self.Mu = None
        self.compute_weights()
        if self.Mu != oldMu:
            raise ValueError(f"Inconsistent Mu {j}")
        if self.Nu != oldNu:
            raise ValueError(f"Inconsistent Nu {j}")
    if self.Us is not None:
        self.gsd()
        oldUs = self.Us
        oldM = self.M
        oldN = self.N
        self.Us = None
        self.gsd()
        if self.Us != oldUs:
            raise ValueError(f"Inconsistent Us {j}")
        if self.M != oldM:
            raise ValueError(f"Inconsistent M {j}")
        if self.N != oldN:
            raise ValueError(f"Inconsistent N {j}")

Lattice.check = check

def test_Lattice():
    seed(0)
    for m, n in ((2,3), (3,3), (4,5)):
        for _ in range(100):
            found = False
            while not found:
                U = Matrix([[ randint(-10,10) for _ in range(m)] for _ in range(n) ])
                found = bool(U.rank() == m)
            k = [ randint(-10,10) for _ in range(m)]
            b = Matrix([ randint(-10,10) for _ in range(n)])
            lat = Lattice(U)
            lat2 = Lattice(U)
            lat2.lll()
            lat2.check()
            lat3 = Lattice(U)
            lat3.hermite()
            lat3.check()
            assert lat == lat2
            a = lat(k)
            assert a in lat
            assert a in lat2
            Us = Matrix(lat2.Us).transpose()
            M = Matrix(lat2.M).transpose()
            assert Us * M == lat2.basis()
            s1 = lat2.svp().norm2()
            s2 = lat2.svp(method = 'enum').norm2()
            assert s1 >= s2
            a0 = lat2.cvp(b, method = 'enum')
            assert a0 in lat
            n2 = (b-a0).norm2()
            for method in ('babai_round', 'babai_plane', 'kannan'):
                a = lat2.cvp(b, method = method)
                if a is not None:
                    assert a in lat
                    assert (b-a).norm2() >= n2


def test_sis():
    seed(0)
    for ring in Zmod(3), Zmod(4), Zmod(5), Zmod(6):
        for _ in range(25):
            A = Matrix([[ring.random() for _ in range(6)] for _ in range(3)])
            xs = sis_search(A)
            assert not A * xs
            x = sis(A)
            assert not A * x
            assert x.norm() == xs.norm()
            b = Matrix([[ring.random()] for _ in range(3)])
            xs = isis_search(A, b)
            if xs is not None:
                assert A * xs == b
            x = isis(A, b)
            if x is not None:
                assert A * x == b
            if xs is not None:
                assert xs.norm() <= x.norm() <= 1.5 * xs.norm()
