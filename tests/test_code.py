import pytest
from random import randint, seed
from kryptools import Zmod, GF2, Poly, Matrix
from kryptools import gen2pchk, Goppa

seed(0)
Z_2 = Zmod(2)
G_Hamming = Matrix([
[1, 0, 0, 0, 1, 1, 0],
[0, 1, 0, 0, 1, 0, 1],
[0, 0, 1, 0, 0, 1, 1],
[0, 0, 0, 1, 1, 1, 1]
], ring = Z_2)


def test_gen2pchk():
    H = gen2pchk(G_Hamming)
    assert not H * G_Hamming.transpose()

num_tests = 10

def test_Goppa():
    for gf, t in [ [GF2(4), 3],  [GF2(6), 4]]:
        found = False
        while not found:
            c = [randint(0,gf.order-1) for _ in range(t)] + [1]
            g = Poly(c, ring = gf)
            found = g.rabin_test()
        alpha = list(gf)

        goppa = Goppa(gf, g, alpha)
        for _ in range(num_tests):
            x = [ randint(0,1) for _ in range(goppa.G.rows) ]
            y = goppa.encode(x)
            for _ in range(t):
                i = randint(0, goppa.G.cols-1)
                y[i] = 1- y[i]
        assert goppa.decode(y) == x
