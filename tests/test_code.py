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

def test_Goppa():
    gf = GF2(4)
    g = Poly([2, 0, 0, 1], ring = gf)
    alpha = list(gf)

    goppa = Goppa(gf, g, alpha)
    for _ in range(10):
        x = [ randint(0,1) for _ in range(goppa.G.rows) ]
        y = goppa.encode(x)
        for _ in range(3):
            i = randint(0, goppa.G.cols-1)
            y[i] = 1- y[i]
    assert goppa.decode(y) == x
