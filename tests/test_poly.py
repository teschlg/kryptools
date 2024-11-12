import pytest
from kryptools import Poly, Zmod, GF2


def test_Poly():
    p = Poly([2, 0, 2])
    q = Poly([1, 2, 3, 4])
    d, m = q.divmod(p)
    assert d * p + m == q

    Z_5 = Zmod(5)
    p.map(Z_5)
    q.map(Z_5)
    d, m = q.divmod(p)
    assert d * p + m == q

    gf = GF2(8, modulus=0b100011011)
    p = Poly([2, 1, 1, 3], modulus=[1, 0, 0, 0, 1], ring=gf)
    q = p.inv()
    assert p * q == Poly([1], modulus=[1, 0, 0, 0, 1], ring=gf)
