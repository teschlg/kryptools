import pytest
from random import randint, seed
from kryptools import GF2, Poly, Zmod
seed(0)

Z_2 = Zmod(2)


def byte2poly(x: int, n: int, poly: list) -> "Poly":
    "Convert an integer to a polynomial."
    return Poly(list(reversed([int(d) for d in str(format(x, "0" + str(n) + "b"))])), ring=Z_2, modulus=poly)


def poly2byte(p: "Poly") -> int:
    "Convert a polynomial to the corresponding integer."
    s = 0
    for i in reversed(p.coeff):
        s = s * 2 + int(i)
    return s


def test_GF2_ops():
    for n, poly in [(4, None), (8, None), (8, 0b100011011), (12, None)]:
        gf = GF2(n, poly=poly)
        poly = [int(d) for d in reversed(bin(gf.poly)[2:])]
        for _ in range(100):
            a = randint(0, gf.order-1)
            b = randint(0, gf.order-1)
            assert int(gf(a) + gf(b)) == poly2byte(byte2poly(a, n, poly) + byte2poly(b, n, poly))
            assert int(gf(a) - gf(b)) == poly2byte(byte2poly(a, n, poly) - byte2poly(b, n, poly))
            assert int(gf(a) * gf(b)) == poly2byte(byte2poly(a, n, poly) * byte2poly(b, n, poly))
            if not (b):
                continue
            assert int(gf(a) / gf(b)) == poly2byte(byte2poly(a, n, poly) / byte2poly(b, n, poly))
        with pytest.raises(ValueError):
            gf(1) / gf(0)
        with pytest.raises(ValueError):
            gf(0)**-1
        for a in range(1, gf.order):
            assert int(gf(a)**-1) == poly2byte(byte2poly(a, n, poly).inv())
