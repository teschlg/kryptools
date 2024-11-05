import pytest
from random import randint, seed
from kryptools import GF2, Poly, Zmod
seed(0)

Z_2 = Zmod(2)
aes = [1, 1, 0, 1, 1, 0, 0, 0, 1]  # x^8 + x^4 + x^3 + x + 1

def PolyAES(c: list) -> "Poly":
    return Poly(c, ring=Z_2, modulus=aes)

def byte2poly(x: int, n: int = 8) -> "Poly":
    "Convert an integer to a polynomial."
    return PolyAES(list(reversed([int(d) for d in str(format(x, "0" + str(n) + "b"))])))

def poly2byte(p: "Poly") -> int:
    "Convert a polynomial to the corresponding integer."
    s = 0
    for i in reversed(p.coeff):
        s = s * 2 + int(i)
    return s


def PolyAES(c: list) -> "Poly":
    return Poly(c, ring=Z_2, modulus=aes)


def test_GF2_ops():
	gf = GF2(8)
	for _ in range(100):
		a = randint(0, gf.order)
		b = randint(0, gf.order)
		assert int(gf(a) + gf(b)) == poly2byte(byte2poly(a) + byte2poly(b))
		assert int(gf(a) - gf(b)) == poly2byte(byte2poly(a) - byte2poly(b))
		assert int(gf(a) * gf(b)) == poly2byte(byte2poly(a) * byte2poly(b))
		if not(b):
			continue
		assert int(gf(a) / gf(b)) == poly2byte(byte2poly(a) / byte2poly(b))
			
	with pytest.raises(ValueError):
		gf(1) / gf(0)
	with pytest.raises(ValueError):
		gf(0)**-1
	for a in range(1, gf.order):
		assert int(gf(a)**-1) == poly2byte(byte2poly(a).inv())
		

