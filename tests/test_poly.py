import pytest
from kryptools import Poly, Zmod


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
