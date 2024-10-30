import pytest
from math import gcd
from kryptools import Zmod, euler_phi

def test_Zmod_ops():
	for n in range(2, 10):
		Z_n = Zmod(n)
		for a in range(n):
			aa = Z_n(a)
			bool(aa) == bool(a)
			for i in range(4):
				assert aa ** i == Z_n(pow(a, i, n))
				if i == 0 or gcd(a, n) == 1:
					assert aa ** -i == Z_n(pow(a, -i, n))
				else:
					with pytest.raises(ValueError):
						aa ** -i
			if gcd(a, n) == 1:
				assert aa ** aa.order() == Z_n(1)
			else:
				with pytest.raises(ValueError):
					aa.order()
			for b in range(n):
				bb = Z_n(b)
				assert aa + bb == Z_n(a + b)
				assert aa - bb == Z_n(a - b)
				assert aa * bb == Z_n(a * b)
				if gcd(b, n) == 1:
					assert aa / bb == Z_n( a * pow(b, -1, n) )

def test_Zmod_methods():
	Z_5 = Zmod(5)
	assert Z_5(3).order() == 4
	assert Z_5(1).is_generator() == False
	assert Z_5(3).is_generator() == True
	assert Z_5.generator(all = True) == Z_5([2, 3])
	assert Z_5.star() == Z_5([1, 2, 3, 4])
	assert [Z_5(i).sharp() for i in range(Z_5.n)] == [0, 1, 2, -2, -1]
	assert [abs(Z_5(i)) for i in range(Z_5.n)] == [0, 1, 2, 2, 1]
	assert str(Z_5(3)) == "3"
	Z_5.short = False
	assert str(Z_5(6)) == "1 (mod 5)"

	Z_6 = Zmod(6)
	assert Z_6(5).order() == 2
	with pytest.raises(ValueError):
		Z_6(3).order()
	assert Z_6(1).is_generator() == False
	assert Z_6(5).is_generator() == True
	assert Z_6.generator(all = True) == Z_6([5])
	assert Z_6.star() == Z_6([1, 5])
	assert str(Z_6(3)) == "3"

def test_Zmod_order():
	for n in range(2,10):
		Z_n = Zmod(n)
		assert Z_n.order() == euler_phi(n)