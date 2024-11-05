import pytest
from random import randint, seed
from kryptools import EC_Weierstrass
seed(0)


def test_EC():
	ec = EC_Weierstrass(239, 3, 1)
	O = ec(None, None)  # point at infinity
	assert O in ec
	P = ec.random()
	assert P in ec
	Q = ec.random()
	assert Q in ec
	assert P + O == P
	assert O + P == P
	assert P - P == O
	assert P + P == 2 * P
	R = O
	for i in range(5):
		assert i * P == R
		R += P
	assert P + Q == Q + P
	assert Q.order() * Q == O
	k = randint(1, Q.order())
	R = k * Q
	assert R.dlog(Q) == k
