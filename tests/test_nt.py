import pytest
from math import gcd
from kryptools import carmichael_lambda, euler_phi, moebius_mu

#https://oeis.org/A000010
OEIS_A000010 = [
		1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8, 12, 10, 22, 8, 20,
		12, 18, 12, 28, 8, 30, 16, 20, 16, 24, 12, 36, 18, 24, 16, 40, 12, 42, 20, 24, 22,
		46, 16, 42, 20, 32, 24, 52, 18, 40, 24, 36, 28, 58, 16, 60, 30, 36, 32, 48, 20,
		66, 32, 44 ]

def test_euler_phi():
	for n in range(1, len(OEIS_A000010)+1):
		assert euler_phi(n) == OEIS_A000010[n - 1]

#https://oeis.org/A002322
OEIS_A002322 = [
		1, 1, 2, 2, 4, 2, 6, 2, 6, 4, 10, 2, 12, 6, 4, 4, 16, 6, 18, 4, 6, 10, 22, 2, 20,
		12, 18, 6, 28, 4, 30, 8, 10, 16, 12, 6, 36, 18, 12, 4, 40, 6, 42, 10, 12, 22, 46,
		4, 42, 20, 16, 12, 52, 18, 20, 6, 18, 28, 58, 4, 60, 30, 6, 16, 12, 10, 66, 16,
		22, 12, 70, 6, 72, 36, 20, 18, 30, 12, 78, 4, 54 ]

def test_carmichael_lambda():
	for n in range(1, len(OEIS_A002322)+1):
		assert carmichael_lambda(n) == OEIS_A002322[n - 1]

#https://oeis.org/A008683
OEIS_A008683 = [
		1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0,
		0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1, -1, 0, 0, 1, -1, 0, 0,
		0, 1, 0, -1, 0, 1, 0, 1, 1, -1, 0, -1, 1, 0, 0, 1, -1, -1, 0, 1, -1, -1, 0, -1, 1,
		0, 0, 1, -1 ]

def test_moebius_mu():
	for n in range(1, len(OEIS_A008683)+1):
		assert moebius_mu(n) == OEIS_A008683[n - 1]



