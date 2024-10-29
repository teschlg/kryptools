import pytest
from random import randint, seed
from kryptools import sieve_eratosthenes
from kryptools import iroot, ilog, perfect_square, perfect_power, prime_power
seed(0)


def test_iroot():
	assert iroot(3, 0) == 0
	assert iroot(3, 1) == 1
	assert iroot(1, 3) == 3
	for _ in range(100):
		n = randint(0, 10**10)
		for k in range(2, 5):
			r = iroot(k, n)
			assert r ** k <= n and (r+1) ** k > n

def test_ilog():
	for n in range(1, 10000):
		for b in range(2,5):
			l = ilog(b, n)
			assert b**l <= n and b**(l+1) > n
	for _ in range(1000):
		n = randint(1, 10**10)
		for b in range(2,5):
			l = ilog(b, n)
			assert b**l <= n and b**(l+1) > n

def test_perfect_power():
	for n in (3* 7, 4 * 3, 4 * 3 * 5 * 7):
		assert perfect_power(n) == None
		for k in sieve_eratosthenes(11):
			assert perfect_power(n**k) ==  (n, k)
			
def test_prime_power():
	for p in sieve_eratosthenes(30):
		for k in range(1,6):
			assert prime_power(p**k) == (p, k)
			