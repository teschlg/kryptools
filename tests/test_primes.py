import pytest
from math import gcd
from kryptools import sieve_eratosthenes, prime_pi, next_prime, previous_prime


#https://oeis.org/A000040
OEIS_A000010 = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
	71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163,
	167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
	269, 271 ]

def test_sieve_eratosthenes():
	for i, p in enumerate(sieve_eratosthenes(OEIS_A000010[-1])):
		assert p == OEIS_A000010[i]
		if i > 0:
			assert previous_prime(p-1) == OEIS_A000010[i-1]
		if i < len(OEIS_A000010) - 1:
			assert next_prime(p+1) == OEIS_A000010[i+1]

#https://oeis.org/A000040
OEIS_A000720 = [ 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9,
	9, 9, 9, 9, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 14, 14,
	15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 17, 17, 18, 18, 18, 18, 18, 18, 19,
	19, 19, 19, 20, 20, 21, 21, 21, 21, 21, 21 ]

def test_prime_pi():
	for n in range(1, len(OEIS_A000720)+1):
		assert prime_pi(n) == OEIS_A000720[n - 1]
