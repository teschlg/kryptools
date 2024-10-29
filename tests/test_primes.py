import pytest
from math import isqrt
from random import randint, seed
from kryptools import sieve_eratosthenes, prime_pi, next_prime, previous_prime, is_prime, is_safeprime, miller_rabin_test, lucas_test
seed(0)

#https://oeis.org/A000040
OEIS_A000010 = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
	71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
	163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
	257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
	359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457,
	461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571,
	577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
	677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
	809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
	919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021,
	1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109,
	1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 ]

def sieve(max: int) -> list:
	is_prime = [True] * (max + 1)  # to begin with, all numbers are potentially prime
	is_prime[0 :: 2] = [False] * (max // 2 + 1)
	is_prime[1] = False
	is_prime[2] = True

	# sieve out the primes starting at 3 in steps of 2 (ignoring even numbers)
	for p in range(3, isqrt(max - 1) + 2, 2):
		if is_prime[p]: # sieve out all multiples; note that numbers p*q with q<p were already sieved out previously
			pp = p * p
			p2 = p + p
			is_prime[pp :: p2] = [False] * ((max - pp) // p2 + 1)
	return is_prime

prime_test = sieve(OEIS_A000010[-1])
assert [ p for p in range(OEIS_A000010[-1]+1) if prime_test[p] ] == OEIS_A000010

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

max_sieve = 25326001
prime_test = sieve(max_sieve)

def test_is_prime():
	for n in range(max_sieve + 1):
		assert is_prime(n) == prime_test[n]

OEIS_A005385 = [ 5, 7, 11, 23, 47, 59, 83, 107, 167, 179, 227, 263, 347, 359, 383, 467,
	479, 503, 563, 587, 719, 839, 863, 887, 983, 1019, 1187, 1283, 1307, 1319, 1367, 1439,
	1487, 1523, 1619, 1823, 1907, 2027, 2039, 2063, 2099, 2207, 2447, 2459, 2579, 2819,
	2879, 2903, 2963 ]

def test_is_safeprime():
	for n in range(OEIS_A000720[-1]+2):
		assert is_safeprime(n) == (n in OEIS_A005385)
	for n in OEIS_A005385:
		assert is_safeprime(n) == True

OEIS_A014233 = [ 2047, 1373653, 25326001, 3215031751, 2152302898747, 3474749660383,
	341550071728321, 341550071728321, 3825123056546413051, 3825123056546413051,
	3825123056546413051, 318665857834031151167461, 3317044064679887385961981 ]

if OEIS_A014233[2] > max_sieve:
	prime_test = sieve(OEIS_A014233[2])

def test_miller_rabin():
	for i, n in enumerate(OEIS_A014233[0:3]):
		bases = OEIS_A000010[:i+1]
		assert miller_rabin_test(n, bases) == True
		if i < 1:
			for m in range(3, n, 2):
				assert miller_rabin_test(m, bases) == prime_test[m]
		else:
			for _ in range(500 * i):
				r = randint(2, n)
				assert miller_rabin_test(m, bases) == prime_test[m]

OEIS_A217255 = [ 5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519,
	75077, 97439, 100127, 113573, 115639, 130139, 155819, 158399, 161027, 162133, 176399,
	176471, 189419, 192509, 197801, 224369, 230691, 231703, 243629, 253259, 268349,
	288919, 313499, 324899 ]

if OEIS_A217255[30] > max_sieve:
	prime_test = sieve(OEIS_A217255[30])

def test_lucas():
	for m in range(2, OEIS_A217255[30]+1):
		if not (lucas_test(m) == prime_test[m]):
			assert m in OEIS_A217255
