"""
Integer factorization: Pollard's p-1 method
"""

from math import gcd, log
from .primes import sieve_eratosthenes

def _pm1_parameters(B1: int, B2: int = None, primes: tuple = None):
    "Precompute parameters for the P-1 method."
    if not B2:
        B2 = 100 * B1
    if not primes:
        primes = sieve_eratosthenes(B2)

    stage_one = 1
    for p in primes:
        if p > B1:
            break
        stage_one *= p**int(log(B1, p))

    stage_two_deltas = []
    for q in primes:
        if q <= B1:
            continue
        if q > B2:
            break
        stage_two_deltas.append(q - p)

    return stage_one, stage_two_deltas


def factor_pm1(n: int, B1: int = 11000, B2: int = 1900000, x: int = 2, pm1_parameters: tuple = None):
    "Factors a number n using Pollard's P-1 method."    

    if pm1_parameters:
        stage_one, stage_two_deltas = pm1_parameters
    else:
        stage_one, stage_two_deltas = _pm1_parameters(B1, B2)

    g = gcd(pow(x, stage_one, n) - 1, n)
    if 1 < g < n:
        return g
    k = 0
    for d in stage_two_deltas:
        saved = {}
        if not d in saved:
            y = pow(x, d, n)
            saved[d] = y
        x = (x * saved[d]) % n
        k += 1
        if k % 20 == 0:
            g = gcd(x-1, n)
            if 1 < g < n:
                return g
    g = gcd(x-1, n)
    if 1 < g < n:
        return g
    return None
