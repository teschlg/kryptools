"""
Integer factorization: Lentra's ECM
"""

from math import gcd, isqrt, log
from random import randint, seed
from .primes import sieve_eratosthenes

# Crandall and Pomerance: Primes (doi=10.1007/0-387-28979-8)
# Algorithm 7.2.7


def dbl(P, c2, p):
    # c2 = c - 2
    t1 = (P[0] + P[1]) % p
    t2 = (P[0] - P[1]) % p
    t3 = P[0] * P[1] % p
    return pow(t1 * t2, 2, p), ((t1 * t1 + c2 * t3) % p) * 4 * t3 % p


def add(P1, P2, P3, p):
    t1 = pow(P1[0] * P2[0] - P1[1] * P2[1], 2, p)
    t2 = pow(P1[0] * P2[1] - P2[0] * P1[1], 2, p)
    return P3[1] * t1 % p, P3[0] * t2 % p


def mult(k, P, c2, p):
    if k == 1:
        return P
    if k == 2:
        return dbl(P, c2, p)
    Q = dbl(P, c2, p)
    R = P
    for i in bin(k)[3:]:
        if i == "1":
            R = add(Q, R, P, p)
            Q = dbl(Q, c2, p)
        else:
            Q = add(R, Q, P, p)
            R = dbl(R, c2, p)
    return R

# Crandall and Pomerance: Primes (doi=10.1007/0-387-28979-8)
# Algorithm 7.4.4 (Inversionless ECM)

def _ecm_parameters(B1: int, B2: int = None, D: int = None, primes: tuple = None):
    "Precompute parameters for the ECM method."

    # Stage-one/two limits must be even
    B1 += B1 & 1
    if not B2:
        B2 = 100 * B1
    else:
        B2 += B2 & 1
    if not D:
        D = isqrt(B2)
    if not primes:
        primes = sieve_eratosthenes(B1 - 1 + ((B2 - B1 + 1) // (2 * D) + 1) * 2 * D)

    stage_one = 1
    for p in primes:
        if p > B1:
            break
        stage_one *= p**int(log(B1, p))

    stage_two_deltas = {}
    r = B1 - 1
    stage_two_deltas[r] = []
    for p in primes:
        if p <= r:
            continue
        if p <= r + 2 * D:
            stage_two_deltas[r].append((p - r) // 2)
        else:
            r += 2 * D
            while p > r + 2 * D:
                stage_two_deltas[r] = []
                r += 2 * D
            stage_two_deltas[r] = [(p - r) // 2]
    return D, stage_one, stage_two_deltas


def factor_ecm(n: int, B1: int = 11000, B2: int = 1900000, curves: int = 74, ecm_parameters: tuple = None):
    "Factors a number n using Lentsta's ECM method."    

    if ecm_parameters:
        curves, D, stage_one, stage_two_deltas = ecm_parameters
    else:
        D, stage_one, stage_two_deltas = _ecm_parameters(B1, B2)

    for _ in range(curves):
        # find a random curve
        sigma = randint(6, n - 1)
        u = (sigma**2 - 5) % n
        v = (4 * sigma) % n
        try:
            c = (pow(v - u, 3, n) * (3 * u + v) * pow(4 * u**3 * v, -1, n) - 2) % n
        except:
            m = gcd(4 * u**3 * v, n)
            return m
        c2 = c - 2 % n
        # and a point on the curve (or its twist)
        Q = (pow(u, 3, n), pow(v, 3, n))
        # Stage one
        Q = mult(stage_one, Q, c2, n)
        if Q[1] == 0:
            continue
        g = gcd(Q[1], n)
        if g > 1:
            return g

        # Stage two
        S = [dbl(Q, c2, n)]
        S.append(dbl(S[0], c2, n))

        for i in range(2, D):
            S.append(add(S[i - 1], S[0], S[i - 2], n))

        beta = []
        for i in range(D):
            beta.append((S[i][0] * S[i][1]) % n)
        g = 1
        T = mult(B1 - 2 * D - 1, Q, c2, n)
        R = mult(B1 - 1, Q, c2, n)
        for r in range(B1 - 1, B2, 2 * D):
            alpha = R[0] * R[1] % n
            for delta in stage_two_deltas[r]:
                g = g * ((R[0] - S[delta - 1][0]) * (R[1] + S[delta - 1][1]) - alpha + beta[delta - 1]) % n
                if g == 0:
                    break
            if g == 0:
                break
            g =gcd(g,n)
            if 1 < g:
                return g
            R, T = add(R, S[D - 1], T, n), R

        g = gcd(g, n)
        if 1 < g < n:
            return g
