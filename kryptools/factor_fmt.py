"""
Integer factorization: Fermat's method
"""

from math import isqrt


def factor_fermat(n: int) -> list:
    """Find factors of n using the method of Fermat."""
    factors = []
    parameters = {11: (12, 6), 23: (12, 0),
                  5: (6, 3), 17: (6, 3),
                  19: (4, 2), 7: (4, 0),
                  1: (2, 1), 13: (2, 1)}
   # Speedup only works if n is neiter a multipl of 2 or 3
    for p in (2, 3):
        while n % p == 0:
            factors.append(p)
            n //= p
    start = isqrt(n - 1) + 1
    step, mod = parameters[n % 24]
    start += (mod - start) % step
    for a in range(start, (n + 9) // 6 + 1, step):
        b = isqrt(a * a - n)
        if b * b == a * a - n:
            factors.append(a - b)
            factors.append(a + b)
            return factors
    factors.append(n)
    return factors
