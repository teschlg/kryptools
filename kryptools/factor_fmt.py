"""
Integer factorization: Fermat's method
"""

from math import isqrt


def factor_fermat(n: int) -> list:
    """Find factors of n using the method of Fermat."""
    factors = []
    # Fermat only works if n has two factors which are either both even or both odd
    while n % 2 == 0:
        factors.append(2)
        n //= 2
    a = isqrt(n - 1) + 1
    step =2
    if n % 3 == 2:  # if n % 3 = 2, then a must be a multiple of 3
        a += 2 - ((a - 1) % 3)
        step = 3
    elif (n % 4 == 1) ^ (a & 1):  # if n % 4 = 1,3 then a must be odd, even, respectively
        a += 1
    while a <= (n + 9) // 6:
        b = isqrt(a * a - n)
        if b * b == a * a - n:
            factors.append(a - b)
            factors.append(a + b)
            return factors
        a += step
    factors.append(n)
    return factors
