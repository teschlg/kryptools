"""
Tools for prime numbers:
    sieve_eratosthenes(B) a tuple of all primes up to including B
    isprime(n) test if n is probably prime
    miller_rabin_test(n, b) Miller-Rabin primality test with base b
"""
from math import isqrt

# Erathostenes

def sieve_eratosthenes(B: int) -> list:
    """ "Returns a list of all primes up to (including) max."""
    B1 = (isqrt(B) -1)//2
    B = (B - 1)//2
    is_prime = [True] * (B + 1)  # to begin with, all numbers are potentially prime
    # sieve out the primes p=2*q+1 starting at 3 in steps of 2 (ignoring even numbers)
    for q in range(1, B1 + 1):
        if is_prime[q]: # sieve out all multiples; numbers p*q with q<p were already sieved out previously
            qq = (q << 1) * (q + 1)
            p = (q << 1) | 1
            is_prime[qq :: p] = [False] * ((B - qq) // p + 1)

    return tuple([2] + [2 * q + 1 for q in range(1, B + 1) if is_prime[q]])

# Primality testing


def miller_rabin_test(n: int, bases: list[int] | int) -> bool:
    """Run a Miller-Rabin test with given bases on n."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    m = (n - 1) // 2
    k = 1
    while m % 2 == 0:
        m //= 2
        k += 1
    if isinstance(bases, int):
        bases = [bases]
    for a in bases:
        b = pow(a, m, n)
        if b == 1 or b == n - 1:
            return True
        for _ in range(1, k):
            b = pow(b, 2, n)
            if b == 1:
                return False
            if b == n - 1:
                return True
    return False

def isprime(n: int) -> bool:
    """Test if an integer n if probable prime."""
    if n < 18446744073709551616:  # https://miller-rabin.appspot.com
        return miller_rabin_test(n, [2, 325, 9375, 28178, 450775, 9780504, 1795265022])
    if n < 3317044064679887385961981:
        return miller_rabin_test(n, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41])
    return miller_rabin_test(n, [2]) and _is_strong_lucas_prp(n)  # Baillie–PSW primality test


def _lucas_sequence(n, D, k):
    """Evaluate a Lucas sequence."""
    # P = 1
    Q = (1 - D)//4
    U = 1
    V = 1
    Qk = Q
    b = k.bit_length()
    while b > 1:
        U = (U * V) % n
        V = (V * V - 2 * Qk) % n
        Qk *= Qk
        b -= 1
        if (k >> (b - 1)) & 1:
            U, V = U + V, V + U * D
            if U & 1:
                U += n
            if V & 1:
                V += n
            U, V = U >> 1, V >> 1
            Qk *= Q
        Qk %= n
    return U % n, V % n, Qk


def _is_strong_lucas_prp(n: int) -> bool:
    """Strong Lucas primality test."""
    from math import gcd
    from .nt import jacobi_symbol

    # remove powers of 2 from n+1 (= k * 2**s)
    k = (n + 1) // 2
    s = 1
    while k % 2 == 0:
        k //= 2
        s += 1
    # gernerate parameters
    D = 5
    while True:
        g = gcd(abs(D), n)
        if 1 < g < n:
            return False
        if jacobi_symbol(D, n) == -1:
            break
        if D > 0:
            D = -D - 2
        else:
            D = -D + 2

    if D == 0:
        return False

    U, V, Qk = _lucas_sequence(n, D, k)

    if U == 0 or V == 0:
        return True
    for _ in range(1, s):
        V = (V*V - 2*Qk) % n
        if V == 0:
            return True
        Qk = pow(Qk, 2, n)
    return False
