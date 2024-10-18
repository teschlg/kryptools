"""
Tools for prime numbers:
    sieve_eratosthenes(B) a tuple of all primes up to including B
    prime_pi(n) number of primes below or equal n
    is_prime(n) test if n is probably prime
    next_prime(n) find the next prime larger or equal n
    previous_prime(n) find the previous prime smaller or equal n
    random_prime(l) find a pseudorandom prime with bit length at least l
    random_strongprime(l) find a pseudorandom strong prime with bit length at least l
    is_safeprime(n) test if n is a safe prime
    random_safeprime(l) find a pseudorandom safe prime with bit length at least l and ord(2)=(p-1)/2
    is_blumprime(n) test if n is a Blum prime
    random_blumprime(l) find a pseudorandom Blum prime with bit length at least l
    miller_rabin_test(n, b) Miller-Rabin primality test with base b
"""
from math import isqrt, gcd, floor
from random import randint
from .nt import jacobi_symbol

# Erathostenes

def sieve_eratosthenes(B: int) -> tuple:
    """Returns a tuple of all primes up to (including) `B`."""
    if B < 2:
        return ()
    if B < 3:
        return [2]
    B= floor(B)
    B1 = (isqrt(B) - 1)//2
    B = (B - 1)//2
    isprime = [True] * (B + 1)  # to begin with, all odd numbers are potentially prime
    # sieve out the primes p=2*q+1 starting at q = 1
    for q in range(1, B1 + 1):
        if isprime[q]: # sieve out all multiples; numbers p*r with r<p were already sieved out previously
            qq = 2 * q * (q + 1)
            p = 2 * q + 1
            isprime[qq :: p] = [False] * ((B - qq) // p + 1)

    return tuple([2] + [2 * q + 1 for q in range(1, B + 1) if isprime[q]])

def prime_pi(x: float) -> int:
    """Prime-counting function pi(x). Returns the number of primes below or equal `x`."""
    return len(sieve_eratosthenes(x))


# Primality testing


def miller_rabin_test(n: int, bases: list[int] | int) -> bool:
    """Run a Miller-Rabin test with given bases on `n`."""
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

def is_prime(n: int) -> bool:
    """Test if an integer `n` if probable prime."""
    if n < 18446744073709551616:  # https://miller-rabin.appspot.com
        return miller_rabin_test(n, [2, 325, 9375, 28178, 450775, 9780504, 1795265022])
    if n < 3317044064679887385961981:
        return miller_rabin_test(n, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41])
    return miller_rabin_test(n, [2]) and _is_strong_lucas_prp(n)  # Baillie–PSW primality test

def next_prime(n: int) -> int:
    """Find the next prime larger or equal `n`."""
    if n < 3:
        return 2
    n |= 1  # make sure n is odd
    while not is_prime(n):
        n += 2
    return n

def previous_prime(n: int) -> int:
    """Find the previous prime smaller or equal `n`."""
    if n < 3:
        return 2
    n += (n & 1) - 1  # make sure n is odd
    while not is_prime(n):
        n -= 2
    return n

def random_prime(l: int) -> int:
    """Find a pseudorandom prime with bit length at least `l`."""
    return next_prime(randint(2 ** (l - 1), 2**l - 1))

def random_strongprime(l: int) -> int:
    """Find a pseudorandom strong prime with bit length at least `l` using Gordon's algorithm."""
    t = random_prime(l)
    s = random_prime(l)
    u = 2 * t
    uu = u * randint(1, 100)
    while not is_prime(uu + 1):
        uu += u
    r = uu + 1
    u = 2 * r * s
    uu = u * randint(1, 100) + 2 * s * pow(s, r - 2, r) - 1
    while not is_prime(uu):
        uu += u
    return t, s, r, uu

def is_safeprime(p: int) -> bool:
    """Tests if a number is a safe prime."""
    return is_prime(p) and is_prime((p - 1) // 2)

def random_safeprime(l: int) -> int:
    """Find a pseudorandom safe prime with bit length at least `l` and `ord(2)=(p-1)/2`."""
    p = randint(2 ** (l - 1), 2**l - 1)
    p = p - (p % 24) + 23
    while not is_safeprime(p):
        p += 24
    return p

def is_blumprime(p: int) -> bool:
    """Tests if a number is a Blum prime."""
    return p % 4 == 3 and is_prime(p)

def random_blumprime(l: int) -> int:
    """Find a pseudorandom Blum prime with bit length at least `l`."""
    p = randint(2 ** (l - 1), 2**l - 1)
    p = p - (p % 4) + 3
    while not is_prime(p):
        p += 4
    return p

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
