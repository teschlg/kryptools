"""
Discrete Log Problem:
    dlog(a, b, n) solve the discrete log problem a^x = b mod n
"""


from math import gcd, log, sqrt
from .nt import crt, order
from .dlp_bsgs import dlog_bsgs
from .dlp_qs import dlog_qs
from .factor import factorint


def dlog_naive(a: int, b: int, n: int, m: int = None) -> int|None:
    """Compute the discrete log_a(b) in Z_n of an element a of order m by exhaustive search."""
    a %= n
    b %= n
    if not m:
        m = n - 1
    aa = 1
    for i in range(m):
        aa = (aa * a) % n
        if aa == b:
            return i + 1
    raise ValueError("DLP not solvable.")

def _dlog_switch(a: int, b: int, n: int, m: int) -> int:
    """Compute the discrete log_a(b) in Z_n of an element a of order m choosing an appropriate method."""
    if m < 1000:
        return dlog_naive(a, b, n, m)
    if log(m) - 6 < 2 * sqrt(log(n) * log(log(n))):  # compare the theoreticaly expected running times ob bsgs and ic; the constant 6 is determined experimentally
        return dlog_bsgs(a, b, n, m)
    return dlog_qs(a, b, n, m)

def _dlog_ph(a: int, b: int, n: int, q: int, k: int) -> int:
    """Compute the discrete log_a(b) in Z_n of an element a of order q^k using Pohlig-Hellman reduction."""
    if k == 1 or q**k < 10000:
        return _dlog_switch(a, b, n, q**k)
    aj = pow(a, q ** (k - 1), n)
    a1 = aj
    bj = pow(b, q ** (k - 1), n)
    xj = _dlog_switch(a1, bj, n, q)
    for j in range(2, k + 1):
        aj = pow(a, q ** (k - j), n)
        bj = pow(b, q ** (k - j), n)
        yj = _dlog_switch(a1, bj * pow(aj, -xj, n) % n, n, q)  # pylint: disable=E1130
        xj = xj + q ** (j - 1) * yj % q**j
    return xj


def dlog(a: int, b: int, n: int, m: int|None = None) -> int:
    """Compute the discrete log_a(b) in Z_n of an element `a` of order `m` using Pohlig-Hellman reduction."""
    a %= n
    b %= n
    assert gcd(a, n) == 1, "a and n must be coprime."
    assert gcd(b, n) == 1, "b and n must be coprime."
    if m:
        mf = factorint(m)
    else:
        m, mf = order(a, n, True)
    if pow(b, m, n) != 1:
        raise ValueError("DLP not solvable.")
    # We first use Pohlig-Hellman to split m into powers of prime factors
    mm = []
    ll = []
    for pj, kj in mf.items():
        aj = pow(a, m // pj**kj, n)
        bj = pow(b, m // pj**kj, n)
        l = _dlog_ph(aj, bj, n, pj, kj)
        mm += [pj**kj]
        ll += [l]
    return crt(ll, mm)
