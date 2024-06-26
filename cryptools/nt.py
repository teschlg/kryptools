"""
Number theory tools:
    lcm(a, b) least common mutiple of a and b
    egcd(a,b) extended Euclidean agorithm
    crt([a1, a2, ...],[m1, m2, ...]) Chinese Remainder Theorem
    cf(Fraction(m,n)) continued fraction expansions
    convergents() convergents of a continued fraction
    sqrt_mod(n, p) square root of n modulo a prime p
    order(a, n) oder of a in the multiplicative group Z_n^*
"""
from math import gcd, prod
from fractions import Fraction
from .factor import factorint

# Euclid and friends

def lcm(a: int, b: int) -> int:
    """Compute the least common multiple of a and b."""
    if b == 0:
        return 0
    if bool(a > 0) != bool(b > 0):
        a = -a
    return (a // gcd(a, b)) * b


def egcd(a: int, b: int) -> (int, int, int):
    """Perform the extended Euclidean agorithm. Returns gcd, x, y such that a x + b y = gcd."""
    r0, r1 = a, b
    x0, x1, y0, y1 = 1, 0, 0, 1
    while r1 != 0:
        q, r = divmod(r0, r1)
        r0, r1 = r1, r
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return r0, x0, y0


# Chinese remainder theorem

def crt(a: list, m: list) -> int:
    """Solve given linear congruences x[j] % m[j] == a[j] using the Chinese Remainder Theorem."""
    l = len(a)
    assert len(m) == l, "The lists of numbers and modules must have equal length."
    M = prod(m)
    Mi = [M // m[i] for i in range(l)]
    MNi = [Mi[i] * pow(Mi[i], -1, m[i]) % M for i in range(l)]
    return sum([a[i] * MNi[i] % M for i in range(l)]) % M

# Continued fractions

# class Fraction:
#     "Rationl number"
#
#     def __init__(self, numerator: int, denominator: int = 1):
#         assert denominator, "Denominator must be nozero."
#         tmp = gcd(denominator, numerator)
#         if tmp != 1:
#             denominator //= tmp
#             numerator //= tmp
#         self.numerator = numerator
#         self.denominator = denominator
#
#     def __repr__(self):
#         if self.denominator == 1:
#             return str(self.numerator)
#         else:
#             return str(self.numerator) + "/" + str(self.denominator)
#
#     def __eq__(self, other):
#         if not isinstance(other, self.__class__):
#             return False
#         return self.denominator == other.denominator and self.numerator == other.numerator
#
#     def __bool__(self):
#         return self.numerator != 0
#
#     def __add__(self, other: "Fraction") -> "Fraction":
#         if isinstance(other, self.__class__):
#             return Rational(self.numerator * other.denominator + other.numerator * self.denominator, self.denominator * other.denominator)
#         if isinstance(other, int):
#             return Rational(self.numerator + other * self.denominator, self.denominator)
#         raise ValueError(f"Cannot add {self} and {other}.")
#
#     def __neg__(self) -> "Fraction":
#         return Rational(- self.numerator, self.denominator)
#
#     def __sub__(self, other: "Fraction") -> "Fraction":
#         if isinstance(other, self.__class__):
#             return Rational(self.numerator * other.denominator - other.numerator * self.denominator, self.denominator * other.denominator)
#         if isinstance(other, int):
#             return Rational(self.numerator - other * self.denominator, self.denominator)
#         raise ValueError(f"Cannot subtract {self} and {other}.")
#
#     def __mul__(self, other: "Fraction") -> "Fraction":
#         if isinstance(other, self.__class__):
#             return Rational(self.numerator * other.numerator, self.denominator * other.denominator)
#         if isinstance(other, int):
#             return Rational(self.numerator * other, self.denominator)
#         raise ValueError(f"Cannot multiply {self} and {other}.")
#
#     def __rmul__(self, other: int) -> "Fraction":
#         return Rational(self.denominator * other, self.numerator)
#
#     def __truediv__(self, other: "Fraction") -> "Fraction":
#         if isinstance(other, self.__class__):
#             return Rational(self.numerator * other.denominator, self.denominator * other.numerator)
#         if isinstance(other, int):
#             return Rational(self.numerator, self.denominator * other)
#         raise ValueError(f"Cannot divide {self} and {other}.")
#
#     def __pow__(self, scalar: int) -> "Fraction":
#         return Rational(self.numerator**scalar, self.denominator**scalar)


def fraction_repr(self):
    if self.denominator == 1:
        return str(self.numerator)
    else:
        return str(self.numerator) + "/" + str(self.denominator)

Fraction.__repr__ = fraction_repr


def cf(x: Fraction) -> list[int]:
    "Compute the continued fraction expansion of the rational number m/n."
    res = []
    m = x.numerator
    n = x.denominator
    while True:
        d, m = divmod(m, n)
        res.append(d)
        if not m:
            return res
        m, n = n, m


def convergents(cont_frac: list[int]) -> float:
    "Compute the convergents of a continued fraction expansion."
    res = []

    def recursion(x: int, m: int, mm: int) -> (int, int):
        return x * m + mm, m

    m, mm = 1, 0
    n, nn = 0, 1
    for x in cont_frac:
        m, mm = recursion(x, m, mm)
        n, nn = recursion(x, n, nn)
        res.append(Fraction(m, n))
    return res


# Legendre symbol

def legendre_symbol(a: int, p: int) -> int:
    """Compute the Legendre symbol of a with respect to the prime p."""
    a = a % p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) == 1:
        return 1
    return -1

def jacobi_symbol(a: int, n: int) -> int:
    """Compute the Jacobi symbol of a with respect to the integer n."""
    # Crandall/Pomerance Algorithm 2.3.5
    a = a % n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            tmp = n % 8
            if tmp == 3 or tmp == 5:
                t *= -1
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t *= -1
        a %= n
    if n == 1:
        return t
    return 0

def sqrt_mod(a: int, p: int) -> list:
    "Compute a square root of a modulo p unsing Cipolla's algorithm."
    a %= p
    if a == 0 or a == 1:
        return a
    if pow(a, (p - 1) // 2, p) != 1:  # Legendre symbol must be equal one
        return None
    if p % 4 == 3:  # easy case
        return pow(a, (p + 1) // 4, p)
    # Cipolla's algorithm
    for c in range(1, p):  # find a field extension
        tmp = (c * c - a) % p
        if tmp == 0:
            return c
        if pow(tmp, (p - 1) // 2, p) == p - 1:
            break
    r = (c * c - a) % p
    i = (p + 1) // 2
    x1, x2 = c, 1  # compute x^i in Z_p(sqrt(r))
    y1, y2 = 1, 0
    while i > 0:
        if i & 1:  # if i is odd, multiply with x
            y1, y2 = (x1 * y1 + x2 * y2 * r) % p, (x1 * y2 + x2 * y1) % p
        x1, x2 = (pow(x1, 2, p) + pow(x2, 2, p) * r) % p, 2 * (x1 * x2) % p  # now square
        i = i >> 1  # i= i/2
    return y1

# Euler phi and Carmichael function

from math import prod


def euler_phi(n: int) -> int:
    """Euler's phi function of n."""
    k = factorint(n)
    return prod([(p - 1) * p ** (k[p] - 1) for p in k])


def carmichael_lambda(n: int) -> int:
    """Carmichael's lambda function of n."""
    k = factorint(n)
    lam_all = []  # values corresponding to the prime factors
    for p in k:
        lam = (p - 1) * p ** (k[p] - 1)
        if p == 2 and k[p] > 2:
            lam = lam // 2
        lam_all += [lam]
    lam = lam_all[0]  # now take the least common multiple of all values
    for l in lam_all[1:]:
        lam = lcm(lam, l)
    return lam

# Order in Z_p^*

def order(a: int, n: int, factor=False) -> int:
    """Compute the order of a in the group Z_n^*."""
    a %= n
    assert a != 0 and gcd(a, n) == 1, f"{a} and {n} are not coprime!"
    factors = dict()  # We compute euler_phi(n) and its factorization in one pass
    for p, k in factorint(n).items():  # first factorize n
        for pm, km in factorint(p - 1).items():  # factor p-1 and add the factors
            if pm in factors:
                factors[pm] += km
            else:
                factors[pm] = km
        if k > 1:  # if the multiplicity of of p is >1, then we need to add p**(k-1)
            if p in factors:
                factors[p] += k - 1
            else:
                factors[p] = k - 1
    order = 1  # compute the group order euler_phi(n) as our current guess
    for p, k in factors.items():
        order *= p**k
    if factor:  # we compute the factorization of the order along the way
        factors_order = {}  # factorization of the order
    for p, k in factors.items():
        i = 0
        for _ in range(k):
            order_try = order // p
            if pow(a, order_try, n) == 1:
                order = order_try
                i += 1
            else:
                break
        if factor and i < k:
            factors_order[p] = k - i
    if factor:
        return order, factors_order
    return order
