"""
Integer factorization: Quadratic sieve
"""

from math import isqrt, gcd, sqrt, log, exp, ceil
from .primes import sieve_eratosthenes
from .nt import legendre_symbol, sqrt_mod

def bytexor(a: bytearray, b: bytes) -> bytes:
    """Xor the bytestring b to the bytearry a."""
    for i in range(len(a)):
        a[i] ^= b[i]


def byteset(a: bytearray, i: int) -> bytes:
    """Set the i'th bit in the bytearray a to 1."""
    i1, i0 = divmod(i, 8)
    a[i1] |= 2**i0


def bytetest(a: bytes, i: int) -> bytes:
    """Test if the i'th bit in the bytestring a is 1."""
    i1, i0 = divmod(i, 8)
    return (a[i1] & 2**i0) != 0


def factor_qs(n: int) -> list:
    """Find factors of n using the quadratic sieve."""
    # first determine the bound B for the factorbase: Choosing B=p^(1/u) Canfield-Erdös-Pomerance gives us
    # the expected running time |B|^2 u^u = u^(u+2) p^(2/u)/log(n). There is no explicit expression for the optimum, hence
    # we use Newton
    u = 2 * sqrt(log(n) / log(log(n)))  # asymptotic value
    for _ in range(3):
        u = (2 * log(n) + u * u * (2 + log(u))) / (
            2 + 3 * u + 2 * u * log(u)
        )  # Newton iteration
    B = int(exp(log(n) / u))
    # B = int(exp(0.5 * sqrt( log(n) * log(log(n)) )*( 1 + 1/log(log(n)) )))

    factorbase = []
    for p in sieve_eratosthenes(B):  # compute the factorbase
        ls = legendre_symbol(n, p)
        if ls == 0:  # we already found a factor;-)
            if n == p:
                return n
            return p
        if ls == 1:  # we only take primes such that n is quadratic residue
            factorbase.append(p)
    lf = len(factorbase)  # length of the factorbase (including -1)
    lfb = ceil((lf+1) / 8)  # the number of bytes we need to store a relation

    factorbase_log = [log(p) for p in factorbase]  # we add these up to test if a number will probably factor
    factorbase_root = [0] * lf # compute the roots of n mod p
    for i, p in enumerate(factorbase):
        r1 = sqrt_mod(n, p)
        assert r1 is not None # only quadratic residues
        r2 = -r1 % p
        if r1 == r2:
            factorbase_root[i] = [ r1 ]
        else:
            factorbase_root[i] = [ r1, -r1 % p ]

    m = isqrt(n - 1) + 1
    d = m**2 - n
    if d == 0: # Our number is a square
        return m
    m2 = 2 * m

    relation_no = -1
    relations = [None] * (lf + 1)  # here we will store the relations in case we found a B-smooth number
    values = [None] * (lf + 1)  # here we will store the values leading to the B-smooth numbers

    # set up the sieve
    sieve_step = 100
    sieve_bound = sieve_step
    # start values for the iterators in the positive/negative direction
    iterator_p = [ [0, 0] for _ in range(lf) ]
    iterator_m = [ [0, 0] for _ in range(lf) ]
    for i, p in enumerate(factorbase):
        for j, r in enumerate(factorbase_root[i]):
            iterator_p[i][j] = (r - m) % p
            iterator_m[i][j] = iterator_p[i][j] - p  # one step in the negative direction

    def do_sieve():
        nonlocal factorbase, factorbase_root, sieve_bound, factors
        for i, p in enumerate(factorbase):
            for r in range(len(factorbase_root[i])):
                j = iterator_p[i][r]
                while j <= sieve_bound:
                    if j in factors:
                        factors[j][0].append(i)
                        factors[j][1] += factorbase_log[i]
                    else:
                        factors[j] = [[i], factorbase_log[i]]
                    j += p
                iterator_p[i][r] = j  # store as start value for the next step
                j = iterator_m[i][r]
                while j > -sieve_bound:
                    if j in factors:
                        factors[j][0].append(i)
                        factors[j][1] += factorbase_log[i]
                    else:
                        factors[j] = [[i], factorbase_log[i]]
                    j -= p
                iterator_m[i][r] = j  # store as start value for the next step

    def process_relation(j: int, relation: bytes):
        nonlocal relation_no, values, relations
        relation_no += 1
        rhs = bytearray(b"\x00") * lfb # construct the k'th row of the identity matrix
        byteset(rhs, relation_no)
        relation += rhs  # extend the relation with this row
        values[relation_no] = j  # store the value which lead to the relation
        # do the Gauss elimination
        index = lf # this will be the index of the first nonzero entry
        # print(f'{j:3}', ' '.join(f'{b:08b}' for b in reversed(relation)))
        for i in range(lf):
            if bytetest(relation, i) and relations[i] is not None:  # make this entry zero if we can (Gauss elimination)
                bytexor(relation, relations[i])
            if bytetest(relation, i) and index == lf:  # is this the index of the first nonzero entry?
                index = i
        # print(f'{j:3}', ' '.join(f'{b:08b}' for b in reversed(relation)))
        if index == lf:  # the new relation is linearly dependent: we have found a linear combination of the 0 vector
            # now we need to determine the factors
            u = 1  # product over all values m - j such that f(m-j) is B-smooth
            v = 1  # sqrt of the product over all f(m-j) which are B-smooth
            w = 1  # save nonsquare terms for the next round
            for i in range(lf):
                if bytetest(relation, lfb * 8 + i):  # select the relations which sum to the 0 vector
                    ui = m + values[i]
                    u = (u * ui) % n
                    vi = ui ** 2 - n
                    d = gcd(w, vi)
                    v = (v * d) % n  # sqrt of the part which is already square
                    w = (w // d) * (vi // d)  # this part is not square yet
            v = v * isqrt(w) % n
            res = gcd(u - v, n)
            if 1 < res < n:
                return res
            relation_no -= 1  # this one did not work, try again
        else:
            relations[index] = relation
        return None

    while True:
        factors = {} # store the index of the primes dividing j

        do_sieve()

        for j in sorted(factors.keys(),key=abs):
            v = d + (m2 + j) * j  # (j + m) ** 2 - n
            if factors[j][1] > 0.49 * log(abs(v)):
                primes = factors[j][0]
                mask = bytearray(b"\x00") * lfb
                if v < 0:
                    byteset(mask, 0)
                    v *= -1
                for i in primes:
                    p = factorbase[i]
                    k = 1  # per our sieve we already know that p devides v
                    v //= p
                    while v % p == 0:  # divide by p as many times as possible
                        k = (k + 1) % 2  # we only want to know if the exponent is even or odd
                        v //= p
                    if k:
                        byteset(mask, i + 1)
                if v == 1:
                    res = process_relation(j, mask)
                    if res:
                        return res
                del(factors[j])
            else:
                del(factors[j])  # the number is unlikely to factor
        sieve_bound += sieve_step  # increase the sieve and continue
