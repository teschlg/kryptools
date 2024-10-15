"""
Factorization of integers:
    factorint(n) factorize the integer n into prime factors
    divisors(n) unsorted list of all divisors of the integer n
"""

from math import isqrt, gcd, prod
from .primes import sieve_eratosthenes, is_prime
from .factor_pm1 import _pm1_parameters, factor_pm1
from .factor_ecm import _ecm_parameters, factor_ecm
#from .factor_qs import factor_qs


# Factoring


def _factor_fermat(n: int, steps: int = 10) -> list:
    "Fermat method"
    parameters = {11: (12, 6), 23: (12, 0),
                  5: (6, 3), 17: (6, 3),
                  19: (4, 2), 7: (4, 0),
                  1: (2, 1), 13: (2, 1)}
    start = isqrt(n - 1) + 1
    step, mod = parameters[n % 24]
    start += (mod - start) % step
    for a in range(start, min(start + steps * step,(n + 9) // 6) + 1, step):
        b = isqrt(a * a - n)
        if b * b == a * a - n:
            return a - b

def factorint(n: int, verbose: int = 0) -> dict:
    "Factor an ineger `n` into prime factors."
    prime_factors = {}
    if not isinstance(n, int):
        raise ValueError("Number to be factored must be an integer!")
    if n == 0:
        return {0: 1}
    if n < 0:
        n *= -1
        prime_factors[-1] = 1

    def add_factors(m: int, mm: tuple) -> None:
        k = remaining_factors[m]
        del remaining_factors[m]

        for m in mm:
            if m in prime_factors:
                prime_factors[m] += k
            elif is_prime(m):
                prime_factors[m] = k
            else:
                if m in remaining_factors:
                    remaining_factors[m] += k
                elif m in new_factors:
                    new_factors[m] += k
                else:
                    new_factors[m] = k

    # trial division
    B = 2500
    factorbase = sieve_eratosthenes(B)
    for p in factorbase:
        k = 0
        while n % p == 0:
            k += 1
            n //= p
        if k:
            prime_factors[p] = k
    if n == 1:
        return prime_factors
    if verbose:
        print("Trial division found:", list(prime_factors))
    if is_prime(n):
        prime_factors[n] = 1
        return prime_factors
    remaining_factors = { n: 1 }

    # https://gitlab.inria.fr/zimmerma/ecm/
    ECM_PARAMETERS = [
        [    11000,        1900000,    74],
        [    50000,       13000000,   221],
        [   250000,      130000000,   453],
        [  1000000,     1000000000,   984],
        [  3000000,     5700000000,  2541],
        [ 11000000,    35000000000,  4949],
        [ 43000000,   240000000000,  8266],
        [110000000,   780000000000, 20158],
        [260000000,  3200000000000, 47173],
        [850000000, 16000000000000, 77666]
    ]

    for parameters in ECM_PARAMETERS:
        B1, B2, num_curves = parameters
        num_curves *= 2
        D = isqrt(B2)
        primes = sieve_eratosthenes(B1 - 1 + ((B2 - B1 + 1) // (2 * D) + 1) * 2 * D)
        pm1_parameters = _pm1_parameters(10 * B1, B2, primes = primes)
        ecm_parameters = tuple([B1, B2, num_curves] + list(_ecm_parameters(B1, B2, D, primes = primes)))

        methods = {_factor_fermat: "fm", factor_pm1: "pm1", factor_ecm: "ecm"}  #, factor_qs: "qs"}
        while remaining_factors:
            new_factors = {}
            for method in [ _factor_fermat, factor_pm1, factor_ecm ]:  # , factor_qs ]:
                factors = list(remaining_factors)
                for m in factors:
                    if verbose > 1:
                        print("Factoring: ",m, "method", methods[method])
                    if method == factor_pm1:
                        tmp = factor_pm1(m, pm1_parameters = pm1_parameters)
                    elif method == factor_ecm:
                        tmp = factor_ecm(m, ecm_parameters = ecm_parameters)
                    else:
                        tmp = method(m)
                    if tmp:
                        tmp2 = m // tmp
                        g = gcd(tmp, tmp2)
                        if g > 1:
                            tmp3 = [ g, g ]
                            for x in (tmp, tmp2):
                                x //= g
                                while x % g == 0:
                                    tmp3.append(g)
                                    x //= g
                                if x > 1:
                                    tmp3.append(x)
                        else:
                            tmp3 = [tmp, tmp2]
                        if verbose > 0: print("Factors found (", methods[method] ,"): ", tmp3)
                        add_factors(m, tmp3)
            if not new_factors:
                break
            for m in new_factors:
                if m in remaining_factors:
                    remaining_factors[m] += new_factors[m]
                else:
                    remaining_factors[m] = new_factors[m]
            if verbose > 1: print("Remaining: ", remaining_factors)

        if len(remaining_factors) == 0:
            return prime_factors

    raise ValueError("Incomplete factorization:", prime_factors, remaining_factors, new_factors)

# Divisors

def divisors(n:int, proper:bool = False) -> int:
    """Returns an unsorted list of all divisors of an integer `n`."""
    facctordict= factorint(n)
    primes = [p for p in facctordict]
    nprimes = len(primes)
    multiplicites = [facctordict[p] for p in primes]

    exponents = [0] * nprimes
    if proper:
        divisors = []
    else:
        divisors = [1]
    i = 0
    while True:
        exponents[i] += 1
        if exponents[i] > multiplicites[i]:
            while exponents[i] > multiplicites[i]:
                exponents[i] = 0
                i += 1
                exponents[i] += 1
            i = 0
        d = prod(primes[i] ** exponents[i] for i in range(nprimes))
        if d == n:
            if proper:
                return divisors
            divisors.append(n)
            return divisors
        divisors.append(d)
