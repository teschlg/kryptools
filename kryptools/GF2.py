"""
Galois field GF(2^n).
"""

from types import new_class
from .conway_polynomials import conway_polynomials2

class GF_2:
    "Represents a point in the Galois field GF(2^n)."
    poly = 0
    power = 0
    order = 0
    print_hex = True
    bitreversed = False

    def __init__(self, x: int):
        if isinstance(x, bytes|bytearray):
            x = int.from_bytes(x, byteorder="big")
        self.x = int(x) % self.__class__.order

    def __repr__(self):
        if self.__class__.print_hex:
            return format(self.x, "0"+str(self.__class__.power // 4)+"x")
        else:
            return self.x

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.x == other.x

    def __bool__(self):
        return bool(self.x)

    def __int__(self):
        return self.x

    def __hash__(self):
        return hash(self.x)

    def __add__(self, other: "GF_2") -> "GF_2":
        if isinstance(other, self.__class__):
            return self.__class__(self.x ^ other.x)
        return NotImplemented

    def __neg__(self) -> "GF_2":
        return self

    def __pos__(self) -> "GF_2":
        return self

    def __sub__(self, other: "GF_2") -> "GF_2":
        if isinstance(other, self.__class__):
            return self.__class__(self.x ^ other.x)
        return NotImplemented

    def __mul__(self, other: "GF_2") -> "GF_2":
        if isinstance(other, int):
            if other % 2:
                return self
            return self.__class__(0)
        if not isinstance(other, self.__class__):
            return NotImplemented
        z = 0
        y = other.x
        if self.__class__.bitreversed:
            for d in bin(other.x)[2:].zfill(self.__class__.power):
                if int(d):
                    z = z ^ y
                if y % 2:
                    y = y >> 1
                    y ^= self.__class__.poly
                else:
                    y = y >> 1
        else:
            for i in range(self.__class__.power):
                z = z << 1
                if z & 0b100000000:
                    z = z ^ self.__class__.poly
                if y & 0b010000000:
                    z = z ^ self.x
                y = y << 1
        return self.__class__(z)

    def __rmul__(self, scalar: int) -> "GF_2":
        if isinstance(scalar, int):
            if scalar % 2:
                return self
            return self.__class__(0)
        return NotImplemented

    def __truediv__(self, other: "GF_2") -> "GF_2":
        if not isinstance(other, self.__class__):
            return NotImplemented
        if not(other.x):
            raise ValueError("Division by zero.")
        return self * other**(self.order -2)

    def __pow__(self, j: int) -> "GF_2":
        if not isinstance(j, int):
            return NotImplemented
        if not(self.x) and j < 0:
            raise ValueError("Division by zero.")
        if self.__class__.bitreversed:
            res = self.__class__(self.__class__.order >> 1)
        else:
            res = self.__class__(1)
        x = self.__class__(self.x)
        j %= self.__class__.order - 1
        while j > 0:
            # If j is odd, multiply with x
            if j & 1:
                res *= x
            # Now square
            j = j >> 1  # j= j/2
            x *= x
        return res

    def __lshift__(self, i: int) -> "GF_2":
        "Cyclic rotation to the left."
        x = self.x << i
        x = (x % self.order) + (x // self.order)
        return(self.__class__(x))

    def __rshift__(self, i: int) -> "GF_2":
        "Cyclic rotation to the right."
        x = self.x >> i
        x += (self.order >> i) * (self.x % 2**i)
        return(self.__class__(x))


def GF2(n: int, poly: int|None = None):
    """
    Create a Galois field GF(2^n).

    Elements of a Galois filed are specified by integers from 0 to 2**n-1, which are
    identified with polynomials whose coefficients are the binary digits of the integer.
    Multiplication is polynomial multiplication modulo a given irreducible polynomial.
    If the polynomial has the n'th bit set, the lowest bit corresponds to the constant coefficient
    (like in AES). Otherwise, the bit order is reversed (like for Ghash). For n = 8 the default
    polynomial is the AES polynomial 1 + x + x^3 + x^4 + x^8. For n = 128 the default polynomial
    is the Ghash polynomial 1 + x + x^2 + x^7 (+ x^128). For ohter n the polynomial has to be given
    explicitely.

    Example:
    
    To define the Galois field GF(2^8) use
    >>> gf=GF(8)
    
    To declare 3 as an element of our Galois field use (elements are displayed as hex numbers)
    >>> gf(3)
    03

    The usual arithmetic operations are supported.
    >>> gf(2)*gf(135)
    15

    """
    if not poly:
        try:
            poly = conway_polynomials2[n]
        except:
            raise ValueError("Unknown power. Please specify a polynomial.")

    name = 'GF2^'+str(n)
    gf = new_class(
        name,
        bases = (GF_2,),
        kwds = {},
    )
    if poly & 2**n == 2**n:
        gf.bitreversed = False
    else:
        gf.bitreversed = True
    gf.power = n
    gf.order = 2**n
    gf.poly = poly
    return gf
