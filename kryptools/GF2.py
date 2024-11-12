"""
Galois field GF(2^n).
"""

from types import new_class
from .poly import Poly
from .Zmod import Zmod
from .conway_polynomials import conway_polynomials2


class GF_2:
    "Represents a point in the Galois field GF(2^n)."
    modulus: int = 0  # integer whose biniary digits are the coefficients of the irreducible modulusnomial
    power: int = 0  # n
    order: int = 0  # 2**n
    print_hex: bool = True
    bitreversed: bool = False
    byteorder: str = "big"  # big|little

    def __init__(self, x: int):
        if isinstance(x, bytes | bytearray):
            x = int.from_bytes(x, byteorder=self.__class__.byteorder)
        self.x = int(x) % self.__class__.order

    def __iter__(self):
        for x in range(self.__class__.order):
            yield self.__class__(x)

    def __repr__(self):
        if self.__class__.print_hex:
            return format(self.x, "0"+str(self.__class__.power // 4)+"x")
        return self.x

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.x == other.x

    def __bool__(self):
        return bool(self.x)

    def __int__(self):
        return self.x

    def __bytes__(self):
        return self.x.to_bytes(self.__class__.power // 8, byteorder=self.__class__.byteorder)

    def poly(self):
        "Convert to a polynomial"
        Z_2 = Zmod(2)
        if self.bitreversed:
            coeff = [Z_2(int(d)) for d in str(
                format(self.x, "0" + str(self.power) + "b"))]
            coeff_modulus = [Z_2(int(d)) for d in str(
                format((self.modulus << 1) + 1, "0" + str(self.power) + "b"))]
        else:
            coeff = reversed([Z_2(int(d)) for d in str(
                format(self.x, "0" + str(self.power) + "b"))])
            coeff_modulus = reversed([Z_2(int(d)) for d in str(
                format(self.modulus, "0" + str(self.power) + "b"))])
        return Poly(list(coeff), modulus=list(coeff_modulus))

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
            for d in bin(self.x)[2:].zfill(self.__class__.power):
                if int(d):
                    z ^= y
                if y % 2:
                    y >>= 1
                    y ^= self.__class__.modulus
                else:
                    y >>= 1
        else:
            bitmap1 = self.__class__.order
            bitmap2 = 1
            for _ in range(self.__class__.power):
                if self.x & bitmap2:
                    z ^= y
                bitmap2 <<= 1
                y <<= 1
                if y & bitmap1:
                    y ^= self.__class__.modulus
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
        if not other.x:
            raise ValueError("Division by zero.")
        return self * other**(self.order - 2)

    def __pow__(self, j: int) -> "GF_2":
        if not isinstance(j, int):
            return NotImplemented
        if not (self.x) and j < 0:
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
        return self.__class__(x)

    def __rshift__(self, i: int) -> "GF_2":
        "Cyclic rotation to the right."
        x = self.x >> i
        x += (self.order >> i) * (self.x % 2**i)
        return self.__class__(x)


def GF2(n: int, modulus: int | None = None):
    """
    Create a Galois field GF(2^n).

    Elements of a Galois filed are specified by integers from 0 to 2**n-1, which are
    identified with polynomials whose coefficients are the binary digits of the integer.
    Multiplication is polynomial multiplication modulo a given irreducible polynomial.
    If the polynomial has the n'th bit set, the lowest bit corresponds to the constant coefficient
    (like in AES). Otherwise, the bit order is reversed (like for Ghash). The default modulus is
    choosen from a list of Conway polynomials.

    Example:

    To define the Galois field GF(2^8) use
    >>> gf=GF2(8)

    To declare 3 as an element of our Galois field use (elements are displayed as hex numbers)
    >>> gf(3)
    03

    The usual arithmetic operations are supported.
    >>> gf(2) * gf(123)
    f6

    """
    if not modulus:
        if n in conway_polynomials2:
            modulus = conway_polynomials2[n]
        else:
            raise ValueError("Unknown power. Please specify a modulus.")

    name = 'GF2^'+str(n)
    gf = new_class(
        name,
        bases=(GF_2,),
        kwds={},
    )
    if modulus & 2**n == 2**n:
        gf.bitreversed = False
    else:
        gf.bitreversed = True
    gf.power = n
    gf.order = 2**n
    gf.modulus = modulus
    return gf
