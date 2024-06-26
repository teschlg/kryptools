"""
Ring of intergers modulo `n`.
"""

class Zmod:
    """
    Ring of intergers modulo `n`.

    Example:
    
    To define a finite Galois field modulo the prime 5 use
    >>> gf=Zmod(5)
    
    To declare 3 as an element of our Galois filed use
    >>> gf(3)
    3 (mod 5)

    The usual arithmetic operations are supported.
    >>> gf(2) + gf(3)
    0 (mod 5)
    """

    def __init__(self, n: int, short: bool = False):
        self.n = n
        self.short = short

    def __call__(self, x: int):
        return ZmodPoint(x, self)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.n == other.n
        return False

    def __contains__(self, other: "ZmodPoint") -> bool:
        return isinstance(other, ZmodPoint) and self.n == other.ring.n


class ZmodPoint:
    "Represents a point in the ring Zmod."

    def __init__(self, x: int, ring: "Zmod"):
        if isinstance(x, self.__class__) and x.ring.n == ring.n:
            self.x = int(x)
        else:
            self.x = int(x) % ring.n
        self.ring = ring

    def __repr__(self):
        if self.ring.short:
            return str(self.x)
        return f"{self.x} (mod {self.ring.n})"

    def __eq__(self, other):
        if not isinstance(other, self.__class__) or self.ring != other.ring:
            return False
        return self.x == other.x

    def __bool__(self):
        return self.x != 0

    def __int__(self):
        return self.x

    def sharp(self):
        "Returns a symmetric (w.r.t. 0) representative."
        tmp = (self.ring.n - 1) // 2
        return (self.x + tmp) % self.ring.n - tmp

    def __hash__(self):
        return hash(self.x)

    def __add__(self, other: "ZmodPoint") -> "ZmodPoint":
        if isinstance(other, self.__class__) and self.ring == other.ring:
            return self.__class__(self.x + other.x, self.ring)
        if isinstance(other, int):
            return self.__class__(self.x + other, self.ring)
        return NotImplemented

    def __radd__(self, scalar: int) -> "ZmodPoint":
        return self.__class__(scalar + self.x, self.ring)

    def __neg__(self) -> "ZmodPoint":
        return self.__class__(-self.x, self.ring)

    def __sub__(self, other: "ZmodPoint") -> "ZmodPoint":
        if isinstance(other, self.__class__) and self.ring == other.ring:
            return self.__class__(self.x - other.x, self.ring)
        if isinstance(other, int):
            return self.__class__(self.x - other, self.ring)
        return NotImplemented

    def __rsub__(self, scalar: int) -> "ZmodPoint":
        return self.__class__(scalar - self.x, self.ring)

    def __mul__(self, other: "ZmodPoint") -> "ZmodPoint":
        if isinstance(other, self.__class__) and self.ring == other.ring:
            return self.__class__(self.x * other.x, self.ring)
        if isinstance(other, int):
            return self.__class__(self.x * other, self.ring)
        return NotImplemented

    def __rmul__(self, scalar: int) -> "ZmodPoint":
        return self.__class__(scalar * self.x, self.ring)

    def __truediv__(self, other: "ZmodPoint") -> "ZmodPoint":
        if not isinstance(other, self.__class__) or self.ring != other.ring:
            return NotImplemented
        return self.__class__(self.x * pow(other.x, -1, self.ring.n), self.ring)

    def __rtruediv__(self, scalar: int) -> "ZmodPoint":
        return self.__class__(scalar * pow(self.x, -1, self.ring.n), self.ring)

    def __pow__(self, scalar: int) -> "ZmodPoint":
        return self.__class__(pow(self.x, scalar, self.ring.n), self.ring)
