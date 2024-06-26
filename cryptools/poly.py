"""
Polynomials
"""

class Poly:
    """
    Represents a polynomial as a list of coefficients.
    
    Example:
    
    To define a polynomial as alist of coefficients use
    >>> Poly([1, 2, 3])
    3 x^2 + 2 x + 1
    """

    def __init__(self, coeff: list, ring = None, modulus: list = None):
        for i in range(len(coeff) - 1, 0, -1):
            if coeff[i]:
                break
            coeff.pop(i)
        self.coeff = coeff
        self.modulus = modulus
        if ring:
            self.map(ring)
        if modulus:
            self.mod(modulus)

    def __getitem__(self, item):
        return self.coeff[item]

    def __repr__(self):
        def prx(i: int):
            if i == 0:
                return ""
            if i == 1:
                return "x"
            return "x^" + str(i)

        if len(self.coeff) == 1:
            return str(int(self.coeff[0]))
        plus = ""
        tmp = ""
        for i in reversed(range(len(self.coeff))):
            s = self.coeff[i]
            if not s:
                continue
            if not s - 1 and i != 0:
                tmp += plus + prx(i)
                plus = " + "
                continue
            try:
                if s == -1 and i != 0:
                    if plus:
                        plus = " "
                    tmp += plus + "- " + prx(i)
                    plus = " + "
                    continue
            except:
                pass
            try:
                if plus and s < 0:
                    tmp += " - " + str(-s) + " " + prx(i)
                    continue
            except:
                pass
            tmp += plus + str(s) + " " + prx(i)
            plus = " + "
        return tmp.strip()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.coeff == other.coeff

    def __bool__(self):
        return bool(self.degree()) or bool(self.coeff[0])

    def degree(self):
        return len(self.coeff) - 1

    def map(self, func):
        self.coeff = list(map(func, self.coeff))

    def __add__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot add {self} and {other}.")
        ls, lo = len(self.coeff), len(other.coeff)
        if ls < lo:
            scoeff = self.coeff + (lo - ls) * [0]
        else:
            scoeff = self.coeff
        if ls > lo:
            ocoeff = other.coeff + (ls - lo) * [0]
        else:
            ocoeff = other.coeff
        modulus = self.modulus
        if not modulus and other.modulus:
            modulus = other.modulus
        return self.__class__([s + o for s, o in zip(scoeff, ocoeff)], modulus=modulus)

    def __neg__(self) -> "Poly":
        return Poly([-s for s in self.coeff], modulus=self.modulus)

    def __sub__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot subtract {self} and {other}.")
        ls, lo = len(self.coeff), len(other.coeff)
        if ls < lo:
            scoeff = self.coeff + (lo - ls) * [0]
        else:
            scoeff = self.coeff
        if ls > lo:
            ocoeff = other.coeff + (ls - lo) * [0]
        else:
            ocoeff = other.coeff
        modulus = self.modulus
        if not modulus and other.modulus:
            modulus = other.modulus
        return self.__class__([s - o for s, o in zip(scoeff, ocoeff)], modulus=modulus)

    def __mul__(self, other: "Poly") -> "Poly":
        if isinstance(other, int):
            return Poly([other * s for s in self.coeff])
        if not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot multiply {self} and {other}.")
        ls, lo = len(self.coeff), len(other.coeff)
        coeff = [0] * (ls + lo - 1)
        for k in range(ls + lo - 1):
            coeff[k] = sum(
                [
                    self.coeff[j] * other.coeff[k - j]
                    for j in range(max(0, k - lo + 1), min(ls, k + 1))
                ]
            )
        modulus = self.modulus
        if not modulus and other.modulus:
            modulus = other.modulus
        return self.__class__(coeff, modulus=modulus)

    def __rmul__(self, other: int) -> "Poly":
        return self.__class__([other * s for s in self.coeff], modulus=self.modulus)

    def __pow__(self, i: int) -> "Poly":
        res = self.__class__([1], modulus=self.modulus)
        if i < 0:
            if not self.modulus:
                raise NotImplementedError(f"Cannot divide.")
            tmp = self.inv()
        else:
            tmp = self
        for _ in range(i):
            res *= tmp
        return res

    def divmod(self, other: "Poly") -> ("Poly", "Poly"):
        if isinstance(other, list):
            other = self.__class__(other)
        elif not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot divide {self} and {other}.")
        if not other:
            raise ValueError(f"{other} must be nonzero.")
        sd, od = self.degree(), other.degree()
        if sd < od:
            return self.__class__([0]), self
        div = [0] * (sd - od + 1)
        lco = other.coeff[-1]
        if bool(lco - 1):
            tmp = 1 / lco
            oth = [c * tmp for c in other.coeff]
            rem = [c * tmp for c in self.coeff]
        else:
            oth = other.coeff
            rem = [c for c in self.coeff]
        for i in range(sd - od + 1):
            tmp = rem[sd - i]
            div[sd - od - i] = tmp
            for j in range(od + 1):
                rem[sd - i - j] -= tmp * oth[od - j]
        if bool(lco - 1):
            rem = [c * lco for c in rem]
        return self.__class__(div, modulus=self.modulus), self.__class__(
            rem, modulus=self.modulus
        )

    def mod(self, other: "Poly") -> None:
        if isinstance(other, list):
            other = self.__class__(other)
        elif not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot divide {self} and {other}.")
        if not other:
            raise NotImplementedError(f"{other} must be nonzero.")
        sd, od = self.degree(), other.degree()
        if sd < od:
            return self
        lco = other.coeff[-1]
        if bool(lco - 1):
            tmp = 1 / lco
            oth = [c * tmp for c in other.coeff]
        else:
            oth = other.coeff
        for i in range(sd - od + 1):
            tmp = self.coeff[sd - i]
            for j in range(od + 1):
                self.coeff[sd - i - j] -= tmp * oth[od - j]
        for i in range(len(self.coeff) - 1, 0, -1):
            if self.coeff[i]:
                break
            self.coeff.pop(i)

    def inv(self, other: "Poly" = None) -> "Poly":
        if not other:
            other = self.modulus
        if isinstance(other, list):
            other = self.__class__(other)
        elif not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot invert {self} modulo {other}.")
        if not other:
            raise NotImplementedError(f"{other} must be nonzero.")
        r0, r1 = other, self
        y0, y1 = self.__class__([0], modulus=self.modulus), self.__class__(
            [1], modulus=self.modulus
        )
        while r1:
            q, r = r0.divmod(r1)
            r0, r1 = r1, r
            y0, y1 = y1, y0 - q * y1
        if r0.degree() != 0:
            raise ValueError(f"{self} is not invertible mod {other}.")
        tmp = 1 / r0[0]
        for i in range(y0.degree()):
            y0.coeff[i] *= tmp
        return y0
