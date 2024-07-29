"""
Polynomials
"""

from numbers import Number

class Poly:
    """
    Represents a polynomial as a list of coefficients.
    
    Example:
    
    To define a polynomial as a list of coefficients use
    >>> Poly([1, 2, 3])
    3 x^2 + 2 x + 1
    """

    def __init__(self, coeff: list, ring = None, modulus: list = None):
        self.coeff = list(coeff)
        for i in range(len(self.coeff) - 1, 0, -1):
            if self.coeff[i]:
                break
            self.coeff.pop(i)
        self.modulus = modulus
        if ring:
            self.map(ring)
        if modulus:
            self.mod(modulus)

    def __call__(self, x):
        return sum(c * x**j for j, c in enumerate(self.coeff))

    def __getitem__(self, item):
        return self.coeff[item]

    def __setitem__(self, item, value):
        self.coeff[item] = value

    def __len__(self):
        return len(self.coeff)

    def __repr__(self):
        def prx(i: int):
            if i == 0:
                return ""
            if i == 1:
                return "x"
            return "x^" + str(i)

        if len(self.coeff) == 1:
            return str(self.coeff[0])
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

    def degree(self) -> int:
        "Return the degree."
        return len(self.coeff) - 1

    def map(self, func):
        "Apply a given function to all coefficients in place."
        self.coeff = list(map(func, self.coeff))

    def applyfunc(self, func) -> "Poly":
        "Apply a function to all coefficients."
        return self.__class__(list(map(func, self.coeff)), modulus=self.modulus)

    def _check_type(self, other):
        return isinstance(other, int) or (isinstance(other, Number) and isinstance(self.coeff[0], Number)) or type(other) == type(self.coeff[0])

    def __add__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            if self._check_type(other):
                tmp = self.coeff[:]
                tmp[0] += other
                return self.__class__(tmp, modulus=self.modulus)
            return NotImplemented
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

    def __radd__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            if self._check_type(other):
                tmp = self.coeff[:]
                tmp[0] += other
                return self.__class__(tmp, modulus=self.modulus)
        return NotImplemented

    def __neg__(self) -> "Poly":
        return Poly([-s for s in self.coeff], modulus=self.modulus)

    def __pos__(self) -> "Poly":
        return self

    def __sub__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            if self._check_type(other):
                tmp = self.coeff[:]
                tmp[0] -= other
                return self.__class__(tmp, modulus=self.modulus)
            return NotImplemented
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
        return self.__class__([s - o for s, o in zip(scoeff, ocoeff)], modulus = modulus)

    def __rsub__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            if self._check_type(other):
                tmp = self.coeff[:]
                tmp[0] -= other
                return self.__class__(tmp, modulus=self.modulus)
        return NotImplemented

    def __mul__(self, other: "Poly") -> "Poly":
        if not isinstance(other, self.__class__):
            if self._check_type(other):
                return Poly([other * s for s in self.coeff], modulus = self.modulus)
            return NotImplemented
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
        return self.__class__(coeff, modulus = modulus)

    def __rmul__(self, other) -> "Poly":
        if self._check_type(other):
            return self.__class__([other * s for s in self.coeff], modulus=self.modulus)
        return NotImplemented

    def __truediv__(self, other) -> "Poly":
        if self._check_type(other):
            return self.__class__([s / other for s in self.coeff], modulus=self.modulus)
        if isinstance(other, self.__class__):
            if not other.modulus:
                raise NotImplementedError("Cannot invert polynomials without modulus.")
            return self * other.inv()
        return NotImplemented

    def __rtruediv__(self, other) -> "Poly":
        if not self.modulus:
            raise NotImplementedError("Cannot invert polynomials without modulus.")
        return other * self.inv()

    def __pow__(self, j: int) -> "Poly":
        if not isinstance(j, int):
            return NotImplemented
        zero = 0 * self.coeff[0]
        one = zero + 1
        res = self.__class__([one], modulus=self.modulus)
        if j < 0:
            if not self.modulus:
                raise NotImplementedError("Cannot divide polynomials without modulus.")
            tmp = self.inv()
            j *= -1
        else:
            tmp = self
        while j > 0:
            # If j is odd, multiply
            if j & 1:
                res *= tmp
            # Now square
            j = j >> 1  # j= j//2
            tmp *= tmp
        return res

    def max(self) -> int:
        "Maximum of the absolute value of all coefficients."
        return max([abs(c) for c in self.coeff])

    def sum(self) -> int:
        "Sum of the absolute value of all coefficients."
        return sum([abs(c) for c in self.coeff])

    def __floordiv__(self, other: "Poly") -> "Poly":
        return self.divmod(other)[0]

    def __mod__(self, other: "Poly") -> "Poly":
        return self.divmod(other)[1]

    def divmod(self, other: "Poly") -> ("Poly", "Poly"):
        "Polynom division with remainder."
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
        "Reduce with respect to a given polynomial."
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
        "Inverse modulo a given polynomial."
        if not other:
            other = self.modulus
        if isinstance(other, list):
            other = self.__class__(other)
        elif not isinstance(other, self.__class__):
            raise NotImplementedError(f"Cannot invert {self} modulo {other}.")
        if not other:
            raise NotImplementedError(f"{other} must be nonzero.")
        zero = 0 * self.coeff[0]
        one = zero +1
        r0, r1 = other, self
        y0, y1 = self.__class__([zero], modulus=self.modulus), self.__class__([one], modulus=self.modulus)
        while r1:
            q, r = r0.divmod(r1)
            r0, r1 = r1, r
            y0, y1 = y1, y0 - q * y1
        if r0.degree() != 0:
            raise ValueError(f"{self} is not invertible mod {other}.")
        tmp = 1 / r0[0]
        for i in range(len(y0)):
            y0.coeff[i] *= tmp
        return y0
