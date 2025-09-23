"""
Linear codes
"""

from math import prod
from .Zmod import Zmod
from .GF2 import GF2
from .poly import Poly
from .la import Matrix


def poly2bits(x: Poly) -> list:
    "Convert a polynomial in GF(2^n)[x] to a list of binary digits."
    out = []
    for c in x.coeff:
        out += c.bits()
    return out

def left_standard_form(G: Matrix) -> (Matrix, Matrix):
    "Compute the left standard form of a generator matrix. Return the standard form and the pertmutation matrix."
    return G.left_standard_form()

def gen2pchk(G: Matrix) -> Matrix:
    "Compute the parity check matrix from a given generator matrix (and vice versa)."
    return G.kernel().transpose()

def hamming_dist(a: iter, b: iter):
    "Hamming distance of two iterables."
    if len(a) != len(b):
        raise ValueError("Both iteralbles must have equal length.")
    return sum(aa != bb for aa, bb in zip(a, b))

class Goppa():
    """
    Binary irreducible Goppa code.

    Example:

    >>> gf = GF2(4) # base field GF(2^4)
    >>> g = Poly([2, 0, 0, 1], ring = gf) # irreducible Polynomial for the code
    >>> alpha = list(gf) # list of points from the base field (must not contain zeros of g)
    >>> goppa = Goppa(gf, g, alpha)
    
    To encode a list of bits (the length must be equal to the dimension k of the code)
    >>> goppa.encode([1, 0, 1, 1])
    [0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1]

    Decoding still works if there are at most deg(g) errors:
    >>> goppa.decode([0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0])
    [1, 0, 1, 1]
    
    """

    def __init__(self, gf: GF2, g: Poly, alpha: list["GF2Point"]):
        if not all(map(lambda x: x in gf, alpha)):
            raise ValueError("All elements of alpha must be from the Galois field!")
        if not g.rabin_test():
            raise ValueError("The polynomial g must be irreducible!")
        self.gf = gf  # underlying Galois field
        self.g = g  # polynomial
        self.alpha = alpha # list of elements from gf
        # Compute the control matrix
        self.Z_2 = Zmod(2)

        tmp =[]
        for x in alpha:
            tmp.append(Poly([-x, 1], ring = gf, modulus = g).inv().bits())
        self.H = Matrix(tmp, ring=self.Z_2).transpose()
        # Compute the generator matrix
        self.G = gen2pchk(self.H)

    def __repr__(self):
        return f"Binary irreducible Goppa [{self.G.cols}, {self.G.rows}] code over GF(2^{self.gf.degree}) with polynomial g(x) = {self.g}."

    def encode(self, x: list[int]):
        "Encode a given list of bits."
        x = Matrix(x, ring = self.Z_2)
        y = self.G.transpose() * x
        y.map(int)
        return list(y)

    def decode(self, y: list[int]):
        "Decode a given list of bits."
        y = Matrix(y, ring = self.Z_2)
        zero = self.gf(0)
        one = zero**0
        if self.H * y: # no code word
            s = zero # syndrom polynomial
            for i, x in enumerate(self.alpha):
                s += int(y[i]) * Poly([-x, one], modulus = self.g).inv()
            v = (s.inv() - Poly([zero, one], modulus = self.g))**(self.gf.order ** self.g.degree() // 2)
            # determine error locator polynomial sigma(x) = a(x^2) + b(x^2) * x
            # find a and b by running EEA until b has the desired degree
            v.modulus = None
            r0, r1 = self.g, v
            y0, y1 = Poly([zero]), Poly([one])
            while 2 * r1.degree() > self.g.degree():
                q, r = r0.divmod(r1)
                r0, r1 = r1, r
                y0, y1 = y1, y0 - q * y1
            a, b = r1, y1
            sigma = a**2 + b**2 * Poly([zero, one])
            # correct the error
            for x in self.alpha:
                if not sigma(x):
                    y[int(x)] += 1
            if self.H * y:
                raise ValueError("Decoding failed!")
        x = self.G.transpose().solve(y)
        x.map(int)
        return list(x)

class CyclicCode():
    """
    Cyclic code.

    Example:

    >>> n = 4 # length of the code
    >>> g = Poly([1, 1], ring = Zmod(2)) # generating Polynomial for the code
    >>> cc = CyclicCode(n, g)
    
    To encode a list of letters (the length must be equal to the dimension k of the code):
    >>> cc.encode([0,1,0], ring = Zmod(2))
    [0, 1, 1, 0]

    Decoding is only implemnted for code words:
    >>> cc.decode([0, 1, 1, 0], ring = Zmod(2))
    [0, 1, 0]
    
    """

    def __init__(self, n:int, g: Poly):
        if n < 2:
            raise ValueError(f"Code lenght {n} must be at least 2!")
        zero = 0 * g[0]
        if hasattr(zero, "ring"):
            order = zero.ring.n  # Zmod
        elif hasattr(zero, "field"):
            order = zero.field.order  # GF2
        elif hasattr(type(zero), "order"):
            order = type(zero).order  # galois
        else:
            raise ValueError("Unknown base field!")
        self.order = order  # order of the base field
        one = zero**0
        modulus = Poly([-one] + [zero] * (n-1) + [one])
        self.modulus = modulus
        if modulus % g:
            raise ValueError(f"The polynomial g must be a factor of x^{n} - 1!")
        self.n = n  # lenght
        self.k = n - g.degree()  # dimension
        self.g = g  # generating polynomial
        self.h = Poly(modulus) // g

    def __repr__(self):
        return f"Cyclic [{self.n}, {self.k}] code over GF({self.order}) with generator polynomial g(x) = {self.g}."

    def encode(self, x: list|Poly, gf = None, systematic_form = False):
        "Encode a given list or polynomial."
        as_list = False
        zero = 0 * self.g[0]
        if isinstance(x, list | tuple):
            as_list = True
            x = Poly(x, ring = gf)
        if x.degree() >= self.k:
            raise ValueError(f"Degree can be at most {self.k-1}.")
        if systematic_form:
            x.coeff = [zero] * self.g.degree() + x.coeff
            y = x - (x % self.g)
        else:
            y = self.g * x
        if as_list:
            return y.coeff + [zero] * (self.n - len(y.coeff))
        return y

    def decode(self, y: list|Poly, gf = None, systematic_form = False):
        "Decode a given list or polynomial."
        as_list = False
        zero = 0 * self.g[0]
        if isinstance(y, list | tuple):
            as_list = True
            y = Poly(y, ring = gf)
        if y.degree() >= self.n:
            raise ValueError(f"Degree can be at most {self.n-1}.")
        x, r = y.divmod(self.g)
        if r:
            y = self.correct(y)
            x, r = y.divmod(self.g)
        if systematic_form:
            x.coeff = y.coeff[self.g.degree():]
        if as_list:
            return x.coeff + [zero] * (self.k - len(x.coeff))
        return x

    def correct(self, y: Poly):
        "Correct a given polynomial."
        raise NotImplementedError("Error correction not implemented.")

    def generator_matrix(self):
        "Return the generator matrix of the code."
        zero = 0 * self.g[0]
        return Matrix([[zero] * j + self.g.coeff + [zero] * (self.k - j - 1 ) for j in range(self.k)])

    def check_matrix(self):
        "Return the parity check matrix of the code."
        zero = 0 * self.g[0]
        tmp = list(reversed(self.h.coeff))
        return Matrix([[zero] * j + tmp + [zero] * (self.n - self.k - j - 1 ) for j in range(self.n-self.k)])

class ReedSolomonCode(CyclicCode):
    """
    Reed-Solomon code.

    Example:

    >>> gf = Zmod(13) # Galois field
    >>> k = 3 # dimension of the code
    >>> rsc = ReedSolomonCode(k, gf)
    
    To encode a list of letters (the length must be equal to the dimension k of the code; conversion to the Galois field is done automatically):
    >>> rsc.encode([1, 2, 3])
    [6, 4, 5, 1, 8, 4, 2, 9, 2, 8, 9, 6]

    To decode use:
    >>> rsc.decode([1, 4, 4, 0, 8, 4, 2, 9, 2, 0, 9, 6])
    [1, 2, 3]
    
    """
    def __init__(self, k:int, gf, alpha = None):  # pylint: disable=W0231
        zero = gf(0)
        one = zero**0
        self.gf = gf
        if hasattr(gf, "n"):
            n = gf.n - 1  # Zmod
        else:
            n = gf.order - 1  # GF2|galois
        self.order = n + 1  # order of the base field
        self.n = n  # lenght
        self.k = k  # dimension
        if k >= n:
            raise ValueError(f"Code dimension k = {k} must be smaller than the code lenght n = {self.n}.")
        self.d = n - k + 1  # minimal distance
        self.modulus = Poly([-one] + [zero] * (n-1) + [one])
        if alpha is None:
            if hasattr(gf, "primitive_element"):
                alpha = gf.primitive_element  # galois
            else:
                alpha = gf.generator()  # Zmod|GF2
        self.alpha = alpha
        self.g = prod([Poly([-alpha**j, one]) for j in range(1,n-k+1)], start = one) # generating polynomial
        self.h = Poly(self.modulus) // self.g  # check polynomial
        assert k == n - self.g.degree()
        assert self.g * self.h == self.modulus

    def __repr__(self):
        return f"Reed-Solomon [{self.n}, {self.k}, {self.d}] code over GF({self.n+1})."

    def encode(self, x: list|Poly) -> list|Poly:  # pylint: disable=W0221
        "Encode a given list or polynomial."
        as_list = False
        zero = 0 * self.g[0]
        if isinstance(x, list | tuple):
            as_list = True
            x = Poly(x, ring = self.gf)
        if x.degree() >= self.k:
            raise ValueError(f"Degree can be at most {self.k-1}.")
        y = [ x(self.alpha**j) for j in range(self.n) ]
        if as_list:
            return y + [zero] * (self.n - len(y))
        return Poly(y)

    def decode(self, y: list|Poly) -> list|Poly:  # pylint: disable=W0221
        "Decode a given list or polynomial."
        as_list = False
        zero = 0 * self.g[0]
        one = zero**0
        if isinstance(y, list | tuple):
            as_list = True
            y = Poly(y, ring = self.gf)
        if y.degree() >= self.n:
            raise ValueError(f"Degree can be at most {self.n-1}.")
        xx = [ -y(self.alpha**(-j)) for j in range(self.n) ]
        x = Poly(xx)
        if x.degree() >= self.k: # No codewort
            t = (self.n - self.k) // 2
            if not t:
                raise ValueError("This code cannot correct any errors.")
            A = Matrix([[ xx[j - l % self.n] for l in range(1, t+1) ] for j in range(self.k+t, self.k+2*t) ])
            b = [ -xx[self.k + t + l ] for l in range(t) ]
            lam = Poly([one] + list(A.solve(b)))
            lamp = lam.derivative()
            s = Poly([x[self.k + j] for j in range(self.n - self.k)])
            omega = s * lam
            omega = Poly(omega.coeff[:self.n - self.k])  # mod x^(n-k)
            # Chien search
            lam = lam.coeff
            alphai = [ self.alpha**i for i in range(len(lam)) ]
            for i in range(self.n):
                if not sum(lam, start = zero):
                    # Fornay
                    ali = self.alpha**i
                    y[i] -= ali**(self.k-1) * omega(ali) / lamp(ali)
                for j in range(len(lam)):  # pylint: disable=C0200
                    lam[j] *= alphai[j]
            x = Poly([ -y(self.alpha**(-j)) for j in range(self.n) ])
        if as_list:
            return x.coeff + [zero] * (self.k - len(x.coeff))
        return x

    def correct(self, y: Poly):
        "Correct a given polynomial."
        x = self.decode(y)
        return self.encode(x)

class BCHCode(CyclicCode):
    """
    Binary narrow sense primitive BCH code.

    Example:

    >>> n = 7 # lenght of the code
    >>> D = 2 # desired distance
    >>> bch = BCHCode(k, D)
    
    To encode a list of letters (the length must be equal to the dimension k of the code; conversion to the Galois field is done automatically):
    >>> bch.encode([0, 1, 0, 1])
    [0, 1, 1, 1, 0, 0, 1]

    To decode use:
    >>> bch.decode([0, 1, 1, 1, 0, 0, 1])
    [0, 1, 0, 1]
    
    """
    def __init__(self, n:int, D: int):
        m = n.bit_length()
        if 2**m - 1 != n:
            raise ValueError(f"The codelenght n={n} plus one must be a power of 2.")
        self.rsc = ReedSolomonCode(n - D + 1, GF2(m))
        # Compute the generator polynomial
        aj = self.rsc.alpha
        g = aj.minpoly()
        for _ in range(2,D):
            aj *= self.rsc.alpha
            mj = aj.minpoly()
            g = g.lcm(mj)
        g.map(GF2(1))
        super().__init__(n, g)

    def __repr__(self):
        return f"BCH [{self.n}, {self.k}] code over GF({self.order}) with generator polynomial g(x) = {self.g}."

    def correct(self, y: Poly):
        "Correct a given polynomial."
        y.map(self.rsc.gf)
        x = self.rsc.decode(y)
        y = self.rsc.encode(x)
        y.map(GF2(1))
        return y
