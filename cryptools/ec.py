"""
Elliptic curves
"""

from math import isqrt, floor, sqrt
from random import randint
from .factor import factorint
from .nt import legendre_symbol, sqrt_mod, crt
from .Zmod import Zmod

class EC_Weierstrass():
    """
    Elliptic curve in Weierstrass normal form y^2 = x^3 + ax + b mod p.
    
    Example:
    
    To define an elliptic curve use
    >>> ec = EC_Weierstrass(239, 3, 1)
    
    To declare a point on the elliptic curve use
    >>> P = ec(43, 222)
    >>> Q = ec(148, 218))
    

    The usual arithmetic operations are supported.
    >>> P + Q
    (195, 41)
    """

    def __init__(self, p: int, a: int, b: int, order: int = None):
        if p < 3:
            raise ValueError(f"{p} must be a prime larger than 2.")
        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p) ) % p == 0:
            raise ValueError(f"Curve is singular!")
        self.p = p
        self.gf = Zmod(p, short = True)
        self.a = self.gf(a % p)
        self.b = self.gf(b % p)
        self.group_order = order
        self.group_order_factors = None
        self.short = False  # display points in short format
        self.hex = False  # display points as hex values in compressed format

    def __call__(self, x: int | str | None = None, y: int | None = None, short: bool = False):
        if isinstance(x, str):
            x = x.replace(" ", "")
            type = x[:2]
            x = x[2:]
            if type == '02' or type == '03':
                short = 1
                x = int(x, 16)
                if type == '02':
                    y = 0
                else:
                    y = 1
            elif type == '04':
                l = len(x)//2
                y = int(x[l:], 16)
                x = int(x[:l], 16)
            else:
                raise ValueError('Unsupported hex type.')
        return ECPoint(x, y, self, short)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.p == other.p and self.a == other.a and self.b == other.b
        return False

    def __contains__(self, P: "ECPoint") -> bool:
        return P.y**2 == P.x**3 + self.a * P.x + self.b

    def info(self):
        format = "Weierstrass curve y^2 = x^3"
        if int(self.a):
            format += f" + {self.a} x"
        if int(self.b):
            format += f" + {self.b}"
        format += f" over Z_{self.p}."
        print(format)

    def add(self, x1, y1, x2, y2):
        if x1 is None:
            return x2, y2
        if x2 is None:
            return x1, y1
        if x1 == x2:
            if y1 == y2:
                return self.dbl(x1, y1)
            return None, None
        s = (y2 - y1) / (x2 - x1)
        x3 = s**2 - x1 - x2
        y3 = s * (x1 - x3) - y1
        return x3, y3

    def dbl(self, x, y):
        if x is None or not y:
            return None, None
        s = (3 * x**2 + self.a) / (2 * y)
        x3 = s**2 - x - x
        y3 = s * (x - x3) - y
        return x3, y3

    def mult(self, j: int, x, y):  # Addition-subtraction ladder
        if j == 0:
            return None, None
        if j < 0:
            y = -y
            j *= -1
        #    xx, yy = None, None
        #    while j > 0:
        #        # If j is odd, add x
        #        if j & 1:
        #            xx, yy = self.add(xx, yy, x, y)
        #        # Now double
        #        j >>= 1  # j= j//2
        #        x, y = self.dbl(x, y)
        xx, yy = x, y
        j3 = 3 * j
        form = "0" + str(j3.bit_length()) + "b"
        for a, b in zip(format(j3, form)[1:-1], format(j, form)[1:-1]):
            xx, yy = self.dbl(xx, yy)
            if a == "1" and b == "0":
                xx, yy = self.add(xx, yy, x, y)
            if a == "0" and b == "1":
                xx, yy = self.add(xx, yy, x, -y)
        return xx, yy

    def random(self):
        j = -1
        while j == -1:
            x = self.gf(randint(0, self.p - 1))
            y2 = int(x**3 + self.a * x + self.b)
            j = legendre_symbol(y2, self.p)
        return ECPoint(x, randint(0, 1), self, short = True)

    def order(self, order: int = None) -> int:
        if order:
            self.group_order = order
        elif not self.group_order:
            if self.p < 230:
                self.group_order = self.order_naive()
            else:
                self.group_order = self.order_shanks_mestre()
        return self.group_order

    def factor_order(self) -> dict:
        if self.group_order_factors:
            return self.group_order_factors
        if not self.group_order:
            self.order()
        self.group_order_factors = factorint(self.group_order)
        return self.group_order_factors

    def order_naive(self) -> int:
        a1 = (int(self.a) + 1) % self.p
        y = int(self.b)
        order = self.p + 1 + legendre_symbol(y, self.p)
        for x in range(self.p - 1):
            y += (3 * x * (x + 1) + a1) % self.p
            order += legendre_symbol(y, self.p)
        return order

    def order_shanks_mestre(self) -> int:
        if self.p < 230:
            return self.order_naive()
        j = 1
        while j != -1:
            g = randint(1, self.p - 1)
            j = legendre_symbol(g, self.p)
        W = floor(sqrt(2 * sqrt(self.p)))
        a, b = int(self.a), int(self.b)
        while True:
            sigma = 0
            while sigma == 0:
                x = randint(0, self.p - 1)
                sigma = legendre_symbol(x**3 + a * x + b, self.p)
            if sigma == 1:
                ec = EC_Weierstrass(self.p, a, b)
            else:
                ec = EC_Weierstrass(
                    self.p, pow(g, 2, self.p) * a, pow(g, 3, self.p) * b
                )
                x = g * x % self.p
            x = ec.gf(x)
            y2 = int(x**3 + ec.a * x + ec.b)
            y = ec.gf(sqrt_mod(y2, self.p))

            A = {}
            xx, yy = ec.mult(ec.p + 1, x, y)
            for j in range(W):
                if xx is None:
                    A[None] = j
                else:
                    A[int(xx)] = j
                xx, yy = ec.add(xx, yy, x, y)
            B = []
            xg, yg = ec.mult(W, x, y)
            xx, yy = None, None
            for j in range(W + 1):
                if xx is None:
                    if xx in A:
                        B.append([A[None], j])
                elif int(xx) in A:
                    B.append([A[int(xx)], j])
                xx, yy = ec.add(xx, yy, xg, yg)
            if len(B) == 1:
                beta, gamma = B[0]
                break
        t = beta + gamma * W
        if ec.mult(ec.p + 1 + t, x, y)[0] is not None:
            t = beta - gamma * W
            assert ec.mult(ec.p + 1 + t, x, y)[0] is None

        return self.p + 1 + sigma * t

class ECPoint:
    "Point on an elliptic curve"
    def __init__(self, x: int, y: int, curve: EC_Weierstrass, short:bool = False):
        self.curve = curve
        if x is None:
            self.x = None
            self.y = None
        else:
            self.x = curve.gf(x)
            if short:
                y2 = int(self.x**3 + curve.a * self.x + curve.b)
                j = legendre_symbol(y2, curve.p)
                if j == -1:
                    raise ValueError("Point not on curve!")
                y1 = sqrt_mod(y2, curve.p)
                if y % 2 == 0:
                    y = y1
                else:
                    y = (curve.p - y1) % curve.p
            self.y = curve.gf(y)
            if not self in curve:
                raise ValueError("Point not on curve!")

    def __repr__(self):
        if self.x is None:
            return "O"
        if self.curve.hex is True:
            if self.curve.short is True:
                pre = '0x02'
                if int(self.y) % 2 == 1:
                    pre = '0x03'
                return pre+f'{int(self.x):x}'
            pre = '0x04'
            return pre+f'{int(self.x):x}{int(self.y):x}'
        if self.curve.short is True:
            return f"({int(self.x)}, {int(self.y) % 2})"
        return f"({int(self.x)}, {int(self.y)})"

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if not self.curve == other.curve:
            return False
        if self.x is None:
            return other.x is None
        if other.x is None:
            return self.x is None
        return self.x == other.x and self.y == other.y

    def __bool__(self):
        return self.x is not None

    def __hash__(self):
        if self.x is None:
            return hash(None)
        return hash((int(self.x), int(self.y)))

    def __add__(self, other: "ECPoint") -> "ECPoint":
        if not isinstance(other, self.__class__) or self.curve != other.curve:
            raise ValueError(f"Cannot add {self} and {other}.")
        x, y = self.curve.add(self.x, self.y, other.x, other.y)
        return ECPoint(x, y, self.curve)

    def __sub__(self, other: "ECPoint") -> "ECPoint":
        return self + other.__neg__()

    def __rmul__(self, scalar: int) -> "ECPoint":
        x, y = self.curve.mult(scalar, self.x, self.y)
        return ECPoint(x, y, self.curve)

    def __neg__(self) -> "ECPoint":
        if self.x is None or not self.y:
            return self
        return ECPoint(self.x, -self.y, self.curve)

    def order(self) -> int:
        """Compute the order of an element."""
        self.curve.factor_order()  # Make sure the factorization is available
        order = self.curve.group_order
        for p, k in self.curve.group_order_factors.items():
            for _ in range(k):
                order_try = order // p
                if self.curve.mult(order_try, self.x, self.y)[0] is None:
                    order = order_try
                else:
                    break
        return order

    def dlog(Q, P: "ECPoint") -> int:
        """Compute the discrete log_P(Q) in EC."""
        m = P.order()
        mf = factorint(m)
        assert m * Q == P.curve(None, None), "DLP not solvable."
        # We first use Pohlig-Hellman to split m into powers of prime factors
        mm = []
        ll = []
        for pj, kj in mf.items():
            Pj = (m // pj**kj) * P
            Qj = (m // pj**kj) * Q
            l = Qj.dlog_ph(Pj, pj, kj)
            if l is None:
                return None
            mm += [pj**kj]
            ll += [l]
        return crt(ll, mm)

    def dlog_ph(Q, P: "ECPoint", q: int, k: int) -> int:
        """Compute the discrete log_P(Q) in EC if P has order q^k using Pohlig-Hellman reduction."""
        if k == 1 or q**k < 10000:
            return Q.dlog_switch(P, q**k)
        Pj = q**(k - 1) * P
        P1 = Pj
        Qj = q**(k - 1) * Q
        xj = Qj.dlog_switch(P1, q)
        for j in range(2, k + 1):
            Pj = q**(k - j) * P
            Qj = q**(k - j) * Q - xj * Pj
            yj = Qj.dlog_switch(P1, q)
            xj = xj + q ** (j - 1) * yj % q**j
        return xj

    def dlog_switch(Q, P: "ECPoint", m: int) -> int:
        """Compute the discrete log_P(Q) in EC if P has order m choosing an appropriate method."""
        if m < 100:
            return Q.dlog_naive(P, m)
        return Q.dlog_bsgs(P, m)

    def dlog_naive(Q, P: "ECPoint", m: int) -> int:
        """Compute the discrete log_P(Q) in EC using an exhaustive search."""
        if not Q.curve == P.curve and not isinstance(Q, P.__class__):
            raise ValueError(f"Points must be on the same curve!")
        j = 0
        xx, yy = None, None
        while xx != Q.x:
            j += 1
            xx, yy = P.curve.add(xx, yy, P.x, P.y)
            if xx is None:
                raise ValueError("DLP not solvabel!")
        if  yy == Q.y:
            return j
        return m - j

    def dlog_bsgs(Q, P: "ECPoint", m: int) -> int:
        """Compute the discrete log_P(Q) in EC if P has order m using Shanks' baby-step-giant-step algorithm."""
        if not Q.curve == P.curve and not isinstance(P, Q.__class__):
            raise ValueError(f"Points must be on the same curve!")
        mm = 1 + isqrt(m - 1)
        m2 = mm//2 + mm % 1  # we use the group symmetry to halve the number of steps
        # initialize baby_steps table
        baby_steps = {}
        baby_step = P
        for j in range(1,m2+1):
            baby_steps[int(baby_step.x)] = j, int(baby_step.y)
            baby_step += P

        # now take the giant steps
        giant_stride = -mm * P
        giant_step = Q
        for l in range(mm+1):
            if giant_step.x == None:
                return l * mm
            if int(giant_step.x) in baby_steps:
                j = baby_steps[int(giant_step.x)][0]
                if int(giant_step.y) != baby_steps[int(giant_step.x)][1]:
                    j *= -1
                return (l * mm + j) % m
            giant_step += giant_stride
        raise ValueError("DLP not solvabel!")
