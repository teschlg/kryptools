"""
Lattice tools 2
"""

from math import prod, floor, inf
from itertools import product
from .la import Matrix, eye
from .lat import Lattice


def hermite_nf(M: Matrix) -> Matrix:
    "Compute the Hermite normal form of a matrix M."
    return M.hnf()


def gram_schmidt(U: Matrix, drop_dependent: bool = True) -> (Matrix, Matrix):
    "Compute the Gram-Schmidt orthogonalization of the column vectors of a matrix M."
    M = U.eye(U.cols)
    Us = U.zeros()
    jj = 0 # offset taking removed vectors into account
    for j in range(0, U.cols):
        tmp = U[:, j]
        for i in range(j - jj):
            M[i, j - jj] = U[:, j].dot(Us[:, i]) / Us[:, i].norm2()
            tmp -= M[i, j - jj] * Us[:, i]
        if not tmp:
            if not drop_dependent:
                raise ValueError("Vectors are linearly dependent.")
            jj += 1
        else:
            Us[:, j - jj] = tmp
    if jj == U.cols:
        raise ValueError("The matrix must be nonzero.")
    if jj:
        Us = Us[:,:-jj]
        M = M[:-jj,:-jj]
    return Us, M


def gram_det(U: Matrix) -> float:
    "Compute the Gram determinant of a matrix."
    Us = gram_schmidt(U, drop_dependent = False)[0]
    return prod([Us[:, i].norm() for i in range(U.rows)])


def hadamard_ratio(M: Matrix) -> float:
    "Compute the Hadamard ratio of a matrix."
    m = M.rows
    return (gram_det(M) / prod([M[:, i].norm() for i in range(m)])) ** (1 / m)


def babai_round_cvp(x: Matrix, U: Matrix) -> Matrix:
    "Babai's rounding algorithm for approximately solving the CVP."
    s = U.inv() * x
    k = s.applyfunc(round)
    return U * k

def babai_round_bnd(U: Matrix) -> float:
    "Bound for Babai's rounding algorithm for solving the CVPwith respect to the sup norm."
    return floor(1 / (2 * max(U.inv()[i, :].norm(1) for i in range(U.rows))))


def babai_plane_cvp(x: Matrix, U: Matrix) -> Matrix:
    "Babai's closest plane algorithm for approximately solving the CVP."
    Us = gram_schmidt(U)[0]
    y = x
    for k in range(U.cols - 1, -1, -1):
        y = y - round(y.dot(Us[:, k]) / Us[:, k].norm2()) * U[:, k]
    return (x - y).applyfunc(round)


def babai_plane_bnd(U: Matrix, p=2) -> float:
    "Bound for Babai's closest plane algorithm for solving the CVP with respect to the Euclidean norm (p=2) or sup norm (p=inf)."
    Us = gram_schmidt(U)[0]
    return float(0.5 * min(Us[:, i].norm(p) for i in range(Us.rows)))

def lll(V: Matrix, delta: float = 0.75, sort: bool = True) -> Matrix:
    "LLL algorithm for lattice reduction."
    lat = Lattice(V)
    lat.lll(delta = delta)
    if sort:
        lat.sort()
    return lat.basis()

def q_ary_lattice(U: Matrix, lll: bool = False) -> Matrix:  # pylint: disable=W0621
    "Create a q-ary lattice and (optinally) LLL reduce the basis."
    if isinstance(U, Matrix) and hasattr(U.matrix[0][0], "ring"):
        q = U.matrix[0][0].ring.n
    else:
        raise ValueError("The matrix does not seem to be over Zmod.")
    V = eye(U.rows, one = q)
    V.append_column(U, ring = int)
    V = V.hnf()
    if lll:
        return globals()['lll'](V)
    return V

def svp(U: Matrix, method: str = 'lll', delta: float = 0.75) -> Matrix:
    "Solve the shortest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q."
    ring = U[0].ring
    lat = Lattice(U)
    x = lat.svp(method = method, delta = delta)
    x.map(ring)
    return x

def cvp(U: Matrix, x: Matrix, method: str = "babai_plane") -> Matrix:
    "Solve the closest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q (LWE)."
    ring = U[0].ring
    lat = Lattice(U)
    x.map(int)
    y = lat.cvp(x, method = method)
    x.map(ring)
    y.map(ring)
    return y

def sis(A: Matrix, method: str = 'all', delta: float = 0.75) -> Matrix:
    "Solve the short integer problem over a finite ring Z_q (SIS)."
    return svp(A.kernel(), method = method, delta = delta)

def isis(A: Matrix, b: Matrix, method: str = 'babai_plane') -> Matrix:
    "Solve the inhomogenous short integer problem over a finite ring Z_q (ISIS)."
    y = A.solve(b)
    if y is None:
        return None
    return y - cvp(A.kernel(), y, method = method)

def svp_search(U: Matrix, m: int = 0, p: int = 2) -> Matrix:
    "Solve the shortest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q."
    ring = U[0].ring
    if m:
        ml, mu = -m, m + 1
    else:
        ml, mu = -((ring.n-1)//2), ring.n//2 + 1
    U = U.applyfunc(int)
    norm_min = inf
    c_min = None
    for j in range(U.cols):
        for c in product(range(ml, mu), repeat = U.cols-j-1):
            for cc in range(1, mu):
                cc = U * Matrix([ 0 ] * j + [ cc ] + list(c))
                cc.map(ring)
                norm = cc.norm(p)
                if norm and norm < norm_min:
                    norm_min = norm
                    c_min = cc
    return c_min

def cvp_search(U: Matrix, b: Matrix, p: int = 2) -> Matrix:
    "Solve the closest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q (LWE)."
    ring = U[0].ring
    U = U.applyfunc(int)
    b = b.applyfunc(int)
    norm_min = inf
    c_min = None
    for c in product(range(-((ring.n-1)//2), ring.n//2 + 1), repeat = U.cols):
        c = U * Matrix(c)
        d = c - b
        d.map(ring)
        norm = d.norm(p)
        if norm < norm_min:
            norm_min = norm
            c.map(ring)
            c_min = c
    return c_min

def sis_search(A: Matrix, p: int = 2) -> Matrix:
    "Solve the short integer problem over a finite ring Z_q (SIS)."
    return svp_search(A.kernel(), p = p)

def isis_search(A: Matrix, b: Matrix, p: int = 2) -> Matrix:
    "Solve the inhomogenous short integer problem over a finite ring Z_q (ISIS)."
    y = A.solve(b)
    if y is None:
        return None
    return y - cvp_search(A.kernel(), y, p = p)
