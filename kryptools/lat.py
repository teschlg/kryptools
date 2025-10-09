"""
Lattice tools
"""

from math import prod, floor, inf
from itertools import product
from fractions import Fraction
from random import choice, sample
from .la import Matrix, eye, zeros


def hermite_nf(M: Matrix) -> Matrix:
    "Compute the Hermite normal form of a matrix M."
    n, m = M.cols, M.rows
    H = M[:, :]
    j = n - 1
    for i in range(m-1,-1,-1):
        j0 = j
        minimum = abs(H[i, j])  # search for the pivot in the present row
        for jj in range(j):
            tmp = abs(H[i, jj])
            if tmp > 0 and (tmp < minimum or minimum == 0):
                minimum = tmp
                j0 = jj
        if minimum == 0:
            continue  # all entries are zero
        if j0 < j:
            H[:, j], H[:, j0] = H[:, j0], H[:, j]  # swap columns, to move the pivot in place
        if H[i, j] < 0:
            H[:, j] *= -1  # make the pivot positive
        jj = j - 1
        while jj >= 0:  # make the row entries left to the pivot zero
            tmp = H[i, jj] // H[i, j]
            H[:, jj] -= tmp * H[:, j]
            if H[i, jj]:
                H[:, j], H[:, jj] = H[:, jj], H[:, j]  # swap columns
            else:
                jj -= 1
        for jj in range(j + 1, n):  # reduce the row entries right to the pivot
            tmp = H[i, jj] // H[i, j]
            H[:, jj] -= tmp * H[:, j]
        j -= 1
        if j < 0:
            break
    while H.cols > 1 and all(not H[i, 0] for i in range(m)):  # remove zero columns
        H = H[:, 1:]
    return H


def norm2(v: Matrix) -> float:
    "Square of the Euclidean norm of a vector v."
    return sum(map(lambda x: x * x, v))


def gram_schmidt(U: Matrix, drop_dependent: bool = True) -> (Matrix, Matrix):
    "Compute the Gram-Schmidt orthogonalization of the column vectors of a matrix M."
    M = U.eye(U.cols)
    Us = U.zeros()
    Us[:, 0] = U[:, 0]
    jj = 0 # offset taking removed vectors into account
    for j in range(1, U.cols):
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
        y = y - round(y.dot(Us[:, k]) / norm2(Us[:, k])) * U[:, k]
    return (x - y).applyfunc(round)


def babai_plane_bnd(U: Matrix, p=2) -> float:
    "Bound for Babai's closest plane algorithm for solving the CVP with respect to the Euclidean norm (p=2) or sup norm (p=inf)."
    Us = gram_schmidt(U)[0]
    return float(0.5 * min(Us[:, i].norm(p) for i in range(Us.rows)))


def lagrange_lr(V: Matrix) -> Matrix:
    "Lagrange lattice reduction."
    if V.cols != 2:
        raise ValueError("Lagrange lattice reduction requires dimension two.")
    v1, v2 = V[:, 0], V[:, 1]
    if norm2(v1) > norm2(v2):
        v1, v2 = v2, v1
    v3 = v2 - round(v1.dot(v2) / norm2(v1)) * v1
    while norm2(v3) < norm2(v1):
        v2, v1 = v1, v3
        v3 = v2 - round(v1.dot(v2) / norm2(v1)) * v1
    return Matrix([list(v1), list(v3)]).transpose()


def lll(V: Matrix, delta: float = 0.75, sort: bool = True) -> Matrix:
    "LLL algorithm for lattice reduction."

    assert 0 < delta <= 1, f"LLL reqires 0 < delta={delta} <= 1"
    j = 1
    U = V[:, :]
    U.map(int)
    Us = U[:, :]
    Us.map(Fraction)
    M = U.zeros(U.cols)
    M.map(Fraction)
    M[0, 0] = norm2(Us[:, 0])  # we store the squared norms on the diagonal
    for l in range(1, U.cols):  # Gram-Schmidt decomposition
        tmp = U[:, l]
        for i in range(l):
            M[i, l] = U[:, l].dot(Us[:, i]) / M[i, i]
            tmp -= M[i, l] * Us[:, i]
        Us[:, l] = tmp
        M[l, l] = norm2(Us[:, l])

    while j < U.cols:
        for i in range(j - 1, -1, -1):  # reduce the weights of the basis vectors
            r = round(M[i, j])
            if r:
                U[:, j] -= r * U[:, i]
                for k in range(j):
                    if k == i:
                        M[k, j] -= r
                    else:
                        M[k, j] -= r * M[k, i]

        newM11 = M[j, j] + M[j - 1, j] ** 2 * M[j - 1, j - 1]
        if (delta * M[j - 1, j - 1] <= newM11):  # Lovasz condition
            j += 1
            continue  # nothing to be done
        # else swap vectors
        U[:, j], U[:, j - 1] = U[:, j - 1], U[:, j]
        # update the Gram-Schmidt decomposition
        oldM11 = M[j - 1, j - 1]
        oldM10 = M[j - 1, j]
        oldM00 = M[j, j]
        oldUs = Us[:, j - 1]
        Us[:, j - 1] = Us[:, j] + M[j - 1, j] * Us[:, j - 1]
        M[j - 1, j - 1] = newM11
        M[j - 1, j] *= oldM11 / M[j - 1, j - 1]
        Us[:, j] = oldUs - M[j - 1, j] * Us[:, j - 1]
        M[j, j] = oldM11 - M[j - 1, j] ** 2 * M[j - 1, j - 1]
        for l in range(j - 1):
            M[l, j], M[l, j - 1] = M[l, j - 1], M[l, j]
        tmp1 = oldM00 / M[j - 1, j - 1]
        tmp2 = oldM10 * oldM11 / M[j - 1, j - 1]
        for l in range(j + 1, U.cols):
            M[j - 1, l], M[j, l] = (
                tmp1 * M[j, l] + tmp2 * M[j - 1, l],
                M[j - 1, l] - oldM10 * M[j, l],
            )
        j = max(j - 1, 1)  # redo the last step

    if sort:  # sort the vectors according to their norm
        tmp = [U[:, j] for j in range(U.cols)]
        tmp.sort(key=norm2)
        for j in range(U.cols):
            U[:, j] = tmp[j]
    return U


def random_unimodular_matrix(n: int, iterations: int = 50, max_val: int = None) -> Matrix:
    "Create a pseudorandom unimodular matrix of dimension n."
    W = zeros(n)
    for i in range(n):
        for j in range(i, n):
            W[i, j] = choice([-1, 1])
    W = W[sample(range(n), n), sample(range(n), n)]
    for _ in range(iterations):
        i, j = sample(range(n), 2)
        tmp = W[i, :] + choice([-1, 1]) * W[j, :]
        if not max_val or max(abs(x) for x in tmp) <= max_val:
            W[i, :] = tmp
        i, j = sample(range(n), 2)
        tmp = W[:, i] + choice([-1, 1]) * W[:, j]
        if not max_val or max(abs(x) for x in tmp) <= max_val:
            W[:, i] = tmp
    return W

def q_ary_lattice(U: Matrix, lll: bool = False) -> Matrix:  # pylint: disable=W0621
    "Create a q-ary lattice and (optinally) LLL reduce the basis."
    if isinstance(U, Matrix) and hasattr(U.matrix[0][0], "ring"):
        q = U.matrix[0][0].ring.n
    else:
        raise ValueError("The matrix does not seem to be over Zmod.")
    V = eye(U.rows, one = q)
    V.append_column(U, ring = int)
    V = hermite_nf(V)
    if lll:
        return globals()['lll'](V)
    return V

def svp_lll(U: Matrix) -> Matrix:
    "Solve the shortest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q."
    gf = U[0].ring
    V = q_ary_lattice(U, lll = True)
    x = V[:, 0]
    x.map(gf)
    return x

def cvp_lll(U: Matrix, x: Matrix) -> Matrix:
    "Solve the closest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q (LWE)."
    gf = U[0].ring
    V = q_ary_lattice(U, lll = True)
    x.map(int)
    y = babai_plane_cvp(x, V)
    x.map(gf)
    return y

def sis_lll(A: Matrix) -> Matrix:
    "Solve the short integer problem over a finite ring Z_q (SIS)."
    return svp_lll(A.kernel())

def isis_lll(A: Matrix, b: Matrix) -> Matrix:
    "Solve the inhomogenous short integer problem over a finite ring Z_q (ISIS)."
    y = A.solve(b)
    return y - cvp_lll(A.kernel(), y)

def svp_search(U: Matrix, m: int = 0, p: int = 2) -> Matrix:
    "Solve the shortest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q."
    gf = U[0].ring
    if m:
        ml, mu = -m, m + 1
    else:
        ml, mu = -((gf.n-1)//2), gf.n//2 + 1
    U.map(int)
    norm_min = inf
    c_min = None
    for j in range(U.cols):
        for c in product(range(ml, mu), repeat = U.cols-j-1):
            for cc in range(1, mu):
                cc = U * Matrix([ 0 ] * j + [ cc ] + list(c))
                cc.map(gf)
                norm = cc.norm(p)
                if norm < norm_min:
                    norm_min = norm
                    c_min = cc
    U.map(gf)
    return c_min

def cvp_search(U: Matrix, b: Matrix, p: int = 2) -> Matrix:
    "Solve the closest vector problem in a q-ary lattice associated with a matrix over a finite ring Z_q (LWE)."
    gf = U[0].ring
    U.map(int)
    b.map(int)
    norm_min = inf
    c_min = None
    for c in product(range(-((gf.n-1)//2), gf.n//2 + 1), repeat = U.cols):
        c = U * Matrix(c)
        d = c - b
        d.map(gf)
        norm = d.norm(p)
        if norm < norm_min:
            norm_min = norm
            c.map(gf)
            c_min = c
    U.map(gf)
    b.map(gf)
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
