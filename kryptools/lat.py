"""
Lattice tools
"""

from math import prod, inf, sqrt
from numbers import Number
from fractions import Fraction
from random import choice, sample
from .la import Matrix, eye, zeros



class Lattice():
    """
    Lattice Object

    Example:

    A lattice is defined with a matrix containing the basis vectors as columns
    >>> lat = Lattice(Matrix([[1, 2], [3, 4]]))
    >>> lat.lll()
    >>> lat  
    [ 1, -1 ]
    [ 1,  1 ]
    """

    def __init__(self, U: Matrix, fraction: bool = True):
        if not isinstance(U, Matrix):
            raise ValueError("The basis must be given as a matrix containing the basis vectors as columns.")
        if hasattr(U.matrix[0][0], "ring"):  # Create a q-ary lattice
            q = U.matrix[0][0].ring.n
            V = eye(U.rows, one = q)
            V.append_column(U, ring = int)
            U = V.hnf()
        self.U = U.transpose().matrix  # the basis
        self.len = U.cols  # the number of basis vectors
        self.dim = U.rows  # the dimension of the space
        self.rank = U.rank()  # the dimension of the lattice
        self.Us = None  # Gram-Schmidt basis
        self.M = None  # Gram-Schmidt weights
        self.N = None  # Gram-Schmidt squared norms
        self.Mu = None  # Weights
        self.Nu = None  # Squared norms
        self.Uhnf = None  # Hermite normal form for checking equality
        if U.is_integer():
            self.fraction = fraction
        else:
            self.fraction = False

    def __repr__(self) -> str:
        return repr(self.basis())

    def _repr_mimebundle_(self, **kwargs):  # pylint: disable=W0613
        return self.basis()._repr_mimebundle_(**kwargs)  # pylint: disable=W0212

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return False
        if self.Uhnf is None:
            self.Uhnf = self.basis().hnf().matrix
        if other.Uhnf is None:
            other.Uhnf = other.basis().hnf().matrix
        return self.Uhnf == other.Uhnf

    def __call__(self, k: "Matrix") -> "Matrix":
        if isinstance(k, list|tuple):
            k = Matrix(k)
        return self.basis() * k

    def __contains__(self, x: "Matrix") -> bool:
        c = self.coordinates(x)
        if c is None:
            return False
        return c.is_integer()

    def coordinates(self, x: "Matrix") -> Matrix:
        "Return the lattice coordinates of a point."
        U = self.basis()
        if self.fraction:
            U.map(Fraction)
        c = U.solve(x)
        if c is None:
            return None
        if c.is_integer():
            c.map(int)
        return c

    def basis(self) -> Matrix:
        "Return the current basis."
        return Matrix(self.U).transpose()

    def dot(self, u: list, v: list) -> Number:
        "Dot product of two vectors u and v."
        return sum(x * y for x, y in zip(u, v))

    def norm2(self, u: list) -> int:
        "Square of the Euclidean norm of a vector u."
        return sum(map(lambda x: x * x, u))

    def delete_vector(self, l: int):
        "Delete a vector from the basis."
        if self.len == 1:
            raise ValueError("Cannot delete all basis vectors!")
        del self.U[l]
        self.len -= 1
        if self.M is not None:
            del self.M[l]
            for row in self.M:
                del row[l]
            del self.N[l]
            self.gsd(l)
        if self.Mu is not None:
            del self.Mu[l]
            for row in self.Mu:
                del row[l]
        if self.Nu is not None:
            del self.Nu[l]

    def compute_weights(self, start: int = 0) -> None:
        "Compute the weights of all basis vectors."
        if self.Mu is None:
            self.Mu = [[ int(i == j) for i in range(self.len) ] for j in range(self.len) ]
            self.Nu = [ 0 ] * self.len  # squared norms
            start = 0
        for l in range(start, self.len):
            self.Nu[l] = self.norm2(self.U[l])
            if self.fraction:
                self.Nu[l] = Fraction(self.Nu[l])
            if not self.Nu[l]:
                self.delete_vector(l)
        for l in range(start, self.len):
            for i in range(start, self.len):
                if i != l:
                    self.Mu[l][i] = self.dot(self.U[l], self.U[i]) / self.Nu[i]

    def weight_reduce(self) -> bool:
        "Pairwise reducte all weights of the basis vectors."
        reduced = False
        j = 1
        while j < self.len:
            for i in range(j - 1, -1, -1):
                r = round(self.Mu[j][i])
                if not r:
                    continue
                reduced = True
                for l in range(self.dim):
                    self.U[j][l] -= r * self.U[i][l]
                for k in range(self.len):
                    if j != k:
                        self.Mu[j][k] -= r * self.Mu[i][k]
                Nuj = self.norm2(self.U[j])
                if not Nuj:
                    self.delete_vector(j)
                    break
                if self.fraction:
                    Nuj = Fraction(Nuj)
                for k in range(self.len):
                    if j != k:
                        self.Mu[k][j] = (self.Nu[j] * self.Mu[k][j] - r * self.Nu[i] * self.Mu[k][i]) / Nuj
                self.Nu[j] = Nuj
            else:
                j += 1
        return reduced

    def hermite(self) -> None:
        "Hermite algorithm for a weight reduced basis."
        self.Us = None
        if self.Mu is None:
            self.compute_weights()
        reduced = True
        while reduced:
            self.sort()
            reduced = self.weight_reduce()

    def hnf(self) -> None:
        "Transforms the basis into Hermite normal form.."
        if self.Uhnf is None:
            self.Uhnf = self.basis().hnf().matrix
            self.Us = None
            self.Mu = None
            self.Nu = None
        self.U = [list(i) for i in zip(*self.Uhnf)]
        self.len = len(self.U)

    def gsd(self, start: int = 0) -> None:
        "Gram-Schmit decomposition."
        if self.Us is None:
            self.Us = [[ 0 for _ in range(self.dim) ] for _ in range(self.len) ]
            self.M = [[ int(i == j) for i in range(self.len) ] for j in range(self.len) ]
            self.N = [ 0 ] * self.len  # squared norms
            start = 0
        for l in range(start, self.len):  # Gram-Schmidt decomposition
            self.Us[l] = self.U[l][:]
            for i in range(l):
                self.M[l][i] = self.dot(self.U[l], self.Us[i]) / self.N[i]
                for j in range(self.dim):
                    self.Us[l][j] -= self.M[l][i] * self.Us[i][j]
            self.N[l] = self.norm2(self.Us[l])
            if self.fraction:
                self.N[l] = Fraction(self.N[l])
            if not self.N[l]:
                raise ValueError("Basis is not linearly independent!")

    def size_reduce(self) -> bool:
        "Size reduction step at all columns."
        reduced = False
        for j in range(1, self.len):
            reduced = self.size_reduce_column(j)
        return reduced

    def size_reduce_column(self, j:int) -> bool:
        "Size reduction step at column j."
        reduced = False
        for i in range(j - 1, -1, -1):  # reduce the weights of the basis vectors
            r = round(self.M[j][i])
            if r:
                reduced = True
                for l in range(self.dim):
                    self.U[j][l] -= r * self.U[i][l]
                for k in range(i+1):
                    self.M[j][k] -= r * self.M[i][k]
        return reduced

    def swap(self, j:int, newN1 = None) -> None:
        "Swap columns j and j-1."
        self.U[j - 1], self.U[j] = self.U[j], self.U[j - 1]
        # update the Gram-Schmidt decomposition
        if newN1 is None:
            newN1 = self.N[j] + self.M[j][j - 1]**2 * self.N[j - 1]
        oldN1 = self.N[j - 1]
        oldM10 = self.M[j][j - 1]
        oldN0 = self.N[j]
        oldUs = self.Us[j - 1][:]
        for l in range(self.dim):
            self.Us[j - 1][l] = self.Us[j][l] + self.M[j][j - 1] * self.Us[j - 1][l]
        self.N[j - 1] = newN1
        self.M[j][j - 1] *= oldN1 / newN1
        for l in range(self.dim):
            self.Us[j][l] = oldUs[l] - self.M[j][j - 1] * self.Us[j - 1][l]
        self.N[j] = oldN1 - self.M[j][j - 1]**2 * self.N[j - 1]
        for l in range(j - 1):
            self.M[j][l], self.M[j - 1][l] = self.M[j - 1][l], self.M[j][l]
        tmp1 = oldN0 / self.N[j - 1]
        tmp2 = oldM10 * oldN1 / self.N[j - 1]
        for l in range(j + 1, self.len):
            self.M[l][j - 1], self.M[l][j] = tmp1 * self.M[l][j] + tmp2 * self.M[l][j - 1],  self.M[l][j - 1] - oldM10 * self.M[l][j]


    def lll(self, delta: float = 0.75, sort = True) -> None:
        "LLL algorithm for lattice reduction."
        self.Mu = None
        self.Nu = None
        if not 0 < delta <= 1:
            raise ValueError(f"LLL reqires 0 < delta={delta} <= 1")
        if self.Us is None:
            self.gsd()
        j = 1
        while j < self.len:
            self.size_reduce_column(j)
            newN1 = self.N[j] + self.M[j][j - 1]**2 * self.N[j - 1]  # new norm N[j-1] if we swap
            if delta * self.N[j - 1] <= newN1:  # Lovasz condition
                j += 1  # move on
                continue
            self.swap(j, newN1) # swap vectors
            j = max(j - 1, 1)  # redo the last step
        if sort:
            self.sort()

    def sort(self) -> bool:
        "Sort according to length."
        if self.Nu is None:
            self.Nu = list(map(self.norm2, self.U))
        sort_idx = list(range(self.len))
        sort_idx.sort(key=lambda j: self.Nu[j])
        # which vectors are unchanged
        for j in range(self.len):
            if sort_idx[j] != j:
                break
        if j == self.len - 1:
            return True
        self.U = [self.U[i] for i in sort_idx]
        self.Nu = [self.Nu[i] for i in sort_idx]
        if self.Us is not None:
            self.gsd(start = j)
        if self.Mu is not None:
            self.Mu = [[self.Mu[i][j] for j in sort_idx] for i in sort_idx]
        return True

    def gram_det(self) -> float:
        "Compute the Gram determinant of the current basis."
        if self.Us is None:
            self.gsd()
        return sqrt(prod(self.N))

    def hadamard_ratio(self) -> float:
        "Compute the Hadamard ratio of the current basis."
        return (self.gram_det() / prod(map(self.norm2, self.U)) ) ** (1 / self.len)

    def svp(self, method: str = 'all', delta: float = 0.75) -> Matrix:
        "Solve the SVP (approximately) using a given method."
        method = method.lower()
        if method not in ('all', 'hermite', 'lll', 'search'):
            raise ValueError("Supported methods are: all, hermite, lll, search")
        if method in ('all', 'hermite'):
            self.hermite()
        if self.len == self.rank and method in ('all', 'lll'):
            self.lll(delta = delta)
        if method == 'search':
            return self.svp_search()
        self.sort()
        return Matrix(self.U[0])

    def svp_search(self) -> (float, Matrix):
        "Solve the SVP exactly using branch & bound."
        # wet try to improve the basis as much as possible
        self.hermite()
        if self.len == self.rank:
            self.lll(delta = 1)
        lam = self.N[0]  # current guess for the shortest length
        k, lam = self.search_bb([], lam , 0)
        return self.basis() * Matrix(k)

    def search_bb(self, k:list, lam: float, C: float, beta: list|None = None, b: Matrix|None = None) -> (list, float):
        "Branch and bound step for SVP and CVP."
        j = self.len - len(k) - 1
        if j == -1:
            if beta is None:
                lam = (self.basis() * Matrix(k)).norm2()
                if not lam:
                    return k, inf
            else:
                lam = (self.basis() * Matrix(k) - b).norm2()
            return k, lam
        cj = sum( (k[l] * self.M[j+l+1][j] for l in range(len(k))))
        if beta is not None:
            cj -= beta[j]
        k_best= None
        normj = self.N[j]
        delta = 0
        bound = sqrt( max(0,lam - C) / normj )
        if beta is None and all( kj == 0 for kj in k):
            alternate = False # test only non-negative values
        else:
            alternate = True # test both signs
        while abs(delta) <= bound + 0.51:
            kj = - round(cj) + delta
            Cj = C + (kj + cj)**2 * normj
            if abs(kj + cj) <= 1.01 * bound and Cj <= 1.01 * lam:
                k_try = [kj] + k
                k_new, lam_new = self.search_bb(k_try, lam, Cj, beta, b)
                if lam_new < lam or (k_best is None and lam_new == lam):
                    lam = lam_new
                    k_best = k_new
                    bound = sqrt( max(0,lam - C) / normj )
            if alternate:
                if delta > 0:
                    delta = -delta
                else:
                    delta = -delta + 1
            else:
                delta += 1
        return k_best, lam

    def cvp(self, x: Matrix, method: str = 'babai_plane', hermite: bool = True, lll: bool = True, delta: float = 0.75) -> Matrix:
        "Solve the CVP (approximately) using a given method."
        if hermite:
            self.hermite()
        if self.len == self.rank and lll:
            self.lll(delta = delta)
        method = method.lower()
        if method == 'lll':
            method = 'kannan'
        if hasattr(self, "cvp_" + method):
            method = getattr(self, "cvp_" + method)
        else:
            raise ValueError("Supported methods are: babai_round, babai_plane, kannan, search")
        return method(x)

    def cvp_babai_round(self, x: Matrix) -> Matrix:
        "Babai's rounding algorithm for approximately solving the CVP."
        c = self.coordinates(x)
        if c is None:
            return None
        return self.basis() * c.applyfunc(round)

    def cvp_babai_plane(self, x: Matrix) -> Matrix:
        "Babai's closest plane algorithm for approximately solving the CVP."
        if self.Us is None:
            self.gsd()
        y = list(x)
        for k in range(self.len - 1, -1, -1):
            m = round( self.dot(y, self.Us[k]) / self.N[k] )
            for l in range(self.dim):
                y[l] -= m * self.U[k][l]
        y = Matrix(y)
        return (x - y).applyfunc(round)

    def cvp_kannan(self, x: Matrix, m: int = 1, delta: float = 0.75) -> Matrix:
        "Kannan's embedding algorithm for approximately solving the CVP."
        V = self.basis()
        V.append_column(x)
        V.append_row(V.zeros(1,V.cols))
        V[-1,-1]= m
        lat = Lattice( V )
        lat.lll(delta = delta)
        e = None
        for i in range(lat.len):
            mm = lat.U[i][-1]
            if mm == m:
                e = Matrix(lat.U[i][:-1])
                break
            if mm == -m:
                e = - Matrix(lat.U[i][:-1])
                break
        if e is None:
            return None
        return self.basis() * self.coordinates(x - e)

    def cvp_search(self, x: Matrix) -> Matrix:
        "Solve the CVP exactly using branch & bound."
        # wet try to improve the basis as much as possible
        self.hermite()
        lam = inf
        if self.len == self.rank:
            self.lll(delta = 1)
            a = self.cvp_kannan(x)
            if a is not None:
                lam = (x - self.cvp_kannan(x)).norm2()  # current guess for the distance
        xi = list(Matrix(self.Us) * x)
        xi = [ x / n for x, n in zip(xi, self.N) ]
        k, lam = self.search_bb([], lam , 0, xi, x)
        return self.basis() * Matrix(k)


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
