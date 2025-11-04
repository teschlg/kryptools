"""
Lattice tools
"""

from math import prod, inf, sqrt
from numbers import Number
from fractions import Fraction
from random import choice, sample
from .la import Matrix, eye, zeros
from .nt import egcd


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
        if isinstance(U, Matrix):
            if hasattr(U.matrix[0][0], "ring"):  # Create a q-ary lattice
                q = U.matrix[0][0].ring.n
                V = eye(U.rows, one = q)
                V.append_column(U, ring = int)
                U = V
            self.U = [row for row in U.transpose().matrix if any(x for x in row)]  # the basis
            if U.is_integer():
                self.fraction = fraction
            else:
                self.fraction = False
        elif isinstance(U, list):
            self.U = U
            self.fraction = fraction
        else:
            raise ValueError("The basis must be given as a matrix containing the basis vectors as columns.")
        self.len = len(self.U)  # the number of basis vectors
        self.dim = len(self.U[0])  # the dimension of the space
        self.myrank = None  # the dimension of the lattice
        self.Us = None  # Gram-Schmidt basis
        self.Us_end = 0  # Gram-Schmidt decomposition is correct below this index
        self.M = None  # Gram-Schmidt weights
        self.N = None  # Gram-Schmidt squared norms
        self.Mu = None  # Weights
        self.Nu = None  # Squared norms
        self.Uhnf = None  # Hermite normal form for checking equality

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
        rank = self.rank()
        if rank < self.len:
            if self.Uhnf is None:
                self.Uhnf = self.basis().hnf().matrix
            U = Matrix(self.Uhnf)
        else:
            U = self.basis()
        if self.fraction:
            U.map(Fraction)
        c = U.solve(x)
        if c is None:
            return False
        return c.is_integer()

    def coordinates(self, x: "Matrix") -> Matrix:
        "Return the lattice coordinates of a point."
        rank = self.rank()
        if rank < self.len:
            raise NotImplementedError("Computation of coordinates requires a linear independent basis.")
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

    def rank(self) -> Matrix:
        "Return the rank."
        if self.myrank is None:
            if self.Us is not None and self.Us_end == self.len - 1:
                self.myrank = sum(bool(i) for i in self.N)
            else:
                self.myrank = Matrix(self.U).rank()
        return self.myrank

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
        if self.Us is not None:
            del self.Us[l]
            del self.M[l]
            for row in self.M:
                del row[l]
            del self.N[l]
            if l < self.Us_end:
                self.Us_end -= 1
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

    def hermite(self, verbose = False) -> None:
        "Hermite algorithm for a weight reduced basis."
        self.Us = None
        if self.Mu is None:
            self.compute_weights()
        if verbose:
            print("Starting Hermite")
        steps = 0
        reduced = True
        while reduced:
            steps += 1
            self.sort()
            reduced = self.weight_reduce()
        if verbose:
            print(f"Total number of steps: {steps}")

    def hnf(self) -> None:
        "Transform the basis into Hermite normal form.."
        if self.Uhnf is None:
            self.Uhnf = self.basis().hnf().matrix
            self.Us = None
            self.Mu = None
            self.Nu = None
        self.U = [list(i) for i in zip(*self.Uhnf)]
        self.len = len(self.U)

    def gsd(self, end: int = -1) -> None:
        "Compute the Gram-Schmit decomposition up to (including) index `end`."
        if self.Us is None:
            self.Us = [[ 0 for _ in range(self.dim) ] for _ in range(self.len) ]
            self.M = [[ int(i==j) for i in range(self.len) ] for j in range(self.len) ]
            self.N = [ 0 ] * self.len  # squared norms
            self.Us_end = 0
        if end < 0:
            end = self.len - 1
        if self.Us_end > end:
            return
        for j in range(self.Us_end, end + 1):  # Gram-Schmidt decomposition
            self.Us[j] = self.U[j][:]
            for i in range(j):
                if self.N[i]:
                    self.M[j][i] = self.dot(self.Us[j], self.Us[i]) / self.N[i]
                for l in range(self.dim):
                    self.Us[j][l] -= self.M[j][i] * self.Us[i][l]
            self.N[j] = self.norm2(self.Us[j])
            if self.fraction:
                self.N[j] = Fraction(self.N[j])
        self.Us_end = end + 1

    def size_reduce(self, j: int) -> bool:
        "Size reduction at index `j`."
        self.gsd(j)
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
        "Swap columns j and j-1 and update the Gram-Schmidt decomposition."
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
        if newN1:
            self.M[j][j - 1] *= oldN1 / newN1
        else:
            self.M[j][j - 1] = 0
        for l in range(self.dim):
            self.Us[j][l] = oldUs[l] - self.M[j][j - 1] * self.Us[j - 1][l]
        self.N[j] = oldN1 - self.M[j][j - 1]**2 * self.N[j - 1]
        for l in range(j - 1):
            self.M[j][l], self.M[j - 1][l] = self.M[j - 1][l], self.M[j][l]
        if self.Us_end <= j + 1:
            return
        if self.N[j - 1]:
            tmp1 = oldN0 / self.N[j - 1]
            tmp2 = oldM10 * oldN1 / self.N[j - 1]
        else:
            tmp1 = 0
            tmp2 = 0
        for l in range(j + 1, self.Us_end):
            self.M[l][j - 1], self.M[l][j] = tmp1 * self.M[l][j] + tmp2 * self.M[l][j - 1],  self.M[l][j - 1] - oldM10 * self.M[l][j]

    def lll(self, start: int = 1, delta: float = 0.75, sort = True, deep = False, verbose = False) -> None:
        "LLL algorithm for lattice reduction."
        if not 0 < delta <= 1:
            raise ValueError(f"LLL reqires 0 < delta={delta} <= 1")
        self.Mu = None
        self.Nu = None
        j = start
        if verbose:
            print(f"Starting LLL at index {j}")
        steps = 0
        while j < self.len:
            steps += 1
            if self.size_reduce(j) and verbose:
                print(f"Reduced basis vector {j}")
            if all(self.U[j][l] == 0 for l in range(self.dim)):
                if verbose:
                    print(f"Deleting vector {j}")
                self.delete_vector(j)
                continue
            if deep:
                newN = self.norm2(self.U[j])  # new norm if we insert at k = 0
                k = 0
                while delta * self.N[k] <= newN and k < j:  # Deep LLL test
                    newN -= self.M[j][k]**2 * self.N[k]  # new norm if we insert at k + 1
                    k += 1
            else:
                newN = self.N[j] + self.M[j][j - 1]**2 * self.N[j - 1]  # new norm N[j-1] if we swap
                k = j
                if delta * self.N[j - 1] > newN:  # Lovasz condition
                    k -= 1  # swap
            if k == j:
                j += 1
                continue
            if k == j-1:
                if verbose:
                    print(f"Swaping basis vectors {j} <-> {j-1}")
                self.swap(j, newN) # swap vectors
            else:
                if verbose:
                    print(f"Deep insertion of basis vector {j} at {k}")
                self.U[k:j+1] = [self.U[j]] + self.U[k:j]
                self.Us_end = k
            j = max(k, 1)
        if verbose:
            print(f"Total number of steps: {steps}")
        if sort:
            self.sort()

    def project(self, start: int, end: int|None = None) -> "Lattice":
        "Project the lattice onto the span of the given range of Gram-Schmidt vectors."
        if end is None:
            end = self.len
        self.gsd(end-1)
        Uss = [ self.Us[start][:] ]
        for j in range(start + 1, end):
            usj =  self.Us[j][:]
            for i in range(start, j):
                for l in range(self.dim):
                    usj[l] += self.M[j][i] * self.Us[i][l]
            Uss.append(usj)
        lat = Lattice(Uss)
        lat.fraction = self.fraction
        if self.myrank == self.len:
            lat.myrank = lat.len
        if self.Us is not None:
            lat.Us = [ self.Us[j] for j in range(start, end) ]
            lat.N = [ self.N[j] for j in range(start, end) ]
            lat.M = [ self.M[j][start:end] for j in range(start, end) ]
            lat.Us_end = lat.len
            if lat.myrank is None:
                lat.myrank = sum(bool(i) for i in lat.N)
        return lat

    def bkz_round(self, blocksize: int = 3, delta = 0.99, verbose = False) -> "Lattice":
        "Perform one full BKZ round with given blocksize."
        changed = False
        if verbose:
            print(f"Starting BKZ round with blocksize {blocksize}.")
        for start in range(self.len):
            end = min(start + blocksize, self.len)
            # project the lattice
            if verbose:
                print(f"Solving SVP in the projected sublattice {start} - {end-1}.")
            lat = self.project(start, end)
            # solve the SVP
            c = lat.svp_enum(coordinates = True, reduce = False)
            uu = lat.basis() * c
            if uu.norm2() == lat.norm2(lat.U[0]):
                if self.size_reduce(start):
                    changed = True
                continue
            # the new vector is shorter
            if verbose:
                print("Inserting the new vector.")
            changed = True
            # lift the solution
            u = Matrix(self.U[start:end]).transpose() * c
            # insert the new short vector
            self.U = self.U[:start] + [ list(u) ] + self.U[start:]
            self.len += 1
            self.Us.append([0] * self.dim)
            for row in self.M:
                row += [ 0 ]
            self.M.append([0] * (self.len -1 ) + [1])
            self.N.append(0)
            self.Us_end = start
            # run LLL
            self.lll(start = start, delta = delta, sort = False, verbose = verbose)
        return changed

    def bkz(self, blocksize: int = 3, maxrounds = inf, delta = 0.99, verbose = False) -> "Lattice":
        "Perform a BKZ reduction of the basis."
        self.lll(delta = delta, deep = True, sort = False, verbose = verbose)
        rounds = 0
        changed = True
        while changed and rounds < maxrounds:
            rounds += 1
            if verbose:
                print(f"Performing BKZ reduction round {rounds}.")
            changed = self.bkz_round(blocksize = blocksize, delta = delta, verbose = verbose)

    def hkz(self, verbose = False) -> "Lattice":
        "Perform a HKZ reduction of the basis."
        self.lll(delta = 1, deep = True, sort = False, verbose = verbose)
        self.bkz_round(blocksize = self.len, delta = 1, verbose = verbose)

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
            self.Us_end = j
        if self.Mu is not None:
            self.Mu = [[self.Mu[i][j] for j in sort_idx] for i in sort_idx]
        return True

    def replace_basis_vector(self, u: Matrix, coordinates: bool = False) -> Matrix:
        "Returns a unimodular matrix which makes the smallest vector parallel to the given one the first basis vector."
        if not coordinates:
            k = self.coordinates(u)
        else:
            k = u
        W = eye(self.len)
        kmin = inf
        jmin = 0
        for j, kj in enumerate(k):
            if abs(kj) > 0 and abs(kj) < kmin:
                kmin = kj
                jmin = j
        if kmin == inf:
            raise ValueError("Cannot use the zero vector as new basis vector.")
        if kmin < 0:
            k[jmin] *= -1
            W[jmin,jmin] = -1
        k[0], k[jmin] = k[jmin], k[0]
        W[:,0], W[:,jmin] = W[:,jmin], W[:,0]
        for j in range(1,self.len):
            if k[j] != 0:
                k[0], a, b = egcd(k[0], k[j])
                W[:,0], W[:,j] = k[0] * W[:,0] + k[j] * W[:,j], b * W[:,0] - a * W[:,j]
        return W

    def gram_det(self) -> float:
        "Compute the Gram determinant of the current basis."
        self.gsd()
        return sqrt(prod(self.N))

    def hadamard_ratio(self) -> float:
        "Compute the Hadamard ratio of the current basis."
        return (self.gram_det() / prod(map(self.norm2, self.U)) ) ** (1 / self.len)

    def svp(self, method: str = 'all', delta: float = 0.75) -> Matrix:
        "Solve the SVP (approximately) using a given method."
        method = method.lower()
        if method not in ('all', 'hermite', 'lll', 'enum'):
            raise ValueError("Supported methods are: all, hermite, lll, enum")
        if method in ('all', 'hermite'):
            self.hermite()
        if method in ('all', 'lll'):
            self.lll(delta = delta)
        if method == 'enum':
            return self.svp_enum()
        self.sort()
        return Matrix(self.U[0])

    def svp_enum(self, coordinates: bool = False, reduce = False) -> Matrix:
        "Solve the SVP exactly using branch & bound."
        # wet try to improve the basis as much as possible
        if reduce:
            self.lll(delta = 1, deep = True, sort = False)
        self.gsd()
        self.Nu = list(map(self.norm2, self.U))
        lam = min(self.Nu)  # current guess for the shortest length
        k, lam = self.enum_bb([], lam , 0)
        if coordinates:
            return Matrix(k)
        return self.basis() * Matrix(k)

    def enum_bb(self, k:list, lam: float, C: float, beta: list|None = None, b: Matrix|None = None) -> (list, float):
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
                k_new, lam_new = self.enum_bb(k_try, lam, Cj, beta, b)
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
        if lll:
            self.lll(delta = delta)
        method = method.lower()
        if method == 'lll':
            method = 'kannan'
        if hasattr(self, "cvp_" + method):
            method = getattr(self, "cvp_" + method)
        else:
            raise ValueError("Supported methods are: babai_round, babai_plane, kannan, enum")
        return method(x)

    def cvp_babai_round(self, x: Matrix) -> Matrix:
        "Babai's rounding algorithm for approximately solving the CVP."
        rank = self.rank()
        if rank < self.dim:
            self.gsd()
            xi = list(Matrix(self.Us) * x)
            xi = [ x / n for x, n in zip(xi, self.N) ]
            x = Matrix(self.Us).transpose() * Matrix(xi) # we approximate the projection
        c = self.coordinates(x)
        if c is None:
            return None
        return self.basis() * c.applyfunc(round)

    def cvp_babai_plane(self, x: Matrix) -> Matrix:
        "Babai's closest plane algorithm for approximately solving the CVP."
        self.gsd()
        y = list(x)
        for k in range(self.len - 1, -1, -1):
            m = round( self.dot(y, self.Us[k]) / self.N[k] )
            for l in range(self.dim):
                y[l] -= m * self.U[k][l]
        y = Matrix(y)
        return (x - y).applyfunc(round)

    def cvp_kannan(self, x: Matrix, m: int = 1, delta: float = 0.99, deep: bool = True) -> Matrix:
        "Kannan's embedding algorithm for approximately solving the CVP."
        V = self.basis()
        V.append_column(x)
        V.append_row(V.zeros(1,V.cols))
        V[-1,-1]= m
        lat = Lattice( V )
        lat.lll(delta = delta, deep = deep)
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

    def cvp_enum(self, x: Matrix) -> Matrix:
        "Solve the CVP exactly using branch & bound."
        # wet try to improve the basis as much as possible
        self.hermite()
        lam = inf
        self.lll(delta = 1, deep = True, sort = False)
        self.gsd()
        a = self.cvp_kannan(x, deep = True)
        if a is not None:
            lam = (x - self.cvp_kannan(x)).norm2()  # current guess for the distance
        xi = list(Matrix(self.Us) * x)
        xi = [ x / n for x, n in zip(xi, self.N) ]
        rank = self.rank()
        if rank < self.dim:
            x = Matrix(self.Us).transpose() * Matrix(xi) # we approximate the projection
        k, lam = self.enum_bb([], lam , 0, xi, x)
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
