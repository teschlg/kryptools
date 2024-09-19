"""
Linear algebra
"""

from math import inf, sqrt, prod
from numbers import Number
from fractions import Fraction

class Matrix:
    """
    Matrix class.

    Example:
    
    To define a matrix use
    >>> Matrix([[1, 2], [3, 4]])
    [1, 2]
    [3, 4]
    """
    def __init__(self, matrix, ring = None):
        if not isinstance(matrix[0], list|tuple):
            matrix = [ [x] for x in matrix ]
        self.matrix = matrix
        self.cols = len(matrix[0])
        self.rows = len(matrix)
        if ring:
            self.map(ring)

    def __repr__(self) -> str:
        out = ["[ "] * self.rows
        for j in range(self.cols):
            tmp = [ "" ] * self.rows
            max_len = 0
            for i in range(self.rows):
                tmp[i] = str(self.matrix[i][j])
                max_len = max(max_len, len(tmp[i]))
            for i in range(self.rows):
                out[i] += " " * (max_len - len(tmp[i])) + tmp[i]
                if j < self.cols - 1:
                    out[i] += ", "
        for i in range(self.rows):
            out[i] += " ]"
        return '\n'.join(out)

    def __len__(self):
        return self.cols * self.rows

    def __getitem__(self, item):
        if isinstance(item, tuple):
            i, j = item
            if isinstance(i, int) and isinstance(j, int):
                return self.matrix[i][j]
            if isinstance(i, int):
                rows = [ i ]
            elif isinstance(i, list):
                rows = i
            else:
                rows = range(self.rows)[i]
            if isinstance(j, int):
                cols = [ j ]
            elif isinstance(j, list):
                cols = i
            else:
                cols = range(self.cols)[j]
            return Matrix([[self.matrix[i][j] for j in cols] for i in rows])
        if isinstance(item, int):
            i, j = divmod(item, self.cols)
            return self.matrix[i][j]
        return Matrix([self.matrix[k // self.cols][k % self.cols] for k in range(self.cols * self.rows)[item]])

    def __setitem__(self, item, value):
        if isinstance(item, tuple):
            i, j = item
            if isinstance(i, int) and isinstance(j, int):
                self.matrix[i][j] = value
                return
            if isinstance(i, int):
                rows = [ i ]
            elif isinstance(i, list):
                rows = i
            else:
                rows = range(self.rows)[i]
            if isinstance(j, int):
                cols = [ j ]
            elif isinstance(j, list):
                cols = i
            else:
                cols = range(self.cols)[j]
            for i, ii in zip(cols,range(len(cols))):
                for j, jj in zip(rows,range(len(rows))):
                    self.matrix[j][i] = value[jj,ii]
            return
        i, j = divmod(item, self.cols)
        self.matrix[i][j] = value

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.matrix == other.matrix

    def map(self, func):
        "Apply a function to all elements in place."
        for row in self.matrix:
            row[:] = map(func, row)

    def applyfunc(self, func):
        "Apply a function to all elements."
        tmp = self[:,:]
        tmp.map(func)
        return tmp

    def norm2(self) -> float:
        "Squared Frobenius/Euclidean norm."
        return sum( sum(x*x for x in row) for row in self.matrix )

    def norm(self, p: int = 2) -> float:
        "p-norm of a matrix regarded as a vector."
        if p == 2:
            return sqrt(self.norm2())
        if p == 1:
            return sum( sum(abs(x) for x in row) for row in self.matrix )
        if p == inf:
            return max( max(abs(x) for x in row) for row in self.matrix )
        tmp = sum( sum(abs(x)**p for x in row) for row in self.matrix )
        return tmp**(1/p)


    def dot(self, other) -> int:
        "Dot product of two vectors."
        if self.rows == 1 and other.rows == 1 and self.cols == other.cols:
            return sum(x * y for x, y in zip(self.matrix[0], other.matrix[0]))
        if self.cols == 1 and other.cols == 1 and self.rows == other.rows:
            return sum(x[0] * y[0] for x, y in zip(self.matrix, other.matrix))
        return NotImplemented

    def transpose(self) -> "Matrix":
        "Transpose of a matrix."
        return Matrix([list(i) for i in zip(*self.matrix)])

    def multiply(self, other) -> "Matrix":
        "Matrix multiplication."
        if not isinstance(other, Matrix):
            return NotImplemented
        if self.cols != other.rows:
            raise NotImplementedError("Matrix dimensions do not match!")
        result = self.zeros(self.rows, other.cols)
        for i in range(self.rows):
            for j in range(other.cols):
                for k in range(other.rows):
                    result.matrix[i][j] += self.matrix[i][k] * other.matrix[k][j]
        if self.rows == 1 and other.cols == 1:
            return result.matrix[0][0]
        return result

    def __add__(self, other) -> "Matrix":
        if isinstance(other, Matrix):
            if other.cols != self.cols or other.rows != self.rows:
                raise NotImplementedError("Matrix dimensions do not match!")
            return Matrix([ [ x1 + y1 for x1, y1 in zip(x,y)] for x, y in zip(self.matrix, other.matrix)])
        return NotImplemented

    def __sub__(self, other) -> "Matrix":
        if isinstance(other, Matrix):
            if other.cols != self.cols or other.rows != self.rows:
                raise NotImplementedError("Matrix dimensions do not match!")
            return Matrix([ [ x1 - y1 for x1, y1 in zip(x,y)] for x, y in zip(self.matrix, other.matrix)])
        return NotImplemented

    def __neg__(self) -> "Matrix":
        return -1 * self

    def __pos__(self) -> "Matrix":
        return self

    def __mul__(self, other) -> "Matrix":
        if isinstance(other, Matrix):
            return self.multiply(other)
        if isinstance(other, Number) or type(other) == type(self.matrix[0][0]):
            return Matrix([ [item * other for item in row] for row in self.matrix ])
        return NotImplemented

    def __rmul__(self, other) -> "Matrix":
        if isinstance(other, Number) or type(other) == type(self.matrix[0][0]):
            return Matrix([ [item * other for item in row] for row in self.matrix ])
        return NotImplemented

    def rref(self) -> "Matrix":
        "Compute the reduced echelon form of a matrix M."
        n, m = self.cols, self.rows
        R = self[:, :]
        i = 0
        for j in range(n):
            if not R[i, j]: # search for a nonzero entry in the present column
                for ii in range(i+1,m):
                    if R[ii, j]:
                        R[i, :], R[ii, :] = R[ii, :], R[i, :]  # swap rows
                        break
                else:
                    continue  # all entries are zero
            if R[i, j] != 1:
                R[i, :] = 1/ R[i, j] * R[i, :]  # make the pivot one
            for ii in range(m):  # remove the column entries above/below the pivot
                if i == ii:
                    continue
                tmp = R[ii, j]
                R[ii, ::] -= tmp * R[i, :]
            i += 1
            if i == m:
                break
        return R

    def det(self) -> int:
        "Compute the determinant of a matrix M."
        if self.rows != self.cols:
            raise ValueError("Matrix must be square!")
        n, m = self.cols, self.rows
        R = self[:, :]
        D = 1
        i = 0
        for j in range(n):
            if not R[i, j]: # search for a nonzero entry in the present column
                for ii in range(i+1,m):
                    if R[ii, j]:
                        D *= -1
                        R[i, :], R[ii, :] = R[ii, :], R[i, :]  # swap rows
                        break
                else:
                    return 0  # all entries are zero
            if R[i, j] != 1:
                D *= R[i, j]
                R[i, :] = 1/ R[i, j] * R[i, :]  # make the pivot one
            for ii in range(i+1,n):  # remove the column entries below the pivot
                if i == ii:
                    continue
                tmp = R[ii, j]
                R[ii, ::] -= tmp * R[i, :]
            i += 1
        return D

    def inv(self) -> "Matrix":
        "Compute the inverse of a square matrix M."
        if self.rows != self.cols:
            raise ValueError("Matrix must be square!")
        n = self.cols
        MM = Matrix([[0 for _ in range(2*n)] for _ in range(n)])
        for i in range(n):
            MM[i,n+i] = 1
        MM[:,0:n] = self
        MM = MM.rref()
        if not prod(MM[i, i] for i in range(n)):
            raise ValueError("Matrix is not invertible!")
        return MM[:,n:]

    def is_unimodular(self) -> bool:
        "Test if the matrix is unimodular."
        if self.rows != self.cols:
            return False
        def is_integer(i):
            if isinstance(i, int) or (isinstance(i, Fraction) and i.denominator == 1):
                return True
            return False
        return all([is_integer(i) for i in self]) and self.det()**2 == 1

    def zeros(self, m: int = None, n: int = None):
        "Returns a zero matrix of the same dimension."
        if not m and not n:
            n, m = self.cols, self.rows
        elif not n:
            n = m
        try:
            zero = 0 * self[0]
        except:
            zero = 0
        return Matrix([[ zero for j in range(n) ] for i in range(m) ])

    def eye(self, m: int = None, n: int = None):
        "Returns an identity matrix of the same dimension."
        def delta(i, j):
            if i == j:
                return one
            return zero
        if not m and not n:
            n, m = self.cols, self.rows
        elif not n:
            n = m
        try:
            zero = 0 * self[0]
        except:
            zero = 0
        one = 1 + zero
        return Matrix([[ delta(i, j) for j in range(n) ] for i in range(m) ])


def zeros(m: int, n: int = None, zero = 0) -> "Matrix":
    "Returns a zero matrix of the given dimension."
    if not n:
        n = m
    return Matrix([[ zero for j in range(n) ] for i in range(m) ])

def eye(m:int, n: int = None) -> "Matrix":
    "Returns an identity matrix of the given dimension."
    def delta(i, j):
        if i == j:
            return 1
        return 0
    if not n:
        n = m
    return Matrix([[ delta(i, j) for j in range(n) ] for i in range(m) ])
