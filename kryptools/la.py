"""
Linear algebra
"""

from math import sqrt, prod

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
            rows = range(self.rows)[i]
            if isinstance(i, int):
                rows = [ rows ]
            cols = range(self.cols)[j]
            if isinstance(j, int):
                cols = [ cols ]
            return Matrix([[self.matrix[i][j] for j in cols] for i in rows])
        i, j = divmod(item, self.cols)
        return self.matrix[i][j]

    def __setitem__(self, item, value):
        if isinstance(item, tuple):
            i, j = item
            if isinstance(i, int) and isinstance(j, int):
                self.matrix[i][j] = value
                return
            rows = range(self.rows)[i]
            if isinstance(i, int):
                rows = [ rows ]
            cols = range(self.cols)[j]
            if isinstance(j, int):
                cols = [ cols ]
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
        "Apply a function to all elements in place"
        for row in self.matrix:
            row[:] = map(func, row)

    def applyfunc(self, func):
        "Apply a function to all elements"
        tmp = self[:,:]
        tmp.map(func)
        return tmp

    def norm2(self) -> float:
        "Squared Frobenius/Euclidean norm."
        return sum( sum(x*x for x in row) for row in self.matrix )

    def norm(self) -> float:
        "Frobenius/Euclidean norm."
        return sqrt(self.norm2())

    def dot(self, other) -> int:
        if self.rows == 1 and other.rows == 1 and self.cols == other.cols:
            return sum(x * y for x, y in zip(self.matrix[0], other.matrix[0]))
        if self.cols == 1 and other.cols == 1 and self.rows == other.rows:
            return sum(x[0] * y[0] for x, y in zip(self.matrix, other.matrix))
        return NotImplemented

    def transpose(self) -> "Matrix":
        return Matrix([list(i) for i in zip(*self.matrix)])

    def multiply(self, other) -> "Matrix":
        if not isinstance(other, Matrix) or self.cols != other.rows:
            return NotImplemented
        result = [[0 for j in range(other.cols)] for i in range(self.rows)]
        for i in range(self.rows):
            for j in range(other.cols):
                for k in range(other.rows):
                    result[i][j] += self.matrix[i][k] * other.matrix[k][j]
        return Matrix(result)

    def __add__(self, other) -> "Matrix":
        if isinstance(other, Matrix) and other.cols == self.cols and other.rows == self.rows:
            return Matrix([ [ x1 + y1 for x1, y1 in zip(x,y)] for x, y in zip(self.matrix, other.matrix)])
        return NotImplemented

    def __sub__(self, other) -> "Matrix":
        if isinstance(other, Matrix) and other.cols == self.cols and other.rows == self.rows:
            return Matrix([ [ x1 - y1 for x1, y1 in zip(x,y)] for x, y in zip(self.matrix, other.matrix)])
        return NotImplemented

    def __neg__(self) -> "Matrix":
        return -1 * self

    def __mul__(self, other) -> "Matrix":
        if isinstance(other, Matrix):
            return self.multiply(other)
        return NotImplemented

    def __rmul__(self, other) -> "Matrix":
        if isinstance(other, Matrix):
            return self.multiply(other)
        return Matrix([ [item * other for item in row] for row in self.matrix ])

    def rref(self) -> "Matrix":
        "Compute the reduced echelon form of a matrix M."
        n, m = self.cols, self.rows
        R = self[:, :]
        i = 0
        for j in range(n):
            if not R[i, j]: # search for am nonzero entry in the present column
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
        if self.rows != self.rows:
            raise ValueError("Matrix must be square!")
        n = self.cols
        R = self[:, :]
        D = 1
        i = 0
        for j in range(n):
            if not R[i, j]: # search for am nonzero entry in the present column
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
        if self.rows != self.rows:
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


def zeros(m:int, n: int) -> "Matrix":
    return Matrix([[ 0 for j in range(n)] for i in range(m) ])

def eye(m:int, n: int) -> "Matrix":
    def delta(i, j):
        if i == j:
            return 1
        return 0
    return Matrix([[ delta(i, j) for j in range(n) ] for i in range(m) ])
