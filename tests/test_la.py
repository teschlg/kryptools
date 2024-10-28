import pytest
from fractions import Fraction
from kryptools import Matrix, Zmod


def test_Matrix():
	M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 12]])
	M.map(Fraction)
	Mi = M.inv()
	assert M * Mi == M.eye()
	
	Z_11 = Zmod(11)
	M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 12]], ring=Z_11)
	Mi = M.inv()
	assert M * Mi == M.eye()