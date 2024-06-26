# kryptools
Gerald Teschl <Gerald.Teschl@univie.ac.at>

This package was written for my course on cryptography. Consequently its intention is
mainly educational, that is, to show how basic algorithms are implemented. In particular,
you are welcome to read (and modify) the source. Any suggestions on how to make the
code more readable or make it better are welcome. However, my main goal is to keep
it simple and readability will be preferred over small speed improvements.

The tools contained are:

* number theory: sqrt modulo primes, crt, continued fractions, etc.
* primes: Sieve of Erathostenes, primality tests
* solvers for discrete logarithms (naive, Pollard rho, Shanks baby step/giant step, index calculus, quadratic sieve)
* integer factorization (Pollard p-1, Lentra's ECM, Dixon, basic quadratic sieve)
* linear algebra: Hermite normal form, Gram-Schmidt
* lattices: Babai rounding/nearest plane, lattice reduction

* Matrix: a class for Matrices (inverse, det, reduced echelon form, etc.)
* Poly: a class for polynomials (division, modulo)
* Zmod: a class for the ring of integers modulo an integer

Documentation can be found in the jupyter notebook
(currently incomplete: use the force, read the source).