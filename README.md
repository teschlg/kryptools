# kryptools
Gerald Teschl <Gerald.Teschl@univie.ac.at>

This package was written for my course on cryptography. Consequently its intention is
mainly educational, that is, to show how basic algorithms are implemented. In particular,
you are welcome to read (and modify) the source. Any suggestions on how to make the
code more readable or make it better are welcome. However, my main goal is to keep
it simple and readability will be preferred over small speed improvements.

It does not require any external libraries.

The tools contained are:

* number theory: sqrt modulo primes, Jacobi/Legendre symbol, Chinese Remainder Theorem, continued fractions, etc.
* primes: Sieve of Erathostenes, primality tests, generation of random safe/strong primes
* solvers for discrete logarithms (naive, Pollard rho, Shanks baby step/giant step, index calculus, quadratic sieve)
* integer factorization (Fermat, Pollard p-1, Pollard rho, Lenstra's ECM, Dixon, basic quadratic sieve)
* elliptic curves (Weierstrass form), group operations, order, discrete logarithms
* linear algebra: Hermite normal form, Gram-Schmidt
* lattices: Hadamard ratio, Babai rounding/nearest plane algorithm, lattice reduction (Lenstra-Lenstra-Lovaz)

* Matrix: a class for Matrices (inverse, det, reduced echelon form, etc.)
* Poly: a class for polynomials (division, modulo)
* Zmod: a class for the ring of integers modulo an integer

Documentation can be found in the jupyter notebook
(mostly done, but might not contain everything: use the force, read the source).