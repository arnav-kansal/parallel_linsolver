# parallel_linsolver
Parallelizing the linear equations solver using MPI + OPENMP.

## solver
`source/gs.c`

Algorithm:

     Starting of with a set of n equations and n unknowns,
		$AX=b$ or, 
		$a_{11}x_1 + a_{12}x_2 + a_{13}x_3 + ... + a_{1n}x_n = b_1$
		$a_{21}x_1 + a_{22}x_2 + a_{23}x_3 + ... + a_{2n}x_n = b_2$
		..
		$a_{n1}x_1 + a_{n2}x_2 + a_{n3}x_3 + ... + a_{nn}x_n = b_n$

Given $a_{ij}$ and $b_1$ to $b_n$. We need to calculate $x_i$ 's.
Using Jacobi iteration to solve this system of equations:
$$x_j = \frac{b_j - \sum_{i=1, i\neq j}^n a_{ij}x_j}{a_{jj}}$$

This is nearly embarrasingly parallel and a combination of MPI and OPENMP has been utilized to solve this in the solver.
