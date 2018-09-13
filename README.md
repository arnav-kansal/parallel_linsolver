# parallel_linsolver
Parallelizing the linear equations solver using MPI + OPENMP.

## solver
`source/gs.c`

Algorithm
Starting of with a set of n equations and n unknowns,
<img src="/tex/d506b65eceb64c6d397701a33eb5fffb.svg?invert_in_darkmode&sanitize=true" align=middle width=56.20989164999998pt height=22.831056599999986pt/> or,

<img src="/tex/14a83d7929b77f3b7d84b266d415d59e.svg?invert_in_darkmode&sanitize=true" align=middle width=290.2777812pt height=22.831056599999986pt/>

<img src="/tex/1ffc031e096dca3bd4e69de7333165ec.svg?invert_in_darkmode&sanitize=true" align=middle width=290.2777812pt height=22.831056599999986pt/>

..
<img src="/tex/010a91dc8d2c2f70203e12951f94d597.svg?invert_in_darkmode&sanitize=true" align=middle width=298.1451692999999pt height=22.831056599999986pt/>

Given <img src="/tex/db0dce2a6a38aedb28d33f6650cb22e8.svg?invert_in_darkmode&sanitize=true" align=middle width=19.44456194999999pt height=14.15524440000002pt/> and <img src="/tex/a7d0e0605a6acafe642d0b54226ac650.svg?invert_in_darkmode&sanitize=true" align=middle width=13.60734374999999pt height=22.831056599999986pt/> to <img src="/tex/935aab151b542081e51a21ca914e3be6.svg?invert_in_darkmode&sanitize=true" align=middle width=15.18082004999999pt height=22.831056599999986pt/>. We need to calculate <img src="/tex/9fc20fb1d3825674c6a279cb0d5ca636.svg?invert_in_darkmode&sanitize=true" align=middle width=14.045887349999989pt height=14.15524440000002pt/> 's.
Using Jacobi iteration to solve this system of equations:
<p align="center"><img src="/tex/2d84ef387dbb9672f5f1d912600ebaef.svg?invert_in_darkmode&sanitize=true" align=middle width=177.83062934999998pt height=42.77822835pt/></p>

This is nearly embarrasingly parallel and a combination of MPI and OPENMP has been utilized to solve this in the solver.
