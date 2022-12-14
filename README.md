# iterative-characteristic-modes

<img src="iterative-image.png" height="100"/>

An example of the iterative, scattering-based algorithm for computing characteristic modes presented in [1].  

## Implementation notes

The example utilizes pre-calculated data from method of moments simulations to reduce both the code complexity and the computational cost of the demonstration.  Pre-calculated data are provided for an infinitely long, perfectly conducting elliptical cylinder with major and minor axes related by an aspect ratio of 4.  All electric fields are polarized parallel to the axis of the cylinder.

Steps listed in the comments of the code correspond to the steps in Algorithm 2 of [1].

The use of these pre-calculated data to emulate calls to a general full-wave solver is clearly delineated within the code.  In more realistic applications of the algorithm, those lines would be replaced with calls to whatever full-wave simulator is best suited for the problem being studied.

Four user settings alter the behavior of the demonstration (see code for further details):
- `kdex` : this selects one of three electrical sizes to show how the number of modes with high modal significance affects convergence of the iterative algorithm
- `undersampling` : this parameter controls undersampling of radiation patterns used to construct the matrix representation of the scattering dyadic
- `fastflag` : this flag precalculates a matrix inverse to accelerate emulation of a full-wave solver 
- `plotting` : this flag enables or disables plotting

By default, the modified Gram-Schmidt procedure is used in place of the Gram-Schmidt process described in Steps 8 and 9 of Algorithm 2 in [1].  This can be disabled by setting `mgs = 0`.

Note that this example does not include any of the time-domain, multi-frequency techniques discussed in Sec. II-B of [1].

## Timing and performance

Live plot updates greatly increase the apparent time of the iterative algorithm. To accurately assess the relative computational cost of calculating characteristic modes using the full scattering dyadic versus the iterative algorithm, disable the plotting and full-wave acceleration features (`fastflag = 0` and `plotting = 0`).  With these settings, along with no undersampling (`undersampling = 1`), the following computational times, numbers of iterations, and numbers of modes with modal significance greater than 0.01 were obtained on a MacBook Air (M1, 2020):

|   | ka = 0.1?? | ka = ?? | ka = 10?? |
| ------------- | ------------- | -- | -- |
| Full calculation time (s) | 34  | 33 | 34 |
| Iterative algorithm time (s) | 1.5  | 4.7 | 21 |
| Iterations | 6 | 12 | 56 |
| \# modes MS>0.01 | 3 | 9  | 51|

## References

[1] Lundgren, Schab, Capek, Gustafsson, and Jelinek, "Iterative calculation of characteristic modes using arbitrary full-wave solvers,'' arXiv preprint arXiv:2209.00097 (2022)

