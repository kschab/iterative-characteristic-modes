# iterative-characteristic-modes

An example of the iterative, scattering-based algorithm for computing characteristic modes presented in [1].  

## Implementation notes

The example utilizes pre-calculated data from method of moments simulations to reduce both the code complexity and the computational cost of the demonstration.  Pre-calculated data are provided for an infinitely long, perfectly conducting elliptical cylinder with major and minor axes related by an aspect ratio of 4.  

The use of these pre-calculated data to emulate calls to a general full-wave solver is clearly delineated within the code.  In more realistic applications of the algorithm, those lines would be replaced with calls to whatever full-wave simulator is best suited for the problem being studied.

Three user settings alter the behavior of the demonstration (see code for further details):
- `kdex` : this selects one of three electrical sizes to show how the number of modes with high modal significance affects convergence of the iterative algorithm
- `undersampling` : this parameter controls undersampling of radiation patterns used to construct the matrix representation of the scattering dyadic
- `fastflag` : this flag precalculates a matrix inverse to accelerate emulation of a full-wave solver 
- `plotting` : this flag enables or disables plotting

## Timing and performance

The computational cost of calculating characteristic modes using the full scattering dyadic versus the iterative algorithm can be assessed by disabling the plotting and full-wave acceleration features (`fastflag = 0` and `plotting = 0`).  With these 

|   | ka = 0.1\pi | ka = \pi |
| ------------- | ------------- | --|
| Content Cell  | Content Cell  | |
| Content Cell  | Content Cell  | |

## References
