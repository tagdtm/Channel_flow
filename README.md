# Channel_flow

## Overview
This code is developed to understand DNS channel flow code with Fortran.
The simulation is for incompressive viscous fluid and density, $\rho$, is assumed constant everywhere.
The mainstream direction is x, and the wall normal direction is y so that z shows spanwise direction.
Since the code is not parallelised, it is slow currently.
### Governing equation
Non-dimensionalise Navier-Stokes equation and continuity equation are employed.
The equations are defined as

$$ \frac{\partial u^ \ast_i}{\partial x^ \ast_i}=0 $$

$$ \frac{Dv^ \ast_i}{Dt^ \ast}=-\frac{\partial p^ \ast}{\partial x^ \ast_i}+\frac{1}{Re_{\tau}}\frac{{\partial}^2{u^ \ast}_i}{\partial x^ \ast _jx^ \ast_j} $$

where left hand side (LHS) represent material derivative. The Einstain notation is used to represent vector equation and symbol * shows that the variable is nondimensionalised.

### Mesh resolution
The each direction of mesh resolution is shown in table 1.

(Still the document is under developing)
Direction|cell number|cell size
---|:---:|---: 
x|64|$y^+$
y|64|$y^+$
z|64|$y^+$

### Boundary Condition
Direction|BC
---|:---:
x|Periodic
y|Non-slip wall
z|Periodic

## The Result
Currently, the result shows velocity, pressure and Q criterion.
The gnuplot script provides these graphs.

<img src="/picture/u_mean.jpg" alt="u_mean" width="350" style="display: block; margin: 0 auto"/>
<img src="/picture/u_rms.jpg" alt="u_rms" width="350" style="display: block; margin: 0 auto"/>

By using additional post-pro software such as paraview supports some visualisation shown below.

<img src="/picture/qq1000_over.png" alt="u_rms" width="350" style="display: block; margin: 0 auto"/>

## Future plan
The code will be parallelised at some point.
Also, I would like to implement immerse boundary condition.
