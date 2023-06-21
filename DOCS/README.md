
### Governing Equations for combustion

This simulation code solves the fully compressible Navier-Stocks equations with energy and species conservation equations are considered in the governing equations.

#### Governing equations in the Cartesian coordinate
In the planar flame simulations, the one-dimensional Navier-Stocks equations in the conservative form in the Cartesian coordinate are written as:

$$
\frac{\partial Q}{\partial t}+\frac{\partial(E-E_\nu)}{\partial x}=S \\
$$

where $t$ and $x$ donate time and space in the Cartesian coordinate, $Q$ is the vector of conserved variables, $E$ is the convective flux vector with convection and pressure terms, $E_\nu$ is the viscous flux vector with viscosity, heat conduction, and diffusion terms, and $S$ is the vector of source terms. These vectors are expressed as follows:

$$
Q=\begin{pmatrix} \rho \\ \rho u \\ \rho e \\
  \rho Y_1 \\ \vdots \\ \rho Y_n \end{pmatrix}, \quad
E=\begin{pmatrix} \rho u \\ \rho u^2+p \\ (\rho e+p)u \\
  \rho uY_1 \\ \vdots \\ \rho uY_n \end{pmatrix}, \quad
E_\nu=\begin{pmatrix} 0 \\ \tau_{xx} \\ q_x \\
  -\rho Y_1V_1 \\ \vdots \\ -\rho Y_nV_n \end{pmatrix}, \quad
S=\begin{pmatrix} 0 \\ 0 \\ 0 \\
  \omega_1 \\ \vdots \\ \omega_n \end{pmatrix}, \quad
$$



$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$



<img src="https://latex.codecogs.com/svg.image?\frac{\partial&space;Q}{\partial&space;t}&plus;\frac{\partial&space;E}{\partial&space;x}=0&space;" />

