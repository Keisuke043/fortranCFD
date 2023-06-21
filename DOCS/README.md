
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
where $ρ$ is the density, $u$ is the flow velocity, $e$ is the total energy per unit mass, $p$ is the pressure, $T$ is the temperature, $Y_k$ is the mass fraction of species $k$, $n$ is the total number of species, $\tau_{xx}$ is the viscous stress, $q_x$ is the heat flux, $V_k$ is the diffusion velocity, $\omega_k$ is the net chemical production rate. In the energy equation, the total energy per unit mass, $e$, is given as:
$$
e=-\frac{p}{\rho}+\frac{u^2}{2}+h \\
h=\sum_{k=1}^nY_kh_k
$$
where $h$ is the mixture total enthalpy and $h_k$ is the enthalpy of species $k$. In addition to Eq. (1), the ideal gas equation of state is considered.
$$
p=\sum_{k=1}^n\rho Y_kR_kT \\
$$
where $R_k$ is the gas constant of species $k$.



$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$



<img src="https://latex.codecogs.com/svg.image?\frac{\partial&space;Q}{\partial&space;t}&plus;\frac{\partial&space;E}{\partial&space;x}=0&space;" />

