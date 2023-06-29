
## Governing Equations for Combustion

This simulation code solves the fully compressible Navier-Stocks equations with energy and species conservation equations are considered in the governing equations.

### Governing equations in the Cartesian coordinate
In the planar flame simulations, the one-dimensional Navier-Stocks equations in the conservative form in the Cartesian coordinate are written as:

$$
\frac{\partial Q}{\partial t}+\frac{\partial(E-E_\nu)}{\partial x}=S \\
$$

where $t$ and $x$ donate time and space in the Cartesian coordinate, $Q$ is the vector of conserved variables, $E$ is the convective flux vector with convection and pressure terms, $E_\nu$ is the viscous flux vector with viscosity, heat conduction, and diffusion terms, and $S$ is the vector of source terms. These vectors are expressed as follows:

$$
Q=\begin{pmatrix} \rho \\ 
                  \rho u \\ 
                  \rho e \\
                  \rho Y_1 \\ 
                  \vdots \\ 
                  \rho Y_n \end{pmatrix}, \quad
E=\begin{pmatrix} \rho u \\ 
                  \rho u^2+p \\ 
                  (\rho e+p)u \\
                  \rho uY_1 \\ 
                  \vdots \\ 
                  \rho uY_n \end{pmatrix}, \quad
E_\nu=\begin{pmatrix} 0 \\ 
                      \tau_{xx} \\ 
                      q_x \\
                      -\rho Y_1V_1 \\ 
                      \vdots \\ 
                      -\rho Y_nV_n \end{pmatrix}, \quad
S=\begin{pmatrix} 0 \\ 
                  0 \\ 
                  0 \\
                  \omega_1 \\ 
                  \vdots \\ 
                  \omega_n \end{pmatrix}
$$

where $ρ$ is the density, $u$ is the flow velocity, $e$ is the total energy per unit mass, $p$ is the pressure, $T$ is the temperature, $Y_k$ is the mass fraction of species $k$, $n$ is the total number of species, $\tau_{xx}$ is the viscous stress, $q_x$ is the heat flux, $V_k$ is the diffusion velocity, $\omega_k$ is the net chemical production rate. In the energy equation, the total energy per unit mass, $e$, is given as:

$$
e=-\frac{p}{\rho}+\frac{u^2}{2}+h
$$
$$
h=\sum_{k=1}^nY_kh_k
$$

where $h$ is the mixture total enthalpy and $h_k$ is the enthalpy of species $k$. In addition to Eq. (1), the ideal gas equation of state is considered.

$$
p=\sum_{k=1}^n\rho Y_kR_kT
$$

where $R_k$ is the gas constant of species $k$.


For the evaluation of the species transport properties, mixture-averaged diffusion is selected in this study.  The ordinary mixture-averaged diffusion velocity, $V_k$, is given in the Curtiss-Hirschfelder approximation [1].

$$
V_k=\frac{1}{X_k}D_{km}\frac{\partial X_k}{\partial x}
$$
$$
D_{km}=\frac{1}{X_k}\frac{\partial X_k}{\partial x}
$$

where $D_{km}$ is the mixture-averaged diffusion coefficient, and $D_{kj}$ is the binary diffusion coefficients of species $k$ into species $j$.
The viscous stress, $\tau_{xx}$, is expressed as:

$$
\tau_{xx}=(2\mu+\lambda)\frac{\partial u}{\partial x}
$$

where $\mu$ is the dynamic viscosity coefficient and $\lambda$ is the second viscosity coefficient. In the assumption of Stokes hypothesis $(\lambda+\frac{2}{3}\mu=0)$ [2], the viscous stress is simplified to

$$
\tau_{xx}=\frac{4}{3}\mu\frac{\partial u}{\partial x}
$$

The heat flux, $q_x$, in the energy equation is expressed as:

$$
q_{x}=u\tau_{xx}+\kappa\frac{\partial T}{\partial x}-\rho\sum_{k=1}^nh_kY_kV_k
$$

where $\kappa$ is the thermal conductivity of the mixture.

### Governing equations in the spherical coordinate
In the spherical flame simulations, the one-dimensional Navier-Stocks equations in the spherical coordinate are written as:

$$
\frac{\partial Q}{\partial t}+\frac{1}{r^2}\frac{\partial (r^2(E-E_\nu))}{\partial r}+\frac{\partial E_p}{\partial r}=S \\
$$

$$
Q=\begin{pmatrix} \rho \\ 
                  \rho u \\ 
                  \rho e \\
                  \rho Y_1 \\ 
                  \vdots \\ 
                  \rho Y_n \end{pmatrix}, \quad
E=\begin{pmatrix} \rho u \\ 
                  \rho u^2 \\ 
                  (\rho e+p)u \\
                  \rho uY_1 \\ 
                  \vdots \\ 
                  \rho uY_n \end{pmatrix}, \quad
E_p=\begin{pmatrix} 0 \\ 
                    p \\ 
                    0 \\
                    0 \\ 
                    \vdots \\ 
                    0 \end{pmatrix}, \quad
E_\nu=\begin{pmatrix} 0 \\ 
                      \tau_{xx} \\ 
                      q_x \\
                      -\rho Y_1V_1 \\ 
                      \vdots \\ 
                      -\rho Y_nV_n \end{pmatrix}, \quad
S=\begin{pmatrix} 0 \\ 
                  0 \\ 
                  0 \\
                  \omega_1 \\ 
                  \vdots \\ 
                  \omega_n \end{pmatrix}
$$

where $r$ donates space in the spherical coordinate.

In the one-dimensional spherical coordinate, the derivative of the convective and viscous fluxes is written as the second term of Eq. (11) because of the change in area is perpendicular to the radial direction. On the other hand, the geometry effect is negligible in the pressure term of the momentum conservation equation due to the isotropic nature of pressure. Transport properties of the spherical configuration are obtained by replacing $x$ with $r$ in Eq. (6)–(10).

[1] J. O. Hirschfelder, C. F. Curtiss, and R. B. Bird, Molecular Theory of Gases and Liquids,
John Wiley and Sons, New York, 1954.

[2] G.G. Stokes, On the theories of the internal friction of fluids in motion, and of the equilibrium and motion of elastic solids. Trans. Camb. Philos. Soc. 8, (1845) 287–319.

---

### Numerical schemes
In the present study, the fluid transport and chemical reactions are solved separately by the time-splitting scheme [3,4]. In this method, Eq. (1) is separated into two equations of partial differential equations (PDEs) of the multicomponent non-reactive flow and ordinary differential equations (ODEs) of chemical reactions. In the case of Cartesian coordinate, the split form of Eq. (1) is given by:

$$
\left\{
\begin{array}{l}
\begin{align}
\frac{\partial Q}{\partial t}+\frac{\partial(E-E_\nu)}{\partial x}=0 \\
\frac{\partial Q}{\partial t}=S \\
\end{align}
\end{array}
\right. \\
\left(R=-\frac{\partial}{\partial x}\left(E(Q)-E_\nu(Q)\right)\right)
$$

In the fluid simulation, the finite volume method (FVM) is applied for the spatial discretization of Eq. (13). In the FVM, the transfer of conserved variables between control volumes is evaluated by the numerical fluxes at the interface of the adjacent control volumes.

To achieve higher-order spatial accuracy, the monotonic upstream-centered scheme for conservation laws (MUSCL) [6] with minmod limiter is used.

The convective numerical fluxes are calculated by the Harten-Lax-van leer contact (HLLC) Riemann solver [5].

The viscous numerical fluxes including the viscosity, heat conduction, and diffusion are calculated by second-order central difference approximation.

After the evaluation of the fluxes, the conserved variables, $Q$, are updated in time with the third-order total variation diminishing (TVD) Runge-Kutta scheme [7].

$$
\left\{
\begin{array}{l}
\begin{align}
{Q}^{(1)} &=Q^n+\Delta t \cdot R\left(Q^n\right) \\
Q^{(2)} &= \frac{3}{4}Q^n+\frac{1}{4}\left(Q^{(1)}+\Delta t\cdot R \left( Q^{(1)} \right) \right) \\
Q^{n+1} &= \frac{1}{3}Q^n+\frac{2}{3}\left(Q^{(2)}+\Delta t\cdot R \left( Q^{(2)} \right) \right) \\
\end{align}
\end{array}
\right. \\
\left(R=-\frac{\partial}{\partial x}\left(E(Q)-E_\nu(Q)\right)\right)
$$

For the time integration of chemical reactions in Eq. (14), the minimum‑error adaptation of a chemical‑kinetic ODE solver (MACKS) [8] is used. This robust Jacobian-free scheme enables a significant reduction in the computational cost even in stiff ODEs of chemical reactions.

[3] G. Strang, On the construction and comparison of difference schemes, SIAM J. Numer. Anal. 5 (3) (1968) 506–517.

[4] R.P. Fedkiw, B. Merriman, S. Osher, High accuracy numerical methods for thermally perfect gas flows with chemistry, J. Comput. Phys. 132 (2) (1997) 175–190.

[5] B.V. Leer, Flux vector splitting for the Euler equations, Conference on Numerical Methods in Fluid Dynamics (1982) 507–512.

[6] E. Toro, M. Spruce, W. Speares, Restoration of the contact surface in the HLL-Riemann solver, Shock Waves 4 (1) (1994) 25–34.

[7] S. Gottlieb, C. Shu, Total variation diminishing Runge-Kutta schemes, Math. Comput. 67 (1998) 73–85.

[8] Y. Morii, E. Shima, Optimization of one-parameter family of integration formulae for solving stiff chemical-kinetic ODEs, Sci. Rep. 10 (2020) 21330.


