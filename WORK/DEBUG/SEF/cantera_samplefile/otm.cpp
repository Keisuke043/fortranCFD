//
// Created by Youhi Morii on 2021/06/21.
// Converted from mfc.cpp to otm.cpp by Akira Tsunoda on 2022/02/28.
//

#include "radcal.hpp"
#include <boost/algorithm/string.hpp>

void radiationHeatLoss::OTMCalcStFlow::readPMACsFile(std::string fileName){
    try {
        YAML::Node config = YAML::LoadFile(fileName);
        PMACs.resize(config.size());
        radSpName.resize(config.size());
        radSpIndex.resize(config.size());

        for (size_t i = 0; i < config.size(); i++){
            radSpName[i] = config[i]["name"].as<std::string>();
            radSpIndex[i] = m_thermo->speciesIndex(radSpName[i]);
            for (size_t j = 0; j < config[i]["params"].size(); j++){
                PMACs[i].push_back(config[i]["params"][j].as<double>());
            }
        }
    } catch (YAML::Exception &e) {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }
}

double radiationHeatLoss::OTMCalcStFlow::OTMCalc(size_t j, double *x){
    // calculation of the mean Planck absorption coefficient
    double radiative_heat_loss = 0;
    double k_P = 0;
    double T04 = pow(init_T, 4);
    double T1 = T(x, j);
    double T2 = pow(T(x, j), 2);
    double T3 = pow(T(x, j), 3);
    double T4 = pow(T(x, j), 4);
    double T5 = pow(T(x, j), 5);

    for (size_t i = 0; i < radSpIndex.size(); i++){
        if (radSpIndex[i] != npos){
            double k_P_i = 0;
            k_P_i = PMACs[i][0]      + PMACs[i][1] * T1 + PMACs[i][2] * T2\
                  + PMACs[i][3] * T3 + PMACs[i][4] * T4 + PMACs[i][5] * T5;

            k_P += k_P_i * X(x, radSpIndex[i], j) * m_press / Cantera::OneAtm;
        }
    }
    // calculation of the radiative heat loss term
    radiative_heat_loss = -4 * k_P * StefanBoltz * ( T4 - T04 );

    return radiative_heat_loss;
}

void radiationHeatLoss::OTMCalcStFlow::evalResidual(double *x, double *rsd, int *diag, double rdt, size_t jmin,
                                                     size_t jmax) {
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the radiative heat loss term
            radiative_heat_loss = OTMCalc(j, x);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U,0)] =
                -(rho_u(x,1) - rho_u(x,0))/m_dz[0]
                -(density(1)*V(x,1) + density(0)*V(x,0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V,0)] = V(x,0);
            if (doEnergy(0)) {
                rsd[index(c_offset_T,0)] = T(x,0);
            } else {
                rsd[index(c_offset_T,0)] = T(x,0) - T_fixed(0);
            }
            rsd[index(c_offset_L,0)] = -rho_u(x,0);

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x,k,0);
                rsd[index(c_offset_Y + k, 0)] =
                    -(m_flux(k,0) + rho_u(x,0)* Y(x,k,0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            rsd[index(c_offset_V,j)]
            = (shear(x,j) - lambda(x,j) - rho_u(x,j)*dVdz(x,j)
               - m_rho[j]*V(x,j)*V(x,j))/m_rho[j]
              - rdt*(V(x,j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x,j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x,j)*dYdz(x,k,j);
                double diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                                / (z(j+1) - z(j-1));
                rsd[index(c_offset_Y + k, j)]
                = (m_wt[k]*(wdot(k,j))
                   - convec - diffus)/m_rho[j]
                  - rdt*(Y(x,k,j) - Y_prev(k,j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {
                setGas(x,j);

                // heat release term
                const vector_fp& h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp& cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5*(m_flux(k,j-1) + m_flux(k,j));
                    sum += wdot(k,j)*h_RT[k];
                    sum2 += flxk*cp_R[k]/m_wt[k];
                }
                sum *= GasConstant * T(x,j);
                double dtdzj = dTdz(x,j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = - m_cp[j]*rho_u(x,j)*dtdzj
                                            - divHeatFlux(x,j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j]*m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt*(T(x,j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x,j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x,j) - lambda(x,j-1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}
