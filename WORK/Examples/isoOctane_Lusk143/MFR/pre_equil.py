
import cantera as ct
import numpy as np
import subprocess


def calc_equil_prod():

    gasEquil = ct.Solution(chemPath_cti)
    gasEquil.TP = (T, ct.one_atm*p)
    gasEquil.set_equivalence_ratio(equivalence_ratio, \
                                   fuel=fuel, oxidizer=oxid)

    eq_params = ['TP','TV','HP','SP','SV','UV']
    eq_param  = eq_params[0]
    gasEquil.equilibrate(eq_param)
    equil_temp   = gasEquil.T
    equil_spname = gasEquil.species_names
    equil_Xi     = gasEquil.X
    outfile = './equil{0}_phi{1}_p{2}_T{3}.txt'.format(eq_param, equivalence_ratio, p, T)
    f = open(outfile, 'w')
    f.write('{}\n'.format(len(equil_spname)))
    for sp_name, Xi in zip(equil_spname, equil_Xi):
        print(sp_name, Xi)
        f.write('{} {}\n'.format(sp_name, Xi))
    f.close()


def calc_0D_prod():
    gas = ct.Solution(chemPath_cti)
    gas.TP = (T, ct.one_atm*p)
    gas.set_equivalence_ratio(equivalence_ratio, \
                              fuel=fuel, oxidizer=oxid)
    r = ct.IdealGasReactor(gas)
    reactorNetwork = ct.ReactorNet([r])
    timeHistory_state = ct.SolutionArray(gas, extra=['t'])
    count = 0
    cur_time = 0.0
    end_time = 1.0
    while(cur_time < end_time):
        cur_time = reactorNetwork.step()
        if count % 10 == 0:
            timeHistory_state.append(r.thermo.state, t=cur_time)
        count += 1

    sim0D_spname = gasEquil.species_names
    sim0D_Xi     = timeHistory_state.X[-1]
    outfile = './sim0D_phi{0}_T{1}_p{2}.txt'.format(equivalence_ratio, T, p)
    f = open(outfile, 'w')
    f.write('{}\n'.format(len(equil_spname)))
    for sp_name, Xi in zip(sim0D_spname, sim0D_Xi):
        print(sp_name, Xi)
        f.write('{} {}\n'.format(sp_name, Xi))
    f.close()


if __name__ == "__main__":

    chemPath_cti = './Lusk143_isoC8H18/chem.cti'
    fuel = 'IC8H18'
    oxid = 'O2:1, N2:3.76'

    p = 1.0
    T = 1300.0
    equivalence_ratio = 1.0

    calc_equil_prod()

