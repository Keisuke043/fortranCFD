
import cantera as ct
import numpy as np
import subprocess


fuel = 'N-C7H16'
oxid = 'O2:1, N2:3.76'

phi = 1.0
temp = 500.0 # K
pres_atm = 10.0 # atm
pres = ct.one_atm*pres_atm
width = 0.05  # m

savedir = 'exthaust_1Dflame'
subprocess.call(['mkdir','-p', savedir])

chemPath_cti = './nHeptaneSoot_Heinz/Detailed/chem.cti'

gas = ct.Solution(chemPath_cti)
gas.TP = (temp, pres)
gas.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxid)

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)

# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
loglevel = 0  # amount of diagnostic output (0 to 8)

f.solve(loglevel=loglevel, auto=False)

print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0]))

spnames = gas.kinetics_species_names
X_burnt = f.X[:,len(f.grid)-1]
T_burnt = f.T[len(f.grid)-1]

outfile = './{0}/exthaust_phi{1}_p{2}_T{3}_Tb{4:.1f}.txt'.format(savedir, phi, pres_atm, temp, T_burnt)
f = open(outfile, 'w')
f.write('{}\n'.format(len(spnames)))
for sp_name, Xi in zip(spnames, X_burnt):
    print(sp_name, Xi)
    f.write('{} {}\n'.format(sp_name, Xi))
f.close()


