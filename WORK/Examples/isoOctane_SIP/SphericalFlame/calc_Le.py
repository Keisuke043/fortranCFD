import cantera as ct
import csv


def write2csv(alpha, Df, Do, phi, isInertConst):

    if isInertConst:
        keyword = "X"
    else:
        keyword = "Z"

    Le = -1
    Le_f = alpha / Df
    Le_o = alpha / Do
    if phi <= 1.0:
        Le = alpha/Df
    else:
        Le = alpha/Do

    with open(output, 'a') as r:
        writer = csv.writer(r)
        writer.writerow([alpha, Df, Do, Le_f, Le_o, Le, phi, keyword])

    return Le

# made by A.Tsunoda
#========================================
mechanism = './C8_SIP-Gr2.0-s2_mech/chem.cti'
tran = ct.Solution(mechanism)
phi = 1.0

fuel = 'iC8H18'
oxid = 'O2:1, N2:3.76'

temp = 300.0 # K
temp = 600.0 # K
pres = 1.0 # atm

isInertConst = False
isInertConst = True

tran.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxid)
tran.transport_model = 'Multi'
#tran()

k = tran.thermal_conductivity
cp = tran.cp_mass
rho = tran.density_mass
thermal_diffusivity = k/cp/rho

idx_fuel = tran.species_index(fuel)
idx_oxid = tran.species_index('O2')

multi_D = tran.mix_diff_coeffs
fuel_diffusivity = multi_D[idx_fuel]
oxid_diffusivity = multi_D[idx_oxid]

print('phi:', phi, 
      'thermal_diffusivity:', thermal_diffusivity, 
      'fuel_diffusivity:', fuel_diffusivity, 
      'oxid_diffusivity:', oxid_diffusivity)

output = "Ledata.csv"
with open(output, 'w') as r:
    writer = csv.writer(r)
    writer.writerow(["alpha [m2/s]", "D(fuel) [m2/s]", "D(O2)[m2/s]", "Le_f", "Le_o", "Le", "phi", "ZorX"])
    Le = write2csv(thermal_diffusivity, fuel_diffusivity, oxid_diffusivity, phi, isInertConst)
    print("phi = {:.3f}: Le = {:.3f}".format(phi, Le))

