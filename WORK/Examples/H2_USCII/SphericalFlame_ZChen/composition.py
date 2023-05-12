
import cantera as ct
import numpy as np

oxid = 'O2:1, N2:3.76'
oxid = 'O2:1, AR:3.76'

chemPath_cti = './USCII_H2CO/chem.cti'
fuel = 'H2'
equivalence_ratio = 2.0
gas = ct.Solution(chemPath_cti)
gas.set_equivalence_ratio(equivalence_ratio, \
                          fuel=fuel, oxidizer=oxid)
print(chemPath_cti)
print(equivalence_ratio, gas.mole_fraction_dict())

phi = 2.0
rs = 0.5
alpha = 0

A = 0.3*phi/(phi+rs)
B = 0.3*rs/(phi+rs)
C = alpha
D = 0.7-alpha
sum_coef = A+B+C+D
print('0.3(phi/(phi+rs)Fule+rs/(phi+rs)O2)+alpha*He+(0.7-alpha)Ar')
print('0.3*phi/(phi+rs) =', A/sum_coef)
print('0.3*rs/(phi+rs) =', B/sum_coef)
print('alpha =', C/sum_coef)
print('0.7-alpha =', D/sum_coef)

# fuel
# oxid
# He
# Ar

