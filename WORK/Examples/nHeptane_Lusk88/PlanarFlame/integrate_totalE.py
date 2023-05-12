
import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt
import subprocess


rig = 200.0 # micro m
tig = 200.0 # micro s
Eig = 1.0 # mJ

rig_m = rig * 1e-6
tig_s = tig * 1e-6
Eig_J = Eig * 1e-3
tig_s = 1.0
Eig_J = 1.0

x = sp.Symbol('x')
r = sp.Symbol('r')
Q_r = Eig_J/((4.0/3.0)*np.pi*rig_m**3*tig_s)*sp.exp(-(np.pi/4.0)*(r/rig_m)**6)

p = 6

# Spherical
alpha_r = 2
k_r = 4.0*np.pi
g_r = math.gamma(1+(alpha_r+1)/p)**int(p/(alpha_r+1))
Q_r = Eig_J / (k_r/(alpha_r+1) * rig_m**(alpha_r+1) * tig_s) * sp.exp(-g_r * (r/rig_m)**p)
print(Q_r.subs([(r, 0)]))

Q_r_0 = Q_r * 4.0*np.pi*r**2
Q_r_int = sp.integrate(Q_r_0, r)
Q_r_total = sp.integrate(Q_r_0, (r, 0.0, sp.oo))
print('Q_r_total, A.Frendi(1990)', float(Q_r_total))

# Spherical, W. Zhang et al. CNF(2012)
Q_x = Eig_J/(np.pi**1.5*rig_m**3*tig_s)*sp.exp(-(r/rig_m)**2)
Q_x_int = sp.integrate(Q_x, r)
Q_x_total = sp.integrate(Q_x, (r, 0.0, sp.oo))
print('Q_x_total, G.Xiao(2022)', float(Q_x_total))

# Planar
alpha_x = 0
k_x = 1
g_x = math.gamma(1+(alpha_x+1)/p)**int(p/(alpha_x+1))
Q_x = Eig_J / (k_x/(alpha_x+1) * rig_m**(alpha_x+1) * tig_s) * sp.exp(-g_x * (r/rig_m)**p)
Q_x_int = sp.integrate(Q_x, r)
Q_x_total = sp.integrate(Q_x, (r, 0.0, sp.oo))
print('Q_x_total, A.Frendi(1990)', float(Q_x_total))

# Planar, G. Xiao et al. (2022)
Q_x = Eig_J/((4.0/3.0)*np.pi*rig_m**3*tig_s)*sp.exp(-(r/rig_m)**8)
Q_x_int = sp.integrate(Q_x, r)
Q_x_total = sp.integrate(Q_x, (r, 0.0, sp.oo))
print('Q_x_total, G.Xiao(2022)', float(Q_x_total))


rm_list = np.linspace(0, 2*rig_m, 100)
rmm_list = np.linspace(0, 2*rig_m*1e+3, 100)
Q_r = Eig_J / (k_r/(alpha_r+1) * rig_m**(alpha_r+1) * tig_s) * np.exp(-g_r * (rm_list/rig_m)**p)
Q_x = Eig_J / (k_x/(alpha_x+1) * rig_m**(alpha_x+1) * tig_s) * np.exp(-g_x * (rm_list/rig_m)**p)

Q_r_total_list = []
V_r = []
for r_r in rm_list:
    Q_r_total_list.append(Q_r_int.subs([(r, r_r)]) - Q_r_int.subs([(r, 0.0)]))
    V_r.append((4.0/3.0)*np.pi*r_r**3)
Q_x_total_list = []
V_x = []
for x_x in rm_list:
    Q_x_total_list.append(Q_x_int.subs([(r, x_x)]) - Q_x_int.subs([(r, 0.0)]))
    V_x.append(x_x)

print(len(rm_list))
print(len(Q_r))
print(len(Q_x))

plt.rcParams["font.family"] = "Arial"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['xtick.minor.width'] = 1.4
plt.rcParams['ytick.minor.width'] = 1.4
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["xtick.minor.size"] = 6
plt.rcParams["ytick.minor.size"] = 6
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['figure.figsize'] = 10, 6

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11,8))
ax0, ax1 = ax[0].twinx(), ax[1].twinx()

ax[0].plot(rmm_list, Q_r, 'r')
ax[1].plot(rmm_list, Q_x, 'r')
ax0.plot(rmm_list, V_r, 'b')
ax1.plot(rmm_list, V_x, 'b')
# ax0.plot(rmm_list, Q_r_total_list, 'b')
# ax1.plot(rmm_list, Q_x_total_list, 'b')

ax[0].set_xlim(0, None)
ax[1].set_xlim(0, None)
ax[0].tick_params(axis='y', labelcolor='r')
ax[1].tick_params(axis='y', labelcolor='r')
ax0.tick_params(axis='y', labelcolor='b')
ax1.tick_params(axis='y', labelcolor='b')
ax[0].set_xlabel('r [mm]')
ax[1].set_xlabel('x [mm]')
ax[0].set_ylabel('E [J/m$^3$s]', c='r')
ax[1].set_ylabel('E [J/m$^3$s]', c='r')
ax0.set_ylabel('Volume [m$^3$]', c='b')
ax1.set_ylabel('Volume [m$^3$]', c='b')
figname = 'totalE.png'
fig.tight_layout()
fig.savefig(figname)
subprocess.call(['imgcat', figname])


