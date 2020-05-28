import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.constants import e, k as kB
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
rc('text', usetex=True)


# Reactant masses in atomic mass units (u).
masses = {'D': 2.014, 'T': 3.016, '3He': 3.016, 
          'O16':15.9949, 'O18': 17.999, 'p': 1.0072, '17F': 17.0020, }

# Energy grid, 1 – 1000 keV, evenly spaced in log-space.
Egrid = np.logspace(0, 4, 100)

def read_xsec(filename,
              collider_particle_mass,
              target_particle_mass,
              energy_scaling_factor=1.e3,
              xs_scaling_factor=1.e-28):
    """Read in cross section from filename and interpolate to energy grid."""

    E, xs = np.genfromtxt(filename, comments='#', unpack=True)

    m1, m2 = collider_particle_mass, target_particle_mass
    E *= m1 / (m1 + m2)

    E = E*energy_scaling_factor
    xs = xs*xs_scaling_factor

    return E,xs


# units are in MeV for energy and mili barns
# https://tendl.web.psi.ch/tendl_2017/proton_file/O/O016/tables/residual/rp009016.tot
O16p_E, O16p_xs = read_xsec(filename='O16_p-17F.txt',
                    collider_particle_mass=masses['O16'],
                    target_particle_mass=masses['p'],
                    xs_scaling_factor=1e-31,
                    energy_scaling_factor=1e6)

# units are in MeV for energy and mili barns
# https://tendl.web.psi.ch/tendl_2019/proton_file/O/O018/tables/residual/rp009019.tot
O18p_E, O18p_xs = read_xsec(filename='O18_p-19F.txt',
                    collider_particle_mass=masses['O18'],
                    target_particle_mass=masses['p'],
                    xs_scaling_factor=1e-31,
                    energy_scaling_factor=1e6)


# # D + T -> α + n
# from endf site? units are in barns and eV
DT_E, DT_xs = read_xsec(filename='D_T_-_a_n.txt',
                  collider_particle_mass=masses['D'],
                  target_particle_mass=masses['T'])


# # D + D -> T + p
# from endf site? units are in barns and eV
DDa_E, DDa_xs = read_xsec(filename='D_D_-_T_p.txt',
                   collider_particle_mass=masses['D'],
                   target_particle_mass=masses['D'])

# # D + D -> 3He + n
# from endf site? units are in barns and eV
DDb_E, DDb_xs = read_xsec(filename='D_D_-_3He_n.txt',
                   collider_particle_mass=masses['D'],
                   target_particle_mass=masses['D'])
# # Total D + D fusion cross section is due to equal contributions from the
# # above two processes.


# D + 3He -> α + p
# from endf site? units are in barns and eV
DHe_E, DHe_xs = read_xsec(filename='D_3He_-_4He_p.txt',
                   collider_particle_mass=masses['D'],
                   target_particle_mass=masses['3He'])

fig, ax = plt.subplots()
ax.loglog(DT_E, DT_xs, lw=2, label='$\mathrm{D,T}$')
ax.loglog(DDa_E, DDa_xs, lw=2, label='$\mathrm{D,D (3He+n)}$')
ax.loglog(DDb_E, DDb_xs, lw=2, label='$\mathrm{D,D (T+p)}$')
ax.loglog(DHe_E, DHe_xs, lw=2, label='$\mathrm{D,^3He}$')
ax.loglog(O18p_E, O18p_xs, lw=2, label='$\mathrm{O^{18},P}$')
ax.loglog(O16p_E, O16p_xs, lw=2, label='$\mathrm{O^{16},P}$')
# ax.grid(True, which='both', ls='-')
ax.set_xlim(1, 1e9)
xticks= np.array([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
ax.set_xticks(xticks)
ax.set_xticklabels(["{:.0e}".format(float(i)) for i in xticks])


ax.set_xlabel('E (Center of Mass) [keV]')

# # A second x-axis for energies as temperatures in millions of K
ax2 = ax.twiny()
ax2.set_xscale('log')
ax2.set_xlim(1,1000)
xticks2 = np.array(1e3 * xticks* (1/(kB/e)))
ax2.set_xticks(xticks * kB/e)
ax2.set_xticklabels(["{:.1e}".format(float(i)) for i in xticks2])
# # ax2.set_xlabel('$T$ /million K')
ax2.set_xlabel('T [K]')


ax.set_ylabel('$\sigma\;[\mathrm{m^2}]$')
ax.set_ylim(1.e-32, 1.e-27)

ax.legend()
# plt.savefig('fusion-xsecs.png')
plt.show()