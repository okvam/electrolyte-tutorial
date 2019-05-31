from scipy.constants import R, Boltzmann

import chemutils.species as spc
from chemutils.simulation import SimulationSystem, methanol
from chemutils.simulation.gromacs import NPT
from chemutils.simulation.cassandra import CassandraMethod

NPT.mdp['nsteps'] = 200000

EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=5000, kappa=20, swap={methanol})
GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=5000, kappa=20, swap={methanol})

t = 298.15
p = spc.methanol.pressure(t)
rho = spc.methanol.density(t) / spc.methanol.molar_mass

N_vap = 100
d_vap = (N_vap * Boltzmann * t / (p))**(1/3)
rho_vap = spc.methanol.pressure(t) / (R * t)

# define initial systems
liq = SimulationSystem({SPCE:rho}, T=t, P=p, size=3.2e-9)
vap = SimulationSystem({SPCE:rho_vap}, T=t, size=d_vap)

liq = NPT(liq)

del liq.constraints[liq.P]

liq, vap = GEMC0(liq, vap)
liq, vap = GEMC1(liq, vap)

print('liquid density: {}'.format(liq.density.mean()))
print('vapor density: {}'.format(vap.density.mean()))
print('vapor pressure: {}'.format(vap.pressure.mean()))
