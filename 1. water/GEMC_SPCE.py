from scipy.constants import R, Boltzmann

import chemutils.species as spc
from chemutils.simulation import SimulationSystem, SPCE
from chemutils.simulation.cassandra import CassandraMethod

NPT = CassandraMethod(method='npt', equilibration=True, sweeps=1000, kappa=1000)
EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=5000, kappa=20, swap={SPCE})
GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=5000, kappa=20, swap={SPCE})

t = 298.15
rho = spc.H2O.density(t) / spc.H2O.molar_mass

N_vap = 100
d_vap = (N_vap * Boltzmann * t / (spc.H2O.pressure(t)))**(1/3)
rho_vap = spc.H2O.pressure(t) / (R * t)

# define initial systems
liq = SimulationSystem({SPCE:rho}, T=t, P=spc.H2O.pressure(t), size=3.2e-9)
vap = SimulationSystem({SPCE:rho_vap}, T=t, size=d_vap)

liq = NPT(liq)

del liq.constraints[liq.P]

liq, vap = GEMC0(liq, vap)
liq, vap = GEMC1(liq, vap)

