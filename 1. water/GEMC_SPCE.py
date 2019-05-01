from scipy.constants import R, Boltzmann

import chemutils.species as spc
from chemutils.simulation import SimulationSystem, SPCE
from chemutils.simulation.cassandra import CassandraMethod

EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=10000, kappa=20, swap={SPCE})
GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=10000, kappa=20, swap={SPCE})

t = 298.15
rho = spc.H2O.density(t) / spc.H2O.molar_mass

N_vap = 50
d_vap = (N_vap * Boltzmann * t / (spc.H2O.pressure(t)))**(1/3)
rho_vap = spc.H2O.pressure(t) / (R * t)

# define initial systems
liq = SimulationSystem({SPCE:rho}, T=t, size=3.2e-9)
vap = SimulationSystem({SPCE:rho_vap}, T=t, size=d_vap)

liq, vap = EQ(liq, vap)
liq, vap = GEMC(liq, vap)

