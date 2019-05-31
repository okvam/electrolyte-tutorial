import chemutils.species as spc
from chemutils.simulation import SimulationSystem, methanol
from chemutils.simulation.cassandra import CassandraMethod

NPT0 = CassandraMethod(method='npt', name='eq', equilibration=True, sweeps=1000, kappa=100)
NPT = CassandraMethod(method='npt', name='npt', sweeps=5000)

# set temperature, pressure, and initial density
t = 298.15
p = 100e3
rho = spc.methanol.density(t) / spc.methanol.molar_mass

# define initial system with size = 3.2 nm
liq = SimulationSystem({methanol:rho}, T=t, P=p, size=3.2e-9)

# do NVT equilibration, then NPT
liq = NPT0(liq)
liq = NPT(liq)

print('liquid density: {}'.format(liq.density.mean()))
