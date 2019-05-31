import chemutils.species as spc
from chemutils.simulation import SimulationSystem, SPCE
from chemutils.simulation.gromacs import NVT, NPT

NVT.mdp['nsteps'] = 200000
NPT.mdp['nsteps'] = 500000

# set temperature, pressure, and initial density
t = 298.15
p = 100e3
rho = spc.H2O.density(t) / spc.H2O.molar_mass

# define initial system with size = 3.2 nm
liq = SimulationSystem({SPCE:rho}, T=t, P=p, size=3.2e-9)

# do NVT equilibration, then NPT
liq = NVT(liq)
liq = NPT(liq)

print('liquid density: {}'.format(liq.density.mean()))
