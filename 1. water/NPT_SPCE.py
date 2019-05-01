import chemutils.species as spc
from chemutils.simulation import SimulationSystem, SPCE
from chemutils.simulation.gromacs.method import NVT, NPT

NVT.mdp['dt'] = 0.0005
NVT.mdp['nsteps'] = 200000
NPT.mdp['nsteps'] = 500000

t = 298.15
p = spc.H2O.pressure(t)
rho = spc.H2O.density(t) / spc.H2O.molar_mass

# define initial system, size = 3.2 nm
liq = SimulationSystem({SPCE:rho}, T=t, P=p, size=3.2e-9)

# do NVT equilibration, then NPT
liq = NVT(liq)
liq = NPT(liq)

