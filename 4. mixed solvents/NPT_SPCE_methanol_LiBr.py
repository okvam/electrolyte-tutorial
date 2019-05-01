import sys
from scipy.constants import R, Boltzmann

import chemutils.species as spc
from chemutils.simulation import SimulationSystem, methanol, SPCE
from chemutils.simulation.species.JC2008.SPCE import sodium, iodide
from chemutils.simulation.cassandra import CassandraMethod
from chemutils.simulation.gromacs.method import NVT, NPT

NVT.mdp['dt'] = 0.0005
NVT.mdp['nsteps'] = 200000
NPT.mdp['nsteps'] = 500000

EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=10000, kappa=20, swap={methanol})
GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=50000, kappa=20, swap={methanol})

def simulate(m):

    t = 298.15

    # liquid box length ~ 4.0 nm
    d_liq = 4.2e-9
    rho_liq = (spc.methanol.density(t) / spc.methanol.molar_mass) * 1/(1 + m * spc.methanol.molar_mass)
    rho_salt = m * rho_liq * spc.methanol.molar_mass

    N_vap = 50
    d_vap = (N_vap * Boltzmann * t / (methanol.pressure(t)))**(1/3)
    rho_vap = methanol.pressure(t) / (R * t)

    # define initial systems
    liq = SimulationSystem({methanol:rho_liq, sodium:rho_salt, iodide:rho_salt}, T=t, P=methanol.pressure(t), size=d_liq)
    vap = SimulationSystem({methanol:rho_vap, sodium:0, iodide:0}, T=t, size=d_vap)

    # do NPT equilibration for liquid phase, then remove pressure constraint
    liq = NVT(liq)
    liq = NPT(liq)

    del liq.constraints[liq.P]

    # do GEMC equilibration to establish weights
    liq, vap = EQ(liq, vap)

    # do GEMC
    liq, vap = GEMC(liq, vap)

if __name__ == "__main__":
    simulate(float(sys.argv[1]))
