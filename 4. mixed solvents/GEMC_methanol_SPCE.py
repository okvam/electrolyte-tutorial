import sys
from scipy.constants import R, Boltzmann

import chemutils.species as spc
from chemutils.simulation import SimulationSystem, methanol, SPCE
from chemutils.simulation.cassandra import CassandraMethod
from chemutils.simulation.gromacs import NVT, NPT

NVT.mdp['dt'] = 0.001
NVT.mdp['nsteps'] = 200000
NPT.mdp['nsteps'] = 500000

EQ = CassandraMethod(method='gemc_npt', equilibration=True, sweeps=10000, kappa=20, swap={methanol, SPCE})
GEMC = CassandraMethod(method='gemc_npt', name='gemc', sweeps=50000, kappa=20, swap={methanol, SPCE})

def simulate(m):

    # Raoult's law: x1 * P1(t) + x2 * P2(t) = 101.325
    t = 353.15
    x1 = (101.325e3 - SPCE.pressure(t)) / (methanol.pressure(t) - SPCE.pressure(t))
    x2 = 1 - x1
    p1, p2 = x1 * methanol.pressure(t), x2 * SPCE.pressure(t)

    print('p1 + p2 = {}'.format(p1 + p2))

    # liquid box length ~ 3.5 nm
    d_liq = 3.5e-9
    rho_liq = 1 / (x1 * spc.methanol.density / spc.methanol.molar_mass + x1 * spc.methanol.density / spc.methanol.molar_mass)

    N_vap = 200
    d_vap = (N_vap * Boltzmann * t / (101.325e3))**(1/3)
    rho1 = p1 / (R * t)
    rho2 = p2 / (R * t)

    # define initial systems
    liq = SimulationSystem({methanol:x1 * rho_liq, SPCE:x2 * rho_liq}, T=t, P=methanol.pressure(t), size=d_liq)
    vap = SimulationSystem({methanol:rho1, SPCE:rho2}, T=t, size=d_vap)

    # do NPT equilibration for liquid phase, then remove pressure constraint
    liq = NVT(liq)
    liq = NPT(liq)

    # do GEMC equilibration to establish weights
    liq, vap = EQ(liq, vap)

    # do GEMC
    liq, vap = GEMC(liq, vap)

    # print mean mole fraction compositions
    print('x1 = {}'.format((liq.molarity[methanol] / (liq.molarity[methanol] + liq.molarity[SPCE])).mean())) 
    print('y1 = {}'.format((vap.molarity[methanol] / (vap.molarity[methanol] + vap.molarity[SPCE])).mean()))

if __name__ == "__main__":
    simulate(float(sys.argv[1]))
