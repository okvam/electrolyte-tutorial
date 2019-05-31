# electrolyte-tutorial
Tutorial for chemical theory and simulation of electrolyte solutions.

## Part 0. Introduction
Electrolytes are compounds which upon contact with a solvent dissociate (from greek _lysis_, loosening) to form electrically charged species (_electro_). The most familiar electrolyte is perhaps the table salt found in your kitchen, but electrolytes are of fundamental important for physical processes ranging from the regulation of biological systems to the battery technology powering our electronic devices. Thermodynamics of electrolyte solutions have a long history of scientific investigation. 

### dependencies

This tutorial assumes you have the simulation packages Gromacs and Cassandra installed and callable with the shell commands 'gmx' and 'cassandra', respectively. For installation instructions, please see their respective home pages. Further, we make use of the python package 'chemutils'. For acess to this package please contact @okvam on github. 

## Part 1. Water
The best studied solvent for electrolyte solutions is water, partly because of its prevalence and partly because of the high solubility of most simple electrolyte species in this solvent. Water alone is a highly complex fluid, showing many unusual thermophysical properties (1).

In the first part of this tutorial we will simulate some properties of water using the SPC/E water model (2). We will answer the following questions;

* How does the liquid density of SPC/E water compare to real water at 298.15K and 1 bar?
* Does SPC/E water have a temperature of maximum liquid density at 1 bar, like real water? 
* What is the saturated vapor pressure of the SPC/E model of water at 298.15K? What happens when temperature increases?

### NPT simulation

We will be using the simulation package Gromacs (3) to study the liquid densities of water. Locate the python scripts for setting up simulations for water in the directory **1. water**, and open the file **NPT_SPCE.py**. First, we load some pre-packaged python modules: 

    import chemutils.species as spc
    from chemutils.simulation import SimulationSystem, SPCE
    from chemutils.simulation.gromacs import NVT, NPT

Then, we adjust some of the built-in simulation methods, reducing the time step for NVT simulation and increasing the length of both NVT and NPT simulation. The grmacs methods NVT and NPT each have the python property <GromacsMethod>.mdp, a dictionary with keys and values corresponding to options for the [.mdp gromacs input file](http://manual.gromacs.org/documentation/2018.6/user-guide/mdp-options.html). 

    NVT.mdp['dt'] = 0.0005
    NVT.mdp['nsteps'] = 200000
    NPT.mdp['nsteps'] = 500000

We also need to decide on the temperature, pressure, and initial molar density of water. We are interested in a temperature of 298.15K and pressure of one bar, equal to 100x10<sup>3</sup> Pa. For initial molar density we use the liquid density of real water, from the pre-packaged _chemutils.species_ module. 

    t = 298.15
    p = 100e3
    rho = spc.H2O.density(t) / spc.H2O.molar_mass

Our simulation system can now be defined. The SimulationSystem class takes a python dictionary (or another SimulationSystem) as first argument. It can also accept T, P, and size as keyword arguments. Here, we set all three, since we want isothermal-isobaric conditions. Then we call each simulation method using the simulation system we just created as argument:

    # define initial system, size = 3.2 nm
    liq = SimulationSystem({SPCE:rho}, T=t, P=p, size=3.2e-9)

    # do NVT equilibration, then NPT
    liq = NVT(liq)
    liq = NPT(liq)

These simulations should take about 1 hour to run on a modern desktop computer. First, a simulation box of 3.2 nm is filled with SPC/E water molecules at a density matching that of real water for 298.15K. Then Gromacs simulates the box for 100 ps in the NVT ensemble. Once this has completed, the final configuration (system coordinates) of the water box is used as the starting configuration for the NPT simulation. This NPT simulations runs for 1 ns. 

### Gibbs ensemble simulation

Having looked at liquid properties at 298.15K and 1 bar, we will now consider vapor-liquid equilibria of water. We will be using the simulation package Cassandra (4) to determine the the saturated vapor pressure for a given temperature. Locate the file **GEMC_SPCE.py**. This file is similar in structure to the one used for liquid density, with some differences shown below. Firstly, we use a different set of simulation methods:

    from chemutils.simulation.cassandra import CassandraMethod

    NPT = CassandraMethod(method='npt', equilibration=True, sweeps=1000, kappa=1000)
    EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=5000, kappa=20, swap={SPCE})
    GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=5000, kappa=20, swap={SPCE})

In this case we build our own methods using the CassandraMethod class. We need three methods: one for generating a reasonable starting configuration, another to equilibrate the vapor-liquid systems, and a third for performing the Gibbs ensemble simulation. Since we are simulating the vapor phase, we also need to decide on the size of this box. We use a particle count of 100 SPCE molecules, with initial box size and density guessed from the properties of real water at 298.15K. 

    N_vap = 50
    d_vap = (N_vap * Boltzmann * t / (spc.H2O.pressure(t)))**(1/3)
    rho_vap = methanol.pressure(t) / (R * t)

We then define our simulation systems. We need two, one for the liquid phase and one for the vapor phase: 

    liq = SimulationSystem({SPCE:rho_liq}, T=t, size=3.2e-9)
    vap = SimulationSystem({SPCE:rho_vap}, T=t, size=d_vap)

Finally we run our simulations. After simulating using the NPT method we delete the pressure constraint, as this is no longer required. Then, we simulate in the Gibbs ensemble, this time passing both simulation systems as arguments: 

    liq = NPT(liq)

    del liq.constraints[liq.P]

    liq, vap = EQ(liq, vap)
    liq, vap = GEMC(liq, vap)

These simulations take a bit longer to run, about 2-3 hours on a modern desktop computer. First, a simulation box of 3.2 nm is filled with SPC/E water molecules at a density matching that of real water at 298.15K. Cassandra simulates the box for 1000 MC sweeps in the NPT ensemble. Once this has completed, the final configuration of the liquid is used together with a (randomly assigned) vapor phase as the starting configurations for Gibbs ensemble simulation. 5000 MC cycles will be used for equilibration, and another 5000 cycles for the final Gibbs ensemble simulation.

After the first of the two simulations have completed (or while its still running, if you are impatient), look at the files **sim.out.box1.prp** and **sim.out.box2.prp** inside the simulation directory. These contain properties of the liquid and vapor systems as a function of simulation cycle. What is the initial behaviour of the two systems? Which properties converge quickly, and which take a long time? How can you be confident the system is stable?

### Further Simulations

Using the scripts provided for Gromacs and Cassandra as a starting point, simulate SPC/E water across a range of temperatures to investigate (i) whether this water model has a density maximum at 1 bar, like real water, and (ii) how the saturated vapor pressure varies as you increase temperature.

### references

1. [Anomalous properties of water, by Martin Chaplin, London South Bank University](http://www1.lsbu.ac.uk/water/water_anomalies.html)
2. [H.J.C. Berendsen, J.R. Grigera, and T.P. Straatsma, J. Phys. Chem. 91, 6269-6271](http://dx.doi.org/10.1021/j100308a038)
3. [Gromacs home page](http://manual.gromacs.org/)
4. [Cassandra home page](https://cassandra.nd.edu/)

## Part 2. Methanol
Methanol is the smallest alcohol, with one of the H-atoms from water replaced with an organic CH<sub>3</sub> group. Many simple electrolytes are reasonably soluble also in methanol, albeit less so than in water. 

In the second part of this tutorial we will simulate some properties of methanol and water + methanol mixtures, using the TraPPE model of methanol (5). We will answer the following questions:

* How does the liquid density of TraPPE methanol compare to real methanol at 298.15K and 1 bar?
* What is the saturated vapor pressure of TraPPE methanol at 298.15K?
* Is the simulated vapor-liquid composition diagram in good agreement with experiment for water + methanol at 1 bar? 

### NPT simulation

We will use similar methodology to simulate the properties of methanol as we did for water. Locate the python scripts for setting up simulations for water in the directory **2. methanol**, and open the file **NPT_methanol.py**. This time, we import methanol instead of SPCE, and we use Cassandra rather than Gromacs for NPT simulation (the results should be the same regardless of which simulation package is used). 

    import chemutils.species as spc
    from chemutils.simulation import SimulationSystem, methanol
    from chemutils.simulation.cassandra import CassandraMethod

We create two Cassandra methods, one for equilibration and one for the main NPT simulation:

    NPT0 = CassandraMethod(method='npt', name='eq', equilibration=True, sweeps=1000, kappa=100)
    NPT = CassandraMethod(method='npt', name='npt', sweeps=5000)

For temperature, pressure, and initial molar density, we use the same approach as for water. We then create our simulation system object and run the simulation. 

    # define initial system with size = 3.2 nm
    liq = SimulationSystem({methanol:rho}, T=t, P=p, size=3.2e-9)

    # do NVT equilibration, then NPT
    liq = NPT0(liq)
    liq = NPT(liq)

After the simulations have completed, we print the average density of the simulation system

    print('liquid density: {}'.format(liq.density.mean()))

The method <SimulationSystem>.density returns a [numpy](http://www.numpy.org) array with instantaneous system densities over the course of the simulation. Numpy arrays come with many useful methods, such as <NumpyArray>.mean(), which computes the geometric mean of all array entries. 

### Gibbs ensemble simulation

As for water, we also want to consider vapor-liquid equilibria of methanol. Locate the file **GEMC_methanol.py**. This file is similar to the one used for water, but this time we will use Gromacs to generate an initial configuration, and then Cassandra to calculate VLE properties. This means we need to import both Gromacs and Cassandra methods:

    from chemutils.simulation.gromacs import NPT
    from chemutils.simulation.cassandra import CassandraMethod

    NPT.mdp['nsteps'] = 200000

    EQ = CassandraMethod(method='gemc', equilibration=True, sweeps=5000, kappa=20, swap={methanol})
    GEMC = CassandraMethod(method='gemc', name='vpd', sweeps=5000, kappa=20, swap={methanol})

We use the NPT Gromacs method, and create our own Cassandra methods as before. The liquid and vapor systems are set up and simulated as for water, with the liquid system initially constrained to 100x10<sup>3</sup> Pa. We then run our simulations. Notice that the liquid system returned from the Gromacs simulation (NPT) is also compatible with Cassandra. The inverse is also true, systems generated by Cassandra can be used as starting points for Gromacs simulations. Finally, we print the liquid density, vapor density, and vapor pressure for this temperature. 

    print('liquid density: {}'.format(liq.density.mean()))
    print('vapor density: {}'.format(vap.density.mean()))
    print('vapor pressure: {}'.format(vap.pressure.mean()))

Saturated vapor and liquid densities as a function of temperature together form the _vapor-liquid coexistence curve_, the fluid phase diagram of the pure fluid up to the critical point. Vapor-liquid coexistence curves are often reported together with the saturated vapor pressure as a function of temperature, and together give a good picture of the phase behavior of the pure fluid. 

### 2-component Gibbs ensemble simulation

Separation processes often take advantage of the differences in composition of fluids between different phases. Perhaps the most well-known fluid separation process is distillation, where lighter components are concentrated in the vapor phase while heavy components are concentrated in the liquid phase. For a given pressure, the _distillation curve_ gives the composition of the liquid and vapor as a function of temperature. For methanol and water at 101.325 kPa, the distillation curve is a simple envelope connecting the boiling temperature of water to that of methanol, shown in Figure 1. For other species such as ethanol, the two curves may intersect, creating an [azeotrope](https://en.wikipedia.org/wiki/Azeotrope).

![methanol distillation curve](figures/methanol.jpg?raw=true "methanol + water") 

_Figure 1: Distillation curve of methanol (1) + water (2) at 101.325 kPa. Points indicate experimental values, line is NRTL correlation. Figure from [Gao et al. 2012](http://dx.doi.org/10.1155/2012/641251)._

Using the SPCE and TraPPE models, we will simulate the vapor-liquid composition of methanol + water, and see how our model compares against experimental values. We guess the initial composition of liquid and vapor phases using Raoult's law, which is usually not too far away. 

### References

5. [B. Chen, J.J. Potoff, and J.I. Siepmann, J. Phys. Chem. B 105, 3093-3104](http://dx.doi.org/10.1021/jp003882x)
6. binary VLE 

## Part 3. Vapor pressure depression
Electrolyte species are for all intents and purposes non-volatile. From Raoult's law we expect the vapor pressure _p_ of a liquid mixture to be proportional to the sum of it components scaled by their liquid mole fraction x<sub>i</sub>, 

P<sub>i</sub> = x<sub>i</sub> * P<sub>i</sub><sup>s</sup>

where P<sub>i</sub><sup>s</sup> is the saturated vapor pressure of the neat solvent. In the case of electrolyte solutions, only the solvent is expected to contribute to the overall vapor pressure of solution - in other words, the electrolyte has a saturated vapor pressure of zero.  (show derivation of osmotic coefficient).

In the third part of this tutorial, we will simulate the vapor pressure depression of water with increasing concentration of LiBr, using the alkali halide model of Joung and Cheatham, 2008. We will compare the resulting osmotic coefficients with those correlated by Pitzer and Mayorga 19XX for real water. 

### references

7. Joung and Cheatham 2008
8. Pitzer and Mayorga 19XX

## Part 4. Mixed solvents
Most applications of electrolytes involve more than one solvent. In some cases, the electrolyte species act to modulate the thermodynamic behaviour of the mixed system. An example of this is purification of alcohols from water, where the addition of electrolyte causes phase separation of solvents which would otherwise be miscible. In other cases, the electrolytes themselves are the focus, such as in conductive electrochemical systems. Regardless of application the thermodynamics of electrolytes and solvents are intimately linked. 

### references 

9. ???