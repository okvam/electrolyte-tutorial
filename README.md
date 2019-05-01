# electrolyte-tutorial
Tutorial for chemical theory and simulation of electrolyte solutions.

## part 0. introduction
Electrolytes are compounds which upon contact with a solvent dissociate (from greek _lysis_, loosening) to form electrically charged species (_electro_). The most familiar electrolyte is perhaps the table salt found in your kitchen, but electrolytes are of fundamental important for physical processes ranging from the regulation of biological systems to the battery technology powering our electronic devices. Thermodynamics of electrolyte solutions have a long history of scientific investigation. 

### dependencies

This tutorial assumes you have the simulation packages Gromacs and Cassandra installed and callable with the shell commands 'gmx' and 'cassandra', respectively. For installation instructions, please see their respective home pages. Further, we make use of the python package 'chemutils'. For acess to this package please contact @okvam on github. 

## part 1. water
The best studied solvent for electrolyte solutions is water, partly because of its prevalence and partly because of the high solubility of most simple electrolyte species in this solvent. Water alone is a highly complex fluid, showing many unusual thermophysical properties.

In the first part of this tutorial we will simulate some properties of water using the SPC/E water model. We will answer the following questions;

* How does the liquid density of SPC/E water compare to real water at 298.15K and 1 bar?
* Does SPC/E water have a temperature of maximum liquid density at 1 bar like real water? 
* What is the saturated vapor pressure of the SPC/E model of water at 298.15K? What happens when temperature increases?

### simulations

We will be using the simulation package Gromacs to study liquid densities and Cassandra for vapor pressures. Locate the python scripts for setting up simulations for water in the directory **1. water**, and open the file **NPT_SPCE.py**. First, we load some pre-packaged python modules. 

    import chemutils.species as spc
    from chemutils.simulation import SimulationSystem, SPCE
    from chemutils.simulation.gromacs.method import NVT, NPT

We adjust some of the built-in simulation methods, reducing the time step for NVT simulation and increasing the length of both NVT and NPT simulation. The python property <GromacsMethod>.mdp is a dictionary with keys and values corresponding to options for the [.mdp](http://manual.gromacs.org/documentation/2018.6/user-guide/mdp-options.html) gromacs input file. 

    NVT.mdp['dt'] = 0.0005
    NVT.mdp['nsteps'] = 200000
    NPT.mdp['nsteps'] = 500000

Then, we set the temperature, pressure, and initial molar density of water. Note that we use physical properties from water in the pre-packaged _chemutils.species_ module. 

    t = 298.15
    p = spc.H2O.pressure(t)
    rho = spc.H2O.density(t) / spc.H2O.molar_mass

The simulation system can now be defined. The SimulationSystem class takes a python dictionary (or another SimulationSystem) as first argument. It can also accept T, P, and size as keyword arguments. Here, we set all three, since we want isothermal-isobaric conditions. Then we call each simulation method using the simulation system as argument.

    # define initial system, size = 3.2 nm
    liq = SimulationSystem({SPCE:rho}, T=t, P=p, size=3.2e-9)

    # do NVT equilibration, then NPT
    liq = NVT(liq)
    liq = NPT(liq)

These simulations should take about 1 hour to run on a modern desktop computer. First, a simulation box of 3.2 nm is filled with SPC/E water molecules at a density matching that of real water for 298.15K. Then Gromacs simulates the box for 100 ps in the NVT ensemble. Once this has completed, the final configuration (system coordinates) of the water box is used as the starting configuration for the NPT simulation. This NPT simulations runs for 1 ns. 

### references

1. Water properties
2. SPCE
3. Gromacs
4. Cassandra

## part 2. methanol
Methanol is the smallest alcohol, with one of the H-atoms from water replaced with an organic CH<sub>3</sub> group. Many simple electrolytes are reasonably soluble also in methanol, albeit less so than in water. 

In the second part of this tutorial we will simulate some properties of methanol and water + methanol mixtures, using the TraPPE model of methanol. We will answer the following questions:

* How does the liquid density of TraPPE methanol compare to real methanol at 298.15K and 1 bar?
* What is the saturated vapor pressure of TraPPE methanol at 298.15K?
* Is the simulated vapor-liquid composition diagram in good agreement with experiment for water + methanol at 1 bar? 

### references

5. Trappe methanol
6. binary VLE 

## part 3. vapor pressure depression
Electrolyte species are for all intents and purposes non-volatile. From Raoult's law we expect the partial vapor pressure P<sub>i</sub> of a component in a liquid mixture to be proportional to its liquid mole fraction x<sub>i</sub>, 

P<sub>i</sub> = x<sub>i</sub> * P<sub>i</sub><sup>s</sup>

where P<sub>i</sub><sup>s</sup> is the saturated vapor pressure of the neat solvent. In the case of electrolytes, (show derivation of osmotic coefficient).

In the third part of this tutorial, we will simulate the vapor pressure depression of water with increasing concentration of LiBr, using the alkali halide model of Joung and Cheatham, 2008. We will compare the resulting osmotic coefficients with those correlated by Pitzer and Mayorga 19XX for real water. 

### references

7. Joung and Cheatham 2008
8. Pitzer and Mayorga 19XX

## part 4. mixed solvents
Most applications of electrolytes involve more than one solvent. In some cases, the electrolyte species act to modulate the thermodynamic behaviour of the mixed system. An example of this is purification of alcohols from water, where the addition of electrolyte causes phase separation of solvents which would otherwise be miscible. In other cases, the electrolytes themselves are the focus, such as in conductive electrochemical systems. Regardless of application the thermodynamics of electrolytes and solvents are intimately linked. 

### references 

9. ???