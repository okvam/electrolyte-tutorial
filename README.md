# electrolyte-tutorial
Tutorial for chemical theory and simulation of electrolyte solutions.

## 0. introduction
Electrolytes are compounds which upon contact with a solvent dissociate (_lysis_) to form electrically charged species (_electro_). The most familiar electrolyte is perhaps the table salt found in your kitchen, but electrolytes are of fundamental important for physical processes ranging from regulation of biological systems to the battery technology powering our electronic devices. The thermodynamics of electrolyte solutions have a long history of scientific investigation. 

## 1. water
The best studied solvent for electrolyte solutions is water, partly because of its prevalence and partly because of the high solubility of most simple electrolyte species in this solvent. Water alone is a highly complex fluid, showing many unusual thermophysical properties. T

In the first point of this tutorial we will simulated some properties of water, using the SPC/E water model. We will answer the following questions;

1. How does the liquid density of SPC/E water compare to real water at 298.15K and 1 bar?
2. Does SPC/E water have a temperature of maximum liquid density? 
3. What is the saturated vapor pressure of the SPC/E model of water at 298.15K?

## 2. methanol
Methanol is the smallest alcohol, replacing one of the H-atoms from water with an organic CH3 group. Many simple electrolytes are reasonably soluble also in methanol, albeit less so than in water. 

In the second point of this tutorial we will simulate some properties of methanol and water + methanol mixtures, using the TraPPE model of methanol. We will answer the following questions:

1. How does the liquid density of TraPPE methanol compare to real methanol at 298.15K and 1 bar?
2. What is the saturated vapor pressure of the TraPPE methanol at 298.15K?
3. Is the simulated isobaric vapor-liquid composition diagram in good agreement with experiment for water + methanol at 1 bar? 

## 3. vapor pressure depression
Electrolyte species are for all intents and purposes non-volatile. From Raoult's law we expect the partial vapor pressure P<sub>i</sub> of a component in a liquid mixture to be proportional to its liquid mole fraction x<sub>i</sub>, 

	P<sub>i</sub> = x<sub>i</sub> * P<sub>i</sub><sup>s</sup>

where P<sub>i</sub><sup>s</sup> is the saturated vapor pressure of the neat solvent. In the case of electrolytes, (show derivation of osmotic coefficient).

In the third point of this tutorial, we will simulate the vapor pressure depression of water with increasing concentration of LiBr, using the alkali halide model of Joung and Cheatham, 2008. We will compare the resulting osmotic coefficients with those correlated by Pitzer and Mayorga 19XX for real water. 

## 4. mixed solvents
Most applications of electrolytes involve more than one solvent. In some cases, the electrolyte species act to modulate the thermodynamic behaviour of the mixed system. An example of this is purification of alcohols from water, where the addition of electrolyte causes phase separation of solvents which would otherwise be miscible. In other cases, the electrolytes themselves are the focus. Regardless of application the thermodynamics of electrolytes and solvents are intimately linked. 

