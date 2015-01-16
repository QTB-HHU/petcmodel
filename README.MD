# Model of the photosynthetic electron transport chain
## Abstract
Photosynthetic eukaryotes house two photosystems with distinct light absorption spectra. Natural fluctuations in light quality and quantity can lead to unbalanced or excess excitation, compromising photosynthetic efficiency and causing photodamage. Consequently, these organisms have acquired several distinct adaptive mechanisms, collectively referred to as non-photochemical quenching (NPQ) of chlorophyll fluorescence, which modulates the organization and function of the photosynthetic apparatus. The ability to monitor NPQ processes fluorometrically has led to substantial progress in elucidating the underlying molecular mechanisms. However, the relative contribution of distinct NPQ mechanisms to variable light conditions in different photosynthetic eukaryotes remains unclear. Here, we present a mathematical model of the dynamic regulation of eukaryotic photosynthesis using ordinary differential equations. We demonstrate that, for Chlamydomonas, our model recapitulates the basic fluorescence features of short-term light acclimation known as state transitions and discuss how the model can be iteratively refined by comparison with physiological experiments to further our understanding of light acclimation in different species.

## Publication
Ebenhoeh, O., Fucile, G., Finazzi, G., Rochaix, J-D & Goldschmidt-Clermont, M. (2014). 'Short-term acclimation of the photosynthetic electron transfer chain to changing light: a mathematical model'. Philosophical Transactions of the Royal Society B: Biological Sciences, vol 369, no. 1640, 20130223. 

Online: [DOI: 10.1098/rstb.2013.0223](http://dx.doi.org/10.1098/rstb.2013.0223)

## Model
We share Python files reproducing the mathematical model described in the paper published by Ebenhoeh et al. 2014, originally developed in MATLAB The model comprises a set of seven coupled ordinary differential equations and captures the temporal evolution of major protein complexes and key products of photosynthetic light reaction for model organism alga Chlamydomonas reinhardtii . The equations are implemented in petcModel.py and their origin is explained in detail in the electronic supplementary material of the publication. 
* lightProtocol.py: includes several light functions necessary for simulations
* parametersPETC.py: stores all parameters
* misc.py: contains two additional methods to calculate lumenal pH
* simulate.py: includes integration methods and methods to calculate steady-state, also based on absolute error tolerance 
* petcResults: additional file helping to reproduce some of the figures included in the publication

## Experimental protocols
* anoxia - simulates fluorescence trace in anoxia induced state transitions
* dld - simulates fluorescence trace in Dark-Light-Dark experiment 
* steadyStateAnalysis - calculates steady state for different combination of light intensity and phosphorylation. Investigates systems behaviour when state transitions switched off

## Dependencies
All Python packages that are needed in order to run the model are listed in the file requirements.txt and can be installed by running 
	pip install -r requirements.txt
It was noticed that under some systems LSODA integrator from the scipy.odeint package is not available. This integrator is necessary for a correct computation. Please check if it is available on your system before running the model.

## Example
For internal debbuging information we produce a number of output on the command line. In order to avoid them, we suggest to suppress them by writing all log information to /dev/null. 

To obtain similar graph as Figure 2 top run:
	python dld.py >/dev/null 2>/dev/null
and to obtain similar graph as Figure 2 bottom run:
	python anoxia.py >/dev/null 2>/dev/null
