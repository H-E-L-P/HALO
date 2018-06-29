

# ![HELP LOGO](https://avatars1.githubusercontent.com/u/7880370?s=75&v=4) HELP Extended halo model

This is a software product of the Herschel Extragalactic Legacy Project ([HELP](http://www.herschel.sussex.ac.uk)).


As part of HELP, a script for constructing and constraining extended halo models in extragalactic surveys has been delivered with a demo set of data and extensive guidance for using the script. The script and demonstration data can be found in the HELP Github repository [https://github.com/H-E-L-P/HALO] ( (https://github.com/H-E-L-P/HALO). The analysis is based on that presented in [Hatfield et al., (2016)](http://adsabs.harvard.edu/abs/2016MNRAS.459.2618H), and the data used in the demo is from the VIDEO survey [(Jarvis et al., 2013)](http://adsabs.harvard.edu/abs/2013MNRAS.428.1281J), a data set included in HELP [DMU VIDEO VISTA](https://github.com/H-E-L-P/dmu_products/tree/master/dmu0/dmu0_VISTA-VIDEO).
The script is written in python 2.7, and the only non-standard module required is *Halomod* (Murray, Power, Robotham, in prep.), which can be installed using pip from github (although alternative clustering codes could be used instead with relative ease).

## Science Applications
The script provided is written in mind for a science goal of simultaneously modelling linear ($\gtrsim$Mpc) and non-linear ($\lesssim$Mpc) galaxy clustering in extragalactic galaxy surveys. It is written in mind for an analysis modelling an angular (as opposed to projected or spatial) correlation function, at redshifts $0.3 \lesssim z \lesssim 3$, on scales $0.001 \lesssim \theta /\mathrm{deg} \lesssim 1$, with a Zheng et al., (2005) parametrisation, but the code is extensively annotated to show how one would change the script to investigate different models and parts of parameter space etc. Similarly the code has a default cosmology, and halo mass function model etc., but points out clearly where this could be changed.

## Summary of Script
The script begins by importing the required python packages, and opening the files containing the clustering data being modelled. It then creates a halo model object, which contains information like the halo mass function and halo bias function etc. Then it performs a Bayesian MCMC fit over the parameters describing how galaxies occupy haloes, producing MCMC chains of samples from the posterior of the parameters.This requires a user given prior on the parameters. It then generates some derived parameters e.g. galaxy bias and satellite fraction, and generates some plots displaying the results. Finally it saves the MCMC chains and figures. The script can take between a day and a week to run on a laptop depending on precise setup.

## Inputs and Outputs

The script requires as input:
* correlation function as a 3xN array of the clustering values, the scale values and error bars (modifications for using a covariance matrix annotated in the script)
* A redshift distribution for the galaxy sample as a 2xM array of redshift values and $p(z)$
* If correcting for integral constraint is required, RR($\theta$) is required as a 2xK array ($\theta$ values and RR values, number of random pairs at that separation)
* A comoving number density of the galaxies and associated uncertainty
* Specifications of how many walkers and how many steps in the MCMC run
* The redshift and angular ranges to consider
 * A prior on the HOD parameters being constrained

The script will output:
 * MCMC chains sampling from the posterior for the HOD parameters
 * Samples from the posterior of derived parameters (in particular galaxy bias and satellite fraction)
 * Print the 16th, 50th and 84th percentiles of each HOD parameter and each derived parameter
 * A corner plot showing the posterior of the HOD parameters
 * A plot showing the correlation function data, and the posterior of the fitted correlation functions
 

-------------------------------------------------------------------------------

**Authors**:  [Seb Oliver](https://www2.physics.ox.ac.uk/contacts/people/hatfield)

 ![HELP LOGO](https://avatars1.githubusercontent.com/u/7880370?s=75&v=4)
 

The Herschel Extragalactic Legacy Project, ([HELP](http://herschel.sussex.ac.uk/)), is a [European Commission Research Executive Agency](https://ec.europa.eu/info/departments/research-executive-agency_en)
funded project under the SP1-Cooperation, Collaborative project, Small or medium-scale focused research project, FP7-SPACE-2013-1 scheme, Grant Agreement
Number 607254.

[Acknowledgements](http://herschel.sussex.ac.uk/acknowledgements)

