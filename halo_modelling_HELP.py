########### HALO MODELLING DEMO FOR THE HERSCHEL EXTRAGALACTIC LEGACY PROJECT (D10.10)
# This is a script to fit a Halo Occupation Model to galaxy clustering measurements (correlation functions) in a deep extragalactic field, in both the linear and non-linear regimes

# It is written presuming fitting a model to an angular correlation function, but could easily be modified for projected or spatial correlation functions

# The analysis in this script is based on that described in Hatfield et al., 2016 (arxiv.org/abs/1511.05476)
# This demo uses a correlation function calculated in the CHFTLS-D1 XMM-LSS 1sq deg field (specifically the 1<z<1.25, M_*>10^10.35 M_sun bin.), using VIDEO-CHFTLS data (Jarvis et al., 2013) that is included in HELP

# It uses the halomod module (github.com/steven-murray/halomod) developed by Steven Murray, who also contributed to the development of this application

# For uses very similar to the demo only the PRELIMINARIES section need be modified to use your own data etc.
# For more advanced uses you may need to edit the HALO MODEL section, and tinker with the rest of the script
# This has been written in mind for an application over z~0.3-3, for 0.001-1deg scales, with comparatively good photometric redshifts and very little cross-bin contamination etc, so some thought may be required for individual applications
# E.g. don't try to fit a z~8 correlation function or 100deg scales and expect to work without checking that the halo model parameters are appropriate!

# This is written assuming you just have error bars on your clustering measurements, but there are annotations on where you would insert a covariance matrix
# This script assumes you have already calculated your correlation function;
# you can either write your own script for this, or alternatively use one of the many packages on github e.g.
# - https://github.com/jcoupon/swot
# - https://github.com/rohinkumar/correlcalc

# The code to generate model correlation functions can take up to ~30s, meaning that the MCMC can take a day or even a week on a labtop to converge, depending on how configured and how many MCMC steps etc.

### About Halomod

# The only `non-standard' module required is halomod (Murray, Power, Robotham, in prep.)
# Most of the dependencies (e.g. python module `hmf') are automatically installed with halomod. It does require a fortran compiler.
# It also requires scipy 1.1.0 or higher

# At time of writing this is written for the `develop' version of halomod
# e.g. install using pip install --upgrade "git+https://github.com/steven-murray/halomod@develop#egg=halomod"
# (although this may change)

# This script could probably be relatively easily adapted for other halo modelling packages e.g. 
# - halomodel as used in Coupon et al., 2015 (https://github.com/jcoupon/halomodel)
# - halotools as used in Hearin et al., 2017 (https://github.com/astropy/halotools)


#### Files needed
# - Correlation function as a function of angular scales, with error bars (can be replaced with covariance matrix) [written to expect an array in log10 degrees, log10 w(theta), and s.d. on log10 w(theta) but you can customise]
# - A redshift distribution [a column of redshifts, and a column of p(z)]
# - RR (number of random pairs as a function of angular scale for the field) if you want to account for integral constraint [written to expect an array in log10 degrees and log10 RR, but you can customise]
# - A measurement of the comoving number density and an associated uncertainty [both in Mpc^-3]. If you only want to model the clustering, you can comment out the appropriate parts
# May need to slightly modify the code below for if your measurements are in radians etc.

# This script uses variables without little-h, but halomod uses units with little-h (Mpc versus Mpc/h etc.), so there are conversions throughout

#### Files produced
# - MCMC chains of sampling over the HOD parameters
# - Sampled chains of derrived parameters (satellite fraction, bias, average halo mass and comoving number density)
# - Triangle plot of HOD params
# - Plot comparing clustering measurements and fits

# Peter Hatfield, June 2018 (peter.hatfield@physics.ox.ac.uk) for HELP

print 'Starting Analysis'


########## Import modules used in analysis

# Import common standard modules used
import random
import time
import numpy as np
import pylab as pb
pb.ion()
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='sans-serif')

# Astronomical specific modules
import halomod as hm # Module that calculates model correlation functions

## Not used in script as is, but could be useful if you need to calculate comoving densities and quantities like that
#import cosmolopy.distance as cd
#import cosmolopy.parameters as cp

# Modules for MCMC/fitting
import corner
import emcee

print 'Modules Successfully Loaded'


#################### PRELIMINARIES - This section has the settings etc. that you need to change for your analysis



# Preliminaries and file locations
version=1 # Version, just so if you are fitting multiple sets of data, testing different settings out etc.
working_directory='/Users/hatfield 1/Documents/new_files/HELP_code/' # Directory where all your files are
correlation_function_file_name='demo_correlation_function.dat' # Name of the correlation function file
redshift_distribution_function_file_name='demo_redshift_distribution.dat' # Name of the redshift distribution file
account_for_integral_constraint=True # If true assumes that needs to correct for integral constraint; if set to false assumes this is already done
# If accounting for integral constraint import RR values
if account_for_integral_constraint:
    RR_file_name='demo_RR_theta_values.dat' # Name of your RR data file


print 'Version '+str(version)

#### If defining own cosmological parameters rather than using Planck...
#little_h_defined=0.7 #Reduced Hubbles constant
#omega_matter=0.3
#omega_baryon=0.045
#n_index=1.0
#sigma_8_value=0.8

# Check that its only fitting on expected range
# Combine the RR into one file


# Halo properties
type_HOD_model="Zheng05" # Name of HOD model used. Zheng05 is the standard(ish) 5-param model, halomod has a few others, but will require slight modification to the script in places

# These are various models for how you believe dark matter halos are arranged in the Universe...
hmf_function='Behroozi' # Halo mass function model (comoving number density of haloes as a function of halo mass), halomod includes a few options
bias_function='Tinker10' # Halo bias function model (bias of haloes as a function of halo mass), halomod includes a few options
concentration_mass_relation='Duffy08' # Concentration-mass relation (concentration parameter of haloes as a function of halo mass), halomod includes a few options


# Redshift and angular scale ranges
z_median=1.12 # Medians of the redshifts range
z_min_value=1 # Min to do redshift integral over
z_max_value=1.25 # Max to do redshift integral over
ang_scale_min=0.001 # Min angular scale interested in (degrees)
ang_scale_max=1.0 # Max angular scale interested in (degrees)
#N.B. if integral constraint already corrected for the ang scales interested in can probably just be the range of your measurements that you are fitting for
# If you are correcting for the integral constraint, the range probably needs to be larger to incorporate all scales that go into the integral constraint


# MCMC properties
MCMC_num_steps=500 # Number of MCMC steps, this should probably be at least 1000...
ndim = 5 # Number of parameters (normally number of params in the hod model, but concievably you might want to leave some fixed etc.)
nwalkers=10 # Number of MCMC Walkers, should be greater than twice number of dimensions
sample_rate=10 # How much of the full sample to use in calculating derrived parameters (e.g. a value of 10 will calculate derrived parameters for each tenth sample), this must be a positive integer factor of ndim*nwalkers

print 'Number of MCMC steps = ' + str(MCMC_num_steps)


### Prior ranges on your HOD parameters
# Halo masses in log10 solar masses
# These are Zheng05 parametrisation
# Mmin is the halo mass required to have a central galaxy
M_min_min=11.2 # Min M_min value for prior
M_min_max=12.8 # Max M_min value for prior
# M1 is the halo mass required for the first satellite
M_1_max=14.1 # Max M_1 value for prior (M_1 assumed to be greater than M_min)
# Alpha is the power law index for number of satellites as a function of halo mass
alpha_max=1.35 # Max alpha value for prior
alpha_min=0.7  #0.2 # Min alpha value for prior
# Sigma is how smooth the cut-off for having a central or not is (related to central galaxy to halo mass ratio scatter)
sigma_min=0.25 # Min sigma value for prior
sigma_max=0.6 # Max sigma value for prior
# M0 is a halo mass below which satellites are forbidden
M_0_min=9.9 # Min M_0 value for prior





## Import data
# Load acf
data = np.genfromtxt(working_directory+correlation_function_file_name) # Load the measurements - take care to do conversion between log and linear, degrees and radians etc.
x = data[:,0]    #  log10 theta in degrees
y = data[:,1]    # log10 w(theta)
yerr_log=data[:,2] # error on log10 w(theta)
#yerr_log=np.maximum(yerr_log,0.05) # Get rid of unrealistic small errors if you want to account for cosmic variance on clustering measurements or something perhaps
yerr=(10**(y+yerr_log))-(10**(y)) # Move out of log space
x=10**x # Move out of log space
y=10**y # Move out of log space

# Comoving number density (change this to your own values)
n_data_gal_den=0.0017316 # This is just calculated as `number of sources' divided by `comoving volume for survey area and redshift range', although one might want to do something more complicated if there is incompleteness etc.
cv_n_error=0.00015305 # This is both Poisson uncertainty, and cosmic variance on the counts. Cosmic variance on the counts here calculated based on Trenti et al., , see also online calculator: http://casa.colorado.edu/~trenti/CosmicVariance.html

# Technically comoving number density is also a cosmology dependent quantity, but typically the cosmic variance will be greater than the uncertainty on number density from cosmology in most contexts


#### If your correlation function is integral constraint corrected, this is not nescessary
#Probably technically marginally more rigorous to have integral constraint as part of model, but I don't think makes a huge different to having a fixed integral constraint from start
# Load RR values
RR_data= np.genfromtxt(working_directory+RR_file_name)
RR=10**RR_data[:,1] # Convert RR out of log10
theta_RR= 10**RR_data[:,0] # Convert RR angular scales out of log10



# Read in N(z) and create function
nz = np.loadtxt(working_directory+redshift_distribution_function_file_name)
redshift_distribution=spline(nz[:,0],nz[:,1]) # The code should correct normalisation for sensible values
# It may warn you that your redshift distribution is not normalised, but it should be able to cope with it


print 'Data successfully loaded and analysis parameters selected'


#################### HALO MODEL - This part can probably be largely left unchanged, but some values might need to be changed to make code run quicker etc.


# Define a halo model
h = hm.AngularCF(z=z_median,hod_model="Zheng05") # Use hm.HaloModel if not interested in angular correlation functions, or hm.ProjectedCF for projected correlation functions
h.update(hod_params={"central":True}) # Normally want this true; specifies that a halo must have a central galaxy before it can have satellites
h.update(hmf_model=hmf_function) # Sets your Halo Mass Function prescription
h.update(bias_model=bias_function) # Sets your halo bias prescription
h.update(concentration_model=concentration_mass_relation) # Set concentration-mass relation
h.update(zmin=z_min_value) # Set range of redshift integral 
h.update(zmax=z_max_value) # Set range of redshift integral 
h.update(znum=50) # Set number of bins for redshift integral 
h.update(p1=redshift_distribution) # Set redshift distribution
h.update(theta_min=ang_scale_min*(np.pi/180.0)) # Set angular scales considered (converted into radians)
h.update(theta_max=ang_scale_max*(np.pi/180.0))# Set angular scales considered (converted into radians)
h.update(theta_num=60) # Number of points in angular space calculated for
h.update(rmin=1e-3) # r range for projection integral # This one can probably be done automatically (units Mpc/h for halomod)
h.update(rmax=320.0) # r range for projection integral (units Mpc/h for halomod)
h.update(rnum=150) # Number of r bins for integral

# Note about the rmin, rmax, rnum etc. These numbers define what ranges and step sizes halomod does the projection interval for
# Typically a higher rnum will be more accurate, but take longer to calculate
# Similarly for a larger range (higher rmax etc.)
# The values here are sensible for the redshift range considered in this example, but you may have to tinker if you want to consider really small scales, or much lower or higher redshifts etc.


h.update(cosmo_model="Planck15") # Sets your cosmology


## If defining own cosmological parameters....
#h.update(sigma_8=sigma_8_value)
#h.update(cosmo_params={"H0":little_h_defined*100, "Om0":omega_matter,"Ob0":omega_baryon})
#h.update(n=n_index)
## You can in principle vary more complicated things like number of neutrino species and so on within the module but this probably requires special care


# Extract little_h value
little_h=h.cosmo.H0.value/100

print 'Halo model sucessfully created'


#################### Do the fitting process

# Define log likelihood, prior and probabilities for model fitting

# Log likelihood
def lnlike(theta, x, y, yerr,n_data_gal_den,cv_n_error):
    
    M_min,M_1,alpha,sig_logm,M_0  = theta # Extract HOD parameters
    
    # Update halomodel
    h.update(hod_params={"M_min":M_min+np.log10(little_h), "M_1":M_1+np.log10(little_h), "alpha":alpha, "sig_logm":sig_logm, "M_0":M_0+np.log10(little_h)}) # Update HOD model with new parameters
    # (The little_h's are because halomod is expecting things in M/h units)
    
    #Create model acf at the data points
    theta_in_degrees=h.theta*(180.0/np.pi)
    # Model acf at halomod theta values
    w=h.angular_corr_gal
    # Spline over the relevant range
    s = spline(theta_in_degrees,w,k=3)
    
    if account_for_integral_constraint:
        # Calculate an estimate of the integral constraint
        C=np.sum(RR*s(theta_RR))/np.sum(RR)
        
    else:
        C=0 # If already corrected for integral constraint just set C=0
        
    # Find model w(theta) over same theta values as data, and take off integral constraint
    model = s(x)-C
    
    # Model comoving number density
    model_num_density=h.mean_gal_den*(little_h**3) # Making little_h correction

    # Calculate the log-liklihood
    inv_sigma2 = 1.0/(yerr**2)
    log_liklihood_clustering=-0.5*np.sum(((y-model)**2)*inv_sigma2) # Log liklihood of clustering
    log_liklihood_counts=-0.5*(((model_num_density-n_data_gal_den)/cv_n_error)**2) # Log liklihood of comoving number density. Set this to 0 if you don't want to fit to number count data
    total_log_liklihood=log_liklihood_clustering+log_liklihood_counts # Add the number counts agreement
    return total_log_liklihood # Log likelihood
    
    # If convariance matrix instead of just error bars can use something like:
    #y_model_difference=np.matrix(y-model)
    #y_model_difference_trans=np.transpose(y_model_difference)
    #log_liklihood_clustering=-0.5*(y_model.dot(cov_matrix_inverse.dot(y_model_trans)))

# Log prior
def lnprior(theta):
    

    M_min,M_1,alpha,sig_logm,M_0 = theta # Extract HOD parameters
    
    # Flat priors over range specified earlier
    if M_min_min < M_min < M_min_max and M_min < M_1 < M_1_max and alpha_min < alpha < alpha_max and sigma_min < sig_logm < sigma_max  and M_0_min < M_0 < M_1: 
        return 0.0 # If in prior range, ignoring normalisation
    return -np.inf # If outside prior range ("log zero = minus infinity")

# Log posterior
def lnprob(theta, x, y, yerr,n_data_gal_den,cv_n_error):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr,n_data_gal_den,cv_n_error) # Log posterior proportional to Log prior + Log liklihood





# Set up initial positions of walkers as random samples from the prior
pos=[] # Empty list
for i in range(nwalkers):
    rand1=random.uniform(M_min_min, M_min_max)
    rand2=random.uniform(rand1, M_1_max)
    rand3=random.uniform(alpha_min, alpha_max)
    rand4=random.uniform(sigma_min, sigma_max)
    rand5=random.uniform(M_min_min, rand2)
    pos.append([rand1,rand2,rand3,rand4,rand5])


print 'Starting MCMC'

# Define MCMC sampler from emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y,yerr,n_data_gal_den,cv_n_error)) 

t1 = time.time() # Time the fitting process
# Run MCMC process (this is the long bit)
sampler.run_mcmc(pos,MCMC_num_steps) 
t2 = time.time()
print 'MCMC Finished'
print 'MCMC took '+str(np.floor((t2-t1)/60))+' minutes' # Prints how long the MCMC fitting took

samples = sampler.chain.reshape((-1, ndim)) # Combine chains into one array (you may wish to remove some burn-in e.g. delete first third of chains or something)
lnprob = sampler.lnprobability.reshape(-1) # Probability values for each point in chain, don't use in this script but may be useful

#################### Now plot stuff, save data etc...

## Triangle plot of the samples
corner.corner(samples, labels=[r'$M_\mathrm{min}$', r'$M_1$', r'$\alpha$',r'$\sigma$',r'$M_0$'],quantiles=[0.16, 0.5, 0.84]) # Show 16th, 50th and 84th percentiles

#Save triangle plot
pb.savefig(working_directory+"hod_fit_triangle_plot_v"+str(version)+".png")


## Calculate derrived parameters



# Set up empty vectors of model correlation functions and derrived parameters
model_omega_distrib=np.zeros(shape=(len(x),nwalkers*MCMC_num_steps/sample_rate)) # Calculate acf for each parameter sampling
bias_effective_distrib=np.zeros(nwalkers*MCMC_num_steps/sample_rate) # Bias for distribution
mass_effective_distrib=np.zeros(nwalkers*MCMC_num_steps/sample_rate) # Effective halo mass for distribution
satellite_fraction_distrib=np.zeros(nwalkers*MCMC_num_steps/sample_rate) # Satellite fraction for distribution
number_density_distrib=np.zeros(nwalkers*MCMC_num_steps/sample_rate) # Galaxy comoving density for distribution

print 'Starting Calculating Derived Parameters'

for i in range(nwalkers*MCMC_num_steps/sample_rate): # Go through some fraction of the samples and calculate derrived parameters
    
    # Update halo model
    h.update(hod_params={"M_min":samples[i*sample_rate,0]+np.log10(little_h), "M_1":samples[i*sample_rate,1]+np.log10(little_h), "alpha":samples[i*sample_rate,2], "sig_logm":samples[i*sample_rate,3], "M_0":samples[i*sample_rate,4]+np.log10(little_h)}) # Update HOD model with new parameters
    
    #Create model acf at the data points
    theta_in_degrees=h.theta*(180.0/np.pi)
    # Model acf at halomod theta values
    w=h.angular_corr_gal
    # Spline over the relevant range
    s = spline(theta_in_degrees,w,k=3)
    
    if account_for_integral_constraint:
        # Calculate an estimate of the integral constraint
        C=np.sum(RR*s(theta_RR))/np.sum(RR)
    else:
        C=0 # If already integral constraint corrected just set C=0
        
    # Find model w(theta) over same theta values as data, and take off integral constraint
    model = s(x)-C
    
    # Save all the derrived parameters
    model_omega_distrib[:,i]=model
    bias_effective_distrib[i]=h.bias_effective # The galaxy bias of the sample
    satellite_fraction_distrib[i]=h.satellite_fraction # The fraction of the galaxies that are satellite galaxies
    mass_effective_distrib[i]=h.mass_effective-np.log10(little_h) #This turns it from M/h to just M
    number_density_distrib[i]=h.mean_gal_den*(little_h**3) #This turns it from h^3/Mpc^3 to just Mpc^-3


print 'Finished Calculating Derived Parameters'


# Print the HOD parameters
print 'M_min$'
print np.percentile(samples[i*sample_rate,0],16), np.percentile(samples[i*sample_rate,0],50), np.percentile(samples[i*sample_rate,0],84)
print 'M_1'
print np.percentile(samples[i*sample_rate,1],16), np.percentile(samples[i*sample_rate,1],50), np.percentile(samples[i*sample_rate,1],84)
print 'alpha'
print np.percentile(samples[i*sample_rate,2],16), np.percentile(samples[i*sample_rate,2],50), np.percentile(samples[i*sample_rate,2],84)
print 'sigma'
print np.percentile(samples[i*sample_rate,3],16), np.percentile(samples[i*sample_rate,3],50), np.percentile(samples[i*sample_rate,3],84)
print 'M_0'
print np.percentile(samples[i*sample_rate,4],16), np.percentile(samples[i*sample_rate,4],50), np.percentile(samples[i*sample_rate,4],84)

# Print the derrived parameters
print 'Satellite Fraction'
print np.percentile(satellite_fraction_distrib[:],16), np.percentile(satellite_fraction_distrib[:],50), np.percentile(satellite_fraction_distrib[:],84)
print 'Bias'
print np.percentile(bias_effective_distrib[:],16), np.percentile(bias_effective_distrib[:],50), np.percentile(bias_effective_distrib[:],84)
print 'Effective Halo Mass'
print np.percentile(mass_effective_distrib[:],16), np.percentile(mass_effective_distrib[:],50), np.percentile(mass_effective_distrib[:],84)
print 'Derrived Number Density'
print np.percentile(number_density_distrib[:],16), np.percentile(number_density_distrib[:],50), np.percentile(number_density_distrib[:],84)
print 'Observed Number Density'
print n_data_gal_den


# Create objects for holding the models
model_best=np.zeros(len(x))
model_lower=np.zeros(len(x))
model_upper=np.zeros(len(x))

# Find percentiles of acf
for i in range(len(x)):
    model_best[i]=np.percentile(model_omega_distrib[i,:],50)
    model_lower[i]=np.percentile(model_omega_distrib[i,:],16)
    model_upper[i]=np.percentile(model_omega_distrib[i,:],84)

# Plot data
plt.figure()
plt.errorbar(x,y,yerr=yerr,linestyle="None",fmt='o',color=(1,0,0), label="Data",elinewidth=2,markersize=6)


# Plot model fits
x_long=10**np.arange(np.log10(ang_scale_min),np.log10(ang_scale_max),0.01)

## Calculate integral constraint for models
# Best
g = spline(x,model_best,k=3)
pb.plot(x_long,g(x_long),'k',linewidth=2.0,label="Model") # Plot HOD model c.f.
# Lower
g = spline(x,model_lower,k=3)
pb.plot(x_long,g(x_long),'k',linewidth=1.0) # Plot HOD model c.f.
# Upper
g = spline(x,model_upper,k=3)
pb.plot(x_long,g(x_long),'k',linewidth=1.0) # Plot HOD model c.f.

# Plot details
pb.xscale("log") # Log scale
pb.yscale("log")
pb.xlim(0.001,0.1) # x-axis limits
pb.ylim(0.01,10**(0.3)) # y-axis limits
pb.legend(loc='lower left',frameon=False, fontsize=11)
pb.xlabel(r'$\theta \mathrm{ (degrees)}$')
pb.ylabel(r'$w(\theta)$')





# Save plot
pb.savefig(working_directory+"demo_model_acf_plot_v"+str(version)+".png")

#### Save MCMC chains etc.

np.savetxt(working_directory+"hod_fit_MCMC_chains_v"+str(version)+".dat",samples) # Save HOD param samples

derrived_parameters=np.transpose([satellite_fraction_distrib,bias_effective_distrib,mass_effective_distrib,number_density_distrib])
np.savetxt(working_directory+"hod_derrived_parameter_samples_v"+str(version)+".dat",derrived_parameters) # Save derrived params

model_acf=np.transpose([model_lower,model_best,model_upper])
np.savetxt(working_directory+"model_acf_percentiles_v"+str(version)+".dat",model_acf) # Save model acfs



print 'MCMC Chains etc. Saved'
