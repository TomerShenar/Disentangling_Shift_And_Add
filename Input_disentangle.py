# shift-and-add & grid disentangling, by Tomer Shenar, with contributions from Matthias Fabry & Julia Bodensteiner
# 21.11.2022, V1.0; feel free to contact at T.Shenar@uva.nl or tomer.shenar@gmail.com for questions/inquires
# Input file

import numpy as np


###############################
##### STAR AND DATA INFO ######
###############################


# Name of object (just for file names)
StarName = 'Test'

# Path to data (folder where all spectra are stored):
ObsPath = 'obs/'

### Type of data format. There are two options:
### OPTION 1: ObsFormat = 'TXT' 
### assumes that the observations are in ascii format, each file containing 2-column tables of wave & normalised flux. 
### In addition, the observation directory MUST contain a file called 'ObsDat.txt', which has the following format:
###   MJD          obsname
###   xxxx          NAME1
###   yyy           NAME2
###   ...           ...
### The paths should be either absolute or relative to the directory in which the script is stored.
### OPTION 2: ObsFormat = 'FITS' 
### The script will look for ALL fits files in the given directory. 
### The script will attempt to retrieve the dates from the fits headers using a user-specified header keyword
### IMPORTANT NOTES:
### 1. It doesn't matter if the dates are "MJD", "JD", "HJD", etc -- important is that the T0 provided by the user matches this!
### 2. For the "fits" option, I include a few built-in functions to read e.g. HERMES, X-SHOOTER, FEROS spectra.... 
### User should feel free to update the reading of the files!

ObsFormat = 'TXT'

#Only important if ObsFormat='FITS'
MJDHeader = 'MJD-OBS'



###############################
##### SYSTEM PROPERTIES #######
###############################


#Number of components; currently possible only 2 components
# up to four components in upcoming version
CompNum = 2


# Orbital parameters 
### P, T0, ecc, omega, and gamma cannot be derived with current version and are assumed by the user
### K1, K2 can be explored via chi2 if required by the user, but initial guesses should be given.
### Important: omega is defined via vr(1) = Gamma1 + K1*(cos(nu) + ecc * cos(omega) )
### If results don't make sense, your best bet is to set omega --> omega + pi
Orbital_Params = {
####### MUST BE FILLED BELOW ALWAYS (inner binary) ######   
    'Period': 473.,
    'T0': 0.,
    'ecc': 0.5,
    'omega': 60.,
    'Gamma': 0.,
    'K1': 87.,
    'K2': 135.  ,
####### Only for triples / quadruples: refers to outer orbit of companion/binary AROUND the inner binary ######       
####### IF outer period very long/negligible, set PeriodOut = 1E10 to neglect long-term motion     
    'PeriodOut': 0.,
    'T0Out' : 0.,
    'eccOut' : 0.,
    'omegaOut' : 0.,
    'KOut' : 0.,
####### Only for Quadruples: refers to orbit of 2nd binary in the system  ######           
    'Period_2' : 0.,
    'T0_2' : 0.,
    'ecc_2' : 0.,
    'omega_2' : 0.,  
    'K3' : 0. ,
    'K4' : 0.}
    

# Vector of light ratios, [l2, l3, l4], i.e. flux_i / sum(flux). Assumed constant throughout range.
lguessVec = [0.3, 0., 0.] 

    
lguess1 = 1. - np.sum(lguessVec)



#Where to measure S2N, only important for  defining continuum and weighting of spectra when co-adding (not critical)
S2Nblue = 4155
S2Nred = 4165



################################
##### Disentangling options ####
################################


# Run grid disentangling? 
# If TRUE: will conduct grid disentangling and derive Ks
# If FALSE: will only peform separation using input K1,K2 
GridDis = False

# Define grid search (only important if GridDis = True). 
# For setting K1, K2, K3, K4 search arrays: Karr = np.arange(IniFacK*K, FinFacK*K, Dense)
# Current version only works for first two columns (K1, K2)
# If DenseKArr[i] = 1, then the search is "1D", i.e. K is fixed to the value specified by the user.
DenseKArr = [15, 15, 1, 1]
IniFacKArr = [0.1, 0.1, .3, .3]
FinFacKArr = [2., 2., 2., 2.]



# Number of iterations
### IMPORTANT NOTES:
### 1. Ideally convergence could be determined via a condition on EPS (See above). However, a suitable condition could not yet be developed
### --> User needs to judge when results are "sufficiently converged", by either comparing the results for different itr numbers, or using
###     options below.
### 2. itrnumlim is the number of iterations per K1,K2 pair; NumItrFinal is the number of iterations for the final separation, 
###    after K1, K2 have been derived / set. Often, itrnumlim < NumItrFinal, for speed, and since individual iterations occur on individual lines.
### 3. See documentation for tips and insights about number of iterations.

itrnumlim =50
NumItrFinal = 1000


# If StrictNegA = True, enforce disentangled spectra to be below continuum except for prespecified regions (given in array).
# Below continuum = ForceNegSigma "sigmas" below continuum. 
# HIGHLY RECOMMENDED for OB-type stars -- otherwise, output often exhibits cosmetic "wings" and continuum offsets.
# For WR stars typically "False" is better.
# "Positive regions" should be regions with expected emission lines etc. For WR stars, 

ForceNegSigma = 2.

StrictNegA = True

#Only relevant if StrictNegA=True
PosLimCondA = np.array([ 
            [3968., 3969.]    
              ])

# Same as StrictNegA for secondary
StrictNegB = True

PosLimCondB = np.array([ 
            [3968., 3969.]
              ])

# Same as StrictNegA for tertiary
StrictNegC = True

PosLimCondC = np.array([ 
            [3968., 3969.]
              ])

# Same as StrictNegA for fourth companion
StrictNegD = True

PosLimCondD = np.array([ 
            [3968., 3969.]
              ])



# Define regions where the solution is allowed to be above continuum (where emission is expected)
PosLimCondNeb = np.array([ 
            [3968., 3971.],
            [4025., 4027.],            
            [4100.5, 4103],
            [4143., 4145],            
            [4339., 4342],
            #[4340., 4345],            
            [4387., 4391.5],       
            #[4385., 4393],             
            #[4335., 4345.],            
            #[4465., 4482.]       
            [4470., 4473.]              
            #[4135., 4150.],            
            #[4228., 4238.],
            #[4330., 4355.], 
              #[4381., 4396.],              
              #[4465., 4485.],  
              #[4840., 4880.], 
              #[4916., 4930.],  [5010., 5024.],
              #[5160., 5180.], [5190., 5210.], [5225., 5240.], [5270., 5290.],
              #[5310, 5325.], [5255., 5370.], [5520., 5540.], [6140., 6160.], [6230., 6260.], [6310., 6330], 
              #[6340., 6390.], [6410., 6460.], [6510, 6520.], [6553., 6574.], [7505, 7525], [7766., 7786.]
              ])



# Plot fits between disentangled spectra, their sum, and the observations at RV extremes.
# Highly recommended for sanity checks and presentation in papers.
# The plot is shown for the K1, K2 pair most closely matching (Velo_plot_usrK1_ext, Velo_plot_usrK2_ext, ...) given by the user.
# Recommended: True
PLOTEXTREMES = True
Velo_plot_usrK1_ext =  Orbital_Params['K1']
Velo_plot_usrK2_ext=  Orbital_Params['K2']
Velo_plot_usrK3_ext=  Orbital_Params['K3']
Velo_plot_usrK4_ext=  Orbital_Params['K4']

# line width and figsize for "Extreme plots"
linewidExt = 7
ExtremesFigSize = (7, 7)




# Plot fits between disentangled spectra, their sum, and each epoch of observation. 
# Useful to examine all data; rather tedious but important for critical systems (e.g., black holes!)
# The plot is shown for the K1, K2 pair most closely matching (Velo_plot_usrK1, Velo_plot_usrK2) given by the user.
# Recommended: False
PLOTFITS = False
Velo_plot_usrK1 = Orbital_Params['K1']
Velo_plot_usrK2 = Orbital_Params['K2']
Velo_plot_usrK3 = Orbital_Params['K3']
Velo_plot_usrK4 = Orbital_Params['K4']

# Plot convergence plot
# If True, will produce converge plot, i.e. EPS vs. itr for each run.
# EPS = SUM(DisSpec[i+1] - DisSpec[i]), where the maximum on all components is taken. 
# Recommended: False
PLOTCONV = False


# Plot disentangled spectra after each "N_Iteration_Plot" iterations; helpful for judging convergence.
# Recommended: False
PLOTITR = False
N_Iteration_Plot = 100


# Type of interpolation in interp1d (see python doc for options);
# 'linear' can lead to artificial increase of S/N due to interpolation
# 'cubic' performs better, but is slower. 
InterKind='linear'         

# Region for fitting parabola of chi2 in index steps from minimum
ParbSize=2



################################
##### Disentangling lines ######
################################

# User chooses in which line/lines the K1,K2 search should occur.
# All that is required is:
# 1. Ranges = [ [l1, l2], [l3, l4], ...]
# 2. Rangestr = 'xxx' -- used below to pick range, but either way, needs to be specified for file-saving purposes.
# For convenience, typical lines (for massive stars) are provided below. 
# USERS: feel free to edit wavelength regions below!!!
# IMPORTANT: the final ranges used/plotted are NOT identical to those provided by the user: the script reduces them to ensure that edge issues are avoided.
# The reduction depends on K1, K2; the user should judge (e.g., using "PLOTEXTREMES") that the lines are well covered and reach continuum at both edges.
# Ideally, disentangled region should be line-dominated (to enhance the signal on chi2), but certainly reach continuum at the edges.


#Rangestr = 'Hdelta'
#Rangestr = 'Hgamma'
#Rangestr = 'Hbeta'
#Rangestr = 'Halpha'
#Rangestr = 'Balmer'
#Rangestr = 'Balmer_noHalpha'
#Rangestr = 'HeI'
Rangestr = 'HeI4472'
#Rangestr = 'HeI4122'
#Rangestr = 'HeI4009'
#Rangestr = 'HeI4026'
#Rangestr = 'HeI4144'
#Rangestr = 'HeII4200'
#Rangestr = 'HeI4388'
#Rangestr = 'HeII4546'
#Rangestr = 'HeI5878'
#Rangestr = '4120Region'
#Rangestr = '4020Region'
#Rangestr = 'IronEmission'
#Rangestr = 'OIII'
#Rangestr = 'OIII8446'
#Rangestr = 'HI8367'
#Rangestr = 'Fe4584'
#Rangestr = 'Fe5168'     
#Rangestr = 'Fe5192'   
#Rangestr = 'Fe5234'   
#Rangestr = 'Fe5275'  
#Rangestr = 'Fe5316'
##Rangestr = 'Fe5362'
#Rangestr = 'AllHeI'
#Rangestr = 'AllHeII'
#Rangestr = 'AllHe'
#Rangestr = 'Indiv'

##### Define ranges corresponding too the strings above.... CHANGE IF NEEDED

RangeHa = [6553, 6570.]
RangeHb = [4840, 4877.]
RangeHg = [4310., 4370.]
RangeHd = [4070, 4140.]
RangeHeI5878 = [5869., 5881.]
RangeHeI4472 = [4457., 4489.]
RangeHeI4144 = [4120., 4170.]
RangeOIIIJulia = [7760., 7785.]
RangeFe4584 = [4580, 4588]
RangeFe4584 = [4580, 4586]
RangeFe5168 = [5162, 5174]
RangeFe5192 = [5190., 5205]
RangeFe5234 = [5230., 5240]
RangeFe5275 = [5268, 5282.]
RangeFe5316 = [5310., 5322.]
RangeFe5362 = [5358., 5367.]
RangeOIII8446 = [8438., 8455.]
RangeHI8367 = [8367., 8500.]
RangeHeI4122 = [4115., 4127.]
RangeHeI4009 = [4003., 4018.]
RangeHeI4026 = [4000., 4050.]
RangeHeI4388 = [4365, 4410.]
RangeHeII4545 = [4515., 4565.]
RangeHeII4200 = [4185., 4215.]

# Define "Ranges" list based on user's choices from above.

if Rangestr == 'Hgamma':
    Ranges = [RangeHg]
elif Rangestr == 'Hbeta':
    Ranges = [RangeHb]    
elif Rangestr == 'Halpha':
    Ranges = [RangeHa]   
elif Rangestr == 'Balmer':
    Ranges = [RangeHd, RangeHg, RangeHb, RangeHa]       
elif Rangestr == 'Balmer_noHalpha':
        Ranges = [RangeHd, RangeHg, RangeHb]    
elif Rangestr == 'HeI':
    Ranges = [RangeHeI4009, RangeHeI4026, RangeHeI4122, RangeHeI4144, RangeHeI4388, RangeHeI4472]      
elif Rangestr == 'Hdelta':
        Ranges = [RangeHd]   
elif Rangestr == 'HeI4472':
    Ranges = [RangeHeI4472]    
elif Rangestr == 'HeI5878':
    Ranges = [RangeHeI5878]    
elif Rangestr == 'HeI4144':
    Ranges = [RangeHeI4144]   
elif Rangestr == 'HeII4200':
    Ranges = [RangeHeII4200]          
elif Rangestr == '4120Region':
    Ranges = [RangeHeI4144]       
elif Rangestr == 'IronEmission':   
    Ranges = [RangeFe5168, RangeFe5192, RangeFe5275, RangeFe5316]           
elif Rangestr == 'Fe4584':     
    Ranges = [RangeFe4584] 
elif Rangestr == 'Fe5168':     
    Ranges = [RangeFe5168] 
elif Rangestr == 'Fe5192':     
    Ranges = [RangeFe5192] 
elif Rangestr == 'Fe5234':     
    Ranges = [RangeFe5234]     
elif Rangestr == 'Fe5275':     
    Ranges = [RangeFe5275]       
elif Rangestr == 'Fe5316':     
    Ranges = [RangeFe5316]     
elif Rangestr == 'Fe5362':     
    Ranges = [RangeFe5362]     
elif Rangestr == 'OIII':
    Ranges = [RangeOIIIJulia]        
elif Rangestr == 'OIII8446':
    Ranges = [RangeOIII8446]    
elif Rangestr == 'HI8367':
    Ranges = [RangeHI8367]        
elif Rangestr == 'HeI4122':
    Ranges = [RangeHeI4122]   
elif Rangestr == 'HeI4009':
    Ranges = [RangeHeI4009]       
elif Rangestr == 'HeI4026':
    Ranges = [RangeHeI4026]        
elif Rangestr == 'HeI4388':
    Ranges = [RangeHeI4388]            
elif Rangestr == 'HeII4546':
    Ranges = [RangeHeII4545]  
elif Rangestr == 'AllHe':
    Ranges = [RangeHeI4026, RangeHeII4200, RangeHeI4472, RangeHeII4545]             
elif Rangestr == 'AllHeI':
    #Ranges = [RangeHeI4026, RangeHeI4388, RangeHeI4472]     
    Ranges = [RangeHeI4026, RangeHeI4144, RangeHeI4388, RangeHeI4472, RangeHeII4545]         
elif Rangestr == 'AllHeII':
    Ranges = [RangeHeII4200, RangeHeII4545]          
elif Rangestr == 'Indiv':
    Ranges = [RangeHeII4200, RangeHeI4472]            



                                
         
################################
##### Fancy options   ##########
################################     


#Clean cosmics?
CleanCos = False


#Renormalise spectra at pre-specified points. 
Renormalise = False
NormPoints = [3961., 4006., 4016., 4038., 4088., 4116., 4129., 4138., 4154., 4195., 4210., 4328., 4362., 4386., 4400., 4462., 4490., 4494., 4530., 4557., 4560]

# Nebular line handling?
NebLines = False


####################
