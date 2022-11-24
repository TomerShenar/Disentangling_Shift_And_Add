# shift-and-add & grid disentangling, by Tomer Shenar, with contributions from Matthias Fabry & Julia Bodensteiner
# 21.11.2022, V1.0; feel free to contact at T.Shenar@uva.nl or tomer.shenar@gmail.com for questions/inquires
# Algorithm and examples in Gonzales & Levato 2006, A&A, 448, 283; Shenar et al. 2021, A&A, 639, 6; Shenar et al. 2022, A&A, 665, 148
# Current version only applicable for binaries. 
# Input: input file and observed spectra.
# Output: chi2 map on K1,K2 plane (if requested) and separated spectra.
# See "input file" for more documentation
# Coming up in upcoming versions: Higher-order multiples and nebular contamination. 



from Input_disentangle import *
from Disentangling.disentangle_functions import * 

##################################################################
####### READING OF DATA --- USER POTENTIALLY NEEDS TO EDIT #######
##################################################################

if ObsFormat=='TXT':
    PhaseFiles = ascii.read(ObsPath + 'ObsDat.txt')
    MJDs = PhaseFiles['MJD']
    specnames = np.array([ObsPath + el for el in PhaseFiles['obsname']])
elif ObsFormat=='FITS':
    specnames = glob.glob(ObsPath+ '/*.fits')
    MJDs = np.array([])
    
    
############################################################
####### Starts code -- don't touch unless necessary! #######
############################################################


# Avoid weird errors:
if PLOTITR and PLOTCONV:
    print("ERROR: avoid having both PLOTITR and PLOTCONV true")
    sys.exit("Exit in disentangle_shift_and_add.py")
if PLOTEXTREMES and PLOTFITS:
    print("ERROR: avoid having both PLOTEXTREMES and PLOTFITS true")
    sys.exit("Exit in disentangle_shift_and_add.py")
if CompNum != 2:
    print("ERROR: Can't handle more than two components with current version.")
    sys.exit("Exit in disentangle_shift_and_add.py")

# vector of light ratios l1, l2, l3, l4
lguessVec = [lguess1] + lguessVec

if lguessVec[1] == 0.:
    print("light ratio of component 2 cannot be zero")
    sys.exit()
if CompNum==3:
    if lguessVec[2] == 0.:
        print("light ratio of component 3 cannot be zero")
        sys.exit()
elif CompNum==4:
    if lguessVec[2] == 0. or lguessVec[3] == 0.:
        print("light ratio of components 3 or 4 cannot be zero")
        sys.exit()
        

S2Ns = []
ObsSpecs = []
# Read spectra and potentially dates
for i, filepath in enumerate(specnames):
# If fits, read dates as well    
    if ObsFormat=='FITS':
        header = fits.getheader(filepath)
        MJDs = np.append(MJDs, header[MJDHeader])
# read_file returns a 2D array of waves vs. flux for given file path. User can edit "read_file" function if needed.
    spec = read_file(filepath)
# Small script for cleaning spectra of cosmics; use at own risk!
    if CleanCos:
        print("Cleaning Cosmics...")
        SpecClean = Cosclean(np.copy(spec))
    else:
        SpecClean = np.copy(spec)
# Can also re-normalise spectra... better avoid
    if Renormalise and GridDis:
        SpecNorm = Normalise(np.copy(SpecClean, points=NormPoints))
    else:
        SpecNorm= np.copy(SpecClean)
    ObsSpecs.append(SpecNorm)
# Compute S2N of spectrum in prespecified range    
    waves = SpecNorm[:,0]
    fluxes = SpecNorm[:,1]
    S2Nrange  = (waves > S2Nblue) * (waves < S2Nred)  
    S2Ns.append(1./np.std(spec[S2Nrange,1]))

S2Ns = np.array(S2Ns)

# Form array of "force negativity" conditions
StrictNeg = [StrictNegA, StrictNegB, StrictNegC, StrictNegD]
PosLimCond = [PosLimCondA, PosLimCondB, PosLimCondC, PosLimCondD, PosLimCondNeb]
    
    
# Compute true anomalies of data    
phisData = (MJDs-Orbital_Params['T0'])/Orbital_Params['Period'] - ((MJDs-Orbital_Params['T0'])/Orbital_Params['Period']).astype(int)
MsData = 2 * np.pi * phisData
EsData =  Kepler(1., MsData, Orbital_Params['ecc'])
eccfac = np.sqrt((1 + Orbital_Params['ecc']) / (1 - Orbital_Params['ecc']))
nusdata = 2. * np.arctan(eccfac * np.tan(0.5 * EsData))



## Determines by how much "negative" spectra can be above 1.

S2Nsmean = np.mean(S2Ns) *np.sqrt(len(S2Ns))


#Poslimall = [ForceNegSigma/(S2Nsmean * lguessVec[0]), ForceNegSigma/(S2Nsmean * lguessVec[1])]
Poslimall = ForceNegSigma/S2Nsmean

#print(Poslimall)
#dsads
        
# weighting for final co-added spectrum
weights = S2Ns**2 / np.sum(S2Ns**2)

# Grid on which the spectra are calculated on (taken here as longest and densenst wavelength grid among spectra):
for i, spec in enumerate(ObsSpecs):
    Delta = np.amin(np.abs(np.diff(spec[:,0])))
    w1 = spec[0,0]
    w2 = spec[-1,0]
    if i==0:
        DeltaMin = Delta
        w1Min = w1
        w2Max = w2
    else:
        DeltaMin = min(DeltaMin, Delta)
        w1Min = min(w1Min, w1)
        w2Max = max(w2Max, w2)
    
wavegridall = np.arange(w1Min, w2Max, DeltaMin)


# Initialize wavelength arrays for individual line disentangling

wavegridDiffCondall = np.array([(wavegridall > Range[0])*(wavegridall < Range[1]) for Range in Ranges])
wavegrid = wavegridall[np.sum(wavegridDiffCondall,axis=0).astype(bool)]
waveRanges = [wavegridall[el.astype(bool)] for el in wavegridDiffCondall]


##Initialize array of components
CompArr = [interp1d(wavegrid, np.zeros(len(wavegrid)),bounds_error=False, fill_value=0.) ] * CompNum


# Initialize search arrays for K1, K2, K3, K4

if DenseKArr[0] == 1:
    K1s = np.array([Orbital_Params['K1'] ])
else:
    K1s = np.linspace(IniFacKArr[0]*Orbital_Params['K1'], FinFacKArr[0]*Orbital_Params['K1'], DenseKArr[0])
if DenseKArr[1] == 1:
    K2s = np.array([Orbital_Params['K2'] ])
else:
    K2s = np.linspace(IniFacKArr[1]*Orbital_Params['K2'], FinFacKArr[1]*Orbital_Params['K2'], DenseKArr[1])
if DenseKArr[2] == 1:
    K3s = np.array([Orbital_Params['K3'] ])
else:
    K3s = np.linspace(IniFacKArr[2]*Orbital_Params['K3'], FinFacKArr[2]*Orbital_Params['K3'], DenseKArr[2])
if DenseKArr[3] == 1:
    K4s = np.array([Orbital_Params['K4'] ])
else:
    K4s = np.linspace(IniFacKArr[3]*Orbital_Params['K4'], FinFacKArr[3]*Orbital_Params['K4'], DenseKArr[3])    



if not os.path.exists('Output/'):
    os.mkdir('Output/')


# Run main disentangling routine
if GridDis:
# Compute RVs for comp1, comp2    
# Compute RVs for comp1, comp2      
   vrads1, vrads2 = v1andv2(nusdata, Orbital_Params)     
   ScalingNeb = np.ones(len(vrads1))
   kcount_extremeplot = np.argmin(np.abs(K1s - Velo_plot_usrK1_ext)) * len(K2s) + np.argmin(np.abs(K2s - Velo_plot_usrK2_ext))
   kcount_usr = np.argmin(np.abs(K1s - Velo_plot_usrK1)) * len(K2s) + np.argmin(np.abs(K2s - Velo_plot_usrK2))
   K1, K2 = Grid_disentangling2D(waveRanges, nusdata, CompArr[1], Orbital_Params, K1s, K2s,  ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall, MJDs, phisData,  specnames, Rangestr, StarName, ScalingNeb, Ini='B', ShowItr=False,   InterKind=InterKind, itrnumlim = itrnumlim,   PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,  PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=kcount_extremeplot, linewidExt=linewidExt, CompNum=CompNum, ParbSize=ParbSize, N_Iteration_Plot=N_Iteration_Plot,  NebLines = NebLines, NebFac=1, kcount_usr=kcount_usr, ExtremesFigSize=ExtremesFigSize)
   Orbital_Params['K1'] = K1
   Orbital_Params['K2'] = K2
else:  
    PLOTEXTREMES = False
    PLOTFITS = False
    K2 = Orbital_Params['K2']
    K1 = Orbital_Params['K1']
    ScalingNeb = np.ones(len(phisData))
    #K2=2*K1

print("K2 found:", K2)
if K2 < 0:
    print('setting K2=2*K1')
    Orbital_Params['K2'] = 2*Orbital_Params['K1']
print("disentangling...., K1, K2:", Orbital_Params['K1'], Orbital_Params['K2'])
vrads1, vrads2 = v1andv2(nusdata, Orbital_Params)
itrnumlim=NumItrFinal
NebSpec = np.zeros(len(wavegridall))
DisSpecVector, redchi2 = disentangle(np.zeros(len(wavegridall)), vrads1, vrads2,  wavegridall, ObsSpecs, weights, StrictNeg, PosLimCond, Poslimall,  nusdata, Orbital_Params, K1s, K2s, MJDs, phisData, specnames, Rangestr, StarName,  ScalingNeb, NebSpec, Resid=False, Reduce=True, ShowItr=True, Once=True,   InterKind=InterKind, itrnumlim=itrnumlim,   PLOTCONV=PLOTCONV, PLOTITR=PLOTITR,  PLOTEXTREMES=PLOTEXTREMES, PLOTFITS=PLOTFITS, kcount_extremeplot=0, linewidExt=linewidExt, CompNum=CompNum, N_Iteration_Plot=N_Iteration_Plot,  NebLines = NebLines, NebFac=1, ExtremesFigSize=ExtremesFigSize)


# These are the final, scaled spectra:
A, B, NebSpec = DisSpecVector

A = (A-1)/ lguessVec[0] + 1.    
B = (B-1)/lguessVec[1] + 1.
if CompNum>=3:
    C = (C-1)/lguessVec[2] + 1.
if CompNum==4:
    D = (D-1)/lguessVec[3] + 1.



plt.plot(wavegridall, A, label='dis A')
plt.plot(wavegridall, B, label='dis B')
if NebLines:
    plt.plot(wavegridall, NebSpec, label='dis Neb')




np.savetxt('Output/ADIS_lguess2_K1K2=' + str(lguessVec[1]) + '_' + str(Orbital_Params['K1']) + '_' +  str(Orbital_Params['K2']) + '.txt', np.c_[wavegridall, A])
np.savetxt('Output/BDIS_lguess2_K1K2=' + str(lguessVec[1]) + '_' + str(Orbital_Params['K1']) + '_' + str(Orbital_Params['K2']) + '.txt', np.c_[wavegridall, B])
if NebLines:
    np.savetxt('Output/NebDIS_lguess2_K1K2=' + str(lguessVec[1]) + '_' + str(Orbital_Params['K1']) + '_' + str(Orbital_Params['K2']) + '.txt', np.c_[wavegridall, NebSpec])

plt.legend()
plt.show()






