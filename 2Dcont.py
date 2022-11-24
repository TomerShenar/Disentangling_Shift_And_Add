# For plotting 1-D and 2-D chi2 maps and compute K1+K2 + error estimate.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure



###########################
######### USER INPUT ######
###########################

# Cut part of the chi2 matrix?
indK1cut1, indK1cut2  = 0, -1
indK2cut1, indK2cut2 = 0, -1


# How many indices next to minimum should parabolas be fit
step1 = 3
step2 = 3

# File to be read
fileK1K2 = 'Output/HeI4472_grid_dis_K1K2.txt'
#fileK1K2 = 'HeII4546_grid_dis_K1K2.txt'
#fileK1K2 = 'HeII4200_grid_dis_K1K2.txt'
#fileK1K2 = 'HeI4388_grid_dis_K1K2.txt'
#fileK1K2 = 'HeI4026_grid_dis_K1K2.txt'
#fileK1K2 = 'HeI4009_grid_dis_K1K2.txt'
#fileK1K2 = 'HeI4144_grid_dis_K1K2.txt'
#fileK1K2 = 'Hdelta_grid_dis_K1K2.txt'
#fileK1K2 = 'Indiv_grid_dis_K1K2.txt'


###########################
#######END USER INPUT #####
###########################


Z = np.loadtxt(fileK1K2)

#MinInds =  np.where(Z == np.min(Z))

try:
    fileK1 = fileK1K2[:-8] + 'K1.txt'
    sigma1 = float((open(fileK1).readlines())[0].split()[-1])
    K1s, chis1 = np.loadtxt(fileK1, unpack=True)
    minind1 = np.argmin(chis1)
except: 
    sigma1=0
    K1s = [0]
    pass
try:
    fileK2 = fileK1K2[:-8] + 'K2.txt'
    sigma2 = float((open(fileK2).readlines())[0].split()[-1])
    K2s, chis2 = np.loadtxt(fileK2, unpack=True)   
    minind2 = np.argmin(chis2)
except:
    K2s = [0]
    sigma2=0
    pass
sigma = (sigma1 + sigma2)/2.




if len(K1s) > 1:
    a,b,c = np.polyfit(K1s[max(minind1-step1, 0): min(minind1+step1, len(K1s)-1)], chis1[max(minind1-step1, 0): min(minind1+step1, len(K1s)-1)] , 2)    
    Chi2K1min = -b/2./a

    K1fine = np.arange(-K1s[-1], K1s[-1]+1, 0.01)

    parb = a*K1fine**2 + b*K1fine + c

    LeftPart = K1fine < Chi2K1min
    K1err = Chi2K1min - K1fine[np.argmin((parb[LeftPart] - sigma1)**2)]

    plt.scatter(K1s, chis1)
    plt.plot(K1fine, parb)
    plt.plot([min(K1s), max(K1s)], [sigma1, sigma1], color='red', label='1-sigma')
    plt.ylabel('reduced chi2')
    plt.xlabel('K1 [km/s]')
    plt.xlim(K1s[0]-1, K1s[-1] + 1)
    plt.ylim(min(chis1)*.95, max(chis1)*1.05)
    plt.text(np.average(K1s), max(chis1), r"$K_1 =$ " + str(round(Chi2K1min,1)) + r'$\pm$' + str(round(K1err,1)) +  " km/s", size=13, horizontalalignment='center')
    plt.show()
    print("K1: ", Chi2K1min, " +- " , K1err)

if len(K2s) > 1:
    a,b,c = np.polyfit(K2s[max(minind2-step2, 0): min(minind2+step2, len(K2s)-1)], chis2[max(minind2-step2, 0): min(minind2+step2, len(K2s)-1)] , 2)    
    Chi2K2min = -b/2./a

    K2fine = np.arange(-K2s[-1], K2s[-1]+1, 0.01)

    parb = a*K2fine**2 + b*K2fine + c

    LeftPart = K2fine < Chi2K2min
    K2err = Chi2K2min - K2fine[np.argmin((parb[LeftPart] - sigma2)**2)]

    plt.scatter(K2s, chis2)
    plt.plot(K2fine, parb)
    plt.plot([min(K2s), max(K2s)], [sigma2, sigma2], color='red', label='1-sigma')
    plt.ylabel('reduced chi2')
    plt.xlabel('K2 [km/s]')
    plt.xlim(K2s[0]-1, K2s[-1] + 1)
    plt.ylim(min(chis2)*.95, max(chis2)*1.05)
    plt.text(np.average(K2s), max(chis2), r"$K_2 =$ " + str(round(Chi2K2min,1)) + r'$\pm$' + str(round(K2err,1)) +  " km/s", size=13, horizontalalignment='center')    
    plt.show()
    print("K2: ", Chi2K2min, " +- " , K2err)




###########################
####### 2D stuff      #####
###########################

if len(K1s) > 1 and len(K2s) > 1:



    #head = open(fileK1K2).readlines()[0].split()
    
    MinChis1, MaxChis1, MinChis2, MaxChis2 = np.amin(chis1), np.amax(chis1), np.amin(chis2), np.amax(chis2)
    
    LimMin = min(MinChis1, MinChis2)
    LimMax = max(MaxChis1, MaxChis2)
    
    Chi2K1min, Chi2K2min =  K1s[minind1], K2s[minind2]


    K1min = minind1

    fig, ax = plt.subplots()

    X, Y = np.meshgrid(K2s[indK2cut1:indK2cut2], K1s[indK1cut1:indK1cut2])

    cs = plt.contourf(X, Y, Z[indK2cut1:indK2cut2, indK1cut1:indK1cut2], 200, cmap='inferno_r',  vmin=LimMin, vmax=LimMax)

    cont = ax.contour(X, Y, Z[indK2cut1:indK2cut2, indK1cut1:indK1cut2], [sigma], linewidths=2, linestyles='dashed', colors='green')

    vertices = (cont.collections[0].get_paths()[0].vertices)
    K2min, K2max = np.amin(vertices[:,0]), np.amax(vertices[:,0])
    K1min, K1max = np.amin(vertices[:,1]), np.amax(vertices[:,1])





    plt.xlabel(r'$K_2 [{\rm km}\,{\rm s}^{-1}]$')
    plt.ylabel(r'$K_1 [{\rm km}\,{\rm s}^{-1}]$')
    cbar = plt.colorbar()
    tickchi = np.arange(0.6, 1., 0.01)
    cbar.ax.set_ylabel(r'$\chi_{\rm reduced}^2$', size=12)
    plt.title(fileK1K2)
    plt.savefig(fileK1K2[-4] + '_2Dmap.pdf')
    plt.show()
