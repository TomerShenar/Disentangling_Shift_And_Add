# Disentangling_Shift_And_Add
21.11.2022, V1.0;  contact: T.Shenar@uva.nl or tomer.shenar@gmail.com

shift-and-add & grid disentangling, written by Tomer Shenar, with contributions from Matthias Fabry & Julia Bodensteiner

Algorithm, examples & tests in Gonzales & Levato 2006, A&A, 448, 283; Shenar et al. 2021, A&A, 639, 6; Shenar et al. 2022, A&A, 665, 148

Current version only applicable for binaries. 

Input: input file and observed spectra.

Output: chi2 map on K1,K2 plane (if requested) and separated (disentangled) component spectra, scaled by user-provided light ratio.

See "input file" for more documentation

Coming up in upcoming version: Higher-order multiples and better handling of nebular contamination. 

#########################
###Quick instructions:###
#########################

To retrieve the files, make a folder in which you want to store the files and type:

git clone https://github.com/shtomer/Disentangling_Shift_And_Add.git

The download includes:
1. "disentangle_shift_and_add.py" is the main code that runs the script. It should not be editted unless you are an experienced user, except perhaps for the first part, where the spectra and dates are read, if needed.
2. "Input_disentangle.py" is where the user inserts multiple options and parameters for the disentangling procedure. This script is documented and should be filled as instructed in the file.
3. "disentangle_functions.py" contains the main functions used by the script. Do not edit unless you know what you're doing.
4. "make_spectra_SB2.py" is a small script that creates mock data of an hypothetical binary with prespecified templates, resolution & S/N (option: nebular contamination)
5. "2DCont.py" reads a "chi2 output" from the "disentangle_shift_and_add.py" output and produces K1,K2 measurements and a 2D color map.
6. "plot.py" is a helpful script to easily plot spectra (ascii/fits)
7. a "Models" directory with a few selected TLUSTY & PoWR models of OB and WR stars, and a HERMES observation of a Be star; for creating mock data.

Instructions:

1. Place the file "disentangle_functions.py" in a directory in which they store typical python libraries in a directory called "Disentangling". For example, on my machine, this file is in: "/Users/tomer/Programs_and_packages/Disentangling/"

2. Add the line PYTHONPATH="/Users/tomer/Programmes_and_packages:$PYTHONPATH" in your .bashrc file in the home directory.

3. Place the "Input_disentangle.py" and "disentangle_shift_and_add.py" in a directory in which you want the analysis to be performed.

4. For the script to work, make sure the following packages are installed (e.g. via pip install):
glob 
os 
numpy 
astropy
scipy
pandas 
matplotlib
sys

5. Prepare a folder that contains all spectra for your object. The data could be in ASCII format or FITS format. If FITS, the script will attempt to read the dates from the .fits files. If ASCII, the script will require a file named ObsDat.txt listing the MJDs and file names in the observational data folder. For example, use the "make_spectra_SB2.py" to create such a folder with mock data.

6. Fill out the "Input_disentangle.py" as per instructions. Importantly, the path to the data and data format should be specified, and the orbital parameters should be appropriate (if the results are to make sense).

7. After setting the relevant options in "Input_disentangle.py", run the script by typing "python disentangle_shift_and_add.py". 

8. If the run was successful, an output directory will be created containing the disentangling spectra stored as: "ADIS_lguess2_K1K2=0.7_87.0_135.0.txt

For the script to run,  the script will need corrects paths to the observations were provided, and if the observations are correctly read.
I tried to make the script relatively flexible (see documentation in script), but it is possible that individual formats will not be read successfuly. The user is encouraged to replace the appropriate reading functions with whatever functions they deem appropriate.

After a successful execution of the script, an output directory will be created containing the disentangled spectra ADIS and BDIS. The names of the files code the K1, K2 pair for which the separation was performed, a few files documenting the chi2 arrays for K1, K2, and a few figures.


Q&A:

Q. What is spectral disentangling?

Spectral disentangling refers to the simultaneous derivation of orbital parameters of a spectroscopic multiple system and the separation of the composite spectra to the component spectra. 
In general, the algorithm receives a set of N spectra, as returns as output the orbital parameters and the individual spectra of each component.
More information, references, and documentation can be found here: http://sail.zpf.fer.hr/fdbinary/

Q. What is the technique used here?

Spectral disentangling can be performed either by minimizing a large set of linear equations, or via an iterative procedure. I refer the reader to http://sail.zpf.fer.hr/fdbinary/ for a description and overview of the various methods.

Here, the iterative shift-and-add algorithm is implemented.

Q. How does shift-and-add work?

The shift-and-add technique uses approximations for the component spectra in the i'th iteration, A[i] and B[i], to compute better approximations in the i+1'th iteration.

Consider a set of continuum-subtracted normalised spectra (i.e., continuum = 0)
Let us assume that each observed spectrum comprises the sum of two stellar spectra that are time-invariant with exception of their orbital Doppler motion. 
Let us further assume that we possess full knowledge of the orbit, i.e. we can compute the radial velocities for each epoch of observation for both stars.
The algorithm is as follows:

1. Start by assuming a flat spectrum for the secondary, B[0] = [0, 0, 0, ...]

2. Compute A[i] by subtracting B[i] from all observed spectra by shifting it to its respective RV, co-add all residuals in the frame of A.

3. Compute B[i+1] by subtracting A[i] from all observed spectra, and co-add residuals in the frame of B.

4. Repeat Nitr times until satisfied with convergence. Obtain final spectra A = A[Nitr] and B = B[Nitr]

In general, however, not all orbital parameters are known. In this case, one can explore the agreement with the data for a given set of orbital parameters.
For example, imagine we perform the separation for a set of orbital parameters O', and obtain A' and B' as disentangled spectra. We can compare the sum of A' and B' at each epoch to the observation, and compute the usual chi2 across all pixles and all observations. In principle, one may expect that the best results (minimum chi2) would be obtained for the correct orbital parameters.

In principle, the chi2 exploration could be performed on all orbital paramreters. However, for feasibility reasons, the script is limited to a K1, K2 exploration (i.e., RV semi-amplitudes). In almost all cases, the period, T0, eccentricity and omega can be well constrained by fitting a single-lined orbit to the primary.

A minimization of chi2(K1, K2) enables the derivation of K1, K2 and their corresponding statistical errors. Once K1, K2 have been derived, the final separation can be performed.

Q. What about the light ratios?

The light ratios determine the final scaling of the disentangled spectra. In principle, this information is not contained (mathematically) in the data; the technique cannot tell whether a star has "intrinsically weak lines" or whether it is diluted. Hence, the user needs to specify the light ratios, which are currently assumed constant throughout the wavelength domain. 
Luckily, the light ratios only impact the final scaling of the spectra, and have no impact on the results of disentangling. They can therefore be separated from the problem of disentangling.

Q. Do I disentangle specific lines or the entire spectrum?

In principle ,this is up to the user. However, it is strongly encouraged to disentangle individual lines or regions of lines. The method becomes less efficient if large chunks of continuum are included. Moreover, by computing K1, K2 for individual lines, one may obtain independent K1,K2 measurements, which can allow the user to compute weighted means for K1 and K2.

Q. How do I pick the number of iterations Nitr?

This is a non-trivial question, and in truth, no satisfactory general condition could still be found. As a rule-of-thumb, when the RV amplitudes are comparable or larger than the typical widths of the line profiles (set by either rotation or resolution), then tens of iterations suffice. However, when the RV amplitudes are smaller, hundreds or even thousands of iterations may be necessary. 
The user is urged to make use of the various convergence plots and tests to study this question. For example, the user may want to compare the appearance of the disentangled spectra each X iterations (which is an option in the script), or they may want to investigate how their derived K1, K2 values differ if the iteration number increases.
The number of iterations further depends on the necessary level of accuracy. For example, for very high S/N data, more iterations could make sense, while for low S/N, the differences become rapidly negligible. Similarly, if the secondary is very faint, then small discrepancies can translate into very large differences in its disentangled spectrum, and more iterations may be needed. 

Q. Is it better to derive K1, K2 with disentangling or standard RV measurement methods?

If possible, it is always better to determine at least K1, and if possible also  K2, via RVs, since it will result in much more accurate measurements. This can be done if the different components show different lines, or if one of them strongly dominates (l1 > ~ 90%) over the other. If this is not the case and the lines are constantly blended, then disentangling is the method of choice. Note, however, that typical errors are at least X10 times larger than with RV measurements due to the large freedom in "shaping" the disentangled spectra. 

Q. Should I enforce the spectra to lie below the continuum?

The s&a technique is known to often cause "emission wings" and "extended troughs" next to broad absorption lines, and the reason for this is well documented in Quintero 2020, AN, 341, 628. 
A very simple yet powerful solution is to enforce the components to lie below the continuum, unless the user suspects that specific regions may exhibit emission. The user is therefore encouraged to make use of this option. In the case of emission line stars such as Be stars or Wolf-Rayet stars, this option should only be used in continuum regions or absorption-line regions.

Q. Can I use different datasets (e.g., different instruments)

Yes; the script should be able to handle this, as long as it can read the observations.


Q. What if the stars are variable (e.g., pulsations)?

Excellent question. Hopefully, I'll have the time to explore this soon enough. However, as long as the variability does not follow the orbital period, one can expect that it will wash away in the process of shifting-and-adding. 

Q. What about triples? Quadruples? 

The method can be extended to multiple components in a straight-forward manner, as has been in fact used on quadruples (e.g., Shenar et al. 2022, A&A, 665, 148). The current script is designed for binaries only, but a follow-up version will include a module for higher-order multiples 

Q. What is the impact of using wrong orbital parameters (specifically K1, K2) on the disentangled spectra?

While it is tough to provide a single answer here, generally it holds: the brighter the component, the more critical your results would depend on its K-value. For example, consider a binary with the primary contribution 95% and the secondary 5%. You will find that differences of ~5% in the value of K1 could translate into quite substantial differences in the disentangled spectrum of the secondary. In contrast, the disentangled spectra would be very weakly dependent on K2 in this case.
An important example is when trying to establish whether or not the companion is a star or a compact object. For SB1 systems, adopting a K1 value 5-10% off the "true" value will result in a disentangled spectrum for the secondary which mimics the appearance of the primary. This may look like a star to the inexperienced user. Hence, test how the disentangled spectrum of the secondary varies are values for K1 are altered.

Q. How can I know my results are sensible?

Test, test, test! Try to simulate spectra corresponding to the case of interest, and see how well (or not) the method is able to retrieve the orbital parameters. This is the only way to know for sure whether you can trust the results. 



