MATLAB/Octave package
Below are definitions of all output/input/defined-in-file objects, and files for different model variants. 
A. Prep files
Prep pop density input from fluscape data:
[Df,xf,yf]=lscanplot(lscan,fluscapeLocations);
If using pixelated transect:
[Ds,Is]=selectSubgrid(Df,xf,yf)
Else just need a grid of pop densities (lscan)
Prep objects for simulation:
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluSub(lscan,R0,stoch,cellSet); (using pixealted fluscape transect)
OR =prepFlu(lscan,R0,stoch); (using grid: Df or lscan)
Generates all objects needed to simulate:
- gamma (recovery rate): defined in the code, and outputted to workspace
- NN: vector of pixel populations (counting down columns)
- n: number of cells
- NNbar: vector of populations by age (length 4n, first n entries are age-group 1 population of each pixel etc.)
- NNrep: 4 copies of NN (same length as NNbar, but values are total population)
- minNind, maxNind: indices of min and max popualtion (or first occurence thereof) - excluding zeros? These are unused at present
- maxN: maximum population of a given pixel - useful for plot axes
- Kbar, K1, Cbar: expanded matrices for FOI (c.f. first paper)
- betaS, betaI, betaD: betas generated for given R_0, allowing for dependency on mobility assumption
- ages0: initial age distribution (used in stochastic model)
Inputs:
- lscan=Df OR Ds: grid of pop densities
- R0: R_0
- stoch: =0 for final size or ODE (solvetype=1/2), =1 if running stochastic comp. model (solvetype=3) - aged populations are rounded to integers
Defined in script:
- aa/alpha/p (a, alpha and p in offsety gravity kernel)
- gamma
- celldist: distance between cells (km)
- a1immobile: =1 to turn off mobility in children (naive option)
- normaliseKernel: =1 for row=normalised K, =0 otherwise
- Cnum/Cdur: define age-mixing matrices (number and duration of contacts). Set C=Cnum/Cnum*Cdur/ones(4) for number/number x duration/age-off. Note: using more or fewer that 4 age-groups requires additional modifications to the code
Other options:
- prepFluAge - for non-uniform distribution of age
- prepFluDE - delta and epsilon (c.f. first paper)
B. Simulation files
[f,g]=finalSizeMulti2sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual);
OR =finalSizeMulti();
Generates attack rates (f=season-specific, g=total immunity)
Additional inputs (other than outputs of prep files):
- solvetype: =1/2/3 for FSC/ODE/SCM (finalSizeMulti only)
- numseed (see below - currently an input in finalSizeMulti)
- isdual: =0/1/2 for S/dual/I-mobility
Defined in script:
- eps: parameter for Prow2 (see below)
- cross: cross immunity
- demog: =1 for aging on, =0 for aging off
- amax: maximum age (cyclic aging)
- tauend: number of consecutive seasons to simulate
- plotTau: option to plot epi curve for a given year (=0 for no plot - useful for debugging)
- tend: time (days) to simulate for each epidemic
- mu: death rate (in addition to demography - unused)
- phi1, phi2: seasonality (phi=phi1-phi2 cos(pi t/180) - currently turned off in "integr8all")
- numseed: seed proportion
- Prow1, Prow2: Theta_1 and Theta_2 (modification required inside of loop for these objects to depend on year)
Other options:
- finalSizeMultiAge - for non-uniform distribution of age (requires prepFluAge)
- [f,g]=MA1D2sub(eps,cross): 1D 2-subtype model (with population 1)
C. Plotting files
FYI - these are the files that are most often messed with for specific occasions.
trimbyk: for square grid only, takes vector of attack rates and dimensions of original grid and removes k-thick boundary (k defined in code), returning new population and attack rate vectors. Called in plotAR, but must be commented out for non-rectangular inputs.
plotAR: Plot for single column of finalSizeMulti output (f or g - same for pandemic case). Plots singe vector of attack rates against log population.
plotMA: Plot for multiple columns of finalSizeMulti outputs (A=f, Acum=g). Plots single-subtype model output as time series, option to plot only a sample of cells. Colormap shows log population density.
plotMAnospace: plot for output of MA1D2sub (Asus=f, Arec=g). Set ar=1 for season-specific attack, ar=0 for susceptibility profiles.
plotMAspace: plot for finalSizeMulti2sub (Asus=f, Arec=g). Only plots season-specific attack.
There are more files entitled "plotSomething.m" - see filenames/descriptions in the comments for details.
X. Other applications
All files starting with "RP" use a modification of the above model that gives mean-field approximations of the behaviour of epidemics on spatially embedded networks. 
