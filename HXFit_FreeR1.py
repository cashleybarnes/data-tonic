#Calculate Fast HX Rates by using lmfit to fit df and t0 globally (used for HX data in reference:https://pubs.acs.org/doi/abs/10.1021/jacs.9b03116)
#Then use scipy curve_fit to fit k and get its error
#CABarnes

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from lmfit import minimize, Minimizer, Parameters, report_fit
from scipy.optimize import curve_fit
import sys, math

pdf4CF = 'HX_freeR1.pdf'

vdList = np.array([0.005,0.025,0.065,0.145,0.255,0.405,1.000]) #variable dephasing delay list in seconds

R1w    = 0.33 #Water R1
T1w    = 1.0/R1w #Water T1
d1     = 5.0 #Recyle Delay
#vd     = 1.0 

#c0     = 0.925

#extract columns to fit using genfromtxt 
tab4R1    = './Exp_HN_R1.txt' #experimental HN R1
resNameR1 = np.genfromtxt(tab4R1,         usecols=(0),         dtype='U'    )
resNumbR1 = np.genfromtxt(tab4R1,         usecols=(1),         dtype='int'  )
R1Data    = np.genfromtxt(tab4R1,         usecols=(4),         dtype='float')

resNameHX = np.genfromtxt('./peakInt.txt',  usecols=(0),         dtype='U'    )
resNumbHX = np.genfromtxt('./peakInt.txt',  usecols=(1),         dtype='int'  )
hxData    = np.genfromtxt('./peakInt.txt',  usecols=range(2,11), dtype='float')

#Determination of the noninverted signal intensity  
wdData    = hxData[:,0:7]
slope     = (hxData[:,8] - hxData[:,7])/(vdList[-1]-vdList[0]) #(Int_noninvert (5ms) - Int_noninvert(1s))/ (1s - 5ms)
wuCorr    = np.kron(slope[:, np.newaxis], vdList-vdList[0])
wuData    = wuCorr + hxData[:, 7, np.newaxis]
hxData    = (wuData - wdData)/hxData[:, 8, np.newaxis] # Corr_Int_noninvert(t) - Int_invert(t)/Int_noninvert(1s)

del slope, wdData, wuCorr, wuData

saveData   = 0 #flag to save the data

if saveData == 1:
   fd = open('hxData.txt', 'w')
   nres, ncol = hxData.shape
   for i in range(nres):
      fd.write('%s %3d' % (resNameHX[i], resNumbHX[i]))
      for j in range(ncol):
         fd.write(' %6.3f' % hxData[i, j])
      fd.write('\n')
   fd.close()
   sys.exit()

if len(resNumbR1)==len(resNumbHX)==len(resNameR1)==len(resNameHX):
   if R1Data.shape[0] == hxData.shape[0]:
      for i in np.arange(resNumbR1.size):
         if (resNameR1[i] != resNameHX[i] or resNumbR1[i] != resNumbHX[i]):
            print('Inconsistent residues between R1 and HX')
else:
   print('Inconsistent list length')
   sys.exit()

#Equation Reference: https://onlinelibrary.wiley.com/doi/full/10.1002/pro.582

def eq4HX(t, df, k, R1, watR1, t0):
    R = (R1+k-watR1)
    return df*(k*np.exp(-(t+t0)*watR1)/R)*(1-np.exp(-(t+t0)*R))

def calcR(params, t, hxData, expR1): #For Calculation of R1 according to the HXdata and starting at the expR1 -- free R1
    nres, ncol = hxData.shape

    df    = params['df'].value #delta f --> fractional steady state water signal 
    t0    = params['t0'].value #correction to make t0 fit from zero
    watR1 = params['watR1'].value #longitudinal spin relaxation rate of water

    rsd   = np.zeros((nres, ncol))

    for i in np.arange(nres): #number of residues
       k         = params['k_%i'  % i].value
       R1        = expR1[i]
       rsd[i, :] = hxData[i, :] - eq4HX(t, df, k, R1, watR1, t0) #residual between fitted and experimental HX rates 

    return rsd.flatten() #flatten to a 1D array for minimize

p = Parameters()
p.add( 'df',    value=1.80,    vary=True,  min=1.00, max=2.500) 
p.add( 't0',    value=0.00125, vary=False, min=0.00, max=0.0025) 
p.add( 'watR1', value=R1w,     vary=False, min=0.20, max=1.000) 

for ires, icol in enumerate(hxData):
    p.add( 'k_%i'  % ires, value=1.0, vary=True, min=0, max=50)

m   = Minimizer(calcR, p, fcn_args=(vdList, hxData, R1Data)) 
mr2 = m.minimize(method='least_squares')

#take globally fit df/t0/R1/R1w from above and fits each residue independently using scipy curve fit

def eq4HX2(x, k, currR1):
   return eq4HX(x, df, k, currR1, watR1, 0.0)

nres = hxData.shape[0]

kex = []
ker = []

R1  = []
R1e = []

df    = mr2.params['df'].value
t0    = mr2.params['t0'].value
watR1 = mr2.params['watR1'].value

tdata = np.array(np.arange(0, 1.1, 0.01))

#print("#ASS   Kex    error      R1")
for i in np.arange(nres):
   fit, cov = curve_fit(eq4HX2, vdList+t0, hxData[i,:], p0=[2.0, 1.0], bounds=[[0.0, 0.0], [50.0, 5.0]])
   kex.append(fit[0])
   ker.append(np.sqrt(cov[0,0]))
   R1.append(fit[1])
   R1e.append(np.sqrt(cov[1,1]))

   print("%s %2d %8.3f %6.3f %8.3f %6.3f" % (resNameHX[i], resNumbHX[i], kex[i], ker[i], R1[i], R1e[i]))

with PdfPages(pdf4CF) as pdf:
   for i in np.arange(nres):
      yfit = eq4HX2(tdata, kex[i], R1[i])
      plt.plot(tdata, yfit, 'g-')
      plt.plot(vdList+t0, hxData[i, :], 'go')
      plt.xlim([0, 1.1])
      plt.ylim([0, 2.0])
      plt.xlabel("time (s)")
      plt.ylabel(r'$\Delta I/I_0$')
      plt.text(0.65, 1.9, 'Residue '+resNameHX[i]+str(resNumbHX[i]))
      plt.text(0.65, 1.8, 't0 = ' + "{:.3f}".format(t0))
      plt.text(0.65, 1.7, r'$\Delta_f$ = ' + "{:.3f}".format(df))
      plt.text(0.65, 1.6, r'$R_1$ = '+"{:.3f}".format(R1[i]) + ' +/- ' + "{:.3f}".format(R1e[i]) + r' ($s^{-1}$)')
      plt.text(0.65, 1.5, r'$R_{1w}$ = '+"{:.3f}".format(watR1))
      plt.text(0.65, 1.4, 'k = '+"{:.3f}".format(kex[i]) + ' +/- ' + "{:.3f}".format(ker[i])  + r' ($s^{-1}$)')
      f = plt.gcf()
      pdf.savefig(f)
      plt.close(f)

