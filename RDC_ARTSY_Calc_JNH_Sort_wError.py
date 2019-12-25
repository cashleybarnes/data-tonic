#Calculate the ARTSY RDCs of Protein
#References:https://pubs.acs.org/doi/abs/10.1021/jacs.9b03116 
#CABarnes
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import math
from matplotlib.backends.backend_pdf import PdfPages

title_font = {'fontname':'Helvetica', 'size':'16', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} 
axis_font = {'fontname':'Helvetica', 'size':'16'}

infile = open ('nlin.tab') #autoFit.tcl output from nmrPipe
contents = infile.read()
infile.close()

outfile = open('RDC_293K.txt', 'w')
#outfile.write ('%s %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n' % ("resN", REF, ATT, Ratio, RDC))

split_lines = contents.split('\n\n')[3].strip().split('\n')

T = 0.01175 #dephasing delay in seconds
JNH = -92 #Isotropic JNH

N1 = 138013.828125 #Noise recorded from spectra

#RDC Calculation from Intensities of interleaved spectra -- Attenuated and Reference
#Reference: https://link.springer.com/article/10.1007%2Fs10858-010-9441-9 

def DNH (T, JNH, Q):
    return (-1/T - JNH + (2/(math.pi*T))*math.asin(Q/2))

def Error (Q, I, N):
    A = 2/(math.pi*T)
    B = I/N
    G = A/B
    C = 1 + Q**2
    D = 4 - Q**2
    return (G*((C/D)**(1/2)))

for row in split_lines:
    x = row.strip().split()[22]
    if x == 'None':
        continue
    else:
        resN = re.sub("[^0-9]", "", x)
        REF = float(row.strip().split()[25])
        Q = float(row.strip().split()[26])
        N = N1/float(row.strip().split()[17])
        Err = Error(Q,REF,N)
        outfile.write('%s %.3f %.3f\n' % (resN, DNH(T,JNH,Q), Err))
        #print('%s %.3f %.3f\n' % (resN, DNH(T,JNH,Q), Err))

outfile.close()

infile2 = open('RDC_293K.txt')
contents2 = infile2.readlines()
infile2.close()

contents2.sort(key=lambda line: int(line.split()[0]))

with open('Sorted_RDC_293K.txt', 'w') as RDC:
    RDC.writelines(contents2)

#Plot and Fit the Correlation Between Two Sets of Data at Two Temps using Pandas, Numpy, and Matplotlib
Columns = ['ASS','RDC','Error']

df1=pd.read_csv("Sorted_RDC_293K.txt", delimiter="\s+", names=Columns)
df2=pd.read_csv("Sorted_RDC_308K.txt", delimiter="\s+", names=Columns)

df1 = df1.merge(df2, on='ASS')

df_comb = df1.values #convert data frame into an array for saving into text file

np.savetxt('Correlation_RDC_293K_308K.txt', df_comb, fmt='%f')

res = df_comb[:,0]
RDC_293K = df_comb[:,1]
Error_293K = df_comb[:,2]
RDC_308K = df_comb[:,3]
Error_308K = df_comb[:,4]

fit_coeff = np.polyfit(RDC_293K,RDC_278K,1) #Fit = m, b --> m*x + b, where x = RDC_293K data
fit_fn = np.poly1d(fit_coeff) #With m and b, take x value and provide estimated y

Figure = 'RDCs_293K_vs_308K.pdf'

with PdfPages(Figure) as pdf:
    plt.errorbar(RDC_293K, RDC_308K, xerr = Error_293K, yerr = Error_308K, ls = 'none', color='black', marker='o', capsize=4, capthick=1, ecolor='black')
    plt.plot(RDC_293K,RDC_308K, 'ko', RDC_293K, fit_fn(RDC_293K), 'k')
    plt.xlabel('RDC at 293K (Hz)', **axis_font)
    plt.ylabel('RDC at 278K (Hz)', **axis_font)
    plt.xticks(np.arange(-2, 40, step=10))
    plt.xlim((-2,40))
    plt.ylim((-2,40))
    #plt.title ('None', **title_font)
    plt.legend(frameon=False)
    fig=plt.gcf()
    pdf.savefig(fig)
    #plt.ion()
    plt.show()

