#Calculate the Normalized i+/-1 Cross Peak Intensities of NOESY spectra
from collections import defaultdict
import math 
import re

#nmrDraw Cross and Diagonal Peak Assignment Tab
infile = open('4DNOESY.tab')
contents = infile.read()
infile.close()

outfile = open('Test.txt', 'w')

split_lines = contents.split('\n\n')[2].strip().split('\n')

for row in split_lines:
    ass = row.strip().split()[38].replace('+1', '.1').replace('-1','.2')
    #ass3 = re.sub("\-", '.', ass2)
    H = float(row.strip().split()[33])
    if ass != 'None':
        outfile.write('%-5s %.f\n' % (ass, H))

outfile.close()

infile2 = open('Test.txt')
contents2 = infile2.readlines()
infile2.close()

contents2.sort(key=lambda line : float(line[1:3]))

cross1 = re.compile("[A-Z]*[0-9]*[.][1]*\Z")
cross2 = re.compile("[A-Z]*[0-9]*[.][2]*\Z")
ass = re.compile("[A-Z]*[0-9]*\Z")

#Create Dictionary to Collect Residue Numbers as keys with Corresponding Intensities of Diagonal and Cross Peaks as values
parsed = defaultdict(lambda: ['', float('inf')]*3)

for row in contents2:
    x = row.split()[0]
    Int = float(row.split()[1])
    resnum = int(x.split('.')[0][1:])
    
    if ass.match(x):
        parsed[resnum][0] = x
        parsed[resnum][1] = Int
     
    elif cross1.match(x):
        parsed[resnum][4] = x
        parsed[resnum][5] = Int
        
    elif cross2.match(x):
        parsed[resnum][2] = x
        parsed[resnum][3] = Int


outfile2 = open('NOESY_Intensity_Calc_RN.txt', 'w')

for i in range(1, 11):
    if i in parsed and i+1 in parsed:
        if (
            parsed[i][1] != float('inf') and
            parsed[i+1][1] != float('inf') and
            parsed[i][3] != float('inf') and
            parsed[i+1][5] != float('inf')
        ):
            v = math.sqrt((parsed[i][3] * parsed[i][5]) / (parsed[i][1] * parsed[i+1][1]))
            outfile2.write ('%s %.3f\n' % (i, v))
            
outfile2.close()
    
