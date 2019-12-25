#Extract Residue Specific Chemical Shifts from nmrDraw Assignment Tabs using Pandas
import pandas as pd

columns = ['INDEX','X_AXIS','Y_AXIS','DX','DY','X_PPM','Y_PPM','X_HZ','Y_HZ','XW','YW','XW_HZ',
           'Y_HZ','X1','X3','Y1','Y3','HEIGHT','DHEIGHT','VOL','PCHI2','TYPE','ASS','CLUSTID','MEMCNT']

df1=pd.read_csv("CS_av900.tab", delimiter="\s+", names=columns, skiprows = [0, 1, 2, 3, 4])

HN = df1[['ASS','X_PPM', 'Y_PPM']]

#Sort Chemical Shifts According to Residue Number while keeping the residue name
HN['sort'] = HN['ASS'].str.extract('(\d+)', expand=False).astype(int)
HN.sort_values('sort',inplace=True, ascending=True)
HN = HN.drop('sort', axis=1)

HN.rename(columns={'ASS': 'ASS', 'X_PPM':'1H_PPM', 'Y_PPM':'15N_PPM'}, inplace=True)

#Correction of Chemical Shifts for TROSY Effects
Add_Sub = 46
N = 91.235
H = 900.274

HN['1H_Corr_ppm'] = HN['1H_PPM'] + Add_Sub/H
HN['15N_Corr_ppm'] = HN['15N_PPM'] - Add_Sub/N

HN_Corr = HN[['ASS', '1H_Corr_ppm', '15N_Corr_ppm']]

#Write Corrected Chemical Shifts to csv file for further analysis
HN_Corr.to_csv('SAH_Chemical_Shift_Table_TrCorr_av900.csv', index=False, float_format='%.2f')
