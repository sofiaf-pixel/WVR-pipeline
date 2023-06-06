import plotElnod as pE
import plotAzscan as pA
import AM_functions as am
import zenith_PWV as zPWV
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
import datetime
import matplotlib.dates as mdates
from math import *
import extract_mod_par as emp
import Read_BICEP_ts as bts
from dateutil import parser



x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()
x_emp = emp.mod_parameters()
x_zen=zPWV.T_zenith()


#WVR options

f = open('analysis_prods/dT150_del_at55.txt','rb')
dTdel_dict=pk.load(f)
f.close()

pl.plot(ax_date, dT_del_list)
pl.show()


sys.exit()

wvr_scan='20200503_170002_skyDip_fast.txt'
days=['10', '20', '30']
months=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
# for month in months:
#     for day in days:
#         os.system('rm -r 2020'+str(month)+day+'_140002_skyDip')
#         os.system('rm -r wvr1_data/2020'+str(month)+day+'_140002_skyDip')
#         os.system('mkdir wvr1_data/2020'+str(month)+day+'_140002_skyDip')
#         os.system('cp /Volumes/LaCie/BICEP/WVR_analysis/wvr1_data_local/2020'+str(month)+day+'_140002_skyDip_fast.txt wvr1_data/2020'+str(month)+day+'_140002_skyDip/2020'+str(month)+day+'_140002_skyDip_fast.txt')






scan_list=[]
for scan in os.listdir('wvr1_data'):
    if scan[:4]=='2020':
        if scan[6:8] in days:
            if scan[9:11]=='14':
                if scan[-6:]=='skyDip':
                    scan_list.append(scan+'_fast.txt')

dT_del_list=[]
ax_date=[]

for wvr_scan in scan_list:
    try:
        dT_del=x_zen.dT_del55(wvr_scan, show_im=0)
        dT_del_list.append(dT_del)
        fn_date=datetime.datetime(int(wvr_scan[:4]), int(wvr_scan[4:6]), int(wvr_scan[6:8]), int(wvr_scan[9:11]))
        ax_date.append(fn_date)
    except Exception as e:
        print(e)

dTdel_dict={'dT_del':dT_del_list, 'date':ax_date}

f = open('analysis_prods/dT150_del_at55.txt','wb')
pk.dump(dTdel_dict, f)
f.close()

pl.plot(ax_date, dT_del_list)
pl.show()
