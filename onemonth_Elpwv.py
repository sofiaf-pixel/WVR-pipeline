import plotElnod as pE
import AM_functions as am
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
import datetime


x_am=am.AM_functions()
x=pE.ReadElnod()

# test_data='20200406_210002_skyDip_fast.txt'
# path_to_test='Trial_folder/20200406_210002_skyDip/'

temp_file='SPole_annual_50.amc'

#x_am.create_am_datfile(test_data, path_to_data=path_to_test, template= temp_file, spline=2, showplots=0)
#x_am.fit_w_am(test_data, path_to_data=path_to_test, template= temp_file, spline=2)
#x_am.fit_w_am(test_data, path_to_data=path_to_test, template=temp_file, spline=0)

#el_list, pwv_list=x_am.plot_am_fit_2(test_data, template= temp_file, spline=2)
#el_list, pwv_list=x_am.plot_am_fit_2(test_data, template= temp_file, spline=0)


month='202007'
trial='Trial_folder/'+month

os.system('rm wvr1_data/Trial_folder/'+month+'*')

if not os.path.exists('wvr1_data/'+trial):
    os.makedirs('wvr1_data/'+trial)


for filelist in os.listdir('wvr1_data/BAElnod_data'):
    print(filelist[0:6])
    if filelist[0:6]==month:
        os.system(f"cp wvr1_data/BAElnod_data/"+filelist+" wvr1_data/"+trial+'/'+filelist)

if not os.path.exists(month+'_elnodpwv.txt'):
    dict_red=x_am.plot_full_folder(folder, picklefn=month+'_elnodpwv.txt', spline=2)



f = open(month+'_elnodpwv.txt','rb')
dict_red = pk.load(f)
f.close()

print(dict_red)

dict_red=dict_red[1]

day_red=dict_red['day']
time_red=dict_red['time']
pwv_red=np.array(dict_red['pwv'])
pwv_err_red=np.array(dict_red['pwv_err'])
print(pwv_red)

for i in range(len(day_red)):
    str_date=day_red[i]
    str_time=time_red[i]

    YY=int(str_date[:4])
    MM=int(str_date[4:6])
    dd=int(str_date[6:8])

    H=int(str_time[:2])
    m=int(str_time[2:4])
    s=int(str_time[4:8])

    date_xaxis_red.append(datetime.datetime(YY, MM, dd, H, m, s, 0))



fig, ax = pl.subplots()
pl.errorbar(date_xaxis_red, pwv_red, yerr=pwv_err_red, fmt='.', markersize=3, c='r', ecolor='k')
#pl.errorbar(bad_pt_x, bad_pt_y, yerr=pwv_err_red[index], c='y', fmt='.', label='fn='+time_red[index])
pl.title('PWV variations for Apr2020')
pl.xlabel('date')
pl.ylabel('pwv[um]')
pl.legend()
pl.text(0.10, 0.90, 'pwv_median[um]='+str(round(np.nanmedian(pwv_red),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
pl.text(0.10, 0.80, 'pwv_std[um]='+str(round(np.nanstd(pwv_red),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
pl.setp(ax.xaxis.get_majorticklabels(), rotation=45)

pl.show()
