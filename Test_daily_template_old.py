import plotElnod as pE
import plotAzscan as pA
import AM_functions as am
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


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()

#x_am.create_template('2020-04-18', '14')
#x_am.create_template('2019-01-03', '15')
#x_am.create_template('2020-04-18', '19')

data_folder='wvr1_data_local/'

try:
    year=str(sys.argv[1])
except:
    print('Error: Please specify date.\nyyyy mm dd\nyyyy mm\nyyyy')
    sys.exit()
    
try:
    month=str(sys.argv[2])
except:
    month='None'
    print('Month not specified.')
try:
    day=str(sys.argv[3])
except:
    day='None'
    print('Day not specified.')
    
file_list=[]



print(year, month, day)

#year='2018'
#month='10'

if day != 'None':
    for flist in os.listdir(data_folder):
        if flist[:4]==year:
            if flist[4:6]==month:
                print(flist)
                if flist[6:8]==day:
                    if flist[-15:]=='skyDip_fast.txt':
                        print(flist)
                        file_list.append(flist)

elif month != 'None':
    for flist in os.listdir(data_folder):
        if flist[:4]==year:
            if flist[4:6]==month:
                if flist[-15:]=='skyDip_fast.txt':
                    print(flist)
                    file_list_mm.append(flist)
            
else:
    for flist in os.listdir(data_folder):
        if flist[:4]==year:
            if flist[-15:]=='skyDip_fast.txt':
                print(flist)
                file_list.append(flist)

                
            
print('file_list=', file_list)
            
for test_file in file_list:
            
    #date='20200418'
    #time='14'
    #time='19'
    
    date=test_file[:8]
    time=test_file[9:11]

    print('date=', date)
    print('time=', time)
    
    #test_data=date+'_'+time+'0002_skyDip_fast.txt'
    #path_to_test=date+'_'+time+'0002_skyDip/'

    test_data=test_file
    path_to_test=''

    print('Creating template file.')
    x_am.create_template(date[:4]+'-'+date[4:6]+'-'+date[6:8], time) 
    
    #temp_file='SPole_annual_50.amc'
    temp_file='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_'+time+'.amc'
    print('using temp_file=', temp_file)

    print('starting fit.')
    #t1=perf_counter()
    #x_am.create_am_datfile(test_data, path_to_data=path_to_test, template= temp_file, spline=2, showplots=0) # not needed-- it is done authomatically inside fit_w_am
    x_am.fit_w_am(test_data, path_to_data=path_to_test, template=temp_file, spline=2)
    #t2=perf_counter()

    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='nscale', template= temp_file, spline=2)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', template= temp_file, spline=2)
