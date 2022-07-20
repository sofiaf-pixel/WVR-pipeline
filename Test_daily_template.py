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
from zipfile import ZipFile
import tarfile


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()

#x_am.create_template('2020-04-18', '14')
#x_am.create_template('2019-01-03', '15')
#x_am.create_template('2020-04-18', '19')

#data_folder='wvr1_data/'
data_folder='../../wvr1_data_local/'
posting_folder='../../Postings/WVR_postings/20210625_SkyDips_Template_Comparison/plots/'

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
                if flist[6:8]==day:
                    print(data_folder+flist)
                    if flist[-7:]=='.tar.gz':
                        tar = tarfile.open(data_folder+flist, "r:gz")
                        tar.extractall(data_folder_tmp)
                        tar.close()
    for flist in os.listdir(data_folder_tmp):
        if flist[-15:]=='skyDip_fast.txt':
            print(flist)
            if flist[:4]==year:
                if flist[4:6]==month:
                    if flist[6:8]==day:
                        fn=data_folder_tmp+flist
                        file_list.append(flist)

elif month != 'None':
    for flist in os.listdir(data_folder):
        if flist[:4]==year:
            if flist[4:6]==month:
                print(flist[-7:])
                print(data_folder+flist)
                if flist[-7:]=='.tar.gz':
                    print(len(flist))
                    if len(flist)==15:
                        try:
                            tar = tarfile.open(data_folder+flist, "r:gz")
                            tar.extractall()
                            tar.close()
                            os.system(f'rm '+data_folder+flist)
                        except:
                            print(data_folder+flist+' failed.')
                    else:
                        try:
                            tar = tarfile.open(data_folder+flist, "r:gz")
                            tar.extractall(data_folder_tmp)
                            tar.close()
                            #os.system(f'rm '+data_folder+flist)
                        except:
                            print(data_folder+flist+' failed.')
    for flist in os.listdir(data_folder_tmp):
        if flist[-15:]=='skyDip_fast.txt':
            if flist[:4]==year:
                if flist[4:6]==month:
                    fn=data_folder_tmp+flist
                    file_list.append(flist)

else:
    for flist in os.listdir(data_folder):
        if flist[:4]==year:
            print(data_folder+flist)
            if flist[-7:]=='.tar.gz':
                tar = tarfile.open(data_folder+flist, "r:gz")
                tar.extractall(data_folder_tmp)
                tar.close()
    for flist in os.listdir(data_folder_tmp):
        if flist[-15:]=='skyDip_fast.txt':
            if flist[:4]==year:
                fn=data_folder_tmp+flist
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
    path_to_test=data_folder_tmp

    print('Creating template file.')
    temp_file='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_'+time+'.amc'
    if not os.path.exists('Templates/'+temp_file):
        x_am.create_template(date[:4]+'-'+date[4:6]+'-'+date[6:8], time)

    standard_temp='SPole_annual_50.amc'
    print('using temp_file=', temp_file)

    print('starting fit.')
    t1=perf_counter()
    #x_am.create_am_datfile(test_data, path_to_data=path_to_test, template= temp_file, spline=2, showplots=0) # not needed-- it is done authomatically inside fit_w_am
    x_am.fit_w_am(test_data, path_to_data=path_to_test, template=standard_temp, spline=2)
    x_am.fit_w_am(test_data, path_to_data=path_to_test, template=temp_file, spline=2)
    t2=perf_counter()

    print('time to fit one SkyDip (2 times) =', t2-t1)

    print('Template:', temp_file)

    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='nscale', pwv_layer='total', template= temp_file, spline=2, pf=posting_folder, pf2='paper_plots/')

    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='total', template= temp_file, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('total:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='strato', template= temp_file, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('strato:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='meso', template= temp_file, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('meso:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='tropo', template= temp_file, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('tropo:', el_list, pwv_list)




    print('Template:', standard_temp)

    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='nscale', pwv_layer='total', template= standard_temp, spline=2, pf=posting_folder, pf2='paper_plots/')

    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='total', template= standard_temp, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('total:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='strato', template= standard_temp, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('strato:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='meso', template= standard_temp, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('meso:', el_list, pwv_list)
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='tropo', template= standard_temp, spline=2, pf=posting_folder, pf2='paper_plots/')
    print('tropo:', el_list, pwv_list)
