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
from scipy import stats
from scipy import signal
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))



#test_file='20190103_150135_scanAz_fast.txt'
test_file='20200418_140135_scanAz_fast.txt'
test_file='20200101_230134_scanAz_fast.txt'
#test_file='20190101_010135_scanAz_fast.txt'

#test_file='20200418_190135_scanAz_fast.txt'
path_to_test='BAElnod_data/'+test_file[:-9]+'/'


path='wvr1_data/BAElnod_data/'
month='202004'


picklefn='l_scale_output.pk'

def fit_lscale(path=path, month=month, pkfn=picklefn):

    l_scale_dict = {'filename': [], 'day':[], 'time':[], 'l_scale':[], 'l_scale_err':[]}

    if not os.path.exists(path):
        print('ERROR. Path does not exist')
    for folderlist in os.listdir(path):
        if folderlist[-6:]=='scanAz':
            print('month=', folderlist[:6])
            if folderlist[:6]==month:
                for filelist in os.listdir(f'{path+folderlist}'):
                    if filelist[-8:]=='fast.txt':
                        print('Fitting '+filelist)
                        print('day=',filelist[4:8])
                        print('time=',filelist[9:11])
                        waz, p, p_err, mod_removed_data, model, calib_data_Gavg=x.read_Az_fast(filelist, pathtofn='BAElnod_data/'+filelist[:-9]+'/', clean_mod=3, clean_method='fit_both', show_mod_plots=0)
                        l_scale_dict['filename'].append(filelist)
                        l_scale_dict['day'].append(filelist[5:9])
                        l_scale_dict['time'].append(filelist[10:12])
                        l_scale_dict['l_scale'].append(p[14])
                        l_scale_dict['l_scale_err'].append(p_err[14])
                        print('l_scale=',p[14])
                        print('l_scale_err=',p_err[14])

                        f=open(pkfn, "wb")
                        pk.dump(l_scale_dict, f)
                        f.close()


def read_lscale(pkfn=picklefn):

        f = open(pkfn,"rb")
        l_scale_dict = pk.load(f)
        f.close()

        fw=open('SP_windspeed_042020.pk', "rb")
        wind_dict = pk.load(fw)
        fw.close()

        t_wind=np.array(wind_dict['t'])
        ws=np.array(wind_dict['ws'])
        wd=np.array(wind_dict['wd'])


        day_list=l_scale_dict['day']
        time_list=l_scale_dict['time']
        filelist=l_scale_dict['filename']
        date_time_str=[fl[6:-20] for fl in filelist]

        time_list_obj=[datetime.datetime.strptime(fl[:-20], '%Y%m%d_%H') for fl in filelist]
        time_list_obj=np.array(time_list_obj)


        ws_ds=[]
        ws_t=[]
        for i_dt in time_list_obj:
            idx=np.where(t_wind == pd.Timestamp(2020, 4, i_dt.day, i_dt.hour, 0, 0))
            idx=idx[0]
            ws_ds.append(ws[idx[0]])
            ws_t.append(pd.Timestamp(2020, 4, i_dt.day, i_dt.hour, 0, 0))


        fl_max=filelist[len(filelist)-1]
        date_max=datetime.datetime.strptime(fl_max[:-20], '%Y%m%d_%H')

        t_ind_match=np.where(t_wind <= date_max)

        l_scale=np.array(l_scale_dict['l_scale'])
        l_scale_err=np.array(l_scale_dict['l_scale_err'])
        avg=np.nanmean(l_scale)
        avg_err=np.nansum([j*j for j in l_scale])/len(l_scale)
        std=np.std(l_scale)

        mask_cut=np.where((l_scale>(avg-std)) & (l_scale<(avg+std)))
        avg_cut=np.nanmean(l_scale[mask_cut])
        avg_err_cut=np.nansum([j*j for j in l_scale_err[mask_cut]])/len(l_scale_err[mask_cut])


        kernel = np.ones(5)
        kernel_gauss=gaussian(np.linspace(0,len(kernel),len(kernel)), len(kernel)/2, len(kernel)/8)
        l_scale_conv = np.convolve(l_scale, kernel_gauss, mode='same')


        ws_y_up=np.mean(ws_ds)+2*np.std(ws_ds)
        ws_y_dn=np.mean(ws_ds)-2*np.std(ws_ds)

        lscale_y_up=np.mean(l_scale)+2*np.std(l_scale)
        lscale_y_dn=np.mean(l_scale)-2*np.std(l_scale)

        ws_ds=np.array(ws_ds)
        ws_t=np.array(ws_t)

        ws_good=np.where((ws_ds<ws_y_up) & (ws_ds>ws_y_dn))
        l_scale_good=np.where((l_scale<lscale_y_up) & (l_scale>lscale_y_dn))

        print(ws_good)
        print(ws_good[0])

        fig, ax1 = pl.subplots()

        color = 'tab:blue'
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Wind Speed[m/s]', color = color)
        ax1.plot(ws_t[ws_good[0]], ws_ds[ws_good[0]], color=color, label='wind speed')
        ax1.set_ylim(ws_y_dn, ws_y_up)
        #ax1.tick_params(axis ='y', labelcolor = color)

        ax2 = ax1.twinx()

        color = 'tab:red'
        ax2.set_ylabel('sin', color = color)
        ax2.plot(ws_t[l_scale_good[0]], l_scale[l_scale_good[0]], color=color, label='l_scale')
        ax2.tick_params(axis ='y', labelcolor = color)
        ax2.set_ylabel('l_scale', color = color)
        ax2.set_ylim(lscale_y_dn, lscale_y_up)
        #fig.suptitle()

        pl.show()





        fig, (ax, ax2) = pl.subplots(2,1)
        col_arr=np.full(len(ws_t), 'r')
        col_arr[l_scale_good]='k'
        #ax.plot(l_scale_conv/len(kernel))
        ax.plot(ws_t, l_scale)
        ax.scatter(ws_t, l_scale, c='r', s=3)
        ax.axhline(y=avg_cut, color='r', linestyle='--', label='lscale_avg='+str(round(avg_cut,2))+'+-'+str(round(avg_err_cut,2)))
        ax.axhline(y=avg-std, color='r', alpha=0.4, linestyle='-', label='lower cut')
        ax.axhline(y=avg+std, color='r', alpha=0.4, linestyle='-', label='upper cut')
        #ax2=ax.twinx()
        #ax.plot(ws[t_ind_match])
        #ax.xticks(rotation=90)
        ax.xaxis.set_major_locator(MultipleLocator(20))


        pl.suptitle(month)
        ax.set_title('l_scale\n')
        #pl.xlabel()
        ax.set_ylabel('l_scale')


        ws_conv=np.convolve(ws[t_ind_match], kernel_gauss, mode='same')

        #ax2.plot(t_wind[t_ind_match], ws_conv/len(kernel))
        ax2.plot(ws_t, ws_ds, c='k', alpha=0.4)
        ax2.scatter(ws_t, ws_ds, c=col_arr, s=3)
        ax2.set_title('wind speed\n')


        pl.legend()
        #pl.savefig('l_scale_onemonth_'+month+'.png')
        pl.show()



        fig, ax= pl.subplots()
        #ax.plot(l_scale_conv/len(kernel))
        ax.plot(ws_t, l_scale)
        ax.scatter(ws_t, l_scale, c='r', s=3)
        ax.axhline(y=avg_cut, color='r', linestyle='--', label='lscale_avg='+str(round(avg_cut,2))+'+-'+str(round(avg_err_cut,2)))
        ax.axhline(y=avg-std, color='r', alpha=0.4, linestyle='-', label='lower cut')
        ax.axhline(y=avg+std, color='r', alpha=0.4, linestyle='-', label='upper cut')
        #ax2=ax.twinx()
        #ax.plot(ws[t_ind_match])
        #ax.xticks(rotation=90)
        #ax.xaxis.set_major_locator()

        pl.suptitle('l_scale\n')
        ax.set_title(month)
        ax.set_ylabel('l_scale')
        pl.legend()
        pl.savefig('l_scale_onemonth_'+month+'.png')
        pl.show()








read_lscale()
