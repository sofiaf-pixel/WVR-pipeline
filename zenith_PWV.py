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
import extract_mod_par as emp
import Read_BICEP_ts as bts
from dateutil import parser


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()
x_emp = emp.mod_parameters()

#WVR options

wvr_scan='20200503_170002_skyDip_fast.txt'



class T_zenith(object):

    def __init__(self, unit=None, verb=True):

        '''


        '''

    def dT_del55(self, wvr_scan, template='SPole_annual_50.amc',show_im=0):

        month = wvr_scan[:6]
        # 
        # if not os.path.exists('wvr1_data/'+wvr_scan[:-9]):
        #     x_emp.unzip(month)
        #
        # if not os.path.exists('wvr1_data/'+wvr_scan[:-9]+'/'+wvr_scan):
        #     x_emp.unzip(wvr_scan[:-9])

        #2. Extract PWv atmogram from Temperatures
        path='am_datfiles/'+template[:-4]+'/'+wvr_scan[:-4]
        fn=pickle_fn=path+'/'+wvr_scan[:-4]+'_fitoutput_corr.txt'

        # if not os.path.exists(fn):
        #     #x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
        #     data=x_am.fit_w_am(wvr_scan)
        # else:
        #     print(fn+' already exists.\n')
        #     answer = input("Do you want to overwrite it? ")
        #     if answer == "y":
        #         #x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
        #         data=x_am.fit_w_am(wvr_scan)
        #     elif answer == "n":
        #         print('Converting PWV to Trj.')
        #     else:
        #         print("Please enter y or n.")
        #
        #
        # #
        #
        # f = open(fn,'rb')
        # data=pk.load(f)
        # f.close()
        #
        # pl.plot(data['El'], data['pwv_tropo'], label='pwv_tropo')
        # pl.plot(data['El'], data['pwv_total'], label='pwv_total')
        # pl.plot(data['El'], data['pwv_los_total'], label='pwv_los_total')
        # pl.legend()
        # pl.show()
        # print(data.keys())

        print(path+'Trj_el_'+wvr_scan[:-9]+'_pk.txt')

        if os.path.exists(path+'/Trj_el_'+wvr_scan[:-9]+'_pk.txt'):
            print('File exists!')
            f=open(path+'/Trj_el_'+wvr_scan[:-9]+'_pk.txt',"rb")
            Trj_dict=pk.load(f)
            f.close()

        else:
            print('File doesnt exists!')
            Trj_dict=x_am.plot_Trj_skydip(wvr_scan)

        Trj=np.array(Trj_dict['Trj'])
        el=np.array(Trj_dict['el'])

        print(np.shape(Trj)[0])

        Trj_30=[]
        Trj_40=[]
        Trj_90=[]
        Trj_150=[]
        Trj_220=[]
        Trj_270=[]

        for i in range(np.shape(Trj)[0]):
            Trj_30.append(Trj[i][0])
            Trj_40.append(Trj[i][1])
            Trj_90.append(Trj[i][2])
            Trj_150.append(Trj[i][3])
            Trj_220.append(Trj[i][4])
            Trj_270.append(Trj[i][5])

        print(np.logical_and(el>54, el<56))

        mask=np.where(np.logical_and(el>54, el<56))

        print(mask)
        Trj_150=np.array(Trj_150)
        Trj_150_55=Trj_150[mask]
        el_55=el[mask]

        T150=np.polyfit(el_55, Trj_150_55, 1)
        T150_fit=np.poly1d(T150)
        T150_f=T150_fit(el)

        # pl.scatter(Trj_dict['el'], Trj_30, s=1, label='30GHz')
        # pl.scatter(Trj_dict['el'], Trj_40, s=1, label='40GHz')
        # pl.scatter(Trj_dict['el'], Trj_90, s=1, label='90GHz')
        pl.scatter(Trj_dict['el'], Trj_150, s=1, label='T_cmb - 150GHz')
        pl.plot(el, T150_f, c='r', label='dT_150/d(el)='+str(round(T150[0],4)))
        # pl.scatter(Trj_dict['el'], Trj_220, s=1, label='220GHz')
        # pl.scatter(Trj_dict['el'], Trj_270, s=1, label='270GHz')

        pl.title(wvr_scan[:-9])

        #pl.xlim(54,56)
        pl.legend()
        pl.savefig(path+'/dT150_del_at55deg.png')
        pl.savefig('plots/'+wvr_scan[:-9]+'dT150_del_at55deg.png')
        if show_im==1:
            pl.show()
        else:
            pl.close()

        return T150[0]
