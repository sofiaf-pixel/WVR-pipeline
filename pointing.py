import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.optimize as opt
import scipy.integrate as integrate
from scipy.integrate import simps
import numpy as np
import pickle as pk
import pylab as pl
import scipy.optimize as sp
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec
from itertools import groupby
from operator import itemgetter
from math import atan, atan2
import math
import Read_BICEP_ts as bts
import plotAzscan as pA
import AM_functions as am
import plotElnod as pE
import pvlib
import pandas
import datetime
from astropy.modeling import models, fitting, functional_models
from scipy.interpolate import griddata

import lmfit
from lmfit.lineshapes import gaussian2d, lorentzian

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


x_am=am.AM_functions()
raz=pA.ReadAzscan()
x_el=pE.ReadElnod()
ets=bts.extract_ts()


class pointing(object):

    def __init__(self, unit=None, verb=True):

        '''


        '''
        self.lat = -89.59464#deg
        self.lon = -44.39#0 #deg
        self.alt = 2843 #mself.lat = -90 #deg


    def read_pointing_model(self, file_list):
        azoffs_list=[]
        eloffs_list=[]
        fn_list=[]

        for fn in file_list:

            f = open('Pointing_Data_'+fn[:4]+'.txt','rb')
            pdata=pk.load(f)
            f.close()

            fn_list.append(fn[:-24])

            azoffs_list.append(pdata['Az_offs'])
            eloffs_list.append(pdata['El_offs'])
            #
            # print('fn=', fn)
            # print('pdata=', pdata)


        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.scatter(fn_list,azoffs_list, s=8, c='blue')
        ax1.plot(fn_list,azoffs_list, c='k', alpha=0.4)
        ax1.axhline(y=np.mean(azoffs_list), ls='--', c='y', label='Az_offs = '+str(round(np.mean(azoffs_list),1)))
        ax1.set_title('Az')
        ax1.set_ylabel('Az offs[deg]')
        ax1.legend()

        ax2.scatter(fn_list,eloffs_list, s=8, c='blue')
        ax2.plot(fn_list,eloffs_list, c='k', alpha=0.4)
        ax2.axhline(y=np.mean(eloffs_list), ls='--', c='y', label='El_offs = '+str(round(np.mean(eloffs_list),1)))
        ax2.set_ylabel('El offs[deg]')
        ax2.set_title('El')
        ax2.legend()
        #pl.xticks(rotation=45)
        pl.suptitle('Pointing Offset')
        pl.close()

        return fn_list, azoffs_list, eloffs_list


    def make_beam_map(self, test_file, pf, sampling='slow', scan_speed='nominal', showplots=0):

        # D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH = x.read_Az_fast(test_file_fast, pathtofn=path_to_test, clean_mod=3, show_mod_plots=1)
        # D, waz_1, calib_data_Gavg, FH = x.read_Az_fast(test_file_fast, pathtofn=path_to_test, clean_mod=0, show_mod_plots=1)
        #

        #
        # path='wvr1_data/pointing_data/'
        # day='20181209'

        # for folderlist in os.listdir(path):
        #     print('folderlist=', folderlist)
        #     if ((folderlist[-7:]=='beamMap') & (folderlist[:8]==day)):
        #         for filelist in os.listdir(f'{path+folderlist}'):
        #             print('filelist=', filelist)
        #             if filelist[-8:]=='slow.txt':
        #                 #try:
        #test_file_slow=filelist

        path_to_test='pointing_data/'+test_file[:-9]+'/'
        # test_file_slow=test_file[:-9]+'_slow.txt'
        # test_file_fast=test_file[:-9]+'_fast.txt'

        if sampling=='slow':
            time, data, D, waz, fs= raz.read_Az_slow(test_file, pathtofn=path_to_test)
        #print('time_slow=', time)

        if sampling=='fast':
            time, D, waz, calib_data_Gavg, FH= raz.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=0) #calib_data_Gavg[:,1+ch] contains ts for ch                                                                                               #calib_data_Gavg[:,0] contains time (?)
            data=calib_data_Gavg

        t_obj=[]

        for string in time:
            yy = int(string[:4])
            mm = int(string[5:7])
            dd = int(string[8:10])
            h = int(string[11:13])
            m = int(string[14:16])
            s = int(string[17:19])
            us = int(string[20:])
            t_obj.append(datetime.datetime(yy, mm, dd, h, m, s))


        #time=FH[:,0]
        # print('time_slow=', time)

        print('t_obj=', t_obj)

        t_panda=pandas.DatetimeIndex(t_obj)

        print(test_file)

        print('t_panda=', t_panda)

        data_frame=pvlib.solarposition.pyephem(t_panda, self.lat, self.lon, altitude=self.alt)#, pressure=101325, temperature=12, horizon='+0:00')
        #columns in data data_frame
        #apparent_elevation, elevation, apparent_azimuth, azimuth, apparent_zenith, zenith.

        sun_az_app = data_frame['apparent_azimuth']
        sun_az = data_frame['azimuth']
        sun_el = data_frame['elevation']
        #
        # pl.plot(sun_el, sun_az, label='sun_az')
        # pl.plot(sun_el, sun_az_app, label='sun_az_app')
        # pl.legend()
        # pl.show()

        # print(sun_az_app/sun_az)
        # pl.plot(sun_az)
        # pl.plot(sun_az_app)
        # pl.show()


        # fig = plt.figure()
        # ax1 = fig.add_subplot(111)
        # #for i in range(4):
        # i=0
        # #pl.scatter(FH[:,17],FH[:,1+4*i], label=i)
        # ax1.scatter(data[:,5]-sun_el, data[:,i+1], s=6, marker='+', label=i)
        # #pl.xlim(0,360)
        # #ax2 = ax1.twiny()
        # #for i in range(4):
        # #pl.scatter(FH[:,17],FH[:,1+4*i], label=i)
        # #ax2.scatter(sun_el, data[:,i+1], s=6, marker='+')
        # #pl.xlim(0,360)
        # pl.legend()
        # pl.show()


        sun_position={'time': [], 'sun_az':[], 'sun_el':[]}

        f = open('sun_position'+test_file[:-9]+'.txt','wb')
        pk.dump(sun_position, f)
        f.close()

        map_el=data[:,5]
        map_az=waz#*np.cos(map_el*np.pi/180.)

        # pl.plot(waz, map_az)
        # pl.show()
        #
        # pl.scatter(map_el, waz, s=10, c='b')
        # pl.scatter(map_el, map_az, marker='+', s=10, c='r')
        # pl.show()

        az_offs=map_az-sun_az #az_offs = wvr_az - real_az --> real_az = wvr_az - az_offs
        el_offs=map_el-sun_el

        az_offs = np.array([az_offs[i] for i in range (len(az_offs))])
        el_offs = np.array([el_offs[i] for i in range (len(el_offs))])
        az_offs=az_offs#*np.cos(map_el*np.pi/180.)

        sun_az = np.array([sun_az[i] for i in range (len(sun_az))])
        sun_el = np.array([sun_el[i] for i in range (len(sun_el))])

        e_max=np.max(map_el)
        e_min=np.min(map_el)

        z = copy(data[:,1]) #for ch0


        fig, ax = plt.subplots(1, 4, figsize=(10,8))
        for i in range (4):
            im = ax[i].scatter(map_el, map_az, s=55, c=data[:,1+i], cmap='plasma')
            #im = ax[i].scatter(map_el-sun_el, map_az-sun_az, s=55, c=data[:,1+i], cmap='plasma')
            #ax[i].set_aspect(15)
            #ax[i].set_ylabel('el_offs[deg]')
            ax[i].set_xlabel('el[deg]')
            cbar = fig.colorbar(im, extend='both', ax=ax[i], orientation="horizontal")
            cbar.set_label('T_sky[K]\nCh%s'%i)

        #ax[3].set_xlabel('az_offs[deg]')
        ax[0].set_ylabel('az[deg]')

        pl.suptitle('Sun Beam map\n\n'+'scanning speed = '+scan_speed+'\n'+test_file[:-9])
        pl.savefig(pf+'/'+test_file[:-9]+'_Sun_Beam_Map_abs_'+sampling+'.png')
        if showplots==1:
            pl.show()
        else:
            pl.close()


        az_e=[]
        az_s=[]

        az_s.append(0)
        for i in range(len(map_az)-1):
            if map_az[i+1]<map_az[i]:
                az_e.append(i)
                az_s.append(i+1)


        D_map = {0:None, 1:None, 2:None, 3:None}
        D_cut = {0:None, 1:None, 2:None, 3:None}

        nscans = len(az_s)-1

        # x_positions=[]
        # x_labels=[]
        # map_el_lab=[]
        # sun_el_lab=[]

        for j in range(4):
            x_positions=[]
            x_labels=[]
            map_el_lab=[]
            sun_el_lab=[]
            tsrc = copy(data[:,j+1])

            y = copy(tsrc)
            # for each 360-scan
            D_smooth = zeros([len(arange(np.min(az_offs), np.max(az_offs))),nscans])
            for i in range(nscans):
                s = az_s[i];e=az_e[i]
                idx=range(s,e)
                yp = np.interp(arange(np.min(az_offs), np.max(az_offs)), az_offs[idx], y[idx], right=nan, left=nan)
                x_positions.append(i)
                x_labels.append(round(el_offs[s+1],1))
                map_el_lab.append(round(map_el[s+1],1))
                sun_el_lab.append(round(sun_el[s+1],1))
                #pl.plot(az_offs[idx], y[idx])
                D_smooth[:,i]=yp
            #pl.show()
            D_map[j]=D_smooth


        D_sun = {'ch0': [], 'ch1':[], 'ch2':[], 'ch3':[]}


        el_struct={'sun_el':[], 'map_el':[], 'el_offs':[]}
        el_struct['sun_el']=sun_el_lab
        el_struct['map_el']=map_el_lab
        el_struct['el_offs']=x_labels

        D_sun['ch0']=D_map[0]
        D_sun['ch1']=D_map[1]
        D_sun['ch2']=D_map[2]
        D_sun['ch3']=D_map[3]

        f = open('beam_maps'+test_file[:-9]+'_'+scan_speed+'.txt','wb')
        pk.dump(D_sun, f)
        f.close()


        scan_i=0
        #el=np.zeros(nscans+1)
        #el_offs=np.zeros(nscans+1)
        # el=[]
        # el_offs=[]
        #
        # scan_list=[]
        # for i in range(len(map_el)-1):
        #     if map_el[i+1]>map_el[i]:
        #         # el[scan_i]=round(map_el[i],1)
        #         # el_offs[scan_i]=round(map_el[i]-sun_el[i],1)
        #         # el.append(round(map_el[i],1))
        #         # el_offs.append(round(map_el[i]-sun_el[i],1))
        #         el.append(round(map_el[i],1))
        #         el_offs.append(round(map_el[i]-sun_el[i],1))
        #         scan_list.append(scan_i)
        #         print(el[scan_i])
        #         scan_i+=1
        #


        #
        # x_positions = scan_list[:-1] # pixel count at label position -- excluding last scan
        # x_labels = el_offs[:-1]# labels you want to see -- excluding last scan

        y_int=arange(np.min(az_offs), np.max(az_offs))
        y_positions = arange(0,len(y_int)) # pixel count at label position

        #y_labels = [round(j+360.,1) for j in y_int] # labels you want to see
        y_labels = [round(j,1) for j in y_int] # labels you want to see

        az_struct={'sun_az':[], 'map_az':[], 'az_offs':[]}
        az_struct['sun_az']=sun_az
        az_struct['map_az']=map_az
        az_struct['az_offs']=y_labels


        fig, ax = plt.subplots(1, 4, figsize=(14,8))
        for i in range (4):
            map=D_map[i]

            im=ax[i].imshow(map, aspect='auto',interpolation='nearest', origin='lower')
            ax[i].set_xticks(x_positions[::3])
            ax[i].set_xticklabels(x_labels[::3], rotation=90)
            ax[i].set_xlabel('el_offs[deg]')
            ax[i].set_yticks(y_positions[::20])
            ax[i].set_yticklabels(y_labels[::20])

            #ax[i].set_xlim(-3,3)
            #ax[i].set_ylim(140,180)
            cbar = fig.colorbar(im, extend='both', ax=ax[i], orientation="horizontal")
            cbar.set_label('T_sky[K]\nCh%s'%i)
        #
        # for i in range(1,4):
        #     ax[i].set_yticks([])
        ax[0].set_ylabel('az_offs[deg]')

        pl.suptitle('Sun Beam map\n\n'+'scanning speed = '+scan_speed+'\n'+test_file[:-9])
        #pl.xticks(arange(0,17), x_labels)
        pl.savefig(pf+'/'+test_file[:-9]+'_Sun_Beam_Map_offs_'+sampling+'.png')
        if showplots==1:
            pl.show()
        else:
            pl.close()


        #zoomed in version
        az_offs_2018=296.8#-360. #from 2018 beam maps
        # print(y_labels)
        # print((np.where(np.array(y_labels)>az_offs_2018)))
        # try:
        #     az_plot_center_px=np.where(np.array(y_labels)>az_offs_2018)[0]
        #     print(az_plot_center_px)
        #     az_lim_dn=az_plot_center_px[0]-10
        #     az_lim_up=az_plot_center_px[0]+10
        # except:
        #     az_plot_center_px=np.where(np.array(y_labels)>az_offs_2018-360.)[0]
        #     print(az_plot_center_px)
        #     az_lim_dn=az_plot_center_px[0]-10
        #     az_lim_up=az_plot_center_px[0]+10



        # fig, ax = plt.subplots(1, 4, figsize=(14,8))
        # for i in range (4):
        #     data=D_map[i]
        #     #data_cut=data[140:180,:]
        #     #data_cut=data[210:250,:]
        #     #data_cut=data[250:275,:]
        #     data_cut=data[az_lim_dn:az_lim_up,:]
        #     D_cut[i]=data_cut
        #     im=ax[i].imshow(data_cut, aspect='auto',interpolation='nearest', origin='lower')
        #     # ax[i].set_xticks(x_positions)
        #     # ax[i].set_xticklabels(x_labels, rotation=90)
        #     #ax[i].set_yticklabels(y_labels[140:180])
        #     #ax[i].set_yticklabels(y_labels[210:250])
        #     #ax[i].set_yticklabels(y_labels[250:275])
        #     ax[i].set_xticks(x_positions[::3])
        #     ax[i].set_xticklabels(x_labels[::3], rotation=90)
        #     #ax[i].set_yticks(y_positions[az_lim_dn:az_lim_up])
        #     ax[i].set_yticklabels(y_labels[az_lim_dn:az_lim_up])
        #     ax[i].set_xlabel('el_offs[deg]')
        #     #ax[i].set_xlim(-3,3)
        #     #ax[i].set_ylim(140,180)
        #     cbar = fig.colorbar(im, extend='both', ax=ax[i], orientation="horizontal")
        #     cbar.set_label('T_sky[K]\nCh%s'%i)
        #
        # #ax[0].set_ylabel('az_offs[deg]')
        # ax[0].set_ylabel('az_offs[deg]')
        # pl.suptitle('Sun Beam map\n\n'+test_file[:-9])
        # pl.savefig(pf+'/'+test_file[:-9]+'_Sun_Beam_Map_offs_zoomedin_'+sampling+'.png')
        # #pl.xticks(arange(0,17), x_labels)
        # if showplots==1:
        #     pl.show()
        # else:
        #     pl.close()


        #return y_labels, x_labels, D_map, D_cut
        return az_struct, el_struct, D_map, D_cut



    def fit_mast(self, az_map_interp, el_map_interp, D_map, test_file, pf):

        path_to_test='pointing_data/'+test_file[:-9]+'/'

        # time, data, D, waz, fs= raz.read_Az_slow(test_file, pathtofn=path_to_test)
        #
        # map_el=data[:,5]
        # map_az=waz

        map_el=el_map_interp
        map_az=az_map_interp

        data=D_map[3]

        def gauss(x, *p):
            A, mu, sigma, offs = p
            return A*np.exp(-(x-mu)**2/(2.*sigma**2))+(offs)

        p0 = [40., 60., 10., 80.]

        g_coeff_list = []
        nscan_list = []
        #el_list = []


        for i in range (len(data[0,:])):
            y=np.array(data[:,i])
            x=np.array(map_az)

            x_masked=x[np.isfinite(y)]
            y_masked=y[np.isfinite(y)]

            pl.scatter(x_masked, y_masked, s=4, c='k', label='data')
            nscan_list.append(i)
            #el_list.append(az_offs_axis[i])
            #try:
            coeff, var_matrix = curve_fit(gauss, x_masked, y_masked, p0=p0, bounds=((0, -400, 0, 0), (1000, 400, 360, 300)))
            x_model=np.linspace(np.min(x_masked), np.max(x_masked), 500)
            g=gauss(x_model, *coeff)
            pl.plot(x_model, g, c='r', label='gaussian model:\nsigma_az[deg] ='+str(round(coeff[2],1)))
            pl.ylabel('T[K]')
            pl.xlabel('Az_offs[deg]')
            #pl.title('A, mu, sigma, offs = '+str(coeff))
            pl.legend()
            g_coeff_list.append(coeff)
            #except:
            #    print('scan failed')
            #    coeff=np.full(len(p0), np.nan)
            #    g_coeff_list.append(coeff)
            pl.title('El = '+str(el_axis[i]))
            pl.savefig(pf+'/'+test_file[:-9]+'_oneDGaussFit_Az_nscan'+str(i)+'_'+sampling+'.png')
            #pl.show()
            pl.close()


    def fit_sun_1D(self, az_axis, el_axis, D_map, test_file, pf):

        Coeff_Map = {'nscan':[], 'el':[], 'ch0': [], 'ch1':[], 'ch2':[], 'ch3':[]}

        #x=y_labels#[140:180]
        g_coeff_list = {0:None, 1:None, 2:None, 3:None}
        nscan_list = {0:None, 1:None, 2:None, 3:None}
        el_list = {0:None, 1:None, 2:None, 3:None}

        #for i in range (len(D_cut)):
        i=3
        #data=D_cut[i]

        data=D_map[i]
        x=az_axis
        data_clean=self.remove_bg(x, data)



        data=data_clean

        def oneD_analysis_Az(x, data):
            def gauss(x, *p):
                A, mu, sigma, offs = p
                return A*np.exp(-(x-mu)**2/(2.*sigma**2))+(offs)

            p0 = [40., np.mean(x), 10., 80.]

            g_coeff_list = []
            nscan_list = []
            #el_list = []


            for i in range (len(data[0,:])):
                y=np.array(data[:,i])
                x=np.array(x)

                x_masked=x[np.isfinite(y)]
                y_masked=y[np.isfinite(y)]


                pl.scatter(x_masked, y_masked, s=4, c='k', label='data')
                nscan_list.append(i)
                #el_list.append(az_offs_axis[i])
                #try:
                coeff, var_matrix = curve_fit(gauss, x_masked, y_masked, p0=p0, bounds=((0, -400, 0, 0), (1000, 400, 360, 300)))
                x_model=np.linspace(np.min(x_masked), np.max(x_masked), 500)
                g=gauss(x_model, *coeff)
                pl.plot(x_model, g, c='r', label='gaussian model:\nsigma_az[deg] ='+str(round(coeff[2],1)))
                pl.ylabel('T[K]')
                pl.xlabel('Az_offs[deg]')
                #pl.title('A, mu, sigma, offs = '+str(coeff))
                pl.legend()
                g_coeff_list.append(coeff)
                #except:
                #    print('scan failed')
                #    coeff=np.full(len(p0), np.nan)
                #    g_coeff_list.append(coeff)
                pl.title('El = '+str(el_axis[i]))
                pl.savefig(pf+'/'+test_file[:-9]+'_oneDGaussFit_Az_nscan'+str(i)+'_'+sampling+'.png')
                pl.close()
                #pl.close()

            return g_coeff_list, nscan_list, el_axis




        def oneD_analysis_El(el_x, Az_amp):
            def gauss(x, *p):
                A, mu, sigma, offs = p
                return A*np.exp(-(x-mu)**2/(2.*sigma**2))+(offs)

            p0 = [60., 0., 2., 0.]

            coeff, var_matrix = curve_fit(gauss, el_x, Az_amp, p0=p0, bounds=((0, -400, 0, 0), (1000, 400, 360, 300)))
            el_x_model=np.linspace(np.min(el_x), np.max(el_x), 500)
            g=gauss(el_x_model, *coeff)
            FWHM=coeff[2]*2.*np.sqrt(np.log(2.))
            pl.figure()
            pl.plot(el_x_model, g, c='r', label='gaussian model:\nsigma_el[deg]='+str(round(coeff[2],1))+'\nFWHM_el='+str(round(FWHM,1)))
            pl.scatter(el_x, Az_amp, s=4, label='Amp')
            pl.axvline(x=coeff[1], ls='--', c='y', label='offset_el='+str(round(coeff[1],1)))
            pl.suptitle('Az scan Gaussian Amplitude')
            #pl.title('offset_el = \n'+str(coeff[1]))
            pl.legend()
            pl.xlabel('El_offs[deg]')
            pl.ylabel('Gaussian Amp[deg]')
            pl.savefig(pf+'/'+test_file[:-9]+'_oneDGaussFit_El.png')
            pl.close()
            #pl.show()
            return coeff, var_matrix



        g_coeff_list[i], nscan_list[i], el_list[i]=oneD_analysis_Az(x,data)
        p=g_coeff_list[i]
        print(g_coeff_list)


        Coeff_Map['ch0']=g_coeff_list[0]
        Coeff_Map['ch1']=g_coeff_list[1]
        Coeff_Map['ch2']=g_coeff_list[2]
        Coeff_Map['ch3']=g_coeff_list[3]

        Coeff_Map['nscan']=nscan_list[i]
        Coeff_Map['el']=el_list[i]


        A=[]
        mu=[]
        sigma=[]
        offs=[]

        for coeff in g_coeff_list[3]:
            A.append(coeff[0])
            mu.append(coeff[1])
            sigma.append(coeff[2])
            offs.append(coeff[3])


        start_good=0
        #start_good=4
        #start_good=7

        end_good=len(A)
        #end_good=len(A)-1

        el_coeff, el_cov=oneD_analysis_El(el_list[i][start_good:end_good], A[start_good:end_good])
        #el_coeff=oneD_analysis_El(el_list[i], A)
        coeff_err=np.sqrt(np.diag(el_cov))
        # print('coeff =', el_coeff)
        # print('coeff_err =', coeff_err)


        FWHM_el=el_coeff[2]*2.*np.sqrt(np.log(2.))
        sigma=np.array(sigma)
        sigma_az_idx=np.where(A==np.max(A))
        sigma_az=np.mean(sigma[start_good:end_good])#sigma[sigma_az_idx]
        FWHM_az=sigma_az*2.*np.sqrt(np.log(2.))

        print('FWHM_az = ', FWHM_az)
        print('FWHM_el = ', FWHM_el)

        print('Az_offs = ', np.mean(mu[start_good:end_good]))
        print('El_offs = ', el_coeff[1])

        f = open('G-Fit_Coeff_'+test_file[:-9]+'.txt','wb')
        pk.dump(Coeff_Map, f)
        f.close()


        point_data = {'FWHM_az':[], 'FWHM_el':[], 'Az_offs': [], 'El_offs':[]}

        point_data['FWHM_az']=FWHM_az
        point_data['FWHM_el']=FWHM_el
        point_data['Az_offs']=np.mean(mu[start_good:end_good])
        point_data['El_offs']=el_coeff[1]

        f = open('Pointing_Data_'+test_file[:4]+'.txt','wb')
        pk.dump(point_data, f)
        f.close()



        az_range=(point_data['Az_offs']-3*FWHM_az, point_data['Az_offs']+3*FWHM_az)


        # except:
        #     print(filelist+' failed.')


        pl.figure()
        pl.scatter(el_list[i][start_good:end_good], A[start_good:end_good], s=5, c='r', label='Amp_az')
        pl.plot(el_list[i][start_good:end_good], A[start_good:end_good], c='k', alpha=0.6)
        pl.xlabel('el_offs[deg]')
        pl.ylabel('A[K]')
        #pl.axhline(y=np.mean(A[start_good:end_good]), c='y', label='Amp_avg = '+str(np.mean(A[start_good:end_good])))
        pl.legend()
        pl.savefig(pf+'/'+test_file[:-9]+'_Gauss_Amp_az.png')
        pl.close()

        pl.figure()
        pl.scatter(el_list[i][start_good:end_good], mu[start_good:end_good], s=5, c='r', label='mu_az')
        pl.plot(el_list[i][start_good:end_good], mu[start_good:end_good], c='k', alpha=0.6)
        pl.axhline(y=np.mean(mu[start_good:end_good]), c='y', label='offset_az = '+str(round(np.mean(mu[start_good:end_good]),1)), ls='--')
        pl.xlabel('el_offs[deg]')
        pl.ylabel('mu[deg]')
        pl.legend()
        pl.savefig(pf+'/'+test_file[:-9]+'_Gauss_mu_az.png')
        pl.close()

        pl.figure()
        pl.scatter(el_list[i][start_good:end_good], sigma[start_good:end_good], s=5, c='r', label='sigma_az')
        pl.plot(el_list[i][start_good:end_good], sigma[start_good:end_good], c='k', alpha=0.6)
        pl.xlabel('el_offs[deg]')
        pl.ylabel('sigma[deg]')
        pl.axhline(y=np.mean(sigma[start_good:end_good]), c='y', label='sigma_avg = '+str(round(np.mean(sigma[start_good:end_good]),1))+'\nFWHM_az='+str(round(FWHM_az,1)), ls='--')
        pl.legend()
        pl.savefig(pf+'/'+test_file[:-9]+'_Gauss_sigma_az.png')
        pl.close()

        pl.figure()
        pl.scatter(el_list[i][start_good:end_good], offs[start_good:end_good], s=5, c='r', label='offs_az')
        pl.plot(el_list[i][start_good:end_good], offs[start_good:end_good], c='k', alpha=0.6)
        pl.xlabel('el_offs[deg]')
        pl.ylabel('offset[K]')
        #pl.axhline(y=np.mean(offs[start_good:end_good]), c='y', label='offs_avg = '+str(np.mean(offs[start_good:end_good])))
        pl.legend()
        pl.savefig(pf+'/'+test_file[:-9]+'_Gauss_offs_az.png')
        pl.close()

        return Coeff_Map, point_data, az_range



    def remove_bg(self, x, data):

        data_clean=np.zeros(np.shape(data))

        def remove_sin(x, y):

            def modulation_sin(x, C1, a1, phi1, a2, phi2, a3, phi3):
                def sin_mod(a, phi, mod):
                    return (a * np.sin(np.radians(mod*x) + phi))
                return C1 + sin_mod(a1, phi1, 1) + sin_mod(a2, phi2, 2)+ sin_mod(a3, phi3, 1)

            def err_modulation_sin(p, x, y):
                return modulation_sin(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6]) - y

            p=[200., 2., 0., 2., 0.,  2., 0.]
            #p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global_modulation, p, args=(theta_az, T_matrix[:,0],  T_matrix[:,1],  T_matrix[:,2],  T_matrix[:,3]), full_output=1)
            res = sp.least_squares(err_modulation_sin, p, args=(x, y))
            #print('res=', res)
            p_best=res.x
            #print('p_best=', p_best)
            cov_x=res.jac

            p_err=np.sqrt(np.diag(cov_x))

            return y-modulation_sin(x, *p_best)

        for i in range (len(data[0,:])):
            y=np.array(data[:,i])
            x=np.array(x)
            x2=arange(len(data[:,i]))

            x_masked=x[np.isfinite(y)]
            y_masked=y[np.isfinite(y)]

            y[np.isfinite(y)]=remove_sin(x_masked, y_masked)
            data_clean[:,i]=y

        return data_clean

#fit 2D gaussian

    def fit_sun_2D(self, az_dict, el_dict, D_map, test_file, pf, scan_speed='nominal', showplots=0):

        az_axis=az_dict['az_offs']
        el_axis=el_dict['el_offs']


        az_abs=az_dict['map_az']
        el_abs=el_dict['map_el']


        sampling=test_file[-8:-4]
        #i=3
        el_center_array=np.zeros(5)
        az_center_array=np.zeros(5)
        sigmaX_deg_array=np.zeros(5)
        sigmaY_deg_array=np.zeros(5)

        D_map_cube=np.zeros(shape=D_map[0].shape + (4,)) #0-4 channels
        D_map_cube_med = np.zeros(shape=D_map[0].shape + (5,)) #0-4 channels + median

        for k in range(4):
            D_map_cube[:,:,k] = D_map[k]
            D_map_cube_med[:,:,k] = D_map[k]

        D_map_cube_med[:,:,4]=np.median(D_map_cube, axis=2)

        # z=D_map[0]
        # slice=np.zeros((len(z[:,10]), 5))

        for i in range(5):
            try:
                if i==4:
                    data=np.median(D_map_cube, axis=2)
                else:
                    data=D_map[i]
                #data=D_map_cube_med[:,:,i]

                az=np.array(az_axis)
                el=np.array(el_axis)


                data_clean=self.remove_bg(az, data)
                x2=arange(len(data[:,10]))

                # az=az[195:215]
                az_peak=np.where(data_clean==np.nanmax(data_clean))
                az_peak=az_peak[0]

                az_labels=az

                try:
                    az_min_idx=az_peak[0]-5
                    az_max_idx=az_peak[0]+5
                    az=az[az_min_idx:az_max_idx]
                except:
                    az_min_idx=0
                    az_max_idx=len(az)-1
                    az=az[az_min_idx:az_max_idx]

                #az=az#[200:300]

                data_clean[~np.isfinite(data_clean)]=np.nanmean(data_clean)
                z=data_clean[az_min_idx:az_max_idx, :]
                #z=data_clean#[200:300]

                if i == 4:
                    title_string='Ch Median'
                    i_string='median'
                else:
                    title_string='Ch '+str(i)
                    i_string=str(i)

                x1=arange(np.shape(data)[1])

                fig, axs = plt.subplots(1, 2, figsize=(10,8))
                axs[0].imshow(data, aspect=1.8)
                axs[0].set_yticks(x2)
                axs[0].set_yticklabels(az_labels)
                axs[0].set_ylim(az_min_idx,az_max_idx)
                axs[0].set_ylabel('az_offs[deg]')
                axs[0].set_xticks(x1[::2])
                axs[0].set_xticklabels(el[::2])
                axs[0].set_xlabel('el_offs[deg]')
                axs[0].set_title('data')

                axs[1].imshow(data_clean, aspect=1.8)
                axs[1].set_yticks(x2)
                axs[1].set_yticklabels(az_labels)
                axs[1].set_ylim(az_min_idx,az_max_idx)
                axs[1].set_xticks(x1[::2])
                axs[1].set_xticklabels(el[::2])
                axs[1].set_xlabel('el_offs[deg]')
                axs[1].set_title('data - bg')
                pl.suptitle(test_file[:-4]+'\nBackground subtraction\n'+title_string+'\nscanning speed = '+scan_speed)
                pl.savefig(pf+'/'+test_file[:-9]+'_bg_subtraction_ch'+i_string+'_'+sampling+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()



                def twoD_GaussianScaledAmp(x_y, xo, yo, sigma_x, sigma_y, amplitude, offset):
                # """Function to fit, returns 2D gaussian function as 1D array"""
                    (x, y)=x_y
                    xo = float(xo)
                    yo = float(yo)
                    g = offset + amplitude*np.exp( - (((x-xo)**2)/(2*sigma_x**2) + ((y-yo)**2)/(2*sigma_y**2)))
                    return g.ravel()

                def getFWHM_GaussianFitScaledAmp(img, az, el):
                    # Get FWHM(x,y) of a blob by 2D gaussian fitting
                    # Parameter:
                    #     img - image as numpy array
                    # Returns:
                    #     FWHMs in pixels, along x and y axes.
                    img_shape=np.shape(img)
                    x_pix = np.linspace(0, img.shape[1], img.shape[1]) #--do the conversion to degrees here before passing it to the fit to not get quantized values
                    y_pix = np.linspace(0, img.shape[0], img.shape[0])


                    xx, yy = np.meshgrid(el, az)

                    #Parameters: xpos, ypos, sigmaX, sigmaY, amp, baseline
                    #initial_guess = (img.shape[1]/2.,img.shape[0]/2.,img.shape[1]/2.,img.shape[0]/2., 1., 0.
                    initial_guess = ((el.max()+el.min())/2.,(az.max()+az.min())/2.,(el.max()-el.min())/2.,(az.max()-az.min())/2., 1., 0.)
                    popt, pcov = opt.curve_fit(twoD_GaussianScaledAmp, (xx, yy), img.ravel(), p0=initial_guess, bounds = ((el.min()-1, az.min()-1, 0., 0., 0., 0.), (el.max()+1, az.max()+1, (el.max()-el.min())+1, (az.max()-az.min())+1, np.inf, np.inf)))

                    model=twoD_GaussianScaledAmp((xx,yy), *popt)
                    model_matrix=model.reshape(img_shape)

                    pl.imshow(model_matrix)
                    pl.close()
                    #return (FWHM_x, FWHM_y)
                    return popt, model_matrix, el, az


                #calling example: img is your image
                popt, model, el, az = getFWHM_GaussianFitScaledAmp(z, az, el)
                xcenter, ycenter, sigmaX, sigmaY, amp, offset = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]

                idx_el=(np.where(el>popt[0]))[0] #check here
                el_center_abs=el_abs[idx_el[0]]


                el_sun=el_dict['sun_el']
                el_center_sun=el_sun[idx_el[0]]


                az_step=(az[5]-az[4])*np.cos(el_center_sun*np.pi/180.)


                el_step=el[5]-el[4]

                if sampling=='slow':
                    el_step=0.5

                ar = el_step/az_step


                # sigmaX_deg=sigmaX*el_step
                # sigmaY_deg=sigmaY*az_step
                #
                # print('sigmaX_deg=', sigmaX_deg)
                # print('sigmaY_deg=', sigmaY_deg)

                sigmaY_proj=sigmaY*np.cos(el_center_sun*np.pi/180.)

                FWHM_x = np.abs(2*sigmaX*np.sqrt(2.*np.log(2)))
                FWHM_y = np.abs(2*sigmaY_proj*np.sqrt(2.*np.log(2))) #to correct for the fact that we are not looking at the horizon

                print('FWHM_el, FWHM_az=', FWHM_x, FWHM_y)

                x = np.linspace(0, len(az), z.shape[1])
                y = np.linspace(0, len(el), z.shape[0])

                if len(el)>len(x):
                    el=el[:len(x)] #check and change
                else:
                    el_new=np.zeros(len(x))
                    el_new[:len(el)]=el
                    el=el_new


                az_sun=az_dict['sun_az']
                # az_center_sun=az_sun[idx_az[0]]
                #
                # print('az_center_sun=', az_center_sun)


                x_positions=[]
                x_labels=[]
                for j in range (len(x)):
                    x_positions.append(x[j])
                    x_labels.append(el[j])


                y_positions=[]
                y_labels=[]
                for j in range (len(y)):
                    y_positions.append(y[j])
                    y_labels.append(az[j])

                az_center=ycenter
                el_center=xcenter

                fig = plt.figure(figsize=(10,10))
                ax = fig.add_subplot(111)
                ax.imshow(z, origin='bottom', cmap='plasma', extent=(np.array(x_positions).min(), np.array(x_positions).max(), np.array(y_positions).min(), np.array(y_positions).max()), aspect=ar)
                #ax.scatter([popt], [el_center])
                pl.suptitle(test_file[:-4]+'\n'+title_string)
                pl.title('Az_offs, El_offs = '+str(round(az_center,2))+', '+str(round(el_center,2)))
                ax.set_xticks(x_positions[::3])
                ax.set_xticklabels(x_labels[::3])
                ax.set_yticks(y_positions)
                ax.set_yticklabels(y_labels)
                ax.set_xlabel('el_offs[deg]')
                ax.set_ylabel('az_offs[deg]')
                levels_plot=[model.max()/8., model.max()/4., model.max()/2., model.max()-(model.max()-model.max()/2.)/2.]
                ax.contour(x, y, model, levels=levels_plot, colors='w')
                #ax.set_xlim(popt[0]-3*popt[2], popt[0]+3*popt[2])
                #ax.set_ylim(popt[1]-3*popt[3], popt[1]+3*popt[3])
                #ax.set_xlim(popt[0]-3*popt[2], popt[0]+3*popt[2])
                #ax.set_ylim(popt[1]-3*popt[3], popt[1]+3*popt[3])
                if i==4:
                    pl.savefig(pf+'/'+test_file[:-9]+'_2DGaussFit_contours_chmedian_'+sampling+'.png')
                else:
                    pl.savefig(pf+'/'+test_file[:-9]+'_2DGaussFit_contours_ch'+str(i)+'_'+sampling+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()


                fig, ax = plt.subplots(1,2, figsize=(10,8))
                ax[0].imshow(z, origin='bottom', cmap='plasma',  extent=(np.array(x_positions).min(), np.array(x_positions).max(), np.array(y_positions).min(), np.array(y_positions).max()), aspect=ar)
                #ax.scatter([popt], [el_center])
                pl.suptitle(test_file[:-4])
                ax[0].set_title('Data')
                ax[0].set_xticks(x_positions[::2])
                ax[0].set_xticklabels(x_labels[::2])
                ax[0].set_yticks(y_positions)
                ax[0].set_yticklabels(y_labels)
                ax[0].set_xlabel('el_offs[deg]')
                ax[0].set_ylabel('az_offs[deg]')
                ax[0].set_ylim(np.min(y_positions), np.max(y_positions))
                #ax[0].set_ylim(ycenter-2*sigmaY_proj, ycenter+2*sigmaY_proj)
                #ax[0].set_xlim(popt[0]-3*popt[2], popt[0]+3*popt[2])
                #ax[0].set_ylim(popt[1]-3*popt[3], popt[1]+3*popt[3])


                ax[1].imshow(model, origin='bottom', cmap='plasma', extent=(np.array(x_positions).min(), np.array(x_positions).max(), np.array(y_positions).min(), np.array(y_positions).max()), aspect=ar)
                #ax.scatter([popt], [el_center])
                pl.suptitle(test_file[:-4])
                ax[1].set_title('Model')
                ax[1].set_xticks(x_positions[::2])
                ax[1].set_xticklabels(x_labels[::2])
                ax[1].set_yticks(y_positions)
                ax[1].set_yticklabels(y_labels)
                ax[1].set_xlabel('el_offs[deg]')
                ax[1].set_ylim(np.min(y_positions), np.max(y_positions))
                #ax[1].set_ylabel('az_offs[deg]')
                #ax[1].set_xlim(popt[0]-3*popt[2], popt[0]+3*popt[2])
                #ax[1].set_ylim(popt[1]-3*popt[3], popt[1]+3*popt[3])
                pl.suptitle('Gauss Fit\n'+title_string+'\nscanning speed = '+scan_speed+'\nz_offs, El_offs = '+str(round(az_center,2))+', '+str(round(el_center,2))+'\nSigma_az, Sigma_el='+str(round(sigmaY_proj,2))+', '+str(round(sigmaX,2)))#'\nFWHM_az, FWHM_el='+str(FWHM_y)+', '+str(FWHM_x))
                pl.savefig(pf+'/'+test_file[:-9]+'_2DGaussFit_ch'+i_string+'_'+sampling+'.png')
                if showplots==1:
                    pl.show()
                else:
                    pl.close()

                el_center_array[i]=el_center
                az_center_array[i]=az_center

                sigmaX_deg_array[i]=sigmaX
                sigmaY_deg_array[i]=sigmaY_proj

            except:
               print('Channel '+str(i)+' failed.')

        return el_center_array[4], az_center_array[4], sigmaX_deg_array[4], sigmaY_deg_array[4]



    def extract_avg_par(self, figfolder, fn, showplots=0):

        sampling=fn[-8:-4]
        yr=fn[14:18]
        scan_speed=fn[-9:-4]

        f = open(fn,'rb')
        D_p=pk.load(f)
        f.close()

        x_t=[]
        x_t_full=[]
        x_t_mask=[]

        print(D_p)

        for dt in D_p['fn']:
            x_t.append(dt[6:-4])

        print(x_t)

        #dd_list=['08', '09', '10']
        dd_list=['06']

        for dd in dd_list:
            for hh in range(00, 24):
                if hh<10:
                    t_string=dd+'_0'+str(hh)
                else:
                    t_string=dd+'_'+str(hh)

                x_t_full.append(t_string)
                if t_string in x_t:
                    x_t_mask.append(True)
                else:
                    x_t_mask.append(False)


        x_t_mask=np.array(x_t_mask)

        print('x_t=', x_t)
        print('x_t_full=', x_t_full)
        print('x_t_mask=', x_t_mask)



        az_offs=D_p['az_offs']
        az_offs_corr=[]

        for az_i in az_offs:
            if az_i<0.:
                az_i_corr=az_i
                while az_i_corr<0.:
                    az_i_corr=360.+az_i_corr
                az_offs_corr.append(az_i_corr)
            elif az_i>360.:
                az_i_corr=az_i
                while az_i_corr>360.:
                    az_i_corr=az_i_corr-360
                az_offs_corr.append(az_i_corr)
            else:
                az_offs_corr.append(az_i)

        #el_offs_final=round(np.mean(D_p['el_offs']),2)

        az_offs_corr_full=np.zeros(len(x_t_full))
        el_offs_corr_full=np.zeros(len(x_t_full))

        az_offs_corr_full[x_t_mask]=az_offs_corr
        az_offs_corr_full[~x_t_mask]=np.nan

        fig, ax = plt.subplots(2,1, figsize=(12,6))
        ax[0].scatter(x_t_full, az_offs_corr_full, s=10, c='b')
        ax[0].plot(x_t_full, az_offs_corr_full, c='b', alpha=0.5)

        #masking zeroes and outliers
        az_offs_corr_finalset=az_offs_corr_full
        # mean_azoffs=np.nanmean(az_offs_corr)
        # std_azoffs=np.nanstd(az_offs_corr)
        az_offs_corr_finalset[np.where(az_offs_corr_finalset==0.)]=np.nan
        az_offs_corr_finalset[np.where(az_offs_corr_finalset<200.)]=np.nan
        # az_offs_corr_full[np.where(az_offs_corr_full<(mean_azoffs-2*std_azoffs))]=np.nan
        # az_offs_corr_full[np.where(az_offs_corr_full>(mean_azoffs+2*std_azoffs))]=np.nan
        # mean_azoffs=np.nanmean(az_offs_corr)
        # std_azoffs=np.nanstd(az_offs_corr)
        # az_offs_corr_full[np.where(az_offs_corr_full<(mean_azoffs-2*std_azoffs))]=np.nan
        # az_offs_corr_full[np.where(az_offs_corr_full>(mean_azoffs+2*std_azoffs))]=np.nan
        az_offs_final=round(np.nanmean(az_offs_corr_finalset),2)
        az_offs_std=np.nanstd(az_offs_corr_finalset)
        print('az_offs_corr =', az_offs_corr)

        ax[0].scatter(x_t_full, az_offs_corr_finalset, marker='*', s=30, c='y')
        ax[0].axhline(az_offs_final, linestyle='--', alpha=0.5, c='k', label='az_avg='+str(az_offs_final))
        ax[0].set_title('Az Offset')
        ax[0].legend()
        ax[0].set_xticks(ticks=x_t_full[::2])
        ax[0].set_xticklabels('')
        ax[0].set_ylim(az_offs_final-3*az_offs_std, az_offs_final+3*az_offs_std)


        el_offs_corr_full[x_t_mask]=D_p['el_offs']
        el_offs_corr_full[~x_t_mask]=np.nan



        ax[1].scatter(x_t_full, el_offs_corr_full, s=10, c='r')
        ax[1].plot(x_t_full, el_offs_corr_full, c='r', alpha=0.5)


        el_offs_corr_finalset=el_offs_corr_full

        #masking zeroes and outliers
        el_offs_corr_finalset[np.where(el_offs_corr_finalset<-2.)]=np.nan
        mean_eloffs=np.nanmean(el_offs_corr_finalset)
        std_eloffs=np.nanstd(el_offs_corr_finalset)
        # el_offs_corr_full[np.where(el_offs_corr_full<(mean_eloffs-2*std_eloffs))]=np.nan
        # el_offs_corr_full[np.where(el_offs_corr_full>(mean_eloffs+2*std_eloffs))]=np.nan
        # mean_eloffs=np.nanmean(el_offs_corr_full)
        # print('el_offs_corr =', el_offs_corr_full)
        # std_eloffs=np.nanstd(el_offs_corr_full)
        # el_offs_corr_full[np.where(el_offs_corr_full<(mean_eloffs-2*std_eloffs))]=np.nan
        # el_offs_corr_full[np.where(el_offs_corr_full>(mean_eloffs+2*std_eloffs))]=np.nan

        def sin_phi(x, phi, C):
            return (C+1 * np.sin(np.radians((360./24.)*x) + phi))

        x_fit=arange(len(x_t_full))

        coeff, var_matrix = curve_fit(sin_phi, x_fit[x_t_mask & np.isfinite(el_offs_corr_full)], el_offs_corr_full[x_t_mask & np.isfinite(el_offs_corr_full)], p0=[0, 0])
        x_model=np.linspace(np.min(x_fit), np.max(x_fit), 200)
        tilt_model=sin_phi(x_model, *coeff)

        el_offs_final = round(coeff[1], 2)
        el_offs_std = np.nanstd(el_offs_corr_finalset)

        ax[1].scatter(x_t_full, el_offs_corr_finalset, marker='*', s=30, c='y')
        ax[1].plot(x_model, tilt_model, c='y', alpha=0.5, label='Expected Modulation due to WVR tilt')
        ax[1].axhline(el_offs_final, linestyle='--', alpha=0.5, c='k',label='el_offs='+str(el_offs_final))
        ax[1].set_ylim(el_offs_final-2*el_offs_std, el_offs_final+2*el_offs_std)
        ax[1].set_title('El Offset')
        ax[1].legend()
        pl.xticks(ticks=x_t_full[::2], labels=x_t_full[::2], rotation=45)
        pl.suptitle('Pointing Model\n'+yr)
        pl.savefig(figfolder+'/Pointing_Model_'+yr+'_'+sampling+'_'+scan_speed+'.png')

        if showplots==1:
            pl.show()
        else:
            pl.close()

        print('sigma_x, sigma_y=', D_p['el_sigma'], D_p['az_sigma'])


        FWHM_x = np.abs(2*np.array(np.array(D_p['el_sigma']))*np.sqrt(2.*np.log(2)))
        FWHM_y = np.abs(2*np.array(np.array(D_p['az_sigma']))*np.sqrt(2.*np.log(2))) #to correct for the fact that we are not looking at the horizon

        print('FWHM_x, FWHM_y=', FWHM_x, FWHM_y)

        FWHM_x_nozeros=FWHM_x
        FWHM_x_nozeros[np.where(FWHM_x_nozeros==0.)]=np.nan

        FWHM_y_nozeros=FWHM_y
        FWHM_y_nozeros[np.where(FWHM_y_nozeros==0.)]=np.nan

        FWHM_x_mean=round(np.nanmean(FWHM_x_nozeros),2)
        FWHM_y_mean=round(np.nanmean(FWHM_y_nozeros),2)

        FWHM_x_full=np.zeros(len(x_t_full))
        FWHM_y_full=np.zeros(len(x_t_full))

        FWHM_x_full[x_t_mask]=FWHM_x_nozeros
        FWHM_x_full[~x_t_mask]=np.nan

        FWHM_y_full[x_t_mask]=FWHM_y_nozeros
        FWHM_y_full[~x_t_mask]=np.nan
        #
        # mean_FWHM_y=np.nanmean(FWHM_y_full)
        # std_FWHM_y=np.nanstd(FWHM_y_full)
        # FWHM_y_full[np.where(FWHM_y_full<(mean_FWHM_y-2*std_FWHM_y))]=np.nan
        # FWHM_y_full[np.where(FWHM_y_full>(mean_FWHM_y+2*std_FWHM_y))]=np.nan
        # mean_FWHM_y=np.nanmean(FWHM_y_full)
        # std_FWHM_y=np.nanstd(FWHM_y_full)
        # FWHM_y_full[np.where(FWHM_y_full<(mean_FWHM_y-2*std_FWHM_y))]=np.nan
        # FWHM_y_full[np.where(FWHM_y_full>(mean_FWHM_y+2*std_FWHM_y))]=np.nan



        fig, ax = plt.subplots(2,1, figsize=(12,6))
        ax[0].scatter(x_t_full, FWHM_y_full, s=10, c='b')
        ax[0].plot(x_t_full, FWHM_y_full, c='b', alpha=0.5)

        FWHM_y_finalset=FWHM_y_full

        FWHM_y_finalset[np.where(FWHM_y_finalset>5.)]=np.nan
        FWHM_y_finalset[np.where(FWHM_y_finalset<0.5)]=np.nan

        FWHM_y_mean=round(np.nanmean(FWHM_y_finalset),2)
        FWHM_y_std=np.nanstd(FWHM_y_finalset)

        ax[0].scatter(x_t_full, FWHM_y_finalset, marker='*', s=30, c='y')

        ax[0].axhline(FWHM_y_mean, linestyle='--', alpha=0.5, c='k', label='az_FWHM='+str(FWHM_y_mean))
        ax[0].set_title('Az FWHM')
        ax[0].set_xticks(ticks=x_t_full[::2])
        ax[0].set_xticklabels('')
        ax[0].legend()
        ax[0].set_ylim(FWHM_y_mean-3*FWHM_y_std, FWHM_y_mean+3*FWHM_y_std)
        ax[0].set_xticks(ticks=x_t_full[::2])

        ax[1].scatter(x_t_full, FWHM_x_full, s=10, c='r')
        ax[1].plot(x_t_full, FWHM_x_full, c='r', alpha=0.5)

        FWHM_x_finalset=FWHM_x_full

        FWHM_x_finalset[np.where(FWHM_x_finalset>5.)]=np.nan
        FWHM_x_finalset[np.where(FWHM_x_finalset<0.5)]=np.nan


        FWHM_x_mean=round(np.nanmean(FWHM_x_finalset),2)
        FWHM_x_std=np.nanstd(FWHM_x_finalset)

        ax[1].scatter(x_t_full, FWHM_x_finalset, marker='*', s=30, c='y')
        ax[1].axhline(FWHM_x_mean, linestyle='--', alpha=0.5, c='k',label='el_FWHM='+str(FWHM_x_mean))
        ax[1].set_title('El FWHM')
        ax[1].legend()
        ax[1].set_ylim(FWHM_x_mean-3*FWHM_x_std, FWHM_x_mean+3*FWHM_x_std)
        pl.xticks(ticks=x_t_full[::2], labels=x_t_full[::2], rotation=45)
        pl.suptitle('Beam Model\n'+yr)
        pl.savefig(figfolder+'/Beam_Model_'+yr+'_'+sampling+'.png')
        if showplots==1:
            pl.show()
        else:
            pl.close()

        param_dict={'az_offs': [], 'el_offs':[], 'az_FWHM':[], 'el_FWHM':[]}

        param_dict['az_offs']=az_offs_final
        param_dict['el_offs']=el_offs_final

        param_dict['az_FWHM']=FWHM_y_mean
        param_dict['el_FWHM']=FWHM_x_mean

        f = open('pointing_parameters_'+yr+'_'+sampling+'_'+scan_speed+'.txt','wb')
        pk.dump(param_dict, f)
        f.close()
