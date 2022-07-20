import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
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
import pointing as pnt
import PWV_to_wind as ptw


#wvr_scan = '20200729_010134_scanAz_fast.txt'

scan_list = ['20200418_140135_scanAz_fast.txt',
            '20200405_090135_scanAz_fast.txt',
            '20200418_150134_scanAz_fast.txt',
            '20200418_190135_scanAz_fast.txt',
            '20200725_210135_scanAz_fast.txt',
            '20200729_010134_scanAz_fast.txt']

scan_list = ['20200418_190135_scanAz_fast.txt']

#posting folder
pf='../../Postings/WVR_postings/20210901_extracting_wind_info_2/plots'

#for atmo dipole removed atmogram
template='SPole_annual_50'

show_plots=1

for wvr_scan in scan_list:

    pk_fn='am_datfiles_Az/'+template+'/'+wvr_scan[:-4]+'/'+wvr_scan[:-4]+'_Dpwv_atmo-dipole_removed.txt'

    x_w=ptw.extract_wind(wvr_scan, pf, show_plots)

    #creating fake sky
    # sky_model = x_w.sim_sky(az_steps=100)
    # x_w.scan_sky_oneturn(sky_model)
    # sim_atmogram=x_w.scan_sky_map(sky_model, v=10*(1./50.), nscans=1000)

    #loading real sky
    D_pwv, waz, pwv = x_w.load_real_atmo(wvr_scan)


    #removing atmospheric dipole
    if not os.path.exists(pk_fn):
        dipole_removed_atmo_dict = x_w.remove_atmo_dipole(D_pwv, pk_fn)
        D_pwv_nodipole=dipole_removed_atmo_dict['mod_removed_data']
    else:
        f = open(pk_fn,'rb')
        dipole_removed_atmo_dict=pk.load(f)
        f.close()
        D_pwv_nodipole=dipole_removed_atmo_dict['mod_removed_data']
        model=dipole_removed_atmo_dict['model']
        fig, axs = pl.subplots(3, 1, figsize=(10,8))
        im0=axs[0].imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
        axs[0].set_title('PWV Atmogram \n(Instrumental dipole and quadrupole removed)')
        cbar = fig.colorbar(im0, extend='both', ax=axs[0])
        cbar.set_label('PWV[um]')
        im1=axs[1].imshow(model, aspect='auto',interpolation='nearest', origin='lower')
        axs[1].set_title('Atmospheric Dipole Model')
        cbar = fig.colorbar(im1, extend='both', ax=axs[1])
        cbar.set_label('PWV[um]')
        im2=axs[2].imshow(D_pwv_nodipole, aspect='auto',interpolation='nearest', origin='lower')
        axs[2].set_title('PWV Atmogram - Atmospheric Dipole Model')
        cbar = fig.colorbar(im2, extend='both', ax=axs[2])
        cbar.set_label('PWV[um]')
        subplots_adjust(hspace=0.25)
        pl.savefig(pf+'/'+wvr_scan[:-4]+'_atmo_data_systematics+atmo_dipole_clean.png')

        if show_plots==1:
            pl.show()
        else:
            pl.close()

    #removing offset from each az scan
    D_pwv_nodipole_nooffs=x_w.remove_offs(D_pwv_nodipole, wvr_scan)


    #extracting wind speed and direction from autocorrelation
    #ws, wd = x_w.find_ws_wd(D_pwv_nodipole_nooffs, wvr_scan, delta_az=20)

    wd_list=[]
    ws_list=[]
    daz_list=[]

    for i in range (30):
        #daz_list.append(5+2*i)
        daz_list.append(i)

    print('daz_list=', daz_list)


    for i in range(len(daz_list)):

        daz=daz_list[i]

        ws, wd = x_w.find_ws_wd(D_pwv_nodipole_nooffs, wvr_scan, delta_az=daz)

        print('wvr_scan=', wvr_scan)
        print('daz=', daz)
        print('ws=', ws)
        print('wd=', wd)

        wd_list.append(wd)
        ws_list.append(ws)

    pl.scatter(daz_list, ws_list, s=6, c='r')
    pl.plot(daz_list, ws_list, c='r', alpha=0.5)
    pl.xlabel('az_step')
    pl.ylabel('wind speed[m/s]')
    pl.title(wvr_scan[:-4])
    pl.savefig(pf+'/'+wvr_scan[:-4]+'_ws_vs_azstep.png')
    pl.show()

    pl.scatter(daz_list, wd_list, s=6, c='b')
    pl.plot(daz_list, wd_list, c='b', alpha=0.6)
    pl.xlabel('az_step')
    pl.ylabel('wind dir[deg]')
    pl.title(wvr_scan[:-4])
    pl.savefig(pf+'/'+wvr_scan[:-4]+'_wd_vs_azstep.png')
    pl.show()


    #sky=x_w.atmo_to_sky(waz, pwv, D_pwv_nodipole, ws, wd)
    #ws_fake, wd_fake = x_w.find_ws_wd(sim_atmogram, 'Simulated Sky.txt', delta_az=20) #the .txt is just to add 4 characters


    #x_w.NOAA_wind_data(wvr_scan)
