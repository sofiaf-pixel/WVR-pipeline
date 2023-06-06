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

#test_file='20200508_180135_scanAz_fast.txt'
#test_file='20200612_150134_scanAz_fast.txt'
#wvr_scan = '20200723_190134_scanAz_fast.txt'
#test_file='20200808_190134_scanAz_fast.txt'

#wvr_scan='20200921_150134_scanAz_fast.txt'
#wvr_scan='20201021_130135_scanAz_fast.txt'

wvr_scan = '20201021_130135_scanAz_fast.txt'

wvr_scan = '20170815_220105_scanAz_fast.txt'


wvr_scan='20170815_120105_scanAz_fast.txt'

template='SPole_annual_50.amc'

month = wvr_scan[:6]


#0. Unzip data
#x_emp.unzip(month)

#1. Extract Modulation Parameters

# if not os.path.exists('mod_param/tilt_angle_fullmonth'+month+'.txt'):
#     x_emp.ModAmptoTilt(month)
# else:
#     print('mod_param/tilt_angle_fullmonth'+month+'.txt already exists.\n')
#     answer = input("Do you want to overwrite it? ")
#     if answer == "y":
#         x_emp.ModAmptoTilt(month)
#     elif answer == "n":
#         print('Using existing tilt model.')
#     else:
#         print("Please enter y or n.")


#2. Extract PWv atmogram from Temperatures
path_to_pwv = 'am_datfiles_Az/SPole_annual_50/'+wvr_scan[:-4]+'/'+wvr_scan[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'
if not os.path.exists(path_to_pwv):
    #x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
    x_am.fit_w_am_Az(wvr_scan, clean_method='fit')
else:
    print(path_to_pwv+' already exists.\n')
    answer = input("Do you want to overwrite it? ")
    if answer == "y":
        # x_am.fit_w_am_Az(wvr_scan, clean_method='import_model')
        x_am.fit_w_am_Az(wvr_scan, clean_method='fit')
    elif answer == "n":
        print('Converting PWV to Trj.')
    else:
        print("Please enter y or n.")



#3. Extract Trj atmograms from PWV
pf='../../../Postings/WVR_postings/20220210_BK_WVR_correlations/plots/'
path_to_Trj='am_datfiles_Az/SPole_annual_50/'+wvr_scan[:-4]+'/spectra/Trj_pk_BAK.txt'

if not os.path.exists(path_to_Trj):
    D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=1, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
else:
    print(path_to_Trj+' already exists.\n')
    answer = input("Do you want to overwrite it? ")
    if answer == "y":
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=1, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    elif answer == "n":
        D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    else:
        print("Please enter y or n.")

wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

sys.exit()

#BK options

bk_tag = '20200723B09_dk203'

rx=270
p3_filt=False
det=1

#4. Extract basic WVR and BK quantites (for later comparison/correlation)

ets=bts.extract_ts(bk_tag, wvr_scan)

wvr_atmo = ets.wvr_struct.D_pwv
az_wvr = ets.wvr_struct.az_real #az_calib
time_ordered_az=ets.wvr_struct.az_wvr
fs = ets.wvr_struct.fs
t_wvr=ets.wvr_struct.tot
wvr_atmo_Trx = wvr_atmo_Trj[str(rx)]

tod = ets.bk_struct
t_bk=np.array(tod.std)
fs_bk = tod.fs


#5. Find good a/b pairs for the given BK tag+configuration

x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = bk_tag+'det_fpu_location_'+str(rx)+'.png')
x_b = x_pair[np.where(det_pol=='b')]
x_a = x_pair[np.where(det_pol=='a')]
y_b = y_pair[np.where(det_pol=='b')]
y_a = y_pair[np.where(det_pol=='a')]

a_det = det_a_list[det]
i_det_b=np.where(x_b==x_a[det])[0]
b_det = det_b_list[i_det_b[0]]


#6.
az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[det], p3=p3_filt, i_det_savefig=det)



#7.
save_fn=pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)
save_txt='BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(det)+'_WVR_'+wvr_scan[:-4]+'_p3_'+str(p3_filt)+'_pk.txt'
ets.BK_WVR_ts_corr(D_diff_bk, D_sum_bk, wvr_scan, save_fn, save_txt)
