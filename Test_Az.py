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



#test_file='20190103_150135_scanAz_fast.txt'
#test_file='20190103_150135'
test_file='20200418_140135_scanAz_fast.txt'
#test_file='20200420_140134_scanAz_fast.txt'
# #test_file='20190101_010135_scanAz_fast.txt'
path_to_test='BAElnod_data/'+test_file[:-9]+'/'

template='SPole_annual_50.amc'

new_file_loc='/Volumes/Data Analysis/BICEP/Postings/WVR_postings/20210406_PWVatmogram_after_tilt_correction/plots/'

test_file_list=['20200418_140135_scanAz_fast.txt', '20200405_090135_scanAz_fast.txt', '20200418_150134_scanAz_fast.txt', ]

# for filename in test_file_list:
#     for clean_method in clean_method_list:
#         path='am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]
#         old_file_loc=path+'/'+filename[:-4]+'_pwvatmo_clean_mod3_clean_method_'+str(clean_method)+'.png'
#         os.system("cp "+old_file_loc+" "+new_file_loc)
#
# sys.exit()


# D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single=x.read_Az_fast(test_file, pathtofn=path_to_test, clean_method='import_model', show_mod_plots=1)
# D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single=x.read_Az_fast(test_file, pathtofn=path_to_test, clean_method='fit', show_mod_plots=1)
# sys.exit()
#


# pickle_fn_temps=path+'/'+filename[:-4]+'_clean_mod1_pickle_temps.txt'
#
# f = open(pickle_fn_temps,'rb')
# pickle_Temps = pk.load(f)
# f.close()
#
# waz=pickle_Temps['wAz']
# scanN=pickle_Temps['scanN']
#
# D_pwv=np.full(np.shape(pickle_Temps['T0']), np.nan)
#
# az,fs =  x.findScans(waz)
#
# pickle_fn_pwv_mod3=path+'/'+filename[:-4]+'_clean_mod1_fitoutput_full.txt'
#
# f = open(pickle_fn_pwv_mod3,'rb')
# fit_output_mod3 = pk.load(f)
# f.close()
#
# pwv_ts_mod3=fit_output_mod3['pwv']
#
#
# D_pwv = x.interpToImage(az, pwv_ts_mod3, fs)
# print('Az=', az)
# print('fs=', fs)
# print('pwv=', pwv_ts_mod3)
#
# fig=pl.figure()
# im=pl.imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
# pl.suptitle('PWV Atmogram\n'+filename[:-4])
# cbar = fig.colorbar(im, extend='both')
# cbar.set_label('PWV[um]')
# #pl.savefig(path+'/'+filename[:-4]+'_pwvatmo'+'_clean_mod3.png')
# #plt.close()
# pl.show()
#








# tilt_par=(1.12, -138.) #tilt,phi
# x_am.fit_w_am_Az(test_file, tilt_correction_par=tilt_par)
#

#
# f2 = open('tilt_angle.txt','wb')
# pk.dump(tilt_deg, f2)
# f2.close()
#
# fig = pl.figure(figsize=(16,18))
#
# pl.scatter(date_xaxis, tilt_deg['tilt0'], c='k', s=3)
# pl.plot(date_xaxis, tilt_deg['tilt0'], c='r', alpha=0.8)
# pl.axhline(y=np.mean(tilt_deg['tilt0']), color='r', linestyle='--', label='Ch0_avg[deg]='+str(round(np.mean(A0_single/dTdEl_0),2)), alpha=0.5)
# pl.scatter(date_xaxis, tilt_deg['tilt1'], c='k', s=3)
# pl.plot(date_xaxis, tilt_deg['tilt1'], c='b', alpha=0.8)
# pl.axhline(y=np.mean(tilt_deg['tilt1']), color='b', linestyle='--', label='Ch1_avg[deg]='+str(round(np.mean(A1_single/dTdEl_1),2)), alpha=0.5)
# pl.scatter(date_xaxis, tilt_deg['tilt2'], c='k', s=3)
# pl.plot(date_xaxis, tilt_deg['tilt2'], c='g', alpha=0.8)
# pl.axhline(y=np.mean(tilt_deg['tilt2']), color='g', linestyle='--', label='Ch2_avg[deg]='+str(round(np.mean(A2_single/dTdEl_2),2)), alpha=0.5)
# pl.scatter(date_xaxis, tilt_deg['tilt3'], c='k', s=3)
# pl.plot(date_xaxis, tilt_deg['tilt3'], c='y', alpha=0.8)
# pl.axhline(y=np.mean(tilt_deg['tilt3']), color='y', linestyle='--', label='Ch3_avg[deg]='+str(round(np.mean(A3_single/dTdEl_3),2)), alpha=0.5)
# pl.legend(loc='upper right')
#
# myFmt = mdates.DateFormatter('%H:%M')
# pl.gca().xaxis.set_major_formatter(myFmt)
# pl.xlim(np.min(date_xaxis)-datetime.timedelta(hours=1), np.max(date_xaxis)+datetime.timedelta(hours=1))
# pl.ylabel('A_tilt[deg]')
# pl.title('Single Modulation Amplitude converted into Tilt Angle\n'+day_str)
# pl.savefig(outpath_day+'/SingMod_tilt.png')
# pl.show()
# #   pl.close()
#
#
# pl.scatter(day_str_list, tilt_deg['tilt0'], c='k', s=3)
# pl.plot(day_str_list, tilt_deg['tilt0'], c='r', alpha=0.5, label='ch0')
# pl.scatter(day_str_list, tilt_deg['tilt1'], c='k', s=3)
# pl.plot(day_str_list, tilt_deg['tilt1'], c='b', alpha=0.5, label='ch1')
# pl.scatter(day_str_list, tilt_deg['tilt2'], c='k', s=3)
# pl.plot(day_str_list, tilt_deg['tilt2'], c='g', alpha=0.5, label='ch2')
# pl.scatter(day_str_list, tilt_deg['tilt3'], c='k', s=3)
# pl.plot(day_str_list, tilt_deg['tilt3'], c='y', alpha=0.5, label='ch3')
# pl.ylabel('tilt[deg]')
# pl.legend()
# pl.suptitle('Tilt Angle from Single Modulation Amp')
# pl.savefig(outpath+'/SingMod_tilt_global.png')
# pl.show()
# #pl.close()

# sys.exit()

#test_file='20200207_010134_scanAz_fast.txt'

#D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single=x.read_Az_fast(test_file)
#p=[a1, a2, a3, a4, phi, C1, C2, C3, C4]

#
# f = open('D_clean_trial.txt','wb')
# pk.dump(D_clean, f)
# f.close()
#
# print(D_clean[0])

#pickle_Temps=x_am.create_am_datfile_Az(test_file, pathtofn=path_to_test, clean_method='import_model')
#
# path='am_datfiles_Az/'+template[:-4]+'/'+test_file[:-4]
# pickle_fn=path+'/'+test_file[:-4]+'_pickle_temps.txt'
#
# f = open(pickle_fn,'rb')
# pickle_Temps = pk.load(f)
# f.close()
#
# print('waz=', pickle_Temps['wAz'])
# print('scanN=', pickle_Temps['scanN'])
#
# pl.plot(pickle_Temps['wAz'], pickle_Temps['scanN'])
# pl.show()
#
# pl.imshow(pickle_Temps['T0'])
# pl.show()

#test_file='20190103_150135_scanAz_fast.txt'
#test_file='20200418_140135_scanAz_fast.txt'
#test_file='20200418_150134_scanAz_fast.txt'
#test_file='20200418_190135_scanAz_fast.txt'
#test_file='20200420_140134_scanAz_fast.txt'
# #test_file='20190101_010135_scanAz_fast.txt'
#test_file='20200405_090135_scanAz_fast.txt'
test_file='20200729_010134_scanAz_fast.txt' #strong scan sync
test_file='20200725_210135_scanAz_fast.txt'
path_to_test='BAElnod_data/'+test_file[:-9]+'/'

#Test_file_list_posting=['20200405_090135_scanAz_fast.txt', '20200418_140135_scanAz_fast.txt', '20200418_150134_scanAz_fast.txt', '20200418_190135_scanAz_fast.txt']
Test_file_list_posting=['20200405_090135_scanAz_fast.txt']


template='SPole_annual_50.amc'

#x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=3, clean_method='fit', show_mod_plots=1)
#x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=3, clean_method='import_model', show_mod_plots=1)

#pickle_Temps=x_am.create_am_datfile_Az(test_file, pathtofn=path_to_test, clean_method='import_model')

#sys.exit()
for test_file in Test_file_list_posting:

    path_to_pickle='am_datfiles_Az/'+template[:-4]+'/'+test_file[:-4]
    pickle_fn_temps_impmod=path_to_pickle+'/'+test_file[:-4]+'_clean_mod3_clean_method_import_model_pickle_temps.txt'
    pickle_fn_temps_fit=path_to_pickle+'/'+test_file[:-4]+'_clean_mod3_clean_method_fit_pickle_temps.txt'

    #to re-extract the clean atmogram
    os.system('rm '+pickle_fn_temps_impmod)
    os.system('rm '+pickle_fn_temps_fit)

    #x_am.fit_w_am_Az(test_file, clean_method='import_model')
    x_am.fit_w_am_Az(test_file, clean_method='fit')

#x_am.fit_w_am_Az(test_file, tilt_correction=0.)
