
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
raz=pA.ReadAzscan()


test_file='20200418_140135_scanAz_fast.txt'

template='SPole_annual_50.amc'

filename=test_file
path='am_datfiles_Az/'+template[:-4]+'/'+filename[:-4]

pickle_fn=path+'/'+filename[:-4]+'_clean_mod2_fitoutput__tilt_corrected.txt'
pickle_fn_mod3=path+'/'+filename[:-4]+'_clean_mod3_fitoutput.txt'

tilt_par=(1.12, -138.)
tilt_par_day=(0.895, -135.) #tilt,phi
(tilt, phi)=tilt_par

if phi<0.:
    phi=360.+phi

dTdEl=np.zeros(4)

dTdEl[0]=0.88
dTdEl[1]=0.63
dTdEl[2]=0.37
dTdEl[3]=0.27


#p=[a1, a2, a3, a4, phi, C1, C2, C3, C4]
D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single=raz.read_Az_fast(filename, dTdEl, tilt_par, clean_mod=3)

f = open(pickle_fn,'rb')
fit_output = pk.load(f)
f.close()

pwv_ts=fit_output['pwv']
az=fit_output['Az']

# El_corr=fit_output['el_correction']
# El_corr=El_corr[0:len(az)]


El_corr=(tilt*(np.sin(np.radians(np.asarray(az)+phi))))

sing_mod_0=p_single[0]*np.sin(np.radians(az)+p_single[4])
sing_mod_1=p_single[1]*np.sin(np.radians(az)+p_single[4])
sing_mod_2=p_single[2]*np.sin(np.radians(az)+p_single[4])
sing_mod_3=p_single[3]*np.sin(np.radians(az)+p_single[4])

print('Sing mod Amp and phase=', p_single[0], np.degrees(p_single[4]))



dT0=dTdEl[0]*tilt*(np.sin(np.radians(np.asarray(az)+phi)))
dT1=dTdEl[1]*tilt*(np.sin(np.radians(np.asarray(az)+phi)))
dT2=dTdEl[2]*tilt*(np.sin(np.radians(np.asarray(az)+phi)))
dT3=dTdEl[3]*tilt*(np.sin(np.radians(np.asarray(az)+phi)))

print('El Corr Amp and phase=', dTdEl[0]*tilt, phi)

ind0=np.where(dT0==np.min(dT0))
ind0=ind0[0]

ind1=np.where(sing_mod_0==np.min(sing_mod_0))
ind1=ind1[0]

print('ind0, ind1=', ind0, ind1)

print('az[ind0[0]]=', az[ind0[0]])
print('az[ind1[0]]=', az[ind1[0]])

dphi=az[ind1[0]]-az[ind0[0]]

print('dphi=', dphi)

dT0_shift=dTdEl[0]*tilt*(np.sin(np.radians(np.asarray(az)+phi-dphi)))
dT1_shift=dTdEl[1]*tilt*(np.sin(np.radians(np.asarray(az)+phi-dphi)))
dT2_shift=dTdEl[2]*tilt*(np.sin(np.radians(np.asarray(az)+phi-dphi)))
dT3_shift=dTdEl[3]*tilt*(np.sin(np.radians(np.asarray(az)+phi-dphi)))

pl.scatter(az, El_corr)
pl.show()

fig, ax = pl.subplots(4,1, sharex=True, figsize=(18,10))

ax[0].scatter(az, dT0, s=3, label='From El correction')
ax[0].scatter(az, sing_mod_0, s=3, label='From Fit')
ax[0].scatter(az, dT0_shift, s=3, label='From El correction - shifted')
ax[0].set_title('Ch0')
pl.legend()
ax[1].scatter(az, dT1, s=3, label='From El correction')
ax[1].scatter(az, sing_mod_1, s=3, label='From Fit')
ax[1].scatter(az, dT1_shift, s=3, label='From El correction - shifted')
ax[1].set_title('Ch1')

ax[2].scatter(az, dT2, s=3, label='From El correction')
ax[2].scatter(az, sing_mod_2, s=3, label='From Fit')
ax[2].scatter(az, dT2_shift, s=3, label='From El correction - shifted')
ax[2].set_title('Ch2')

ax[3].scatter(az, dT3, s=3, label='From El correction')
ax[3].scatter(az, sing_mod_3, s=3, label='From Fit')
ax[3].scatter(az, dT3_shift, s=3, label='From El correction - shifted')
ax[3].set_title('Ch3')
pl.legend()


pl.suptitle('Single Mod Models in comparison')
pl.show()
