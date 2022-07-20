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


x_am=am.AM_functions()
x=pA.ReadAzscan()


#test_file='20190103_150135_scanAz_fast.txt'
test_file='20200418_140135_scanAz_fast.txt'
template='SPole_annual_50.amc'

# D_clean, waz, mod_removed_data=x.read_Az_fast(test_file)
#
# f = open('az_Fourier_trial.txt','wb')
# pk.dump(waz, f)
# f.close()

f = open('az_Fourier_trial.txt','rb')
waz=pk.load(f)
f.close()

az_trial=waz[0:100]

A1=7.
A2=6.
phi1=0.4
phi2=2.9
sing_mod=A1*np.sin(az_trial+phi1)
double_mod=A2*np.cos(2.*az_trial+phi2)
f_trial=sing_mod+double_mod

phase_rad0, amp0, ffourier0, fcos0=x.TestAngFourier(az_trial, f_trial, k=0)
phase_rad1, amp1, ffourier1, fcos1=x.TestAngFourier(az_trial, f_trial, k=1)
phase_rad2, amp2, ffourier2, fcos2=x.TestAngFourier(az_trial, f_trial, k=2)

phase_rad0=round(phase_rad0,2)
phase_rad1=round(phase_rad1,2)
phase_rad2=round(phase_rad2,2)
amp0=round(amp0,2)
amp1=round(amp1,2)
amp2=round(amp2,2)


pl.plot(az_trial, f_trial, label='Data')
pl.plot(az_trial, ffourier1+ffourier2, ls='--', label='Reconstructed function\navg='+str(np.mean(f_trial)))
pl.plot(az_trial, fcos0, ls='--', label='Fourier_k0\nA0='+str(amp0)+'\nphase0[rad]='+str(phase_rad0))

pl.plot(az_trial, sing_mod, label='Data Single Mod')
pl.plot(az_trial, ffourier1, ls='--', label='Fourier_k1\nA1='+str(amp1)+'\nphase1[rad]='+str(phase_rad1))

pl.plot(az_trial, double_mod, label='Data Double Mod')
pl.plot(az_trial, ffourier2, ls='--', label='Fourier_k1\nA2='+str(amp2)+'\nphase2[rad]='+str(phase_rad2))


pl.legend()
pl.suptitle('Raw Test Data vs Data Reconstructed w Fourier')
pl.title('Original Parameters:\nA1='+str(A1)+' - phi1='+str(phi1)+'\nA2='+str(A2)+' - phi2='+str(phi2))
pl.show()
