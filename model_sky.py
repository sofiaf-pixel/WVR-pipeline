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
import scipy.signal as ss
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
import random



raz=pA.ReadAzscan()
ets=bts.extract_ts()
x_am=am.AM_functions()

#in which range I should choose x0 and y0
# x_min=np.sin(0)=0
# x_max=101*0.05+np.sin(2*pi)=5
#
# y_min=0
# y_max=1

nscans=111

atmogram=np.zeros((361,nscans))
dt=np.linspace(0,1,361)

#simple model
#2 spots, wind just on x direction


x0=3 ; y0=.5 # signal centre
x1=1.5; y1=-.3


def model_2points(t, x0, y0, x1, y1, R):
    for k in range(nscans):
        dt = t[k]
        x=(k+dt)*0.05 +np.sin(2.*np.pi*dt+R)
        y=np.cos(2.*np.pi*dt+R)
        ro=(x-x0)**2+(y-y0)**2
        r1=(x-x1)**2+(y-y1)**2
        z=np.exp(-ro*2)+np.exp(-r1*4)
        atmogram[:,k]=z
    return atmogram


#slighty more complex
#3 spots, wind direction in a arbitrary direction R

xo=3 ; yo=.5 # signal centre
x1=1.5; y1=-.7
x2=4 ; y2=0.


def model_3points_R(t, x0, y0, x1, y1, x2, y2, R):
    for k in range(nscans):
        dt = t[k]
        x=(k+dt)*0.08 +np.sin(2.*np.pi*dt+R)
        y=np.cos(2.*np.pi*dt+R)
        ro=(x-x0)**2+(y-y0)**2
        r1=(x-x1)**2+(y-y1)**2
        r2=(x-x2)**2+(y-y2)**2
        z=np.exp(-ro*2)+np.exp(-r1*4)+np.exp(-r2*6)
        atmogram[:,k]=z


#even more complex
#10 spots in random position, wind direction in a arbitrary direction R

n_bubbles=100

x0=np.zeros(n_bubbles)
y0=np.zeros(n_bubbles)
r0=np.zeros((360,n_bubbles))

for i in range (n_bubbles):
    x0[i] = 20.*random.uniform(0, 1)
    y0[i] = 20.*random.uniform(0, 1)
    print('('+str(x0[i])+', '+str(y0[i])+')')





#R=-np.pi/3.
R=2.*np.pi*random.uniform(0, 1)

for k in range(100):
    x=(k+dt)*0.08 +np.sin(2.*np.pi*dt+R)
    y=np.cos(2.*np.pi*dt+R)
    z=np.zeros(np.shape(x))
    for i in range(n_bubbles):
        r0[:,i]=(x-x0[i])**2+(y-y0[i])**2
        print('r0_i=', r0[i])
        z=z+np.exp(-r0[:,i]*2)
    atmogram[:,k]=z




#loading real data to compare to the model
#first loading just the x axis

wvr_scan = '20200418_140135_scanAz_fast.txt'
path_to_test='BAElnod_data/'+wvr_scan[:-9]+'/'

D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH =raz.read_Az_fast(wvr_scan, pathtofn=path_to_test, clean_mod=3, clean_method='import_model')
waz, az_pwv, fs, idx, pwv, D_pwv = raz.read_pwvatmo(wvr_scan, show=1)



TIMEWVR=np.array(FH[:,0])
az=FH[:,18]

dt_wvr=TIMEWVR[10]-TIMEWVR[9]
daz_wvr=az[10]-az[9]
samp_rate=1./dt_wvr
tot_time=TIMEWVR[len(TIMEWVR)-1]-TIMEWVR[0]

t_scans = []
t_flat = []

nscans=len(fs.s)

for k in range(nscans):
    t_scans.append(TIMEWVR[fs.s[k]:fs.e[k]]-TIMEWVR[fs.s[k]])
    t_flat.append(TIMEWVR[fs.s[k]:fs.e[k]])

t_flat = np.concatenate(t_flat)


atmo=model_2points(t_scans, x0, y0, x1, y1)

fig=pl.figure()
im=pl.imshow(atmo, aspect='auto',interpolation='nearest', origin='lower')
pl.suptitle('PWV Atmogram simulation\n')
cbar = fig.colorbar(im, extend='both')
pl.show()




x=np.zeros(len(waz))
y=np.zeros(len(waz))

for i in range (len(waz)):
    waz_i=waz[i]
    t_i=t_flat[i]
    (scannum_i, az_i) = divmod(waz_i,360)
    y[i]=(1./np.tan(math.radians(el)))*np.sin(math.radians(az_i))+ w_y*t_i
    x[i]=(1./np.tan(math.radians(el)))*np.cos(math.radians(az_i))+ w_x*t_i


pl.scatter(x, y, s=4, c=pwv)#, cmap='plasma')
pl.colorbar()
pl.suptitle('PWV atmogram projected')
pl.title(wvr_scan[:-4])
pl.show()




#now I do the real fit

#defining function to remove dipole and offset

def modulation(x, a, phi, C):#to removee atmospheric (non instrumental) dipole
    return C + (a * np.sin(np.radians(x) + phi))

def err_modulation(p_j, x, y):
    return modulation(x, p_j[0], p_j[1], p_j[2]) - y

def remove_dipole_free(theta_az, AzScan): #mod=1 for Single Mod - mod=2 for Double Mod
    model=np.zeros(len(AzScan))
    mod_removed_data=np.zeros(len(AzScan))
    p=[5, 0, np.nanmean(AzScan)]
    res = sp.least_squares(err_modulation, p, bounds=((0,-np.pi, 0), (np.inf, np.pi, np.inf)), args=(theta_az, AzScan))
    popt=res.x
    cov_x=res.jac
    p_err=np.sqrt(np.diag(cov_x))
    model = modulation(theta_az, *popt)
    mod_removed_data = AzScan - model
    return popt, p_err, mod_removed_data, model


#extracting mod and offset removed data

AzScan_matrix=[]
x_az= np.arange(0,361)

D_pwv_nodipole=np.full(np.shape(D_pwv), np.nan)

for i in range (len(D_pwv[0,:])-1):
    print('i=',i)
    y_data=D_pwv[:,i] #D_pwv is alseady clean from single and double instrumental mod
    y_data_all=np.full(np.shape(y_data), np.nan)
    y_data_finite=y_data[np.isfinite(y_data)]
    x_az_cut=x_az[np.isfinite(y_data)]
    p_best, p_err, mod_removed_data, model = remove_dipole_free(x_az_cut, y_data_finite)
    y_data_all[np.isfinite(y_data)]=mod_removed_data
    D_pwv_nodipole[:,i]=y_data_all





def err_model_2points(p, x, y):
    t=arange(0,361)
    y_interp=np.interp(t, x, y)
    model=model_2points(x, p[0], p[1], p[2], p[3], p[4])
    model=model[np.isfinite(y)]
    y=y[np.isfinite(y)]
    return model.flatten() - y.flatten()



data_to_fit=D_pwv_nodipole

p_0=[x0, y0, x1, y1, R] #[x0, y0, x1, y1]


bounds_inf=[-5., -1., -5., -1., -2*np.pi]
bounds_sup=[5., 1., 5., 1., 2*np.pi]

res = sp.least_squares(err_model_2points, p_0, bounds=(bounds_inf, bounds_sup), args=(t_scans, D_pwv_nodipole))

print('res=', res)
p_best=res.x
print('p_best=', p_best)
cov_x=res.jac
p_err=np.sqrt(np.diag(cov_x))

amp_out=p_best[0]






def err_model_3points(p, x, y):
    model=model_3points_R(x, p[0], p[1], p[2], p[3], p[4], p[5], p[6])
    model=model[np.isfinite(y)]
    y=y[np.isfinite(y)]
    return model.flatten() - y.flatten()



data_to_fit=D_pwv_nodipole

p_0=[x0, y0, x1, y1, x2, y2, -np.pi/3.] #[x0, y0, x1, y1]


bounds_inf=[-5., -1., -5., -1., -5., -1., -np.pi]
bounds_sup=[5., 1., 5., 1., 5., 1., np.pi]

res = sp.least_squares(err_model_3points, p_0, bounds=(bounds_inf, bounds_sup), args=(t_scans, D_pwv_nodipole))

print('res=', res)
p_best=res.x
print('p_best=', p_best)
cov_x=res.jac
p_err=np.sqrt(np.diag(cov_x))

amp_out=p_best[0]
