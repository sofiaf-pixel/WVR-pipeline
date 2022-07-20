import plotElnod as pE
import AM_functions as am
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time


x_am=am.AM_functions()
x=pE.ReadElnod()


#Produces scatter plots of one channel against the other

T0=[]
T1=[]
T2=[]
T3=[]
el=[]

filelist=[]
filelistcut=[]
for file in os.listdir('wvr1_data/BAElnod_data/'):
    if file[-6:]=='skyDip':
        if file[:8]=='20200410':
            print(file)
            filepath='BAElnod_data/'+file+'/'+file+'_fast.txt'
            print(filepath)
            try:
                filelist.append(filepath)
                filelistcut.append(file[:15])
                FH_fast=x.read_elnod_fast(filepath)
                T0_file=FH_fast[:,1]
                T1_file=FH_fast[:,2]
                T2_file=FH_fast[:,3]
                T3_file=FH_fast[:,4]
                el_file=FH_fast[:,0]

                el_file=el_file[np.where(el_file>89.)]
                T0_file=T0_file[np.where(el_file>89.)]
                T1_file=T1_file[np.where(el_file>89.)]
                T2_file=T2_file[np.where(el_file>89.)]
                T3_file=T3_file[np.where(el_file>89.)]

                T0.append(T0_file)
                T1.append(T1_file)
                T2.append(T2_file)
                T3.append(T3_file)
                el.append(el_file)

            except:
                print('File '+filepath+' does not exist.')

T0 = np.concatenate(T0, axis=0)
T1 = np.concatenate(T1, axis=0)
T2 = np.concatenate(T2, axis=0)
T3 = np.concatenate(T3, axis=0)
el = np.concatenate(el, axis=0)

print('T0=', T0)
#print('T0[0]=', T0)

T0_smooth=np.linspace(T0.min(), T0.max(), 500)

T1_coef=np.polyfit(T0, T1, 3)
T1_fit=np.poly1d(T1_coef)
T1_new=T1_fit(T0_smooth)
print('T1_coef=', T1_coef)

T2_coef=np.polyfit(T0, T2, 3)
T2_fit=np.poly1d(T2_coef)
T2_new=T2_fit(T0_smooth)
print('T2_coef=', T2_coef)

T3_coef=np.polyfit(T0, T3, 3)
T3_fit=np.poly1d(T3_coef)
T3_new=T3_fit(T0_smooth)
print('T3_coef=', T3_coef)

print(T0)
print(np.shape(T0))

fig, axes = pl.subplots(3, 1)
axes[0].scatter(T0, T1, s=3, c='r')
axes[0].plot(T0_smooth, T1_new, c='k', label='3rd order poly fit\np0='+str(round(T1_coef[len(T1_coef)-1],3))+'\np1='+str(round(T1_coef[len(T1_coef)-2],3))+'\np2='+str(round(T1_coef[len(T1_coef)-3],3))+'\np3='+str(round(T1_coef[len(T1_coef)-4],3)))
axes[0].legend(loc='upper right')
axes[0].set_title('Channel 1')
axes[0].set_ylabel('T_Ch1[K]')

axes[1].scatter(T0, T2, s=3, c='r')
axes[1].plot(T0_smooth, T2_new, c='k', label='3rd order poly fit\np0='+str(round(T2_coef[len(T1_coef)-1],3))+'\np1='+str(round(T2_coef[len(T1_coef)-2],3))+'\np2='+str(round(T2_coef[len(T1_coef)-3],3))+'\np3='+str(round(T2_coef[len(T1_coef)-4],3)))
axes[1].legend(loc='upper right')
axes[1].set_title('Channel 2')
axes[1].set_ylabel('T_Ch2[K]')

axes[2].scatter(T0, T3, s=3, c='r')
axes[2].plot(T0_smooth, T3_new, c='k', label='3rd order poly fit\np0='+str(round(T3_coef[len(T1_coef)-1],3))+'\np1='+str(round(T3_coef[len(T1_coef)-2],3))+'\np2='+str(round(T3_coef[len(T1_coef)-3],3))+'\np3='+str(round(T3_coef[len(T1_coef)-4],3)))
axes[2].legend(loc='upper right')
axes[2].set_title('Channel 3')
axes[2].set_ylabel('T_Ch3[K]')
axes[2].set_xlabel('T_Ch0[K]')

pl.savefig('Output/ZenithT_scatterplots')

pl.suptitle('Zenith Temperatures Scatter plots\n One Day of data (2020-04-10)')
pl.show()
