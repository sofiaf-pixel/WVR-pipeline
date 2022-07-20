T0=[]
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

# FH_fast=x.read_elnod_fast('20190103_150002_skyDip_fast.txt')
# FH_slow=x.read_elnod_slow('20190103_150002_skyDip_slow.txt')
#
# fig, axes = pl.subplots(2, 2)
# axes[0, 0].plot(FH_fast[:,0], FH_fast[:,1], label="Fast")
# axes[0, 0].plot(FH_slow[:,0], FH_slow[:,1], label="Slow")
# axes[0, 0].legend(loc='upper right')
# axes[0, 0].set_title('Ch0')
#
# axes[0, 1].plot(FH_fast[:,0], FH_fast[:,2], label="Fast")
# axes[0, 1].plot(FH_slow[:,0], FH_slow[:,2], label="Slow")
# axes[0, 1].legend(loc='upper right')
# axes[0, 1].set_title('Ch1')
#
# axes[1, 0].plot(FH_fast[:,0], FH_fast[:,3], label="Fast")
# axes[1, 0].plot(FH_slow[:,0], FH_slow[:,3], label="Slow")
# axes[1, 0].legend(loc='upper right')
# axes[1, 0].set_title('Ch2')
#
# axes[1, 1].plot(FH_fast[:,0], FH_fast[:,4], label="Fast")
# axes[1, 1].plot(FH_slow[:,0], FH_slow[:,4], label="Slow")
# axes[1, 1].legend(loc='upper right')
# axes[1, 1].set_title('Ch3')
#
# pl.show()
#
# FH, fit_par=x.fit_elnod_vsel('20190103_150002_skyDip_fast.txt')


# x.fit_elnod_vsfreq('20190103_150002_skyDip_fast.txt', el_range=(80.,90.))

#x_am.fit_w_am('20190103_150002_skyDip_fast.txt', el_range=(0.,90.))

t1_start = perf_counter()

#filename='20190103_150002_skyDip_fast.txt'
#filename='20190101_010002_skyDip_fast.txt'
#filename='20190102_030002_skyDip_fast.txt'
#filename='20200402_020002_skyDip_fast.txt'
#temp= 'SPole_annual_50.amc'
#temp='SPole_MAM_25.amc'

#filelist=['20190103_150002_skyDip_fast.txt', '20190101_010002_skyDip_fast.txt', '20190102_030002_skyDip_fast.txt']
# filelist=[]
# filelistcut=[]
# for file in os.listdir('wvr1_data'):
#      if file[-15:]=='skyDip_fast.txt':
#          filelist.append(file)
#          filelistcut.append(file[:8])
# #
# zenith_all=np.zeros((4, len(filelist)))
# zenith_corr_all=np.zeros((4, len(filelist)))
# #fit_output=x_am.create_am_datfile(filename, path_to_data='onedayofSkyDips/20200402/', spline=1)
# # fit_output=x_am.create_am_datfile(filename, spline=1)
# x=list(range(len(filelist)))
# for i in range (len(filelist)):
#     filename=filelist[i]
#     zenith_corr_all[:, i]=x_am.create_am_datfile(filename, spline=1)
#     zenith_all[:,i]= 90. - zenith_corr_all[:, i]
#     print(zenith_all[:, i])
# #
# #
# fig, axes = pl.subplots(4, 1)
# axes[0].scatter(x, zenith_corr_all[0, :], s=3, color='k')
# axes[0].axhline(y=np.mean(zenith_corr_all[0, :]), color='r', linestyle='--', label='Zenith Correction = '+str(round(np.mean(zenith_corr_all[0, :]),2)))
# axes[0].set_title('Channel 0')
# axes[0].legend()
#
# axes[1].scatter(x, zenith_corr_all[1, :], s=3, color='k')
# axes[1].axhline(y=np.mean(zenith_corr_all[1, :]), color='r', linestyle='--', label='Zenith Correction = '+str(round(np.mean(zenith_corr_all[1, :]),2)))
# axes[1].set_title('Channel 1')
# axes[1].legend()
#
# axes[2].scatter(x, zenith_corr_all[2, :], s=3, color='k')
# axes[2].axhline(y=np.mean(zenith_corr_all[2, :]), color='r', linestyle='--', label='Zenith Correction= '+str(round(np.mean(zenith_corr_all[2, :]),2)))
# axes[2].set_title('Channel 2')
# axes[2].legend()
#
# axes[3].scatter(x, zenith_corr_all[3, :], s=3, color='k')
# axes[3].axhline(y=np.mean(zenith_corr_all[3, :]), color='r', linestyle='--', label='Zenith Correction= '+str(round(np.mean(zenith_corr_all[3, :]),2)))
# axes[3].set_title('Channel 3')
# axes[3].legend()
# pl.show()
#
#
# fig, axes = pl.subplots(4, 1)
# axes[0].scatter(x, zenith_all[0, :], s=3, color='k')
# axes[0].axhline(y=np.mean(zenith_all[0, :]), color='r', linestyle='--', label='Zenith = '+str(round(np.mean(zenith_all[0, :]),2)))
# axes[0].set_title('Channel 0')
# axes[0].legend()
#
# axes[1].scatter(x, zenith_all[1, :], s=3, color='k')
# axes[1].axhline(y=np.mean(zenith_all[1, :]), color='r', linestyle='--', label='Zenith = '+str(round(np.mean(zenith_all[1, :]),2)))
# axes[1].set_title('Channel 1')
# axes[1].legend()
#
# axes[2].scatter(x, zenith_all[2, :], s=3, color='k')
# axes[2].axhline(y=np.mean(zenith_all[2, :]), color='r', linestyle='--', label='Zenith = '+str(round(np.mean(zenith_all[2, :]),2)))
# axes[2].set_title('Channel 2')
# axes[2].legend()
#
# axes[3].scatter(x, zenith_all[3, :], s=3, color='k')
# axes[3].axhline(y=np.mean(zenith_all[3, :]), color='r', linestyle='--', label='Zenith = '+str(round(np.mean(zenith_all[3, :]),2)))
# axes[3].set_title('Channel 3')
# axes[3].legend()
#
# pl.show()
#
#
# print(filelistcut)

#
#
#
#
# x=list(range(4*len(filelist)))
# pl.scatter(x, zenith_all)
# pl.show()
#
#
# print(filelistcut)
# print('zenith', zenith_all)
#
#
#
#
#
#
#
# # el_list_red=np.array(fit_output_red['El'])
# # pwv_list_red=np.array(fit_output_red['pwv'])
# #
# # fit_output=x_am.create_am_datfile(filename, spline=0)
# # fit_output=x_am.fit_w_am(filename, template=temp, spline=0)
# # el_list=np.array(fit_output['El'])
# # pwv_list=np.array(fit_output['pwv'])
# #
# #
# #
# # #el_list, pwv_list=x_am.plot_am_fit(filename, spline=0)
# # #el_list_red, pwv_list_red=x_am.plot_am_fit(filename, spline=1)
# #
# # print(el_list, pwv_list)
# #
# # pl.scatter(el_list, pwv_list, s=5, color='black')
# # pl.scatter(el_list_red, pwv_list_red, s=15, marker='*', color='red')
# # pl.xlabel('El[deg]')
# # pl.ylabel('pwv[um]')
# # pl.text(np.max(el_list)-50.,  np.max(pwv_list)-80., 'pwv_mean[um]={:.3f}'.format(np.mean(pwv_list))+'±{:.3f}'.format(np.std(pwv_list)), fontsize=10)
# # pl.text(np.max(el_list)-50.,  np.max(pwv_list)-88., 'pwv_mean_redfit[um]={:.3f}'.format(np.mean(pwv_list_red))+'±{:.3f}'.format(np.std(pwv_list_red)), fontsize=10, color='red')
# # pl.axhline(y=np.mean(pwv_list), color='k', linestyle='--')
# # pl.axhline(y=np.mean(pwv_list_red), color='red', linestyle='-')
# # pl.suptitle(filename[:-4])
# # pl.title(template[:-4])
# # pl.savefig('fit_output_plots/'+template[:-4]+'/'+filename[:-4]+'_fitoutput_pwv_scatter_redfitcomparison.png')
# # pl.show()
# # pl.close()
#
#
# t2_end = perf_counter()
#
# #x_am.plot_am_fit('20190103_150002_skyDip_fast.txt')
# print('t1_start=', t1_start)
# print('t2_end=', t2_end)
# print("Time to fit 1 Elnod:", t2_end-t1_start)
# #print("fit_output=", fit_output_red)

#beta=[fit_par[1,0],fit_par[2,0],fit_par[3,0],fit_par[4,0]]
#print(beta)

#x_am.plot_out_file('h2o_test.out')
#x.extract_pwv(FH, beta)
# template_list=['SPole_MAM_5.amc', 'SPole_DJF_25.amc', 'SPole_DJF_50.amc', 'SPole_DJF_75.amc', 'SPole_DJF_95.amc', 'SPole_MAM_25.amc', 'SPole_MAM_50.amc', 'SPole_MAM_75.amc', 'SPole_MAM_95.amc']
# for template_file in template_list:
#     x_am.plot_full_folder('fitwdifftemplate', picklefn='trial.txt', template=template_file)


#x_am.plot_full_folder('BAElnod_data', picklefn='pwv_april.txt')

#x_am.plot_Apr_pwv('pwv_april.txt')


# filelist=[
#     '20200421_050002_skyDip_fast.txt',
#     '20200419_130003_skyDip_fast.txt',
#     '20200418_140002_skyDip_fast.txt',
#     '20200412_060003_skyDip_fast.txt',
#     '20200411_070002_skyDip_fast.txt',
#     '20200406_210002_skyDip_fast.txt'
#     ]
#
#
# pwv, pwv_err, Nscale=x_am.plot_am_fit('20200406_210002_skyDip_fast.txt')
#
#
# print(pwv)



os.system(f"rm wvr1_data/onedayofSkyDips/20200418/*calib*")
os.system(f"rm wvr1_data/onedayofSkyDips/20200418/*interp*")
os.system(f"rm wvr1_data/onedayofSkyDips/20200402/*calib*")
os.system(f"rm wvr1_data/onedayofSkyDips/20200402/*interp*")
os.system(f"rm wvr1_data/onedayofSkyDips/20190101/*calib*")
os.system(f"rm wvr1_data/onedayofSkyDips/20190101/*interp*")
os.system(f"rm wvr1_data/onedayofSkyDips/20190103/*calib*")
os.system(f"rm wvr1_data/onedayofSkyDips/20190103/*interp*")



T0=np.zeros((4,17048))
T1=np.zeros((4,17048))
T2=np.zeros((4,17048))
T3=np.zeros((4,17048))

T0_all=[]
T1_all=[]
T2_all=[]
T3_all=[]

col=['b','r','g','y']

day_name=["" for x in range(4)]
day_ind=0
j=0
for folderlist in os.listdir('wvr1_data/onedayofSkyDips/'):
    for fn in os.listdir('wvr1_data/onedayofSkyDips/'+folderlist+'/'):
        print(fn)
        data=x.read_elnod('onedayofSkyDips/'+folderlist+'/'+fn)
        for i in range (len(data[:,5])):
            if (data[i,5]>89.9):
                if data[i,2]<10:
                    print(day_ind,j)
                    print(data[i,1],data[i,2],data[i,3],data[i,4])
                print(folderlist)
                day_name[day_ind]=folderlist
                T0[day_ind,j]=data[i,1]
                T1[day_ind,j]=data[i,2]
                T2[day_ind,j]=data[i,3]
                T3[day_ind,j]=data[i,4]

                T0_all.append(data[i,1])
                T1_all.append(data[i,2])
                T2_all.append(data[i,3])
                T3_all.append(data[i,4])

                j+=1

    day_ind+=1

def polyfunc(x_data,c,b,a):
    return a+b*x_data+c*x_data*x_data

fig, axes = pl.subplots(3, 1)
#p01=np.polyfit(T0_all, T1_all, deg=2)
p01, pcov=spo.curve_fit(polyfunc, T0_all, T1_all, p0=[0.001,0.7,0.001])
print(p01)
Ty=np.poly1d(p01)

T0_smooth=np.linspace(np.min(T0_all), np.max(T0_all), 500)

for day in range (4):
    print(day, day_name[day])
    axes[0].scatter(T0[day,:], T1[day,:], c=col[day], s=4, label=day_name[day])
axes[0].scatter(T0_smooth, Ty(T0_smooth), c='k', s=2, label='fit p0='+str(round(p01[2],2))+'\n'+'fit p1='+str(round(p01[1],2))+'\n'+'fit p2='+str(round(p01[0],3)))
axes[0].legend(loc='upper right')
axes[0].set_xlabel('Ch0 T[K]')
axes[0].set_ylabel('Ch1 T[K]')
axes[0].set_xlim([65,190])
axes[0].set_ylim([0,np.asarray(T1_all).max()+10.])

#p02=np.polyfit(T0_all, T2_all, deg=2)
p02, pcov=spo.curve_fit(polyfunc, T0_all, T2_all, p0=[0.001,0.4,0.001])
Ty=np.poly1d(p02)
for day in range (4):
    axes[1].scatter(T0[day,:], T2[day,:], c=col[day], s=4, label=day_name[day])
axes[1].scatter(T0_smooth, Ty(T0_smooth), c='k',s=2, label='fit p0='+str(round(p02[2],2))+'\n'+'fit p1='+str(round(p02[1],2))+'\n'+'fit p2='+str(round(p02[0],3)))
axes[1].legend(loc='upper right')
axes[1].set_xlabel('Ch0 T[K]')
axes[1].set_ylabel('Ch2 T[K]')
axes[1].set_xlim([65,190])
axes[1].set_ylim([0,np.asarray(T2_all).max()+10.])

 #p03=np.polyfit(T0_all, T3_all, deg=2)
p03, pcov=spo.curve_fit(polyfunc, T0_all, T3_all, p0=[0.001,0.25,0.001])
Ty=np.poly1d(p03)
for day in range (4):
    axes[2].scatter(T0[day,:], T3[day,:], c=col[day], s=4, label=day_name[day])
axes[2].scatter(T0_smooth, Ty(T0_smooth), c='k', s=2, label='fit p0='+str(round(p03[2],2))+'\n'+'fit p1='+str(round(p03[1],2))+'\n'+'fit p2='+str(round(p03[0],3)))
axes[2].legend(loc='upper right')
axes[2].set_xlabel('Ch0 T[K]')
axes[2].set_ylabel('Ch3 T[K]')
axes[2].set_xlim([65,190])
axes[2].set_ylim([0,np.asarray(T3_all).max()+10.])

fig.suptitle('All Days SkyDips Zenith Temperatures')#+folderlist)
fig.savefig('Output/All_days_ZenithTemperatures_2ndordfit.png')

pl.show()
