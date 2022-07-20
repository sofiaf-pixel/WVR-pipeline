import plotElnod as pE
import AM_functions as am
import pickle as pk
import pylab as pl
import os, sys
import numpy as np
import scipy.optimize as spo
from time import perf_counter
import time
import datetime


x_am=am.AM_functions()
x=pE.ReadElnod()

#Fig3
# test_data='20190102_030002_skyDip_fast.txt'
# path_to_test='20190102_030002_skyDip/'
# x.read_elnod_fast(path_to_test+test_data)
#
# sys.exit()

#Fig8 -> T vs time fit for the 4 channels

test_file1='20200424_170002_skyDip_fast.txt'

#FH, fit_par=x.fit_elnod_vsel(test_file)


#Fig5 -> PWV vs El WVR moving up & down

#test_data='20200402_060002_skyDip_fast.txt'
#path_to_test='onedayofSkyDips/20200402/'

#test_data='20200418_180002_skyDip_fast.txt'
#path_to_test='onedayofSkyDips/20200418/'

#test_data='20190103_100002_skyDip_fast.txt'
#path_to_test='onedayofSkyDips/20190103/'

#test_data='20190103_090002_skyDip_fast.txt'
#path_to_test='onedayofSkyDips/20190103/'


test_data='20200418_140002_skyDip_fast.txt'
path_to_test='BAElnod_data/20200418_140002_skyDip/'

#test_data='20200419_130003_skyDip_fast.txt'
#path_to_test='Trial_folder/20200419_130003_skyDip/'

#test_data='20200424_170002_skyDip_fast.txt'
#path_to_test='Trial_folder/20200424_170002_skyDip/'

#test_data='20190103_090002_skyDip_fast.txt'
#path_to_test='BAElnod_data/20190103/'

#test_data='20200406_210002_skyDip_fast.txt'
#path_to_test='BAElnod_data/20200406_210002_skyDip/'

#test_data='20200114_100002_skyDip_fast.txt'
#path_to_test='BAElnod_data/'+test_data[:-9]+'/'

temp_file='SPole_annual_50.amc'

#=['SPole_annual_50.amc', 'SPole_annual_25.amc', 'SPole_annual_75.amc', 'SPole_MAM_50.amc', 'SPole_MAM_25.amc', 'SPole_MAM_75.amc']
temp_list=['SPole_annual_50.amc']

for temp_file in temp_list:
    #try:
    print('starting fit.')
    t1=perf_counter()
    #x_am.create_am_datfile(test_data, path_to_data=path_to_test, template= temp_file, spline=2, showplots=0)
    x_am.fit_w_am(test_data, path_to_data=path_to_test, template= temp_file, spline=2)
    #x_am.fit_w_am(test_data, path_to_data=path_to_test, template=temp_file, spline=0)
    t2=perf_counter()
    el_list, pwv_list=x_am.plot_am_fit_2(test_data, template= temp_file, spline=2)
    #el_list, pwv_list=x_am.plot_am_fit_2(test_data, template= 'SPole_MAM_50.amc', spline=2)
    #el_list, pwv_list=x_am.plot_am_fit_2(test_data, template= temp_file, spline=0)#wo correction
    #except:
    #    print('Template '+temp_file+' failed.')

sys.exit()


print('fit time=', t2-t1)

#folder_name='BAElnod_data/bad_zenith/'
folder_name='BAElnod_data/'
path='wvr1_data/'+folder_name
pickle_fn=path+'zenith_data_pk.txt'


# zenith_eff=[]
# zenith_eff_err=[]
# day=[]
# time=[]
# zenith_data = {'date_xaxis': [], 'zenith_eff':[], 'zenith_eff_err':[]}
# for folderlist in os.listdir(path):
#     if folderlist[-6:]=='skyDip':
#         path_to_test=folder_name+folderlist+'/'
#         for filelist in os.listdir(f'{path+folderlist}'):
#             if filelist[-8:]=='fast.txt':
#                 print('Fitting file '+path+folderlist+filelist+'.')
#                 test_data=filelist
#                 try:
#                     print('day=', filelist[0:8])
#                     print('time=', filelist[9:15])
# #                     # x_am.fit_w_am(filelist, path_to_data=folder_name+folderlist+'/', spline=0)
# #                     # el_list, pwv_list=x_am.plot_am_fit_2(filelist, spline=0)
# #                     #
# #                     # x_am.fit_w_am(filelist, path_to_data=folder_name+folderlist+'/', spline=1)
# #                     # el_list, pwv_list=x_am.plot_am_fit_2(filelist, spline=1)
# #                     #
# #                     # x_am.fit_w_am(filelist, path_to_data=folder_name+folderlist+'/', spline=2)
#                     z, err=x_am.create_am_datfile(test_data, path_to_data=path_to_test, template=temp_file, spline=2, showplots=0, write_dat=0)
#                     zenith_eff.append(z)
#                     zenith_eff_err.append(err)
#                     day.append(filelist[0:8])
#                     time.append(filelist[9:15])
# #                     #el_list, pwv_list=x_am.plot_am_fit_2(filelist, spline=2)
#                 except:
#                     print(filelist+' Failed.')
# date_xaxis=[]
# for i in range(len(day)):
#     str_date=day[i]
#     str_time=time[i]
#
#     YY=int(str_date[:4])
#     MM=int(str_date[4:6])
#     dd=int(str_date[6:8])
#
#     H=int(str_time[:2])
#     m=int(str_time[2:4])
#     s=int(str_time[4:8])
#
#     date_xaxis.append(datetime.datetime(YY, MM, dd, H, m, s, 0))
#     print(i)
#     print(date_xaxis[i])
#
#     zenith_data['date_xaxis']=date_xaxis
#     zenith_data['zenith_eff']=zenith_eff
#     zenith_data['zenith_eff_err']=zenith_eff_err
#
#     print('writing on file '+ pickle_fn)
#     f = open(pickle_fn,'wb')
#     pk.dump(zenith_data, f)
#     f.close()

# f = open(pickle_fn,'rb')
# zenith_data = pk.load(f)
# f.close()
#
# #sys.exit()
#
# dates=np.array(zenith_data['date_xaxis'])
# zenith_eff=np.array(zenith_data['zenith_eff'])
# err=np.array(zenith_data['zenith_eff_err'])
# bad_ind=np.where(zenith_eff<89.)
# bad_days=dates[bad_ind]
# bad_zenith=zenith_eff[bad_ind]
# bad_zenith_err=err[bad_ind]
#
# print('bad_ind=', bad_ind)
# print('bad_days=', bad_days)
#
# pl.errorbar(zenith_data['date_xaxis'], zenith_data['zenith_eff'], yerr=zenith_data['zenith_eff_err'], fmt='.', markersize=3, c='r', ecolor='k')
# pl.axhline(y=np.mean(zenith_data['zenith_eff']), linestyle='--', label='Z_eff_avg[deg]='+str(round(np.mean(zenith_data['zenith_eff']),2))+'\nzenith_eff_std[deg]='+str(round(np.std(zenith_data['zenith_eff']),2)))
# pl.ylabel('Zenith Angle[deg]')
# pl.suptitle('Effective Zenith Scatter Plot')
# pl.legend()
# pl.show()
#
# pl.errorbar(bad_days, bad_zenith, yerr=bad_zenith_err, fmt='.', markersize=3, c='r', ecolor='k')
# pl.suptitle('Effective Zenith Scatter Plot for Bad days')
# pl.legend()
# pl.show()


#Fig10 -> Apr2020 PWV variations


real='BAElnod_data'
month='202007'
trial='Trial_folder/'+month

os.system('rm wvr1_data/Trial_folder/'+month+'*')

if not os.path.exists('wvr1_data/'+trial):
    os.makedirs('wvr1_data/'+trial)

#dict_red=x_am.plot_full_folder(real, spline=1)
# f=open('AprTrial_pwv.txt',"wb")
# pk.dump(dict, f)
# f.close()

folder=trial

for filelist in os.listdir('wvr1_data/BAElnod_data'):
    print(filelist[0:6])
    if filelist[0:6]==month:
        os.system(f"cp wvr1_data/BAElnod_data/"+filelist+" wvr1_data/"+trial+'/'+filelist)

dict_red=x_am.plot_full_folder(folder, picklefn=month+'_elnodpwv.txt', spline=2)

print(dict_red)
#f_red=open('Apr_pwv_red_afterzcorr.txt',"wb")
#pk.dump(dict, f_red)
#f_red.close()


#f = open('AprTrial_pwv.txt','rb')
#dict = pk.load(f)
#f.close()

# day=dict['day']
# time=dict['time']
# pwv=dict['pwv']
# pwv_err=dict['pwv_err']
#
# date_xaxis=[]
#
# for i in range(len(day)):
#     str_date=day[i]
#     str_time=time[i]
#
#     YY=int(str_date[:4])
#     MM=int(str_date[4:6])
#     dd=int(str_date[6:8])
#
#     H=int(str_time[:2])
#     m=int(str_time[2:4])
#     s=int(str_time[4:8])
#
#     date_xaxis.append(datetime.datetime(YY, MM, dd, H, m, s, 0))
#     print(i)
#     print(date_xaxis[i])


# f = open('Apr_pwv_red.txt','rb')
# dict_red = pk.load(f)
# f.close()

#dict_red=dict_red[1]

day_red=dict_red['day']
time_red=dict_red['time']
pwv_red=np.array(dict_red['pwv'])
pwv_err_red=np.array(dict_red['pwv_err'])
print(pwv_red)


# index=np.where(pwv_err_red==np.nanmax(pwv_err_red))
# print('Index:', index)
# index=(index[0])[0]

# print('Fitting the incriminated file.')
# huge_err_fn=day_red[index]+'_'+time_red[index]+'_skyDip_fast.txt'
# folderlist=day_red[index]+'_'+time_red[index]+'_skyDip'
# #x_am.fit_w_am(huge_err_fn, path_to_data=real+'/'+folderlist+'/', spline=1)
# #x_am.create_am_datfile(test_data, path_to_data=path_to_test, spline=1)
# el_list, pwv_list=x_am.plot_am_fit_2(test_data, spline=1)
#
# bad_pt_x=np.zeros(1)
# bad_pt_y=['None']
# bad_pt_y[0]=pwv_red[index]
#
# date_xaxis_red=[]
for i in range(len(day_red)):
    str_date=day_red[i]
    str_time=time_red[i]

    YY=int(str_date[:4])
    MM=int(str_date[4:6])
    dd=int(str_date[6:8])

    H=int(str_time[:2])
    m=int(str_time[2:4])
    s=int(str_time[4:8])

    date_xaxis_red.append(datetime.datetime(YY, MM, dd, H, m, s, 0))

# print('Huge Errorbar Skydip is:')
# print('pwv_err_red:', pwv_err_red)
# print('pwv_err_red Max:', np.nanmax(pwv_err_red))

# print('Index:', index)
# print('Day:', day_red[index])
# print('Time:', time_red[index])
# print('Error:', pwv_err_red[index])
#
# print('date_xaxis_red[index]:', date_xaxis_red[index])

fig, ax = pl.subplots()
pl.errorbar(date_xaxis_red, pwv_red, yerr=pwv_err_red, fmt='.', markersize=3, c='r', ecolor='k')
#pl.errorbar(bad_pt_x, bad_pt_y, yerr=pwv_err_red[index], c='y', fmt='.', label='fn='+time_red[index])
pl.title('PWV variations for Apr2020')
pl.xlabel('date')
pl.ylabel('pwv[um]')
pl.legend()
pl.text(0.10, 0.90, 'pwv_median[um]='+str(round(np.nanmedian(pwv_red),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
pl.text(0.10, 0.80, 'pwv_std[um]='+str(round(np.nanstd(pwv_red),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
pl.setp(ax.xaxis.get_majorticklabels(), rotation=45)

pl.show()
