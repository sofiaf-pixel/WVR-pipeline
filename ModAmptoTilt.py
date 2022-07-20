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
from matplotlib.dates import DateFormatter
from scipy.signal import butter, filtfilt

def butter_lowpass(normal_cutoff, order=5):
    #nyq = 0.5 * fs
    #normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, order=5):
    b, a = butter_lowpass(cutoff, order=order)
    y = filtfilt(b, a, data)
    return y


#pl.rcParams.update({'font.size': 22})
savefig=0

x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()


#scatter plot of mod parameters vs times
path='../../wvr1_data_local/'
outpath='output_plots/single_mod_tilt/'
outhpath_singleslabfit='output_plots/singleslabfit/'

if not os.path.exists(outpath):
    os.makedirs(outpath)

#day_str='20200401'
#day_str='20200410'
#day_str='20200425'
#day_str='20200429'

#day_str_list=['20200403', '20200408', '20200412', '20200415', '20200418', '20200420', '20200423', '20200427']
day_str_list=['20200410', '20200418']
day_str_list_full=[]


#month='202007'
# month='201907'
#month_list=['202004', '202001', '201909',  '201904']
month_list=['202010', '202011', '202012']

for month in month_list:

    # filelist='20200418_140002_skyDip_fast.txt'
    # folderlist='20200418_140002_skyDip'
    # z_offs, z_offs_err, delta=x_am.create_am_datfile(filelist, path_to_data=folder_name+folderlist+'/', template='SPole_annual_50.amc', spline=2, showplots=1, write_dat=0)
    #
    # sys.exit()

    for d in range(1, 30):
        if d<10:
            day_str=month+'0'+str(d)
        else:
            day_str=month+str(d)
        day_str_list_full.append(day_str)

    #print('day_str_list_full=', day_str_list_full)

    convert_to_tilt=1

    #

    tilt_deg={'tilt0':[], 'tilt1':[], 'tilt2':[], 'tilt3':[], 'dT0':[], 'dT1':[], 'dT2':[], 'dT3':[], 'x_axis_onemonth':[]}
    mod_data={'A_double0':[], 'A_double1':[], 'A_double2':[], 'A_double3':[], 'phi_double':[], 'tilt0':[], 'tilt1':[], 'tilt2':[], 'tilt3':[], 'phi_single':[], 'dT0':[], 'dT1':[], 'dT2':[], 'dT3':[], 'x_axis_onemonth':[]}
    mod_par={'single':[], 'double':[]}
    mod_par['single']=np.zeros(5)
    mod_par['double']=np.zeros(5)

    x_axis_onemonth=[]
    A0_single=[]
    A1_single=[]
    A2_single=[]
    A3_single=[]
    phi_single=[]
    C0_single=[]
    C1_single=[]
    C2_single=[]
    C3_single=[]

    dTdEl_0=[]
    dTdEl_1=[]
    dTdEl_2=[]
    dTdEl_3=[]

    z_offs_fit=[]
    z_offs_fit_err=[]

    A0_double=[]
    A1_double=[]
    A2_double=[]
    A3_double=[]
    phi_double=[]
    C0_double=[]
    C1_double=[]
    C2_double=[]
    C3_double=[]


    for day_str in day_str_list_full:
    #
        pickle_fn='../mod_data/'+day_str+'_mod_parameters_constrained_pk.txt'
    #     #
        #To plot the pickle file

        f = open(pickle_fn,'rb')
        params_data = pk.load(f)
        f.close()

        p1=params_data['p_single']
        p1_err=params_data['p_single_err']
        p2=params_data['p_double']
        p2_err=params_data['p_double_err']
        date_xaxis=params_data['date_xaxis']
        print('date_xaxis=', date_xaxis)


        for i in range(len(date_xaxis)):

            #print('date_xaxis=', date_xaxis)


            params_ifile_single=p1[i]
            params_ifile_double=p2[i]

            x_axis_onemonth.append(date_xaxis[i])

            A0_single.append(params_ifile_single[0])
            A1_single.append(params_ifile_single[1])
            A2_single.append(params_ifile_single[2])
            A3_single.append(params_ifile_single[3])

            A0_double.append(params_ifile_double[0])
            A1_double.append(params_ifile_double[1])
            A2_double.append(params_ifile_double[2])
            A3_double.append(params_ifile_double[3])


        # if params_ifile_single[4] < 0.:
        #     params_ifile_single[4] = (2.*np.pi)-params_ifile_single[4]
        #     params_ifile_single[4]=params_ifile_single[4] % (2.*np.pi)
    #
            phi_single.append(degrees(params_ifile_single[4])) #params_ifile_single[4] is in radians
            C0_single.append(params_ifile_single[5])
            C1_single.append(params_ifile_single[6])
            C2_single.append(params_ifile_single[7])
            C3_single.append(params_ifile_single[8])



            phi_double.append(degrees(params_ifile_double[4])) #params_ifile_double[4] is in radians
            C0_double.append(params_ifile_double[5])
            C1_double.append(params_ifile_double[6])
            C2_double.append(params_ifile_double[7])
            C3_double.append(params_ifile_double[8])

            t=date_xaxis[i]
            time=t.strftime ('%H')

            #print('time=', time)


            if convert_to_tilt==1:
                path='../../wvr1_data_local/'
                template='SPole_annual_50.amc'
                print(path+day_str+'_'+time+'0002_skyDip_fast.txt')
                #print(path+day_str+'_'+time+'0002_skyDip'+'/'+day_str+'_'+time+'0002_skyDip_fast.txt')
                if os.path.exists(path+day_str+'_'+time+'0002_skyDip_fast.txt'):
                    filelist=day_str+'_'+time+'0002_skyDip_fast.txt'
                    folderlist=day_str+'_'+time+'0002_skyDip'
                    print('Using Elnod file '+filelist)
                    #print('filelist=', path+folderlist+'/'+filelist)
                    #print('Using new delta values:\n'+str(delta))
                elif os.path.exists(path+day_str+'_'+time+'0003_skyDip_fast.txt'):
                    filelist=day_str+'_'+time+'0003_skyDip_fast.txt'
                    folderlist=day_str+'_'+time+'0003_skyDip'
                    print('Using Elnod file '+filelist)
                    #print('filelist=', path+folderlist+'/'+filelist)
                    #print('Using new delta values:\n'+str(delta))
                else:
                    print('Elnod file not found. Using previous delta values:\n'+str(delta))

                path_to_deltas='am_datfiles/'+template[:-4]+'/'+filelist[:-4]
                pickle_deltas=path_to_deltas+'/'+filelist[:-4]+'_pickle_deltas.txt'
                pickle_z=path_to_deltas+'/'+filelist[:-4]+'_pickle_z_offs.txt'

                if os.path.exists(pickle_deltas):
                    if os.path.exists(pickle_z):
                        print('Reading from file '+pickle_deltas)
                        f1 = open(pickle_deltas,'rb')
                        delta = pk.load(f1)
                        f1.close()
                        print('Reading from file '+pickle_z)
                        f2 = open(pickle_z,'rb')
                        z_offs_elnod = pk.load(f2)
                        f2.close()

                        z_offs=z_offs_elnod['z_offs']
                        z_offs_err=z_offs_elnod['z_offs_err']
                else:
                    print('Extracting z_offs, z_offs_err, delta.')
                    z_offs, z_offs_err, delta=x_am.create_am_datfile(filelist, path_to_data=folder_name+folderlist+'/', template='SPole_annual_50.amc', spline=2, showplots=0, write_dat=0)


                #print('delta[0]=', delta[0])
                #print('delta[1]=', delta[1])
                #print('delta[2]=', delta[2])
                #print('delta[3]=', delta[3])

                dTdEl_0.append(delta[0])
                dTdEl_1.append(delta[1])
                dTdEl_2.append(delta[2])
                dTdEl_3.append(delta[3])

                z_offs_fit.append(z_offs)
                z_offs_fit_err.append(z_offs_err)

    pl.scatter(x_axis_onemonth, phi_single)
    pl.title('phi single raw', fontsize='xx-large')
    if savefig==1:
        pl.savefig(outpath+'/PhiSingleRaw_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/PhiSingleRaw.png')
    pl.show()

    my_figsize=(18,10)

    #plot single mod parameters
    fig, axes = pl.subplots(3,1, sharex=True, figsize=my_figsize)

    axes[0].scatter(x_axis_onemonth, A0_single, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A0_single, c='r', alpha=0.5, label='Ch0')
    axes[0].scatter(x_axis_onemonth, A1_single, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A1_single, c='b', alpha=0.5, label='Ch1')
    axes[0].scatter(x_axis_onemonth, A2_single, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A2_single, c='g', alpha=0.5, label='Ch2')
    axes[0].scatter(x_axis_onemonth, A3_single, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A3_single, c='y', alpha=0.5, label='Ch3')
    axes[0].legend(loc='upper right', fontsize='xx-large')

    #myFmt = mdates.DateFormatter('%H:%M')
    #axes[0].xaxis.set_major_formatter(myFmt)
    axes[0].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[0].set_ylabel('A_singlemod[K]', fontsize='xx-large')
    axes[0].set_title('Single Modulation Amplitude (A)', fontsize='xx-large')

    axes[1].scatter(x_axis_onemonth, C0_single, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C0_single, c='r', alpha=0.5, label='Ch0')
    axes[1].scatter(x_axis_onemonth, C1_single, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C1_single, c='b', alpha=0.5, label='Ch1')
    axes[1].scatter(x_axis_onemonth, C2_single, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C2_single, c='g', alpha=0.5, label='Ch2')
    axes[1].scatter(x_axis_onemonth, C3_single, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C3_single, c='y', alpha=0.5, label='Ch3')
    axes[1].legend(loc='upper right', fontsize='xx-large')

    #axes[1].xaxis.set_major_formatter(myFmt)
    axes[1].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[1].set_ylabel('C_singlemod[K]', fontsize='xx-large')
    axes[1].set_title('Single modulation offset (C)', fontsize='xx-large')

    axes[2].scatter(x_axis_onemonth, np.array(A0_single)/np.array(C0_single), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A0_single)/np.array(C0_single), c='r', alpha=0.5, label='Ch0 - avg='+str(round(np.mean(np.array(A0_single)/np.array(C0_single)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A1_single)/np.array(C1_single), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A1_single)/np.array(C1_single), c='b', alpha=0.5, label='Ch1- avg='+str(round(np.mean(np.array(A1_single)/np.array(C1_single)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A2_single)/np.array(C2_single), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A2_single)/np.array(C2_single), c='g', alpha=0.5, label='Ch2- avg='+str(round(np.mean(np.array(A2_single)/np.array(C2_single)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A3_single)/np.array(C3_single), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A3_single)/np.array(C3_single), c='y', alpha=0.5, label='Ch3- avg='+str(round(np.mean(np.array(A3_single)/np.array(C3_single)),2)))
    axes[2].legend(loc='upper right', fontsize='xx-large')
    #     #
    #axes[2].xaxis.set_major_formatter(myFmt)
    axes[2].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[2].set_ylabel('A/C_singlemod', fontsize='xx-large')
    axes[2].set_title('Percentual Variation of Single modulation Amplitude (A/C)', fontsize='xx-large')
    #     #
    pl.setp(axes[0].get_xticklabels(), visible=False)
    pl.setp(axes[1].get_xticklabels(), visible=False)
    pl.setp(axes[2].get_xticklabels(), Fontsize='xx-large')
    pl.setp(axes[0].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[1].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[2].get_yticklabels(), Fontsize='xx-large')
    #     #
    pl.suptitle('Single Modulation Amplitude Data', fontsize='xx-large')#\n'+str(x_axis_onemonth[0].date()))
    if savefig==1:
        pl.savefig(outpath+'/SingMod_amp_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_amp.png')

    axes[0].set_ylim(0,2)
    axes[2].set_ylim(0,0.02)
    pl.setp(axes[0].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[1].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[2].get_yticklabels(), Fontsize='xx-large')
    if savefig==1:
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_amp_zoomedin.png')

    pl.show()
    #pl.close()
    #     #
    #     #



    phi_single = np.array(phi_single) #% 360.

    nwraps, phi_single_bad=divmod(phi_single,360)

    nsigma_cut=100

    #print('nwraps=', nwraps)
    #print('phi_single=', phi_single)
    #print('phi_single_bad=', phi_single_bad)

    #unwrapping the positive part
    ones=np.full(len(phi_single[np.where(nwraps>0)]),1)
    phi_single[np.where(nwraps>0)]=(360.*(nwraps[np.where(nwraps>0)]+ones))-phi_single[np.where(nwraps>0)]
    #unwrapping the negative part
    ones=np.full(len(phi_single[np.where(nwraps<-1)]),1)
    #print('phi_single[np.where(nwraps<-1)]=', phi_single[np.where(nwraps<-1)])
    #print('nwraps[np.where(nwraps<-1)]=', nwraps[np.where(nwraps<-1)])
    #print('360.*(nwraps[np.where(nwraps<-1)]+ones)=', 360.*(nwraps[np.where(nwraps<-1)]+ones))
    #print('phi_single[np.where(nwraps<-1)]+360.*(nwraps[np.where(nwraps<-1)]+ones)=', phi_single[np.where(nwraps<-1)]-360.*(nwraps[np.where(nwraps<-1)]+ones))
    phi_single[np.where(nwraps<-1)]=phi_single[np.where(nwraps<-1)]-360.*(nwraps[np.where(nwraps<-1)]+ones)
    #print('phi_single[np.where(phi_single>0.)]=', phi_single[np.where(phi_single>0.)])
    phi_single[np.where(phi_single>0.)]=phi_single[np.where(phi_single>0.)]-360.

    ax=pl.figure(figsize=(18,10))
    pl.scatter(x_axis_onemonth, phi_single, c='k', s=3)
    pl.plot(x_axis_onemonth, phi_single, c='c', alpha=0.5)

    phi_cut=phi_single

    phi_std=np.std(phi_cut)
    phi_mean=np.mean(phi_cut)

    # phi_cut[np.where(phi_cut>(phi_mean+(nsigma_cut*phi_std)))]=np.nan
    # phi_cut[np.where(phi_cut<(phi_mean-(nsigma_cut*phi_std)))]=np.nan

    pl.axhline(y=np.nanmean(phi_cut), color='k', linestyle='--', label='phi_avg[deg]='+str(round(np.nanmean(phi_cut)))+'\nphi_std[deg]='+str(round(np.nanstd(phi_cut))), alpha=0.5)
    pl.axhline(y=phi_mean+(nsigma_cut*phi_std), color='r', linestyle='--', label='upper cut', alpha=0.5)
    pl.axhline(y=phi_mean-(nsigma_cut*phi_std), color='b', linestyle='--', label='lower cut', alpha=0.5)
    #myFmt = mdates.DateFormatter('%H:%M')
    #pl.gca().xaxis.set_major_formatter(myFmt)
    pl.xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    #pl.ylim(-180.,+180.)
    pl.ylabel('phi_singlemod [deg]', fontsize='xx-large')
    pl.legend(fontsize='xx-large')
    pl.suptitle('Single Modulation phase', fontsize='xx-large')#\n'+str(x_axis_onemonth[0].date()))
    pl.xticks(fontsize='xx-large')
    pl.yticks(fontsize='xx-large')
    if savefig==1:
        pl.savefig(outpath+'/SingleMod_phase_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_phase.png')
    pl.show()

    mod_data['phi_single']=phi_cut

    #plot double mod parameters
    fig, axes = pl.subplots(3,1, sharex=True, figsize=my_figsize)
    #     #
    axes[0].scatter(x_axis_onemonth, A0_double, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A0_double, c='r', alpha=0.5, label='Ch0')
    axes[0].scatter(x_axis_onemonth, A1_double, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A1_double, c='b', alpha=0.5, label='Ch1')
    axes[0].scatter(x_axis_onemonth, A2_double, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A2_double, c='g', alpha=0.5, label='Ch2')
    axes[0].scatter(x_axis_onemonth, A3_double, c='k', s=3)
    axes[0].plot(x_axis_onemonth, A3_double, c='y', alpha=0.5, label='Ch3')
    axes[0].legend(loc='upper right', fontsize='xx-large')
    #     #
    #myFmt = mdates.DateFormatter('%H:%M')
    #axes[0].xaxis.set_major_formatter(myFmt)
    axes[0].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[0].set_ylabel('A_doublemod[K]', fontsize='xx-large')
    axes[0].set_title('Double Modulation Amplitude (A)', fontsize='xx-large')
    #     #
    axes[1].scatter(x_axis_onemonth, C0_double, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C0_double, c='r', alpha=0.5, label='Ch0')
    axes[1].scatter(x_axis_onemonth, C1_double, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C1_double, c='b', alpha=0.5, label='Ch1')
    axes[1].scatter(x_axis_onemonth, C2_double, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C2_double, c='g', alpha=0.5, label='Ch2')
    axes[1].scatter(x_axis_onemonth, C3_double, c='k', s=3)
    axes[1].plot(x_axis_onemonth, C3_double, c='y', alpha=0.5, label='Ch3')
    axes[1].legend(loc='upper right', fontsize='xx-large')
    #     #
    #axes[1].xaxis.set_major_formatter(myFmt)
    axes[1].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[1].set_ylabel('C_doublemod[K]', fontsize='xx-large')
    axes[1].set_title('Double modulation offset (C)', fontsize='xx-large')
    #     #
    axes[2].scatter(x_axis_onemonth, np.array(A0_double)/np.array(C0_double), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A0_double)/np.array(C0_double), c='r', alpha=0.5, label='Ch0 - avg='+str(round(np.mean(np.array(A0_double)/np.array(C0_double)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A1_double)/np.array(C1_double), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A1_double)/np.array(C1_double), c='b', alpha=0.5, label='Ch1 - avg='+str(round(np.mean(np.array(A1_double)/np.array(C1_double)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A2_double)/np.array(C2_double), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A2_double)/np.array(C2_double), c='g', alpha=0.5, label='Ch2 - avg='+str(round(np.mean(np.array(A2_double)/np.array(C2_double)),2)))
    axes[2].scatter(x_axis_onemonth, np.array(A3_double)/np.array(C3_double), c='k', s=3)
    axes[2].plot(x_axis_onemonth, np.array(A3_double)/np.array(C3_double), c='y', alpha=0.5, label='Ch3 - avg='+str(round(np.mean(np.array(A3_double)/np.array(C3_double)),2)))
    axes[2].legend(loc='upper right', fontsize='xx-large')
    #     #
    #axes[2].xaxis.set_major_formatter(myFmt)
    axes[2].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[2].set_ylabel('A/C_doublemod', fontsize='xx-large')
    axes[2].set_title('Percentual Variation of double modulation Amplitude (A/C)', fontsize='xx-large')
    #     #
    pl.setp(axes[0].get_xticklabels(), visible=False)
    pl.setp(axes[1].get_xticklabels(), visible=False)
    pl.setp(axes[2].get_xticklabels(), Fontsize='xx-large')
    pl.setp(axes[0].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[1].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[2].get_yticklabels(), Fontsize='xx-large')

    #     #
    pl.suptitle('Double Modulation Amplitude Data', fontsize='xx-large')#\n'+str(x_axis_onemonth[0].date()))
    if savefig==1:
        pl.savefig(outpath+'/DoubleMod_amp_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/DoubleMod_amp.png')

    axes[0].set_ylim(0.15,0.75)
    pl.setp(axes[0].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[1].get_yticklabels(), Fontsize='xx-large')
    pl.setp(axes[2].get_yticklabels(), Fontsize='xx-large')

    if savefig==1:
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/DoubleMod_amp_zoomedin.png')

    pl.show()
    #pl.close()
    #     #
    #     #


    mod_data['A0_double']=A0_double
    mod_data['A1_double']=A1_double
    mod_data['A2_double']=A2_double
    mod_data['A3_double']=A3_double

    phi_double = np.array(phi_double)

    #print('phi_double[np.where(phi_double>0.)]=', phi_double[np.where(phi_double>0.)])
    phi_double[np.where(phi_double>0.)]=phi_double[np.where(phi_double>0.)]-360. #is it right? Double mod phase wraps at 180 deg

    ax=pl.figure(figsize=(18,10))
    pl.scatter(x_axis_onemonth, phi_double, c='k', s=3)
    pl.plot(x_axis_onemonth, phi_double, c='c', alpha=0.5)

    phi_cut_2=np.array(phi_double)

    phi_std=np.std(phi_cut_2)
    phi_mean=np.mean(phi_cut_2)

    phase_tilt_ind=np.where(phi_cut_2>(phi_mean+(nsigma_cut*phi_std)))
    #print('phase_tilt_ind=', phase_tilt_ind)

    phi_cut_2[np.where(phi_cut_2>(phi_mean+(nsigma_cut*phi_std)))]=np.nan
    phi_cut_2[np.where(phi_cut_2<(phi_mean-(nsigma_cut*phi_std)))]=np.nan

    # phi_std=np.nanstd(phi_cut_2)
    # phi_mean=np.nanmean(phi_cut_2)
    #
    # phi_cut_2[np.where(phi_cut_2>(phi_mean+phi_std))]=np.nan
    # phi_cut_2[np.where(phi_cut_2<(phi_mean-phi_std))]=np.nan

    phi_cut_2=phi_cut_2[np.where(phi_cut_2!=np.nan)]

    mod_data['phi_double']=phi_cut_2

    pl.axhline(y=np.nanmean(phi_cut_2), color='k', linestyle='--', label='phi_avg[deg]='+str(round(np.nanmean(phi_cut_2)))+'\nphi_std[deg]='+str(round(np.nanstd(phi_cut_2))), alpha=0.5)
    pl.axhline(y=phi_mean+(nsigma_cut*phi_std), color='r', linestyle='--', label='upper cut', alpha=0.5)
    pl.axhline(y=phi_mean-(nsigma_cut*phi_std), color='b', linestyle='--', label='lower cut', alpha=0.5)
    #myFmt = mdates.DateFormatter('%H:%M')
    #pl.gca().xaxis.set_major_formatter(myFmt)
    pl.xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    #pl.ylim(-180.,+180.)
    pl.ylabel('phi_doublemod [deg]', fontsize='xx-large')
    pl.legend(fontsize='xx-large')
    pl.suptitle('Double Modulation phase', fontsize='xx-large')#\n'+str(x_axis_onemonth[0].date()))
    pl.xticks(fontsize='xx-large')
    pl.yticks(fontsize='xx-large')
    if savefig==1:
        pl.savefig(outpath+'/DoubleMod_phase_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/DoubleMod_phase.png')
    pl.show()

    A0_single=np.array(A0_single)
    A1_single=np.array(A1_single)
    A2_single=np.array(A2_single)
    A3_single=np.array(A3_single)


    #print('Mean Tilt0 Before nan')
    #print(np.nanmean(A0_single/dTdEl_0))

    A0_single[phase_tilt_ind]=np.nan
    A1_single[phase_tilt_ind]=np.nan
    A2_single[phase_tilt_ind]=np.nan
    A3_single[phase_tilt_ind]=np.nan

    dTdEl_0=np.array(dTdEl_0)
    dTdEl_1=np.array(dTdEl_1)
    dTdEl_2=np.array(dTdEl_2)
    dTdEl_3=np.array(dTdEl_3)

    dTdEl_0[phase_tilt_ind]=np.nan
    dTdEl_1[phase_tilt_ind]=np.nan
    dTdEl_2[phase_tilt_ind]=np.nan
    dTdEl_3[phase_tilt_ind]=np.nan

    #print('Mean Tilt0 After nan')
    #print(np.nanmean(A0_single/dTdEl_0))

    if convert_to_tilt==1:

        tilt_deg['dT0']=dTdEl_0
        tilt_deg['dT1']=dTdEl_1
        tilt_deg['dT2']=dTdEl_2
        tilt_deg['dT3']=dTdEl_3

        tilt_deg['tilt0']=A0_single/dTdEl_0
        tilt_deg['tilt1']=A1_single/dTdEl_1
        tilt_deg['tilt2']=A2_single/dTdEl_2
        tilt_deg['tilt3']=A3_single/dTdEl_3

        mod_data['dT0']=dTdEl_0
        mod_data['dT1']=dTdEl_1
        mod_data['dT2']=dTdEl_2
        mod_data['dT3']=dTdEl_3

        mod_data['tilt0']=A0_single/dTdEl_0
        mod_data['tilt1']=A1_single/dTdEl_1
        mod_data['tilt2']=A2_single/dTdEl_2
        mod_data['tilt3']=A3_single/dTdEl_3

        tilt_deg['x_axis_onemonth']=x_axis_onemonth
        mod_data['x_axis_onemonth']=x_axis_onemonth

        f2 = open('tilt_angle_fullmonth'+month+'.txt','wb')
        pk.dump(tilt_deg, f2)
        f2.close()

        fig = pl.figure(figsize=my_figsize)

        pl.scatter(x_axis_onemonth, A0_single/dTdEl_0, c='k', s=3)
        pl.plot(x_axis_onemonth, A0_single/dTdEl_0, c='r', alpha=0.8)
        pl.axhline(y=np.nanmean(A0_single/dTdEl_0), color='r', linestyle='--', label='Ch0_avg[deg]='+str(round(np.nanmean(A0_single/dTdEl_0),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, A1_single/dTdEl_1, c='k', s=3)
        pl.plot(x_axis_onemonth, A1_single/dTdEl_1, c='b', alpha=0.8)
        pl.axhline(y=np.nanmean(A1_single/dTdEl_1), color='b', linestyle='--', label='Ch1_avg[deg]='+str(round(np.nanmean(A1_single/dTdEl_1),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, A2_single/dTdEl_2, c='k', s=3)
        pl.plot(x_axis_onemonth, A2_single/dTdEl_2, c='g', alpha=0.8)
        pl.axhline(y=np.nanmean(A2_single/dTdEl_2), color='g', linestyle='--', label='Ch2_avg[deg]='+str(round(np.nanmean(A2_single/dTdEl_2),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, A3_single/dTdEl_3, c='k', s=3)
        pl.plot(x_axis_onemonth, A3_single/dTdEl_3, c='y', alpha=0.8)
        pl.axhline(y=np.nanmean(A3_single/dTdEl_3), color='y', linestyle='--', label='Ch3_avg[deg]='+str(round(np.nanmean(A3_single/dTdEl_3),2)), alpha=0.5)
        pl.legend(loc='upper right', fontsize='xx-large')

        myFmt = mdates.DateFormatter('%D')
        pl.gca().xaxis.set_major_formatter(myFmt)
        pl.xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
        pl.ylabel('A_tilt[deg]', fontsize='xx-large')
        pl.title('Single Modulation Amplitude converted into Tilt Angle', fontsize='xx-large')
        pl.xticks(fontsize='xx-large')
        pl.yticks(fontsize='xx-large')
        if savefig==1:
            pl.savefig(outpath+'/SingleMod_tilt_'+month+'.png')
            pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_Tilt.png')
        pl.ylim(0.4,2.)
        if savefig==1:
            pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_Tilt_zoomedin.png')
        pl.xticks(fontsize='xx-large')
        pl.yticks(fontsize='xx-large')
        pl.show()
    #   pl.close()


        normal_cutoff=0.06

        tilt_smooth_0 = butter_lowpass_filtfilt(A0_single/dTdEl_0, normal_cutoff)
        tilt_smooth_1 = butter_lowpass_filtfilt(A1_single/dTdEl_1, normal_cutoff)
        tilt_smooth_2 = butter_lowpass_filtfilt(A2_single/dTdEl_2, normal_cutoff)
        tilt_smooth_3 = butter_lowpass_filtfilt(A3_single/dTdEl_3, normal_cutoff)

        fig = pl.figure(figsize=my_figsize)

        pl.scatter(x_axis_onemonth, tilt_smooth_0, c='k', s=3)
        pl.plot(x_axis_onemonth, tilt_smooth_0, c='r', alpha=0.8)
        pl.axhline(y=np.nanmean(A0_single/dTdEl_0), color='r', linestyle='--', label='Ch0_avg[deg]='+str(round(np.nanmean(A0_single/dTdEl_0),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, tilt_smooth_1, c='k', s=3)
        pl.plot(x_axis_onemonth, tilt_smooth_1, c='b', alpha=0.8)
        pl.axhline(y=np.nanmean(A1_single/dTdEl_1), color='b', linestyle='--', label='Ch1_avg[deg]='+str(round(np.nanmean(A1_single/dTdEl_1),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, tilt_smooth_2, c='k', s=3)
        pl.plot(x_axis_onemonth, tilt_smooth_2, c='g', alpha=0.8)
        pl.axhline(y=np.nanmean(A2_single/dTdEl_2), color='g', linestyle='--', label='Ch2_avg[deg]='+str(round(np.nanmean(A2_single/dTdEl_2),2)), alpha=0.5)
        pl.scatter(x_axis_onemonth, tilt_smooth_3, c='k', s=3)
        pl.plot(x_axis_onemonth, tilt_smooth_3, c='y', alpha=0.8)
        pl.axhline(y=np.nanmean(A3_single/dTdEl_3), color='y', linestyle='--', label='Ch3_avg[deg]='+str(round(np.nanmean(A3_single/dTdEl_3),2)), alpha=0.5)
        pl.legend(loc='upper right', fontsize='xx-large')

        myFmt = mdates.DateFormatter('%D')
        pl.gca().xaxis.set_major_formatter(myFmt)
        pl.xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
        pl.ylabel('A_tilt[deg]', fontsize='xx-large')
        pl.title('Single Modulation Amplitude converted into Tilt Angle', fontsize='xx-large')
        pl.xticks(fontsize='xx-large')
        pl.yticks(fontsize='xx-large')
        if savefig==1:
            pl.savefig(outpath+'/SingleMod_tilt_'+month+'.png')
            pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_Tilt_smooth.png')
        pl.ylim(0.4,2.)
        if savefig==1:
            pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_Tilt_smooth_zoomedin.png')
        pl.xticks(fontsize='xx-large')
        pl.yticks(fontsize='xx-large')
        pl.show()


        #plot conversion
        fig, axes = pl.subplots(3,1, sharex=True, figsize=my_figsize)
        #     #
        axes[0].scatter(x_axis_onemonth, A0_single, c='k', s=3)
        axes[0].plot(x_axis_onemonth, A0_single, c='r', alpha=0.5, label='Ch0')
        axes[0].scatter(x_axis_onemonth, A1_single, c='k', s=3)
        axes[0].plot(x_axis_onemonth, A1_single, c='b', alpha=0.5, label='Ch1')
        axes[0].scatter(x_axis_onemonth, A2_single, c='k', s=3)
        axes[0].plot(x_axis_onemonth, A2_single, c='g', alpha=0.5, label='Ch2')
        axes[0].scatter(x_axis_onemonth, A3_single, c='k', s=3)
        axes[0].plot(x_axis_onemonth, A3_single, c='y', alpha=0.5, label='Ch3')
        axes[0].legend(loc='upper right', fontsize='xx-large')
        #     #
        #myFmt = mdates.DateFormatter('%H:%M')
        #axes[0].xaxis.set_major_formatter(myFmt)
        axes[0].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
        axes[0].set_ylabel('A_singlemod[K]', fontsize='xx-large')
        axes[0].set_title('Single Modulation Amplitude (A)', fontsize='xx-large')
        #     #
        axes[1].scatter(x_axis_onemonth, dTdEl_0, c='k', s=3)
        axes[1].plot(x_axis_onemonth, dTdEl_0, c='r', alpha=0.5, label='Ch0_avg[K/deg]='+str(round(np.nanmean(dTdEl_0),2)))
        axes[1].scatter(x_axis_onemonth, dTdEl_1, c='k', s=3)
        axes[1].plot(x_axis_onemonth, dTdEl_1, c='b', alpha=0.5, label='Ch1_avg[K/deg]='+str(round(np.nanmean(dTdEl_1),2)))
        axes[1].scatter(x_axis_onemonth, dTdEl_2, c='k', s=3)
        axes[1].plot(x_axis_onemonth, dTdEl_2, c='g', alpha=0.5, label='Ch2_avg[K/deg]='+str(round(np.nanmean(dTdEl_2),2)))
        axes[1].scatter(x_axis_onemonth, dTdEl_3, c='k', s=3)
        axes[1].plot(x_axis_onemonth, dTdEl_3, c='y', alpha=0.5, label='Ch3_avg[K/deg]='+str(round(np.nanmean(dTdEl_3),2)))
        axes[1].legend(loc='upper right', fontsize='xx-large')
        #     #
        #axes[1].xaxis.set_major_formatter(myFmt)
        axes[1].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
        axes[1].set_ylabel('ΔT/Δ(el)[K/deg]', fontsize='xx-large')
        axes[1].set_title('dT/dEl from Elnod', fontsize='xx-large')
        #     #
        axes[2].scatter(x_axis_onemonth, A0_single/dTdEl_0, c='k', s=3)
        axes[2].plot(x_axis_onemonth, A0_single/dTdEl_0, c='r', alpha=0.5)
        axes[2].axhline(y=np.nanmean(A0_single/dTdEl_0), color='r', linestyle='--', label='Ch0_avg[deg]='+str(round(np.nanmean(A0_single/dTdEl_0),2)), alpha=0.5)
        axes[2].scatter(x_axis_onemonth, A1_single/dTdEl_1, c='k', s=3)
        axes[2].plot(x_axis_onemonth, A1_single/dTdEl_1, c='b', alpha=0.5)
        axes[2].axhline(y=np.nanmean(A1_single/dTdEl_1), color='r', linestyle='--', label='Ch1_avg[deg]='+str(round(np.nanmean(A1_single/dTdEl_1),2)), alpha=0.5)
        axes[2].scatter(x_axis_onemonth, A2_single/dTdEl_2, c='k', s=3)
        axes[2].plot(x_axis_onemonth, A2_single/dTdEl_2, c='g', alpha=0.5)
        axes[2].axhline(y=np.nanmean(A2_single/dTdEl_2), color='r', linestyle='--', label='Ch2_avg[deg]='+str(round(np.nanmean(A2_single/dTdEl_2),2)), alpha=0.5)
        axes[2].scatter(x_axis_onemonth, A3_single/dTdEl_3, c='k', s=3)
        axes[2].plot(x_axis_onemonth, A3_single/dTdEl_3, c='y', alpha=0.5)
        axes[2].axhline(y=np.nanmean(A3_single/dTdEl_3), color='r', linestyle='--', label='Ch3_avg[deg]='+str(round(np.nanmean(A3_single/dTdEl_3),2)), alpha=0.5)
        axes[2].legend(loc='upper right', fontsize='xx-large')
        #     #
        #axes[2].xaxis.set_major_formatter(myFmt)
        axes[2].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
        axes[2].set_ylabel('A_tilt[deg]', fontsize='xx-large')
        axes[2].set_title('Tilt Angle', fontsize='xx-large')
        #     #
        pl.setp(axes[0].get_xticklabels(), visible=False)
        pl.setp(axes[1].get_xticklabels(), visible=False)
        pl.setp(axes[2].get_xticklabels(), Fontsize='xx-large')
        pl.setp(axes[0].get_yticklabels(), Fontsize='xx-large')
        pl.setp(axes[1].get_yticklabels(), Fontsize='xx-large')
        pl.setp(axes[2].get_yticklabels(), Fontsize='xx-large')
        #     #
        pl.suptitle('Tilt Angle from Single Mod Amp data', fontsize='xx-large')#\n'+str(x_axis_onemonth[0].date()))
        if savefig==1:
            pl.savefig(outpath+'/SingleMod_Amp_and_Tilt_'+month+'.png')
            pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/SingleMod_Amp_and_Tilt.png')
        pl.show()

    mod_par['single'][0]=np.nanmean(A0_single/dTdEl_0)
    mod_par['single'][1]=np.nanmean(A1_single/dTdEl_1)
    mod_par['single'][2]=np.nanmean(A2_single/dTdEl_2)
    mod_par['single'][3]=np.nanmean(A3_single/dTdEl_3)


    tilt0_El=(A0_single/dTdEl_0)*np.sin(phi_cut*np.pi/180.)
    tilt1_El=(A1_single/dTdEl_1)*np.sin(phi_cut*np.pi/180.)
    tilt2_El=(A2_single/dTdEl_2)*np.sin(phi_cut*np.pi/180.)
    tilt3_El=(A3_single/dTdEl_3)*np.sin(phi_cut*np.pi/180.)


    fig, axes = pl.subplots(2,1, sharex=True, figsize=my_figsize)
    axes[0].errorbar(x_axis_onemonth, z_offs_fit, yerr=z_offs_fit_err, fmt='.', markersize=3, c='k')
    axes[0].axhline(y=np.nanmean(z_offs_fit), color='r', linestyle='--', label='avg[deg]='+str(round(np.nanmean(z_offs_fit),2)), alpha=0.5)
    axes[0].legend(loc='upper right', fontsize='xx-large')
    axes[0].set_ylabel('zenith_tilt[deg]')
    axes[0].set_title('Tilt angle from Elnods')
    axes[1].scatter(x_axis_onemonth, tilt0_El, c='r', s=3)
    axes[1].axhline(y=np.nanmean(tilt0_El), color='r', linestyle='--', label='Ch0_avg[deg]='+str(round(np.nanmean(tilt0_El),2)), alpha=0.5)
    axes[1].scatter(x_axis_onemonth, tilt1_El, c='b', s=3)
    axes[1].axhline(y=np.nanmean(tilt1_El), color='r', linestyle='--', label='Ch1_avg[deg]='+str(round(np.nanmean(tilt1_El),2)), alpha=0.5)
    axes[1].scatter(x_axis_onemonth, tilt2_El, c='g', s=3)
    axes[1].axhline(y=np.nanmean(tilt2_El), color='r', linestyle='--', label='Ch2_avg[deg]='+str(round(np.nanmean(tilt2_El),2)), alpha=0.5)
    axes[1].axhline(y=np.nanmean(tilt3_El), color='r', linestyle='--', label='Ch3_avg[deg]='+str(round(np.nanmean(tilt3_El),2)), alpha=0.5)
    axes[1].set_ylabel('A_tilt*sin(phase)[deg]')
    axes[1].legend(loc='upper right', fontsize='xx-large')
    axes[1].set_xlim(np.min(x_axis_onemonth)-datetime.timedelta(hours=1), np.max(x_axis_onemonth)+datetime.timedelta(hours=1))
    axes[1].set_title('Tilt angle from AzScans')
    if savefig==1:
        pl.savefig(outhpath_singleslabfit+'/Tilt_Elnod_compared_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/Tilt_Elnod_compared.png')
    pl.close()
    #
    # x_axis_onemonth_cut=np.zeros(len(x_axis_onemonth))
    # for i in range (len(x_axis_onemonth)):
    #     date=x_axis_onemonth[i]
    #     x_axis_onemonth_cut[i]=date[5:]
    #

    fig, ax = pl.subplots(figsize=my_figsize)
    pl.errorbar(x_axis_onemonth, z_offs_fit, yerr=z_offs_fit_err, fmt='.', markersize=3, c='k')
    pl.axhline(y=np.nanmean(z_offs_fit), color='r', linestyle='--', label='avg[deg]='+str(round(np.nanmean(z_offs_fit),2)), alpha=0.5)
    pl.legend(loc='upper right', fontsize='xx-large')
    pl.ylabel('zenith_tilt[deg]', fontsize='xx-large')
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    pl.suptitle('Tilt angle from Elnods', fontsize='xx-large')
    pl.title(month[:4]+'-'+month[4:], fontsize='xx-large')
    pl.tick_params(axis='y', labelsize=18)
    pl.tick_params(axis='x', labelsize=18, rotation=45)

    if savefig==1:
        pl.savefig(outhpath_singleslabfit+'/Tilt_Elnod_'+month+'.png')
        pl.savefig('../Postings/WVR_postings/20210215_Single_Mod_Amp_to_Tilt/plots/single_mod_tilt/Tilt_Elnod.png')
    pl.close()

    f = open('modulation_data_'+month+'.txt','wb')
    pk.dump(mod_data, f)
    f.close()

    mod_par['single'][4]=np.nanmean(mod_data['phi_single'])

    mod_par['double'][0]=np.mean(mod_data['A0_double'])
    mod_par['double'][1]=np.mean(mod_data['A1_double'])
    mod_par['double'][2]=np.mean(mod_data['A2_double'])
    mod_par['double'][3]=np.mean(mod_data['A3_double'])

    mod_par['double'][4]=np.nanmean(mod_data['phi_double'])

    A=mod_par['double']
    print(A)

    #
    #
    f = open('modulation_parameters_'+month+'.txt','wb')
    pk.dump(mod_par, f)
    f.close()

    if not (os.path.exists('doubmemod_data/')):
        os.system('mkdir doubmemod_data')

    double_pkfn='doubmemod_data/doublemod_'+month+'.txt'

    doublemod_dict={'date':[], 'A0_double':[], 'A1_double':[], 'A2_double':[], 'A3_double':[], 'offset0':[], 'offset1':[], 'offset2':[], 'offset3':[]}
    doublemod_dict['date']=mod_data['x_axis_onemonth']

    doublemod_dict['A0_double']=mod_data['A0_double']
    doublemod_dict['A1_double']=mod_data['A1_double']
    doublemod_dict['A2_double']=mod_data['A2_double']
    doublemod_dict['A3_double']=mod_data['A3_double']

    doublemod_dict['offset0']=C0_double
    doublemod_dict['offset1']=C1_double
    doublemod_dict['offset2']=C2_double
    doublemod_dict['offset3']=C3_double

    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.plot(doublemod_dict['date'], doublemod_dict['A0_double'])
    pl.title('Double Mod Amp')
    pl.show()

    pl.plot(doublemod_dict['date'], doublemod_dict['offset0'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset1'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset2'])
    pl.plot(doublemod_dict['date'], doublemod_dict['offset3'])
    pl.title('Atm T')
    pl.show()

    print(doublemod_dict)

    f = open(double_pkfn,'wb')
    pk.dump(doublemod_dict, f)
    f.close()
