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
from zipfile import ZipFile
import tarfile
import math


x_am=am.AM_functions()
x=pA.ReadAzscan()
x_el=pE.ReadElnod()


class make_plots(object):

    def __init__(self, unit=None, verb=True):

        '''

        '''


    def extract_PWV_skydip(self, year='None', month='None', day='None', time='14', out_folder='paper_plots/', posting_folder='None'):

        #data_folder='wvr1_data/'
        data_folder='../../wvr1_data_local/'
        data_folder_tmp='../../wvr1_data_local/wvr1_data_tmp/'
        # posting_folder='../../Postings/WVR_postings/20210625_SkyDips_Template_Comparison/plots/'

        if year=='None':
            print('Error: Please specify date (at least year).\nyyyy mm dd\nyyyy mm\nyyyy')
            sys.exit()

        file_list=[]

        print(year, month, day)

        #year='2018'
        #month='10'

        if day != 'None':
            for flist in os.listdir(data_folder):
                if flist[:4]==year:
                    if flist[4:6]==month:
                        if flist[6:8]==day:
                            print(data_folder+flist)
                            if flist[-7:]=='.tar.gz':
                                tar = tarfile.open(data_folder+flist, "r:gz")
                                tar.extractall(data_folder)
                                tar.close()
            for flist in os.listdir(data_folder):
                if flist[-15:]=='skyDip_fast.txt':
                    print(flist)
                    if flist[:4]==year:
                        if flist[4:6]==month:
                            if flist[6:8]==day:
                                time_fn=flist[9:11]
                                if time_fn==time:
                                    fn=data_folder+flist
                                    file_list.append(flist)

        elif month != 'None':
            for flist in os.listdir(data_folder):
                if flist[:4]==year:
                    if flist[4:6]==month:
                        if flist[-7:]=='.tar.gz':
                            if not os.path.exists(flist[:-7]+'_fast.txt'):
                                try:
                                    tar = tarfile.open(data_folder+flist, "r:gz")
                                    tar.extractall(data_folder)
                                    tar.close()
                                    #os.system(f'rm '+data_folder+flist)
                                except:
                                    print(data_folder+flist+' failed.')
            for flist in os.listdir(data_folder):
                if flist[-15:]=='skyDip_fast.txt':
                    if flist[:4]==year:
                        if flist[4:6]==month:
                            time_fn=flist[9:11]
                            if time_fn==time:
                                fn=data_folder+flist
                                file_list.append(flist)

        else:
            #uncomment below, eventually
            # for flist in os.listdir(data_folder):
            #     if flist[:4]==year:
            #         print(data_folder+flist)
            #         if flist[-7:]=='.tar.gz':
            #             if not os.path.exists(flist[:-7]+'_fast.txt'):
            #                 try:
            #                     tar = tarfile.open(data_folder+flist, "r:gz")
            #                     tar.extractall(data_folder)
            #                     tar.close()
            #                     #os.system(f'rm '+data_folder+flist)
            #                 except:
            #                     print(data_folder+flist+' failed.')
            for flist in os.listdir(data_folder):
                if flist[-15:]=='skyDip_fast.txt':
                    if flist[:4]==year:
                        time_fn=flist[9:11]
                        if time_fn==time:
                            fn=data_folder+flist
                            file_list.append(flist)

        print('file_list=', file_list)

        for test_file in file_list:

            #date='20200418'
            #time='14'
            #time='19'

            date=test_file[:8]
            time=test_file[9:11]

            print('date=', date)
            print('time=', time)

            #test_data=date+'_'+time+'0002_skyDip_fast.txt'
            #path_to_test=date+'_'+time+'0002_skyDip/'

            test_data=test_file
            path_to_test=data_folder

            print('Creating template file.')
            temp_file='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_'+time+'.amc'
            if not os.path.exists('Templates/'+temp_file):
                x_am.create_template(date[:4]+'-'+date[4:6]+'-'+date[6:8], time)

            standard_temp='SPole_annual_50.amc'
            print('using temp_file=', temp_file)

            print('starting fit.')
            t1=perf_counter()
            #x_am.create_am_datfile(test_data, path_to_data=path_to_test, template= temp_file, spline=2, showplots=0) # not needed-- it is done authomatically inside fit_w_am
            x_am.fit_w_am(test_data, path_to_data=path_to_test, template=standard_temp, spline=2)
            x_am.fit_w_am(test_data, path_to_data=path_to_test, template=temp_file, spline=2)
            t2=perf_counter()

            print('time to fit one SkyDip (2 times) =', t2-t1)

            print('Template:', temp_file)

            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='nscale', pwv_layer='total', template= temp_file, spline=2, pf=out_folder, pf2=posting_folder)

            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='total', template= temp_file, spline=2, pf=out_folder, pf2=posting_folder)
            print('total:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='strato', template= temp_file, spline=2, pf=out_folder, pf2=posting_folder)
            print('strato:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='meso', template= temp_file, spline=2, pf=out_folder, pf2=posting_folder)
            print('meso:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='tropo', template= temp_file, spline=2, pf=out_folder, pf2=posting_folder)
            print('tropo:', el_list, pwv_list)



            print('Template:', standard_temp)

            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='nscale', pwv_layer='total', template= standard_temp, spline=2, pf=out_folder, pf2=posting_folder)

            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='total', template= standard_temp, spline=2, pf=out_folder, pf2=posting_folder)
            print('total:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='strato', template= standard_temp, spline=2, pf=out_folder, pf2=posting_folder)
            print('strato:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='meso', template= standard_temp, spline=2, pf=out_folder, pf2=posting_folder)
            print('meso:', el_list, pwv_list)
            el_list, pwv_list=x_am.plot_am_fit_2(test_data, var='pwv', pwv_layer='tropo', template= standard_temp, spline=2, pf=out_folder, pf2=posting_folder)
            print('tropo:', el_list, pwv_list)



    def make_ch_scatterplots(self, date, out_folder='paper_plots/', posting_folder='None'):

        #date='20200410'

        T0=[]
        T1=[]
        T2=[]
        T3=[]
        el=[]

        filelist=[]
        filelistcut=[]
        for file in os.listdir('../../wvr1_data_local/'):
            #if file[:8]==date:
            if file[:len(date)]==date:
                if file[-7:]=='.tar.gz':
                    filepath='../../wvr1_data_local/'+file
                    if not os.path.exists(filepath[:-7]+'_fast.txt'):
                        print('unzipping '+file)
                        try:
                            tar = tarfile.open('../../wvr1_data_local/'+file, "r:gz")
                            tar.extractall('../../wvr1_data_local/')
                            tar.close()
                        except Exception as e:
                            print('can\'t unzip '+file)
                            print(e)
        for file in os.listdir('../../wvr1_data_local/'):
            if file[:len(date)]==date:
                if file[-15:]=='skyDip_fast.txt':
                    #if file[:8]==date:
                    filepath='../../wvr1_data_local/'+file
                    try:
                        filelist.append(filepath)
                        filelistcut.append(file[:15])
                        fast_data=x_el.read_elnod_fast(filepath)
                        T0_file=fast_data[:,1]
                        T1_file=fast_data[:,2]
                        T2_file=fast_data[:,3]
                        T3_file=fast_data[:,4]
                        el_file=fast_data[:,0]

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
        #
        T0 = np.concatenate(T0, axis=0)
        T1 = np.concatenate(T1, axis=0)
        T2 = np.concatenate(T2, axis=0)
        T3 = np.concatenate(T3, axis=0)
        el = np.concatenate(el, axis=0)

        T0_smooth=np.linspace(T0.min(), T0.max(), 500)

        T1_coef=np.polyfit(T0, T1, 1) #original order = 3
        T1_fit=np.poly1d(T1_coef)
        T1_new=T1_fit(T0_smooth)

        T2_coef=np.polyfit(T0, T2, 1)
        T2_fit=np.poly1d(T2_coef)
        T2_new=T2_fit(T0_smooth)

        T3_coef=np.polyfit(T0, T3, 1)
        T3_fit=np.poly1d(T3_coef)
        T3_new=T3_fit(T0_smooth)


        fig, axes = pl.subplots(3, 1, figsize=(12, 10))
        axes[0].scatter(T0, T1, s=6, c='r')
        axes[0].plot(T0_smooth, T1_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T1_coef[0],2)))
        axes[0].legend(loc='upper right')
        axes[0].set_title('Channel 1')
        axes[0].set_ylabel('T_Ch1[K]')

        axes[1].scatter(T0, T2, s=6, c='darkorange')
        axes[1].plot(T0_smooth, T2_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T2_coef[0],2)))
        axes[1].legend(loc='upper right')
        axes[1].set_title('Channel 2')
        axes[1].set_ylabel('T_Ch2[K]')

        axes[2].scatter(T0, T3, s=6, c='gold')
        axes[2].plot(T0_smooth, T3_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T3_coef[0],2)))
        axes[2].legend(loc='upper right')
        axes[2].set_title('Channel 3')
        axes[2].set_ylabel('T_Ch3[K]')
        axes[2].set_xlabel('T_Ch0[K]')

        pl.suptitle('Zenith Temperatures\n'+date)

        pl.savefig(out_folder+'ZenithT_scatterplot_xch0_'+date+'.png')

        if posting_folder != 'None':
            pl.savefig(posting_folder+'ZenithT_scatterplot_xch0_'+date+'.png')


        #pl.show()
        pl.close()


        T1_smooth=np.linspace(T1.min(), T1.max(), 500)

        T21_coef=np.polyfit(T1, T2, 1)
        T21_fit=np.poly1d(T21_coef)
        T21_new=T21_fit(T1_smooth)

        T31_coef=np.polyfit(T1, T3, 1)
        T31_fit=np.poly1d(T31_coef)
        T31_new=T31_fit(T1_smooth)


        fig, axes = pl.subplots(2, 1, figsize=(12, 10))

        axes[0].scatter(T1, T2, s=6, c='darkorange')
        axes[0].plot(T1_smooth, T21_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T21_coef[0],2)))
        axes[0].legend(loc='upper right')
        axes[0].set_title('Channel 2')
        axes[0].set_ylabel('T_Ch2[K]')

        axes[1].scatter(T1, T3, s=6, c='gold')
        axes[1].plot(T1_smooth, T31_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T31_coef[0],2)))
        axes[1].legend(loc='upper right')
        axes[1].set_title('Channel 3')
        axes[1].set_ylabel('T_Ch3[K]')
        axes[1].set_xlabel('T_Ch1[K]')

        pl.suptitle('Zenith Temperatures\n'+date)

        pl.savefig(out_folder+'ZenithT_scatterplot_xch1_'+date+'.png')

        if posting_folder != 'None':
            pl.savefig(posting_folder+'ZenithT_scatterplot_xch1_'+date+'.png')

        #pl.show()
        pl.close()


        T2_smooth=np.linspace(T2.min(), T2.max(), 500)


        T32_coef=np.polyfit(T2, T3, 1)
        T32_fit=np.poly1d(T32_coef)
        T32_new=T32_fit(T2_smooth)


        fig = pl.figure(figsize=(12, 10))

        pl.scatter(T2, T3, s=6, c='gold')
        pl.plot(T2_smooth, T32_new, c='k', alpha=0.6, label='Lin fit\np0='+str(round(T32_coef[0],2)))
        pl.legend(loc='upper right')
        pl.title('Channel 3')
        pl.ylabel('T_Ch3[K]')
        pl.xlabel('T_Ch2[K]')

        pl.suptitle('Zenith Temperatures\n'+date)

        pl.savefig(out_folder+'ZenithT_scatterplot_xch2_'+date+'.png')

        if posting_folder != 'None':
            pl.savefig(posting_folder+'ZenithT_scatterplot_xch2_'+date+'.png')

        #pl.show()
        pl.close()

        return T1_coef, T2_coef, T3_coef, T21_coef, T31_coef, T32_coef


    def make_ch_scatterplots_fit(self, date_list, pk_fn, out_folder='paper_plots/', posting_folder='None', write_pk=1):

        if write_pk==1:

            T_zenith_LinFit = {'T1_vs_T0_coef': [], 'T2_vs_T0_coef':[], 'T3_vs_T0_coef':[],
            'T2_vs_T1_coef': [], 'T3_vs_T1_coef':[], 'T3_vs_T2_coef':[], 'date_list':[]}

            for date_i in date_list:
                try:
                    T10_coef, T20_coef, T30_coef, T21_coef, T31_coef, T32_coef=self.make_ch_scatterplots(date=date_i, posting_folder=posting_folder)

                    T_zenith_LinFit['date_list'].append(date_i)
                    T_zenith_LinFit['T1_vs_T0_coef'].append(T10_coef[0])
                    T_zenith_LinFit['T2_vs_T0_coef'].append(T20_coef[0])
                    T_zenith_LinFit['T3_vs_T0_coef'].append(T30_coef[0])
                    T_zenith_LinFit['T2_vs_T1_coef'].append(T21_coef[0])
                    T_zenith_LinFit['T3_vs_T1_coef'].append(T31_coef[0])
                    T_zenith_LinFit['T3_vs_T2_coef'].append(T32_coef[0])

                except:
                    print(date_i+' failed.')


            f=open(pk_fn,"wb")
            pk.dump(T_zenith_LinFit, f)
            f.close()


        else:
            f=open(pk_fn,"rb")
            T_zenith_LinFit=pk.load(f)
            f.close()

        fig, axes = pl.subplots(3, 1, figsize=(12, 10))
        axes[0].scatter(T_zenith_LinFit['date_list'], T_zenith_LinFit['T1_vs_T0_coef'], s=6, c='r')
        axes[0].plot(T_zenith_LinFit['date_list'], T_zenith_LinFit['T1_vs_T0_coef'], c='k', alpha=0.6)
        axes[0].legend(loc='upper right')
        axes[0].set_title('Ch 1 vs Ch 0 slope')
        #axes[0].set_ylabel('T_Ch1[K]')

        axes[1].scatter(T_zenith_LinFit['date_list'], T_zenith_LinFit['T2_vs_T0_coef'], s=6, c='darkorange')
        axes[1].plot(T_zenith_LinFit['date_list'], T_zenith_LinFit['T2_vs_T0_coef'], c='k', alpha=0.6)
        axes[1].legend(loc='upper right')
        axes[1].set_title('Ch 2 vs Ch 0 slope')
        #axes[1].set_ylabel('T_Ch2[K]')

        axes[2].scatter(T_zenith_LinFit['date_list'], T_zenith_LinFit['T3_vs_T0_coef'], s=6, c='gold')
        axes[2].plot(T_zenith_LinFit['date_list'], T_zenith_LinFit['T3_vs_T0_coef'], c='k', alpha=0.6)
        axes[2].legend(loc='upper right')
        axes[2].set_title('Ch 3 vs Ch 0 slope')
        #axes[2].set_ylabel('T_Ch3[K]')
        #axes[2].set_xlabel('T_Ch0[K]')

        pl.savefig(out_folder+'ZenithT_fitcoeff_'+date_list[0]+'-'+date_list[len(date_list)-1]+'.png')
        if posting_folder != 'None':
            pl.savefig(posting_folder+'ZenithT_fitcoeff_'+date_list[0]+'-'+date_list[len(date_list)-1]+'.png')

        pl.suptitle('Zenith Temperatures\nLin Fit Coeff\n'+date_list[0]+'-'+date_list[len(date_list)-1])
        #pl.show()
        pl.close()



    def plot_WVR_ps(self, test_file, out_folder='paper_plots/', posting_folder='None'):

        #path_to_test='BAElnod_data/'+test_file[:-9]+'/'
        path_to_test=''

        #The data order  in FH is: TIMEWVR CH0C CH0A CH0H CH0B CH1C CH1A CH1H CH1B CH2C CH2A CH2H CH2B CH3C CH3A CH3H CH3B EL AZ
        D_clean, waz, mod_removed_data, p_double, p_err_double, p_single, p_err_single, calib_data_Gavg, FH =x.read_Az_fast(test_file, pathtofn=path_to_test, clean_mod=3, clean_method='import_model', show_mod_plots=0)

        TIMEWVR=np.array(FH[:,0])
        CH0C=np.array(FH[:,1])
        CH0A=np.array(FH[:,2])
        CH0H=np.array(FH[:,3])
        CH0B=FH[:,4]
        T0=np.array(FH[:,0])
        az=FH[:,18] #unwrapped

        dt_wvr=TIMEWVR[10]-TIMEWVR[9]
        daz_wvr=az[10]-az[9]
        samp_rate=1./dt_wvr
        tot_time=TIMEWVR[len(TIMEWVR)-1]-TIMEWVR[0]
        data_len=len(CH0A)

        time_ordered_time=TIMEWVR#[mask]
        time_ordered_az=az#[mask]

        time_ordered_time=calib_data_Gavg[:,0]

        samp_freq=1/(calib_data_Gavg[1,0]-calib_data_Gavg[0,0])

        time_ordered_az=time_ordered_time*(360./30.)

        pickle_fn='am_datfiles_Az/SPole_annual_50/'+test_file[:-4]+'/'+test_file[:-4]+'_clean_mod3_method_import_model_fitoutput.txt'

        waz, az_wvr, az_real, fs, idx, pwv_full_ts, D_pwv = x.read_pwvatmo(test_file, show=0)

        #Extract Time Domain PS
        time_ordered_time=np.array([i*30. for i in D_pwv[1,:]])
        dt_wvr=30.
        time_ordered_az=time_ordered_time*(360./30.)

        ps_matrix=[]

        date=test_file[:8]
        time=test_file[9:11]
        temp_file='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_'+time+'.amc'

        print(test_file[:-20])

        test_file_skydip=test_file[:-20]+'0002_skyDip_fast.txt'

        el_list, pwv_list=x_am.plot_am_fit_2(test_file_skydip, var='pwv', pwv_layer='tropo', template=temp_file, spline=2)
        el_55_idx=np.where(el_list>=55.)[0]
        el_55_idx=el_55_idx[0]

        rms_skydip=np.nanstd(pwv_list[el_55_idx:])
        print(rms_skydip)

        pl.figure(figsize=(10, 6))

        for a in range (2, len(D_pwv[:,10])-2):

            pwv_ts= D_pwv[a,:]
            mask = np.isfinite(pwv_ts)
            pwv_ts_masked = pwv_ts[mask]

            time_ordered_az_masked = time_ordered_az[mask]
            time_ordered_time_masked = time_ordered_time[mask]
            ps_pwv = x.extract_PowerSpectrum(time_ordered_time_masked, dt_wvr, time_ordered_az_masked, daz_wvr, pwv_ts_masked)

            ps_matrix.append(ps_pwv.t.ps)

            pl.plot(ps_pwv.t.freq, ps_pwv.t.ps, c='k', alpha=0.2)
            pl.scatter(ps_pwv.t.freq, ps_pwv.t.ps)#, label='az='+str(a))


        freq_avg=ps_pwv.t.freq

        ps_matrix_array=np.full((len(ps_matrix[1]), len(ps_matrix)), np.nan)

        ps_avg=np.zeros(len(freq_avg))

        for i in range (len(ps_matrix)-1):
                try:
                    ps_matrix_array[:,i]=ps_matrix[i]
                except:
                    print(str(i)+' failed.')

        for k in range (len(freq_avg)-1):
            ps_avg[k]=np.nanmean(ps_matrix_array[k,:])


        # pl.axvline(x=1./dt_wvr, c='r', label='sampling rate')
        # pl.axvline(x=1./(50*60), c='g', label='Az scan lenght')


        pl.plot(freq_avg[:-1], ps_avg[:-1], c='r', alpha=0.8, label= 'Az average')

        #pl.axvline(x=1./60., c='y', label='1 minute') --> Nyquist frequency = max frequency
        pl.axvline(x=1./(5*60.), c='b', ls='--', alpha=0.6, label='5 minutes')
        pl.axvline(x=1./(10*60.), c='y', ls='--', alpha=0.6, label='10 minutes')

        pl.axhline(y=(rms_skydip*rms_skydip), c='c', ls='-.', label='Expected high-frequency noise level')


        pl.suptitle('PWV Power Spectrum')
        pl.title(test_file[:-4])
        pl.xlabel('1/t[s^(-1)]')
        #pl.xlabel('1/az[deg^(-1)]')
        pl.ylabel('PS(PWV)[um^2]')
        pl.ylim(0.1, 10000000)
        pl.legend()
        pl.loglog()
        pl.savefig(posting_folder+'/'+test_file[:-4]+'_PWV_PS_time.png')
        pl.savefig('paper_plots/'+test_file[:-4]+'_PWV_PS_time.png')

        pl.close()


        #Extract Space Domain PS
        #
        # f = open(pickle_fn,'rb')
        # fit_output = pk.load(f)
        # f.close()
        # pwv_ts=fit_output['pwv']
        # az_pk=fit_output['Az'] #wrapped



        #
        # pwv_ts_masked[~mask]=np.nanmean(pwv_ts_masked)
        # time_ordered_az_masked[~mask]=np.nanmean(time_ordered_az_masked)
        # time_ordered_time_masked[~mask]=np.nanmean(time_ordered_time_masked)

        #
        # print('masked values =', len(pwv_ts_masked[~mask]))

        time_ordered_az = az_real

        pl.figure(figsize=(10, 6))
        rms_i=[]
        daz=time_ordered_az[10]-time_ordered_az[9]
        for scan_i in range(110):

            #pwv_ts_onescan= D_pwv[:,scan_i]
            pwv_ts_onescan = pwv_full_ts[fs.s[scan_i]:fs.e[scan_i]]

            mask = np.isfinite(pwv_ts_onescan)

            pwv_ts_onescan_masked = pwv_ts_onescan[mask]

            time_ordered_az_onescan = time_ordered_az[fs.s[scan_i]:fs.e[scan_i]]

            rms_i.append(np.nanstd(pwv_ts_onescan))
            # pl.scatter(time_ordered_az_onescan, pwv_ts_onescan, label='rms='+str(rms))
            # pl.legend()
            # pl.show()
            time_ordered_az_onescan_masked = time_ordered_az_onescan[mask]

            #(scannum, wrapped_az) = divmod(time_ordered_az_masked_oneturn,360)

            #time_ordered_time below is wrong but I don't care because I just want the space domain

            print('daz=', daz)
            ps_pwv=x.extract_PowerSpectrum(time_ordered_time_masked, dt_wvr, time_ordered_az_onescan_masked, daz, pwv_ts_onescan_masked)


            # Kol_lowerlim=1./360.
            # Kol_upperlim=1./(5)
            #
            # i_kol_lower=np.where(ps_pwv.az.freq>=Kol_lowerlim)[0]
            # i_kol_upper=np.where(ps_pwv.az.freq>=Kol_upperlim)[0]
            #
            # print(i_kol_lower)
            # print(i_kol_lower[0])
            #
            # freq_kol=ps_pwv.az.freq[i_kol_lower[0]:i_kol_upper[0]]
            # ps_kol=ps_pwv.az.ps[i_kol_lower[0]:i_kol_upper[0]]
            #
            #
            # freq_kol_log=[math.log(freq_kol_i, 10) for freq_kol_i in freq_kol]
            # ps_kol_log=[math.log(ps_kol_i, 10) for ps_kol_i in ps_kol]
            #
            # PS_lin_coef=np.polyfit(freq_kol_log, ps_kol_log, 1)
            # PS_lin_fit=np.poly1d(PS_lin_coef)
            # ps_kol_model=PS_lin_fit(freq_kol_log)
            #
            # print(PS_lin_coef)

            pl.scatter(ps_pwv.az.freq, ps_pwv.az.ps)
            pl.plot(ps_pwv.az.freq, ps_pwv.az.ps, c='blue', alpha=0.4)
            # pl.plot(freq_kol_log, ps_kol_log, c='blue', label='pwv')
            # pl.plot(freq_kol_log, ps_kol_model, c='green', ls='--', label='sp_idx='+str(round(PS_lin_coef[0],2)))

            rms=np.nanmean(rms_i)

        pl.axvline(x=1./(2.7), c='y', label='FWHM[Az]=2.7deg')
        pl.axhline(y=(rms*rms), c='c', ls='-.', label='Expected noise level')
        pl.suptitle('Power Spectrum - PWV')
        pl.title(test_file[:-4])
        #pl.xlabel('1/t[s^(-1)]')
        pl.xlabel('1/az[deg^(-1)]')
        pl.ylabel('PS(PWV)[um^2]')
        pl.legend()
        pl.loglog()
        pl.savefig('paper_plots/'+test_file[:-4]+'_PWV_PS_space.png')
        pl.savefig(posting_folder+'/'+test_file[:-4]+'_PWV_PS_space.png')
        #pl.show()
        pl.close()


    def extract_zenith_PWV(self, date='2020', out_folder='paper_plots/', posting_folder='None'):

        T0=[]
        T1=[]
        T2=[]
        T3=[]
        el=[]

        filelist=[]
        filelistcut=[]


        skydip_time_list=['01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23']

        for file in os.listdir('../../wvr1_data_local/'):
            if file[:len(date)]==date:
                time_fn=file[9:11]
                print('time_fn=', time_fn)
                if time_fn in skydip_time_list:
                    if file[-7:]=='.tar.gz':
                        filepath='../../wvr1_data_local/'+file
                        if not os.path.exists(filepath[:-7]+'_fast.txt'):
                            print('unzipping '+file)
                            try:
                                tar = tarfile.open('../../wvr1_data_local/'+file, "r:gz")
                                tar.extractall('../../wvr1_data_local/')
                                tar.close()
                            except Exception as e:
                                print('can\'t unzip '+file)
                                print(e)


        zenith_pwv = {'filename': [], 'date':[], 'time':[], 'pwv_zenith':[], 'real_el':[]}


        for file in os.listdir('../../wvr1_data_local/'):
            if file[:len(date)]==date:
                time_fn=file[9:11]
                if time_fn in skydip_time_list:
                    if file[-15:]=='skyDip_fast.txt':
                        t1=perf_counter()
                        print('Adding file ', file)
                        filepath='../../wvr1_data_local/'+file
                        fit_output=x_am.fit_w_am_zenith(file, path_to_data='../../wvr1_data_local/', template= 'SPole_annual_50.amc', path_to_temp='Templates/')
                        print(fit_output)

                        zenith_pwv['filename'].append(file)
                        zenith_pwv['date'].append(file[:8])
                        zenith_pwv['time'].append(time_fn)
                        zenith_pwv['pwv_zenith'].append(np.mean(fit_output['pwv_tropo']))
                        zenith_pwv['real_el'].append(np.mean(fit_output['El']))

                        t2=perf_counter()

                        print('time for one file [s]:', t2-t1)

                        pl.scatter(zenith_pwv['filename'], zenith_pwv['pwv_zenith'])
                        pl.xticks(ticks=[])
                        pl.title('Last Point:'+file[:8]+'_'+time_fn+'\nTime for fit :'+str(t2-t1)+'s')
                        pl.savefig('paper_plots/zenith_pwv_'+date+'.png')
                        pl.close()




        f = open('zenith_pwv_'+date+'.txt','wb')
        pk.dump(zenith_pwv, f)
        f.close()


    def extract_monthly_PS(self, target_date_list):

        TODs_array=[]
        TOTs_array=[]

        for i in range(len(target_date_list)):

            target_date=target_date_list[i]

            f = open('zenith_pwv_'+target_date+'.txt','rb')
            zenith_pwv=pk.load(f)
            f.close()

            print(zenith_pwv.keys())

            fn=np.array(zenith_pwv['filename'])
            date=np.array(zenith_pwv['date'])
            time=np.array(zenith_pwv['time'])

            pwv=np.array(zenith_pwv['pwv_zenith'])
            real_el=np.array(zenith_pwv['real_el'])


            pl.scatter(real_el, pwv)
            pl.xlabel('Zenith[deg] = SkyDip Max El')
            pl.ylabel('PWV[um]')
            pl.close()

            pl.scatter(fn, real_el)
            pl.ylabel('Zenith[deg] = SkyDip Max El')
            pl.close()

            xticks_pos=[]
            xticks_lab=[]
            time_ordered_time=[]
            count=0
            h_count=0
            for flist in fn:
                xticks_pos.append(count)
                xticks_lab.append(flist[6:8])
                time_ordered_time.append(h_count)
                count+=1
                h_count+=2

            pl.figure(figsize=(14,8))
            pl.scatter(fn, pwv, c='r')
            pl.plot(fn, pwv, c='k', alpha=0.6)
            pl.ylabel('PWV Zenith [um]')
            pl.xticks(ticks=xticks_pos[::24], labels=xticks_lab[::24])
            pl.suptitle('PWV Monthly Variations')
            pl.title(target_date)
            pl.savefig('paper_plots/'+target_date+'_PWV_Monthly_variations.png')
            #pl.show()
            pl.close()


            nan_location=np.isnan(pwv)
            print('How many nans:', len(pwv[nan_location]))
            pwv[nan_location]=0.

            TODs_array.append(pwv)
            TOTs_array.append(time_ordered_time)

        TODs_array=np.array(TODs_array)
        TOTs_array=np.array(TOTs_array)


        def ps_time(time_ordered_time, time_ordered_data):

            #Power Spectrum in time domain
            dt=time_ordered_time[10]-time_ordered_time[9]

            time = len(time_ordered_time)
            df=1./dt

            ps    = np.square(np.abs(np.fft.fft(time_ordered_data)/time))#*np.hanning(len(time_ordered_data)))))
            freq_ = np.fft.fftfreq(time,dt)

            idx_pos=np.where(freq_>=0)

            return freq_[idx_pos], ps[idx_pos]


        color_cycle = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499']
        pl.figure(figsize=(10, 6))


        for i in range(len(TODs_array)):

            time_ordered_time=TOTs_array[i]
            pwv=TODs_array[i]

            freq, ps = ps_time(time_ordered_time, pwv)

            pl.plot(freq, ps, c='k', alpha=0.2)
            pl.scatter(freq, ps, c=color_cycle[i], label=target_date_list[i])


        pl.axvline(x=1./24, c=color_cycle[i+5], ls='--', label='1 day')
        pl.axvline(x=1./(24*7), c='y', ls='--', label='1 week')

        pl.suptitle('PWV Power Spectrum')
        pl.title('Monthly Data')
        pl.xlabel('1/t[h^(-1)]')
        pl.ylabel('PS(PWV)[um^2]')
        #pl.ylim(0.1, 10000000)
        pl.legend()
        pl.loglog()
        pl.savefig('paper_plots/'+target_date+'_PWV_Monthly_PS_time.png')

        pl.show()



    def extract_PWV_atmo(self, date, template='SPole_annual_50.amc', out_folder='paper_plots', posting_folder='None'):

        data_folder='../../wvr1_data_local/'

        day_list=['05', '15', '25']
        fn_list=[]
        for day_i in day_list:
            if os.path.exists('../../wvr1_data_local/'+date+day_i+'_140135_scanAz_fast.txt'):
                fn_list.append(date+day_i+'_140135_scanAz_fast.txt')
            else:
                fn_list.append(date+day_i+'_140134_scanAz_fast.txt')

        print('fn_list=',fn_list)

        # for flist in os.listdir(data_folder):
        #     if flist[-15:]=='scanAz_fast.txt':
        #         if flist[:len(date)]==date:
        #             time=flist[9:11]
        #             if time == '14': #to make just one per day
        #                 print('file:', flist)
        #                 path_to_test='../../wvr1_data_local/'+flist[:-9]+'/'
        #                 path_to_pickle='am_datfiles_Az/'+template[:-4]+'/'+flist[:-4]
        #                 pickle_fn_temps_impmod=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_import_model_pickle_temps.txt'
        #                 #pickle_fn_temps_fit=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_fit_pickle_temps.txt'
        #
        #                 print('Starting fit.')
        #                 x_am.fit_w_am_Az(flist, clean_method='import_model', out_path=out_folder, template=template)
        #                 #x_am.fit_w_am_Az(flist, clean_method='fit')
        for flist in fn_list:

            print('file:', flist)
            path_to_test='../../wvr1_data_local/'+flist[:-9]+'/'
            path_to_pickle='am_datfiles_Az/'+template[:-4]+'/'+flist[:-4]
            pickle_fn_temps_impmod=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_import_model_pickle_temps.txt'
            #pickle_fn_temps_fit=path_to_pickle+'/'+flist[:-4]+'_clean_mod3_clean_method_fit_pickle_temps.txt'

            print('Starting fit.')
            x_am.fit_w_am_Az(flist, clean_method='import_model', out_path=out_folder, template=template)
            #x_am.fit_w_am_Az(flist, clean_method='fit')



    def extract_PWV_skydip_slopes(self, date_list, out_folder='paper_plots/', posting_folder='None'):

        slopes={'date':[], 'SPole_annual_LowEl':[], 'SPole_annual_HighEl':[], 'MERRA_LowEl':[], 'MERRA_HighEl':[], 'PWV_avg':[], 'PWV_avg_MERRA':[]}

        failed_list=[]

        for date in date_list:

            try:

                try:
                    filename=date+'_140002_skyDip_fast.txt'
                    el_list, pwv_list=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template='SPole_annual_50.amc', spline=2, show=0)
                    merra_temp='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_14.amc'
                    el_list_merra, pwv_list_merra=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template=merra_temp, spline=2, show=0)
                except:
                    filename=date+'_140003_skyDip_fast.txt'
                    el_list, pwv_list=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template='SPole_annual_50.amc', spline=2, show=0)
                    merra_temp='MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_14.amc'
                    el_list_merra, pwv_list_merra=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template=merra_temp, spline=2, show=0)


                #SPole_annual_50
                i_low=np.where(el_list<=50)[0]
                i_high=np.where(el_list>50)[0]

                coef_low=np.polyfit(el_list[i_low], pwv_list[i_low], 1)
                fit_low=np.poly1d(coef_low)
                pwv_fit_low=fit_low(el_list)
                print(coef_low)

                coef_high=np.polyfit(el_list[i_high], pwv_list[i_high], 1)
                fit_high=np.poly1d(coef_high)
                pwv_fit_high=fit_high(el_list)

                pl.figure(figsize=(14,10))
                pl.plot(el_list, pwv_list, c='k', alpha=0.7)
                pl.scatter(el_list[i_low], pwv_list[i_low], c='y')
                pl.scatter(el_list[i_high], pwv_list[i_high], c='r')
                pl.plot(el_list, pwv_fit_low, c='y', ls='--', label='LowEL-annual='+str(round(coef_low[0],2)))
                pl.plot(el_list, pwv_fit_high, c='r', ls='--', label='HighEl-annual='+str(round(coef_high[0],2)))


                #MERRA
                i_low_m=np.where(el_list_merra<=55)[0]
                i_high_m=np.where(el_list_merra>45)[0]

                coef_low_merra=np.polyfit(el_list_merra[i_low_m], pwv_list_merra[i_low_m], 1)
                fit_low_merra=np.poly1d(coef_low_merra)
                pwv_fit_low_merra=fit_low_merra(el_list_merra)
                print(coef_low_merra)

                coef_high_merra=np.polyfit(el_list_merra[i_high_m], pwv_list_merra[i_high_m], 1)
                fit_high_merra=np.poly1d(coef_high_merra)
                pwv_fit_high_merra=fit_high(el_list_merra)

                pl.plot(el_list_merra, pwv_list_merra, c='k', alpha=0.7)
                pl.scatter(el_list_merra[i_low_m], pwv_list_merra[i_low_m], c='c')
                pl.scatter(el_list_merra[i_high_m], pwv_list_merra[i_high_m], c='orange')
                pl.plot(el_list_merra, pwv_fit_low_merra, c='c', ls='--', label='LowEL-MERRA='+str(round(coef_low_merra[0],2)))
                pl.plot(el_list_merra, pwv_fit_high_merra, c='orange', ls='--', label='HighEl-MERRA='+str(round(coef_high_merra[0],2)))
                pl.legend(title='Slopes')
                pl.savefig(out_folder+'SkyDips_PWV_slopes_'+date+'.png')
                pl.close()


                slopes['date'].append(date)
                slopes['SPole_annual_LowEl'].append(coef_low[0])
                slopes['SPole_annual_HighEl'].append(coef_high[0])
                slopes['MERRA_LowEl'].append(coef_low_merra[0])
                slopes['MERRA_HighEl'].append(coef_high_merra[0])
                slopes['PWV_avg'].append(np.nanmean(pwv_list))
                slopes['PWV_avg_MERRA'].append(np.nanmean(pwv_list_merra))

            except:

                failed_list.append(date)


        f = open('SkyDips_PWV_slopes_2020.txt','wb')
        pk.dump(slopes, f)
        f.close()

        print('failed_list=', failed_list)
