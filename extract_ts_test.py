import Read_BICEP_ts as bts
import plotAzscan as pA
import numpy as np
import pylab as pl
import os, sys
import pickle as pk
import AM_functions as am
from scipy.interpolate import interp1d
from scipy import interpolate
from datetime import *
from dateutil import parser
import scipy

raz=pA.ReadAzscan()

class struct(object):
    pass





x_am=am.AM_functions()

pf='../../../Postings/WVR_postings/20220110_BK_WVR_comparison_v2/plots/'

rx=210
det=20
rx_list = [210, 270]
det_list = [0,5,10,15,20,25,30]
res_150 = 1.5 #deg
rx_res = []
#rx_list = [210]

f = open('pointing_parameters_2018_fast.txt','rb')
point_par = pk.load(f)
f.close()
az_offs = point_par['az_offs']


def find_rx_detlist(tod):

    # while (np.min(az_bk)<0.):
    #     az_bk = [toaz_i+360. for toaz_i in az_bk]
    #     print('adding 360')

    rx_list = [30, 40, 210, 270]

    i_det_list_allrx={'30':[], '40':[], '210':[], '270':[]}

    i_det_list_allrx['30'] = [j for j in range (0, len(tod.ind.rgl030a))]
    i_det_list_allrx['40'] = [j for j in range (0, len(tod.ind.rgl040a))]
    i_det_list_allrx['210'] = [j for j in range (0, len(tod.ind.rgl210a))]
    i_det_list_allrx['270'] = [j for j in range (0, len(tod.ind.rgl270a))]

    return i_det_list_allrx


def pwv_to_Trj_maps(wvr_scan, az_bk, D_sum_bk, fs_bk):


    nscans_bk = len(D_sum_bk[10,:])
    az_bk = az_bk[int(fs_bk.sf[0]):int(fs_bk.ef[nscans_bk-1])]
    az_bk_shifted = [(az_bk_i + az_offs) for az_bk_i in az_bk]
    while (np.nanmin(az_bk_shifted)<0.):
        print(np.nanmin(az_bk_shifted))
        az_bk_shifted = [toaz_i+360. for toaz_i in az_bk_shifted]
        print('+ 360')
    while (np.nanmax(az_bk_shifted)>360.):
        print(np.nanmax(az_bk_shifted))
        az_bk_shifted = [toaz_i-360. for toaz_i in az_bk_shifted]
        print('- 360')

    #if not os.path.exists('am_datfiles_Az/SPole_annual_50/'+wvr_scan[:-4]+'/spectra/Trj_pk_BAK.txt'):
    D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, Az_min=np.nanmin(az_bk_shifted), Az_max=np.nanmax(az_bk_shifted), remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    # else:
    #     D30, D40, D90, D150, D210, D270, newaz = x_am.plot_Trj_az(wvr_scan, Az_min=np.nanmin(az_bk_shifted), Az_max=np.nanmax(az_bk_shifted), remake_post_fig=1, rewrite_txt=0, out_pk='Trj_pk_BAK.txt', posting_folder=pf)
    #


    wvr_atmo_Trj = {'30':D30, '40':D40, '90':D90, '150':D150, '210':D210, '270':D270}

    return wvr_atmo_Trj


def find_planck_to_rj_coeff():
    h=6.626176e-27 #erg*s
    k=1.380662e-16 #erg/K
    c=2.99792458e8 #m/s
    T_cmb=2.725 #K

    alpha = h/(k*T_cmb)

    freq_GHz=np.linspace(0, 300, 100000)
    freq_Hz=freq_GHz*(10**9)

    x = alpha*freq_Hz

    f_v = (x**2.)*np.exp(x)/((np.exp(x)-1)**2.)

    pl.figure(figsize=(12,8))
    pl.plot(freq_Hz, f_v)
    pl.title('T_rj/T_p')

    i210=np.where(freq_Hz>=210e9)[0][0]
    i270=np.where(freq_Hz>=270e9)[0][0]

    f_210=f_v[i210]
    f_270=f_v[i270]
    pl.axvline(x=270e9, c='r', label='270GHz')
    pl.axhline(y=f_270, c='r', label='f_270 ='+str(round(f_270,2)))
    pl.axvline(x=210e9, c='y', label='210GHz')
    pl.axhline(y=f_210, c='y', label='f_210 ='+str(round(f_210,2)))
    pl.legend()
    pl.savefig(pf+'PlanckCMBu_to_RJu_coeff.png')

    pl.close()

    f_v = {'210':f_210, '270':f_270}

    return f_v



def p3_filter(Df, x_diff, extent_wvr, xticks_wvr_mask, xlabels_wvr_mask, rx):
    Df_p3=np.full(np.shape(Df), np.nan)
    nscans=len(Df[10, :])
    for i in range(nscans-1):
        ts=Df[:, i]
        ts_p3_atmo=np.full(len(ts), np.nan)
        az=np.arange(len(ts))
        mask = np.isfinite(ts)
        ts_p3_coeff = np.polyfit(az[mask], ts[mask], deg=3)#, rcond=None, full=False, w=None, cov=False)
        ts_p3_f = np.poly1d(ts_p3_coeff)
        ts_p3 = ts_p3_f(az[mask])
        ts_p3_atmo[mask] = ts_p3
        Df_p3[:, i] = ts_p3_atmo

    fig,(ax1, ax2) = pl.subplots(2,1, figsize=(12,8))
    pos1=ax1.imshow(Df_p3, aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
    ax1.set_title('p3')
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('PWV[um]')
    ax1.set_ylabel('Az[deg]')
    ax1.set_ylim(np.min(x_diff), np.max(x_diff))
    ax1.set_xticks(xticks_wvr_mask[::6])
    ax1.set_xticklabels(xlabels_wvr_mask[::6])

    pos2=ax2.imshow(Df-Df_p3, aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
    ax2.set_title('Trj210 - p3 [residuals]')
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('PWV[um]')
    ax2.set_ylabel('Az[deg]')
    ax2.set_ylim(np.min(x_diff), np.max(x_diff))
    ax2.set_xticks(xticks_wvr_mask[::6])
    ax2.set_xticklabels(xlabels_wvr_mask[::6])
    pl.suptitle(wvr_scan[:8]+'\nTrj - '+str(rx)+'GHz', fontsize=14)

    pl.savefig(pf+wvr_scan[:-4]+'_Trj'+str(rx)+'_p3polyfit.png')

    return Df_p3




def BK_TOD(ets, rx, p3_filt, az_star, wvr_scan):

    tod = ets.bk_struct
    fs_bk = tod.fs
    calib_az = ets.wvr_struct.az_real #az_calib

    #defining time axis

    time_UTC=[]
    time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)
    time_UTC=[parser.parse(time_str) for time_str in time_UTC_str]
    fs_wvr = ets.wvr_struct.fs
    full_time_wvr=[time_UTC[int(fs_s_i)] for fs_s_i in fs_wvr.s]
    xlabels_wvr=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time_wvr])
    full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr[:-1]])

    nscans_wvr = len(fs_wvr.s)
    xticks_wvr = np.arange(nscans_wvr)

    t_bk=np.array(tod.std)
    full_time=[t_bk[int(fs_s_i)] for fs_s_i in fs_bk.sf]
    xlabels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time])
    nscans_bk = len(fs_bk.sf)
    xticks = np.arange(nscans_bk)
    full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])

    t_min_wvr=datetime.time(full_time_wvr[0])
    t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
    t_min_bk=datetime.time(full_time[0])
    t_max_bk=datetime.time(full_time[len(full_time)-1])

    bk_time_mask = np.where((full_time_dt>=t_min_wvr) & (full_time_dt<=t_max_wvr))[0]
    wvr_time_mask = np.where((full_time_wvr_dt<=t_max_bk) & (full_time_wvr_dt>=t_min_bk))[0]

    extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]

    #
    # if (t_min_bk <= t_min_wvr):
    #
    #     bk_time_mask = np.where((full_time_dt>=t_min_wvr) & (full_time_dt<=t_max_wvr))[0]
    #     wvr_time_mask = np.where((full_time_wvr_dt<=t_max_bk) & (full_time_wvr_dt>=t_min_bk))[0]
    #
    #     extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]
    #
    # else:
    #
    #     bk_time_mask = np.where(full_time_dt<=t_max_wvr)[0]
    #     wvr_time_mask = np.where(full_time_wvr_dt>=t_min_bk)[0]
    #
    #     extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]

    xlabels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time])
    xlabels_wvr=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time_wvr])
    xlabels=np.array(xlabels)
    xlabels_wvr=np.array(xlabels_wvr)
    xlabels_mask=xlabels[bk_time_mask]
    xticks_wvr_mask=xticks_wvr[wvr_time_mask]

    xticks_mask=xticks[bk_time_mask]
    xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

    #WVR tod
    #just to extract az_bk and D for shape

    x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')

    x_b = x_pair[np.where(det_pol=='b')]
    x_a = x_pair[np.where(det_pol=='a')]
    y_b = y_pair[np.where(det_pol=='b')]
    y_a = y_pair[np.where(det_pol=='a')]


    i_det=10 #doesn't matter

    a_det = det_a_list[i_det]
    i_det_b=np.where(x_b==x_a[i_det])[0]
    b_det = det_b_list[i_det_b[0]]

    az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[i_det], p3=p3_filt, i_det_savefig=i_det)
    #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=10, p3=p3_filt)

    x_sum, D_sum = raz.interpToImage_BK(az_bk, T_rx_psum_bk, fs_bk)
    x_diff, D_diff = raz.interpToImage_BK(az_bk, T_rx_pdiff_bk, fs_bk)

    wvr_atmo_Trj = pwv_to_Trj_maps(wvr_scan, az_bk, D_sum_bk, fs_bk)
    wvr_atmo_f = wvr_atmo_Trj[str(rx)]
    wvr_atmo_f_p3 = p3_filter(wvr_atmo_f, x_diff, extent_wvr, xticks_wvr_mask, xlabels_wvr_mask, rx)

    x_az_atmo=np.arange(len(wvr_atmo_f[:,1]))

    wvr_az_mask = np.where((x_az_atmo>=(np.nanmin(x_diff)+az_offs)) & (x_az_atmo<=(np.nanmax(x_diff)+az_offs)))[0]
    x_az_atmo_cut=x_az_atmo[wvr_az_mask]
    x_az_atmo_cal=[(xaz_i-az_offs) for xaz_i in x_az_atmo_cut]

    #bringing x_az_atmo_cal in the wanted range
    while (np.nanmin(x_az_atmo_cal)<0.):
        x_az_atmo_cal = np.array([toaz_i+360. for toaz_i in x_az_atmo_cal])
        print('+ 360')

    while (np.nanmax(x_az_atmo_cal)>360.):
        print(np.nanmax(x_az_atmo_cal))
        x_az_atmo_cal = np.array([toaz_i-360. for toaz_i in x_az_atmo_cal])
        print('- 360')

    wvr_atmo_matchbk_f = wvr_atmo_f
    wvr_atmo_matchbk_f = wvr_atmo_matchbk_f[:,wvr_time_mask]
    wvr_atmo_matchbk_f = wvr_atmo_matchbk_f[wvr_az_mask,:]

    wvr_atmo_matchbk_f_p3 = wvr_atmo_f_p3
    wvr_atmo_matchbk_f_p3 = wvr_atmo_matchbk_f_p3[:,wvr_time_mask]
    wvr_atmo_matchbk_f_p3 = wvr_atmo_matchbk_f_p3[wvr_az_mask,:]

    nscans_wvr=len(wvr_atmo_matchbk_f_p3[10,:])
    tod_wvr_f = np.zeros(nscans_wvr)
    tod_wvr_f_p3 = np.zeros(nscans_wvr)
    tod_wvr_scans=[]
    tod_wvr_scans_p3=[]

    mask_wvr_azstar=np.where(x_az_atmo_cal>az_star)[0]
    idx_wvr_azstar=mask_wvr_azstar[0]

    for iscan in range(len(wvr_atmo_matchbk_f[10,:])):
            scan_tod = wvr_atmo_matchbk_f[:,iscan]
            tod_wvr_scans.append(scan_tod)
            scan_tod_az = scan_tod[idx_wvr_azstar]

            tod_wvr_f[iscan] = scan_tod_az

            scan_tod_p3 = wvr_atmo_matchbk_f_p3[:,iscan]
            scan_tod_az_p3 = scan_tod_p3[idx_wvr_azstar]
            tod_wvr_scans_p3.append(scan_tod_p3)

            tod_wvr_f_p3[iscan] = scan_tod_az_p3

    #BK tod
    # i_det_list_allrx=find_rx_detlist(tod)
    # i_det_list = i_det_list_allrx[str(rx)]

    tod_bk_sum_list=[]
    tod_bk_diff_list=[]

    # print('number of det pairs = ', len(i_det_list))
    #
    # i_det_list_cut = i_det_list[::10]

    i_det_list=[]

    #for i_det_star in i_det_list_cut:
    for i_det_star in range (len(det_a_list)):

        i_det_list.append(i_det_star)

        print('p3_filt = ', p3_filt)

        a_det = det_a_list[i_det_star]
        i_det_b=np.where(x_b==x_a[i_det_star])[0]
        b_det = det_b_list[i_det_b[0]]
        az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[i_det_star], p3=p3_filt, i_det_savefig=i_det_star)
        #az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=i_det_star, p3=p3_filt)

        D_diff_bk = D_diff_bk[:,bk_time_mask]
        D_sum_bk = D_sum_bk[:,bk_time_mask]

        while (np.min(az_bk)<0.):
            az_bk = [toaz_i+360. for toaz_i in az_bk]
            print('adding 360')

        az_bk = np.array(az_bk)

        nscans_bk = len(D_sum_bk[10,:])
        tod_bk_sum = np.zeros(nscans_bk)
        tod_bk_diff = np.zeros(nscans_bk)

        tod_bk_sum_scans = []
        tod_bk_diff_scans = []

        for iscan in range(nscans_bk-1):
            az_bk_onescan = az_bk[int(fs_bk.sf[iscan]):int(fs_bk.ef[iscan])]
            x_az_atmo_bk=np.arange(np.min(az_bk_onescan), np.max(az_bk_onescan))
            mask_bk=np.where(x_az_atmo_bk>az_star)[0]

            idx_bk=mask_bk[0]

            scan_tod_bk_sum = D_sum_bk[:,iscan]
            scan_tod_bk_diff = D_diff_bk[:,iscan]
            tod_bk_sum_scans.append(scan_tod_bk_sum)
            tod_bk_diff_scans.append(scan_tod_bk_diff)
            scan_tod_az_bk_sum = scan_tod_bk_sum[idx_bk]
            scan_tod_az_bk_diff = scan_tod_bk_diff[idx_bk]
            tod_bk_sum[iscan] = scan_tod_az_bk_sum
            tod_bk_diff[iscan] = scan_tod_az_bk_diff

        tod_bk_sum_list.append(tod_bk_sum)
        tod_bk_diff_list.append(tod_bk_diff)

    return tod_bk_sum_list, tod_bk_diff_list, i_det_list, tod_wvr_f, tod_wvr_f_p3, full_time_dt, bk_time_mask, full_time_wvr_dt, wvr_time_mask, tod_bk_sum_scans, tod_bk_diff_scans, tod_wvr_scans, tod_wvr_scans_p3


def plot_TOD(rx, scale_f, p3_filt, Az_tod, bk_tag, std_threshold=2):

    f = open('TOD_'+str(rx)+'GHz_'+'p3'+str(p3_filt)+'_'+bk_tag+'_Az_star_'+str(Az_tod)+'.txt','rb')
    TOD_dict = pk.load(f)
    f.close()

    #trial_fn = 'TOD_270GHz_p3True_20200418B01_dk203_Az_star_250.txt'

    bk_tod_sum = TOD_dict['BK_tod_sum']
    bk_tod_diff = TOD_dict['BK_tod_diff']
    bk_detlist = TOD_dict['BK_detlist']
    wvr_tod = TOD_dict['WVR_tod']
    wvr_tod_p3 = TOD_dict['WVR_tod_p3']
    time_wvr = TOD_dict['WVR_time_axis']
    time_bk = TOD_dict['BK_time_axis']
    t_mask_wvr = TOD_dict['wvr_time_mask']
    t_mask_bk = TOD_dict['bk_time_mask']

    alpha=0.87 #dt_scan_bk/dt_scan_wvr

    # wvr_tod = tod_wvr_f[t_mask_wvr]
    # wvr_tod_p3 = tod_wvr_f_p3[t_mask_wvr]

    t_wvr = time_wvr[t_mask_wvr]
    t_bk = time_bk[t_mask_bk]
    t_bk_ticks=np.arange(len(t_bk))
    t_wvr_ticks=np.arange(len(t_wvr))*(1./alpha)
    t_wvr_labels=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in time_wvr[t_mask_wvr]])

    print(TOD_dict.keys())



    i_lab=0
    pl.figure(figsize=(12,8))
    for tod_bk_i in bk_tod_sum:
        tod_bk_rescaled = tod_bk_i - np.nanmean(tod_bk_i[:-1])
        if np.nanstd(tod_bk_rescaled)<std_threshold:
            if i_lab == 0:
                pl.plot(t_bk_ticks[:-1], tod_bk_rescaled[:-1], c='orange', label='bk - pair sum')
                i_lab=1
            else:
                pl.plot(t_bk_ticks[:-1], tod_bk_rescaled[:-1], c='orange')
            #pl.xticks(t_bk_ticks, t_bk)

    i_lab=0
    for tod_bk_i in bk_tod_diff:
        tod_bk_rescaled = tod_bk_i - np.nanmean(tod_bk_i[:-1])
        if np.nanstd(tod_bk_rescaled)<std_threshold:
            if i_lab == 0:
                pl.plot(t_bk_ticks[:-1], tod_bk_rescaled[:-1], c='cyan', label='bk - pair diff')
                i_lab=1
            else:
                pl.plot(t_bk_ticks[:-1], tod_bk_rescaled[:-1], c='cyan')
            #pl.xticks(t_bk_ticks, t_bk)


    tod_wvr_f_rescaled = wvr_tod - np.nanmean(wvr_tod)
    tod_wvr_f_rescaled_p3 = wvr_tod_p3 - np.nanmean(wvr_tod_p3)
    pl.plot(t_wvr_ticks, tod_wvr_f_rescaled*scale_f, linewidth = 2, c='r', label='wvr T[K_cmb]')
    pl.plot(t_wvr_ticks, tod_wvr_f_rescaled_p3*scale_f, linewidth = 2, c='blue', label='wvr T_p3[K_cmb]')
    pl.xticks(t_wvr_ticks[::20], t_wvr_labels[::20])
    pl.legend(fontsize=16)
    #pl.ylim(-1,2)
    pl.xticks(fontsize=14)
    pl.yticks(fontsize=14)
    #pl.xlabel('scan_n', fontsize=14)
    pl.ylabel('(T-T_scan_avg)[k]', fontsize=14)
    pl.suptitle(str(rx)+' GHz - TOD compared\nAz = '+str(Az_tod)+' deg\n'+'p3_filt = '+str(p3_filt), fontsize=16)
    pl.savefig(pf+'TOD_'+str(rx)+'GHz_p3'+str(p3_filt)+'_'+bk_tag+'_Az_star_'+str(az_star)+'.png')
    #pl.show()
    pl.close()


def pl_atmo_pwv(correlate_bk = False):


    wvr_scan1='20200415_140134_scanAz_fast.txt'
    bk_tag1='20200415D03_dk293'

    wvr_scan2='20200418_190135_scanAz_fast.txt'
    bk_tag2='20200418B01_dk203'


    wvr_scan_list=[wvr_scan1, wvr_scan2]
    bk_tag_list=[bk_tag1, bk_tag2]


    for i_tag in range(2):

        wvr_scan=wvr_scan_list[i_tag]
        bk_tag=bk_tag_list[i_tag]

        ets=bts.extract_ts(bk_tag, wvr_scan)

        tod = ets.bk_struct
        fs_bk = tod.fs

        for rx in rx_list:
            x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png')
            x_b = x_pair[np.where(det_pol=='b')]
            x_a = x_pair[np.where(det_pol=='a')]
            y_b = y_pair[np.where(det_pol=='b')]
            y_a = y_pair[np.where(det_pol=='a')]
            for i_det in range(len(det_a_list)):
                if (i_det < 20):
                    print('i_det=', i_det)
                    a_det = det_a_list[i_det]
                    i_det_b=np.where(x_b==x_a[i_det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    pl.figure(figsize=(10,6))
                    pl.scatter(x,y, c='k')
                    pl.scatter(x_a[i_det],y_a[i_det], s=150, marker='o', c='r', label='det A')
                    pl.scatter(x_b[i_det_b],y_b[i_det_b], s=150, marker='*', c='y', label='det B')
                    pl.xlabel('Az[deg]')
                    pl.ylabel('El[deg]')
                    pl.suptitle('FPU Map')
                    pl.legend(title='Selected Pair')
                    pl.title('Az offs = '+str(round(x_a[i_det],2))+' - El Offs = '+str(round(y_a[i_det],2)))
                    pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(i_det)+'.png')
                    pl.close()

                    for p3_filt in [True, False]:

                        f = open('pointing_parameters_2018_fast.txt','rb')
                        point_par = pk.load(f)
                        f.close()

                        wvr_atmo = ets.wvr_struct.D_pwv
                        az_wvr = ets.wvr_struct.az_real #az_calib
                        time_ordered_az=ets.wvr_struct.az_wvr
                        fs = ets.wvr_struct.fs
                        t_wvr=ets.wvr_struct.tot

                        print('i_det = ', (a_det, b_det))

                        az_bk, T_rx_psum_bk, T_rx_pdiff_bk, D_sum_bk, D_diff_bk = ets.pl_tod_atmo(ets.bk_tag, tod, rx, i_det=(a_det, b_det), az_offs_det=x_a[i_det], p3=p3_filt, i_det_savefig=i_det)

                        t_bk=np.array(tod.std)

                        waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan)#, az_lim=(np.min(x_az), np.max(x_az)))
                        time_UTC_str, D_bad, waz_bad, calib_data_Gavg_bad, FH_bad = raz.read_Az_fast(wvr_scan, clean_mod=0)

                        time_UTC=[]
                        time_UTC=[parser.parse(time_str) for time_str in time_UTC_str]
                        # for time_str in time_UTC_str:
                        #     time_UTC.append(parser.parse(time_str))

                        fs_wvr=fs
                        full_time_wvr=[time_UTC[int(fs_s_i)] for fs_s_i in fs_wvr.s]
                        xlabels_wvr=np.array(["{:d}:{:02d}".format(time_i.hour, time_i.minute) for time_i in full_time_wvr])

                        #pl.plot(fs_wvr.s, xlabels_wvr)

                        nscans_wvr = len(fs_wvr.s)
                        xticks_wvr = np.arange(nscans_wvr)
                        #pl.plot(xticks_wvr[::10], xlabels_wvr[::10])


                        x_sum, D_sum = raz.interpToImage_BK(az_bk, T_rx_psum_bk, fs_bk)
                        x_diff, D_diff = raz.interpToImage_BK(az_bk, T_rx_pdiff_bk, fs_bk)

                        fs_bk=tod.fs
                        full_time=[t_bk[int(fs_s_i)] for fs_s_i in fs_bk.sf]
                        xlabels=np.array(["{:d}:{:d}".format(time_i.hour, time_i.minute) for time_i in full_time])
                        #pl.plot(fs_bk.sf, xlabels)

                        nscans_bk = len(fs_bk.sf)
                        xticks = np.arange(nscans_bk)
                        #pl.plot(xticks[::10], xlabels[::10])

                        full_time_dt = np.array([datetime.time(ft_i) for ft_i in full_time])
                        full_time_wvr_dt = np.array([datetime.time(ft_i) for ft_i in full_time_wvr])

                        t_min_wvr=datetime.time(full_time_wvr[0])
                        t_max_wvr=datetime.time(full_time_wvr[len(full_time_wvr)-1])
                        t_min_bk=datetime.time(full_time[0])
                        t_max_bk=datetime.time(full_time[len(full_time)-1])

                        if (t_min_bk <= t_min_wvr):

                            print(t_min_wvr)
                            print(t_max_bk)

                            bk_time_mask = np.where(full_time_dt>=t_min_wvr)[0]
                            wvr_time_mask = np.where(full_time_wvr_dt<=t_max_bk)[0]

                            extent_wvr = [0, len(wvr_time_mask), np.min(calib_az), np.max(calib_az)]

                        else:

                            print(t_min_bk)
                            print(t_max_wvr)

                            bk_time_mask = np.where(full_time_dt<=t_max_wvr)[0]
                            wvr_time_mask = np.where(full_time_wvr_dt>=t_min_bk)[0]

                            extent_wvr = [wvr_time_mask[0], wvr_time_mask[len(wvr_time_mask)-1], np.min(calib_az), np.max(calib_az)]


                        fig,(ax1, ax2, ax3) = pl.subplots(3,1, figsize=(12,8))

                        xlabels=np.array(xlabels)
                        xlabels_wvr=np.array(xlabels_wvr)

                        xticks_mask=xticks[bk_time_mask]
                        xlabels_mask=xlabels[bk_time_mask]
                        xticks_wvr_mask=xticks_wvr[wvr_time_mask]
                        xlabels_wvr_mask=xlabels_wvr[wvr_time_mask]

                        pos1=ax1.imshow(D_sum_bk[:,bk_time_mask], aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_sum), np.max(x_sum)], origin='lower')
                        ax1.set_title('BK'+str(rx)+' Pair Sum')
                        cbar1 = pl.colorbar(pos1, ax=ax1)
                        cbar1.set_label('T[K]')
                        ax1.set_ylabel('Az[deg]')
                        ax1.set_xticks(xticks_mask[::10])
                        #ax1.set_xticklabels(xlabels_mask[::10])
                        ax1.set_xticklabels([])
                        #pl.yticks(fontsize=fs_ticks)

                        pos2=ax2.imshow(D_diff_bk[:,bk_time_mask], aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_diff), np.max(x_diff)], origin='lower')
                        ax2.set_title('BK'+str(rx)+' Pair Diff')
                        cbar2 = pl.colorbar(pos2, ax=ax2)
                        cbar2.set_label('T[K]')
                        ax2.set_ylabel('Az[deg]')
                        ax2.set_xticks(xticks_mask[::10])
                        #ax2.set_xticklabels(xlabels_mask[::10])
                        ax2.set_xticklabels([])
                        #pl.yticks(fontsize=fs_ticks)

                        pos3=ax3.imshow(wvr_atmo[:,wvr_time_mask], aspect='auto', interpolation='nearest', extent=extent_wvr, origin='lower')
                        ax3.set_title('PWV')
                        cbar3 = pl.colorbar(pos3, ax=ax3)
                        cbar3.set_label('PWV[um]')

                        #ax3.set_xlabel('Nscan')
                        ax3.set_ylabel('Az[deg]')
                        ax3.set_ylim(np.min(x_diff), np.max(x_diff))
                        ax3.set_xticks(xticks_wvr_mask[::6])
                        ax3.set_xticklabels(xlabels_wvr_mask[::6])
                        ax3.set_title('WVR - '+wvr_scan[:-16])

                        pl.suptitle(bk_tag+'\ngcp_idx[a-b]='+str(a_det+1)+'-'+str(b_det+1), fontsize=14)

                        pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(i_det)+'_p3_'+str(p3_filt)+'_WVR_'+wvr_scan[:-4]+'.png')
                        pl.savefig(pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(i_det)+'_p3_'+str(p3_filt)+'.png')

                        pl.close()

                        fn_corr = pf+'BK_'+bk_tag+'_rx_'+str(rx)+'_idet_'+str(i_det)+'_p3_'+str(p3_filt)+'_summ-diff_correlation.png'

                        if correlate_bk == True:
                            correlate_sumdiff(D_sum_bk, D_diff_bk, rx, a_det, b_det, fn_corr)



                    #
                    #
                    # #extracting TODs
                    # print('Saving TODs.')
                    #
                    # wvr_atmo_Trj = pwv_to_Trj_maps(wvr_scan, az_bk, D_sum_bk, tod.fs)
                    #
                    # az_star=250
                    #
                    # tod_bk_sum_list, tod_bk_diff_list, i_det_list, tod_wvr_f = BK_TOD(ets, az_star, D_sum_bk, x_sum, D_diff_bk, x_diff, bk_time_mask, wvr_atmo, calib_az, wvr_time_mask, wvr_atmo_Trj, xlabels_mask, xlabels_wvr_mask, p3_filt)
                    #
                    # TOD_dict = {'BK_tod_sum': tod_bk_sum_list, 'BK_tod_diff': tod_bk_diff_list, 'BK_detlist': i_det_list, 'BK_time_axis': full_time_dt[bk_time_mask],'WVR_tod': tod_wvr_f, 'wvr_time_axis': full_time_wvr_dt[wvr_time_mask]}
                    #
                    #
                    # f = open('TOD_'+str(rx)+'GHz_p3'+str(p3_filt)+'_'+bk_tag+'_Az_star_'+str(az_star)+'.txt','wb')
                    # pk.dump(TOD_dict, f)
                    # f.close()
                    #
                    # print('\n\n\nOne is done.\n\n\n')




def correlate_sumdiff(map_sum, map_diff, rx, a_shift, b_shift, fn_corr):

    corr = scipy.signal.correlate2d(map_sum, map_diff)
    fig,(ax) = pl.subplots(1,1, figsize=(12,8))
    #pos=ax.imshow(corr, aspect='auto', interpolation='nearest', extent=[bk_time_mask[0], bk_time_mask[len(bk_time_mask)-1], np.min(x_sum), np.max(x_sum)], origin='lower')
    pos=ax.imshow(corr, aspect='auto', interpolation='nearest', origin='lower')
    ax.set_title('BK Pair Sum - Pair Diff Correlation\nRx '+ str(rx)+'\ngcp_idx[a-b]='+str(a_shift+1)+'-'+str(b_shift+1))
    cbar = pl.colorbar(pos, ax=ax)
    #cbar.set_label('T[K]')
    #ax.set_ylabel('Az[deg]')
    #ax.set_ylim(np.min(x_diff), np.max(x_diff))
    #ax.set_xticks(xticks_mask[::6])
    #ax.set_xticklabels(xlabels_mask[::6])
    pl.savefig(fn_corr)

    pl.show()



wvr_scan1='20200415_140134_scanAz_fast.txt'
bk_tag1='20200415D03_dk293'

wvr_scan2='20200418_190135_scanAz_fast.txt'
bk_tag2='20200418B01_dk203'


wvr_scan_list=[wvr_scan1]#, wvr_scan1]
bk_tag_list=[bk_tag1]#, bk_tag1]


#to make BK-WVR atmograms
# for i_tag in range(1):
#     wvr_scan=wvr_scan_list[i_tag]
#     bk_tag=bk_tag_list[i_tag]
#     ets=bts.extract_ts(bk_tag, wvr_scan)
#     for rx in rx_list:
#         print('Rx=', rx)
#         # fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png'
#         # det_azoffs, det_eloffs = ets.det_fpu_location(rx, fn_save)
#         for p3_filt in [True, False]:
#             print('p3_filt=',p3_filt)
# pl_atmo_pwv()
#
# sys.exit()

#to plot TODs
#az_star=250
az_star_list1=[180]#, 200, 220, 240]
az_star_list2=[200, 220, 240, 250, 280]
f_v = find_planck_to_rj_coeff()

for i_tag in range(1):

    wvr_scan=wvr_scan_list[i_tag]
    bk_tag=bk_tag_list[i_tag]

    print(wvr_scan)

    ets=bts.extract_ts(bk_tag, wvr_scan)

    for az_star in az_star_list1:
        print('Az=',az_star)
        for rx in rx_list:
            print('Rx=', rx)
            scale_f = 1./f_v[str(rx)]
            for p3_filt in [True, False]:
                print('p3_filt=',p3_filt)
                tod_bk_sum_list, tod_bk_diff_list, i_det_list, tod_wvr_f,  tod_wvr_f_p3, full_time_bk, bk_time_mask, full_time_wvr, wvr_time_mask, tod_bk_sum_scans, tod_bk_diff_scans, tod_wvr_scans, tod_wvr_scans_p3 = BK_TOD(ets, rx, p3_filt, az_star, wvr_scan)

                TOD_dict = {'BK_tod_sum': tod_bk_sum_list, 'BK_tod_diff': tod_bk_diff_list, 'BK_detlist': i_det_list, 'BK_time_axis': full_time_bk, 'bk_time_mask': bk_time_mask, 'WVR_tod': tod_wvr_f, 'WVR_tod_p3': tod_wvr_f_p3, 'WVR_time_axis': full_time_wvr, 'wvr_time_mask': wvr_time_mask,
                            'BK_tod_sum_scans': tod_bk_sum_scans, 'BK_tod_diff_scans': tod_bk_diff_scans, 'WVR_tod_scans': tod_wvr_scans, 'WVR_tod_p3_scans': tod_wvr_scans_p3}

                f = open('TOD_'+str(rx)+'GHz_p3'+str(p3_filt)+'_'+bk_tag+'_Az_star_'+str(az_star)+'.txt','wb')
                pk.dump(TOD_dict, f)
                f.close()

                print('\n\n\nOne is done.\n\n\n')

                plot_TOD(rx, scale_f, p3_filt, az_star, bk_tag, std_threshold=3)
