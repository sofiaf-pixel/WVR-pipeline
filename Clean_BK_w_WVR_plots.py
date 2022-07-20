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
import scipy.stats

raz=pA.ReadAzscan()
x_am=am.AM_functions()


class struct(object):
    pass


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



def scale(im, nR, nC):
    number_rows = len(im)     # source number of rows
    number_columns = len(im[0])  # source number of columns
    return [[ im[int(number_rows * r / nR)][int(number_columns * c / nC)]
                 for c in range(nC)] for r in range(nR)]

def make_maps(D_folder, fn, bk_tag, wvr_scan, rx, det, p3_filt, gs_filt):


    f = open(D_folder+fn,'rb')
    D_dict= pk.load(f)
    f.close()
    #
    # fn_p3false='BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtFalse_groundsub_'+str(gs_filt)+'.txt'
    #
    # f = open(D_folder+fn_p3false,'rb')
    # D_dict_p3false= pk.load(f)
    # f.close()

    # a_det = D_dict_p3false['a_det']
    # b_det = D_dict_p3false['b_det']
    # D_sum = D_dict_p3false['D_sum']
    # D_diff = D_dict_p3false['D_diff']
    # D_polA = D_dict_p3false['D_polA']
    # D_polB = D_dict_p3false['D_polB']
    # wvr_atmo = D_dict_p3false['wvr_atmo']

    a_det = D_dict['a_det']
    b_det = D_dict['b_det']
    D_sum = D_dict['D_sum']
    D_diff = D_dict['D_diff']
    D_polA = D_dict['D_polA']
    D_polB = D_dict['D_polB']
    wvr_atmo = D_dict['wvr_atmo']


    #plot0: Maps psum vs WVR
    #p3 False
    fig,(ax1, ax4) = pl.subplots(2,1, figsize=(12,8))
    std_sum = np.nanstd(D_sum.D_sum_bk_tcut)

    pos1=ax1.imshow(D_sum.D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum\nstd = '+str(round(std_sum,3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    ax1.set_xticks(D_sum.xticks_mask[::10])
    ax1.set_xticklabels([])

    std_wvr = np.nanstd(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_p0)

    pos4=ax4.imshow(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_p0 , aspect='auto', interpolation='nearest', extent=wvr_atmo.extent_wvr_cut, origin='lower')
    ax4.set_title('WVR T_CMB '+str(rx)+'\nstd = '+ str(round(std_wvr,3)))
    cbar4 = pl.colorbar(pos4, ax=ax4)
    cbar4.set_label('T[K]')
    ax4.set_ylabel('Az[deg]')
    ax4.set_xticks(wvr_atmo.xticks_wvr_mask[::6])
    ax4.set_xticklabels(wvr_atmo.xlabels_wvr_mask[::6])
    pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4])

    pl.savefig(pf+'PairSum-WVR_comparison_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filt'+str(p3_filt)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    # #p3 True
    # fig,(ax1, ax4) = pl.subplots(2,1, figsize=(12,8))
    # p3_filt = True
    # std_sum_p3 = np.nanstd(D_sum_p3.D_sum_bk_tcut)
    #
    # pos1=ax1.imshow(D_sum_p3.D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_sum_p3.bk_time_mask[0], D_sum_p3.bk_time_mask[len(D_sum_p3.bk_time_mask)-1], np.min(D_sum_p3.x_sum), np.max(D_sum_p3.x_sum)], origin='lower')
    # ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum\nstd = '+str(round(std_sum_p3,3)))
    # cbar1 = pl.colorbar(pos1, ax=ax1)
    # cbar1.set_label('T[K]')
    # ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum_p3.xticks_mask[::10])
    # ax1.set_xticklabels([])
    #
    # pos4=ax4.imshow(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_p0 , aspect='auto', interpolation='nearest', extent=wvr_atmo.extent_wvr_cut, origin='lower')
    # ax4.set_title('WVR T_CMB '+str(rx)+'\nstd = '+ str(round(std_wvr,3)))
    # cbar4 = pl.colorbar(pos4, ax=ax4)
    # cbar4.set_label('T[K]')
    # ax4.set_ylabel('Az[deg]')
    # ax4.set_xticks(wvr_atmo.xticks_wvr_mask[::6])
    # ax4.set_xticklabels(wvr_atmo.xlabels_wvr_mask[::6])
    # pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4])
    # pl.savefig(pf+'PairSum-WVR_comparison_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filt'+str(p3_filt)+'_groundsub_'+str(gs_filt)+'.png')
    # pl.close()
    #




def make_maps_clean(D_folder, fn, bk_tag, wvr_scan, rx, det, gs_filt):

    fp3 = open(D_folder+fn,'rb')
    D_dict_p3true= pk.load(fp3)
    fp3.close()

    fn_p3false='BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtFalse_groundsub_'+str(gs_filt)+'.txt'

    f = open(D_folder+fn_p3false,'rb')
    D_dict_p3false= pk.load(f)
    f.close()

    a_det = D_dict_p3false['a_det']
    b_det = D_dict_p3false['b_det']
    D_sum = D_dict_p3false['D_sum']
    D_diff = D_dict_p3false['D_diff']
    D_polA = D_dict_p3false['D_polA']
    D_polB = D_dict_p3false['D_polB']
    wvr_atmo = D_dict_p3false['wvr_atmo']


    a_det_p3 = D_dict_p3true['a_det']
    b_det_p3 = D_dict_p3true['b_det']
    D_sum_p3 = D_dict_p3true['D_sum']
    D_diff_p3 = D_dict_p3true['D_diff']
    D_polA_p3 = D_dict_p3true['D_polA']
    D_polB_p3 = D_dict_p3true['D_polB']
    wvr_atmo_p3 = D_dict_p3true['wvr_atmo']



    wvr_atmo_Tplanck_tcut_azcut_resized = np.array(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_resized)
    D_sum_bk_tcut = np.array(D_sum.D_sum_bk_tcut)
    alpha, offset = np.polyfit(wvr_atmo_Tplanck_tcut_azcut_resized.ravel(), D_sum_bk_tcut.ravel(), 1)

    #alpha = 1.

    pl.scatter(wvr_atmo_Tplanck_tcut_azcut_resized.ravel(), D_sum_bk_tcut.ravel(), c='k', alpha=0.6)
    pl.plot(wvr_atmo_Tplanck_tcut_azcut_resized.ravel(), offset + alpha * wvr_atmo_Tplanck_tcut_azcut_resized.ravel(), c='r', label='alpha = '+str(alpha))
    pl.legend()
    pl.savefig(pf+'alpha_fit_diagnosticplot_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()

    #plot1: Maps psum_p3 vs psum_res
    fig,(ax1, ax2) = pl.subplots(2,1, figsize=(12,8))

    std_sum_p3 = np.nanstd(D_sum_p3.D_sum_bk_tcut)

    pos1=ax1.imshow(D_sum_p3.D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_sum_p3.bk_time_mask[0], D_sum_p3.bk_time_mask[len(D_sum_p3.bk_time_mask)-1], np.min(D_sum_p3.x_sum), np.max(D_sum_p3.x_sum)], origin='lower')
    ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum - p3 filtered\nstd = '+str(round(std_sum_p3,3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    ax1.set_xticks(D_sum_p3.xticks_mask[::10])
    ax1.set_xticklabels([])

    nan_mask = ~np.isnan(D_sum.D_sum_bk_tcut)
    D_sum_res = np.full(np.shape(D_sum.D_sum_bk_tcut), np.nan)
    D_sum_res[nan_mask] = D_sum.D_sum_bk_tcut[nan_mask]-(alpha*wvr_atmo_Tplanck_tcut_azcut_resized[nan_mask])
    std_sum_clean = np.nanstd(D_sum_res)

    pos2=ax2.imshow(D_sum_res, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    ax2.set_title('BK'+str(rx)+' WVR-Cleaned Pair Sum\nstd = '+str(round(std_sum_clean,3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4])
    pl.savefig(pf+'PairSum_p3-PairSum_wvrcleaned_comparison_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    #plot2: Maps psum vs psum_res
    fig,(ax1, ax2) = pl.subplots(2,1, figsize=(12,8))

    std_sum = np.nanstd(D_sum.D_sum_bk_tcut)

    pos1=ax1.imshow(D_sum.D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    ax1.set_xticks(D_sum.xticks_mask[::10])
    ax1.set_xticklabels([])

    pos2=ax2.imshow(D_sum_res, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    ax2.set_title('BK'+str(rx)+' WVR-Cleaned Pair Sum\nstd = '+str(round(std_sum_clean,3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4])
    pl.savefig(pf+'PairSum-PairSum_wvrcleaned_comparison_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    #plot3: Maps pdiff vs psum_res
    fig,(ax1, ax2) = pl.subplots(2,1, figsize=(12,8))

    std_diff = np.nanstd(D_diff.D_diff_bk_tcut)

    pos1=ax1.imshow(D_diff.D_diff_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_diff.bk_time_mask[0], D_diff.bk_time_mask[len(D_diff.bk_time_mask)-1], np.min(D_diff.x_diff), np.max(D_diff.x_diff)], origin='lower')
    ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Diff\nstd = '+str(round(std_diff,3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    ax1.set_xticks(D_diff.xticks_mask[::10])
    ax1.set_xticklabels([])

    pos2=ax2.imshow(D_sum_res, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    ax2.set_title('BK'+str(rx)+' WVR-Cleaned Pair Sum\nstd = '+str(round(std_sum_clean,3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4])
    pl.savefig(pf+'PairDiff-PairSum_wvrcleaned_comparison_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()

    return std_sum_clean, std_sum, std_diff, std_sum_p3




def make_ts(D_folder, fn, bk_tag, wvr_scan, rx, det, gs_filt):

    fp3 = open(D_folder+fn,'rb')
    D_dict_p3true= pk.load(fp3)
    fp3.close()

    fn_p3false='BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtFalse_groundsub_'+str(gs_filt)+'.txt'

    f = open(D_folder+fn_p3false,'rb')
    D_dict_p3false= pk.load(f)
    f.close()

    a_det = D_dict_p3false['a_det']
    b_det = D_dict_p3false['b_det']
    D_sum = D_dict_p3false['D_sum']
    D_diff = D_dict_p3false['D_diff']
    D_polA = D_dict_p3false['D_polA']
    D_polB = D_dict_p3false['D_polB']
    wvr_atmo = D_dict_p3false['wvr_atmo']


    a_det_p3 = D_dict_p3true['a_det']
    b_det_p3 = D_dict_p3true['b_det']
    D_sum_p3 = D_dict_p3true['D_sum']
    D_diff_p3 = D_dict_p3true['D_diff']
    D_polA_p3 = D_dict_p3true['D_polA']
    D_polB_p3 = D_dict_p3true['D_polB']
    wvr_atmo_p3 = D_dict_p3true['wvr_atmo']


    #for the moment I am just doing this for non-p3 filtered data

    std_sum = np.nanstd(D_sum.D_sum_bk_tcut)
    std_diff = np.nanstd(D_diff.D_diff_bk_tcut)
    std_wvr = np.nanstd(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_p0)

    wvr_atmo_Tplanck_tcut_azcut_resized = np.array(wvr_atmo.wvr_atmo_Tplanck_tcut_azcut_resized)
    D_sum_bk_tcut = np.array(D_sum.D_sum_bk_tcut)


    wvr_atmo_Tplanck_tcut_azcut_resized_fit = wvr_atmo_Tplanck_tcut_azcut_resized.ravel()
    D_sum_bk_tcut_fit = D_sum_bk_tcut.ravel()

    nan_mask = (~np.isnan(wvr_atmo_Tplanck_tcut_azcut_resized_fit)) & (~np.isnan(D_sum_bk_tcut_fit))

    alpha, offset = np.polyfit(wvr_atmo_Tplanck_tcut_azcut_resized_fit[nan_mask], D_sum_bk_tcut_fit[nan_mask], 1)

    nan_mask = ~np.isnan(D_sum.D_sum_bk_tcut)
    D_sum_res = np.full(np.shape(D_sum.D_sum_bk_tcut), np.nan)
    D_sum_res[nan_mask] = D_sum.D_sum_bk_tcut[nan_mask]-(alpha*wvr_atmo_Tplanck_tcut_azcut_resized[nan_mask])
    std_sum_clean = np.nanstd(D_sum_res)


    diff_ts = []
    sum_ts = []
    pwv_ts = []
    sum_m_pwv_ts = []

    for i in range(np.shape(D_diff.D_diff_bk_tcut)[1]):
        diff_ts.append(np.nanmean(D_diff.D_diff_bk_tcut[:,i]))
        sum_ts.append(np.nanmean(D_sum.D_sum_bk_tcut[:,i]))
        pwv_ts.append(np.nanmean(wvr_atmo_Tplanck_tcut_azcut_resized[:,i]))
        sum_m_pwv_ts.append(np.nanmean((D_sum.D_sum_bk_tcut-(alpha*np.array(wvr_atmo_Tplanck_tcut_azcut_resized)))[:,i]))


    pl.figure(figsize=(10,8))
    pl.plot(diff_ts, label='pdiff - std='+str(round(np.nanstd(diff_ts),3)))
    pl.plot(sum_ts, label='psum - std='+str(round(np.nanstd(sum_ts),3)))
    pl.plot(alpha*np.array(pwv_ts), label=str(round(alpha,2))+' * wvr - std='+str(round(np.nanstd(alpha*np.array(pwv_ts)),3)))
    pl.plot(sum_m_pwv_ts, label='psum-wvr - std='+str(round(np.nanstd(sum_m_pwv_ts),3)))
    pl.xlabel('scanN')
    pl.ylabel('T_cmb[K]')
    pl.legend()
    pl.suptitle('BK '+bk_tag+' - WVR '+wvr_scan[:-4]+'\nGCP idx a/b = '+str(a_det)+'/'+str(b_det))
    pl.savefig(pf+'psum-pdiff_az_average_ts_date_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtFalse_groundsub_'+str(gs_filt)+'.png')
    pl.close()







def make_std_histo(bk_tag, wvr_scan, alpha = 0):

    #alpha=0 --> fit for alpha
    #alpha=1 --> no scaling factor alpha

    D_folder = 'bk_wvr_clean_data/'
    std_fn = D_folder+'BK-WVR_maps_std_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'.txt'

    f = open(std_fn,'rb')
    std_dict = pk.load(f)
    f.close()

    std_psum_clean_210=std_dict['std_psum_clean_210']
    std_psum_210=std_dict['std_psum_210']
    std_pdiff_210=std_dict['std_pdiff_210']
    std_psum_p3_210=std_dict['std_psum_p3_210']

    std_psum_clean_270=std_dict['std_psum_clean_270']
    std_psum_270=std_dict['std_psum_270']
    std_pdiff_270=std_dict['std_pdiff_270']
    std_psum_p3_270=std_dict['std_psum_p3_270']

    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_psum_clean_210, bins=np.linspace(0, 0.2, 200), histtype='step', linewidth=2, color='k', facecolor='c', fill=True)
    # pl.suptitle('Rx210 PSum wvr subtracted std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_psum_clean_210)))
    # pl.savefig(pf+'PairSum_wvrcleaned_std_date_'+bk_tag[:8]+'_rx_210.png')
    # pl.close()
    #
    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_psum_210, bins=np.linspace(0, 0.4, 400), histtype='step', linewidth=2, color='k', facecolor='r', fill=True)
    # pl.suptitle('Rx210 Psum std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_psum_210)))
    # pl.savefig(pf+'PairSum_std_date_'+bk_tag[:8]+'_rx_210.png')
    # pl.close()
    #
    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_pdiff_210, bins=np.linspace(0, 0.1, 100), histtype='step', linewidth=2, color='k', facecolor='blue', fill=True)
    # pl.suptitle('Rx210 Pdiff std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_pdiff_210)))
    # pl.savefig(pf+'PairDiff_std_date_'+bk_tag[:8]+'_rx_210.png')
    # pl.close()
    #
    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_psum_p3_210, bins=np.linspace(0, 0.1, 100), histtype='step', linewidth=2, color='k', facecolor='orange', fill=True)
    # pl.suptitle('Rx210 Psum p3 filtered std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_psum_p3_210)))
    # pl.savefig(pf+'PairSum_p3_std_date_'+bk_tag[:8]+'_rx_210.png')
    # pl.close()

    std_min210=0.
    std_max210=1.

    std_min270=0.
    std_max270=1.2


    pl.figure(figsize=(8,6))
    n, x, _ = pl.hist(std_psum_clean_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='c', facecolor='c', fill=True, label='psum_wvrclean')
    n, x, _ = pl.hist(std_psum_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='r', facecolor='r', fill=True, label='psum')
    n, x, _ = pl.hist(std_pdiff_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='blue', facecolor='blue', fill=True, label='pdiff')
    pl.legend()
    pl.suptitle('Rx210')
    pl.xlabel('atmogram std')
    #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
    pl.savefig(pf+'PairSum_wvrcleaned_std_date_'+bk_tag[:8]+'_rx_210_alpha'+str(alpha)+'.png')
    pl.close()

    pl.figure(figsize=(8,6))
    n, x, _ = pl.hist(std_psum_clean_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='c', facecolor='c', fill=True, label='psum_wvrclean')
    n, x, _ = pl.hist(std_psum_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='r', facecolor='r', fill=True, label='psum')
    n, x, _ = pl.hist(std_pdiff_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='blue', facecolor='blue', fill=True, label='pdiff')
    pl.suptitle('Rx270')
    pl.xlabel('atmogram std')
    pl.legend()
    #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
    pl.savefig(pf+'PairSum_wvrcleaned_std_date_'+bk_tag[:8]+'_rx_270_alpha'+str(alpha)+'.png')
    pl.close()


    pl.figure(figsize=(8,6))
    n, x, _ = pl.hist(std_psum_clean_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='c', facecolor='c', fill=True, label='psum_wvrclean')
    n, x, _ = pl.hist(std_psum_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='r', facecolor='r', fill=True, label='psum')
    n, x, _ = pl.hist(std_pdiff_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='blue', facecolor='blue', fill=True, label='pdiff')
    n, x, _ = pl.hist(std_psum_p3_210, bins=np.linspace(std_min210, std_max210, 100), histtype='step', alpha = 0.5, linewidth=2, color='orange', facecolor='orange', fill=True, label='psum_p3')
    pl.legend()
    pl.suptitle('Rx210')
    pl.xlabel('atmogram std')
    #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
    pl.savefig(pf+'PairSum_wvrcleaned_std_date_wp3_'+bk_tag[:8]+'_rx_210_alpha'+str(alpha)+'.png')
    pl.close()

    pl.figure(figsize=(8,6))
    n, x, _ = pl.hist(std_psum_clean_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='c', facecolor='c', fill=True, label='psum_wvrclean')
    n, x, _ = pl.hist(std_psum_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='r', facecolor='r', fill=True, label='psum')
    n, x, _ = pl.hist(std_pdiff_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='blue', facecolor='blue', fill=True, label='pdiff')
    n, x, _ = pl.hist(std_psum_p3_270, bins=np.linspace(std_min270, std_max270, 100), histtype='step', alpha = 0.5, linewidth=2, color='orange', facecolor='orange', fill=True, label='psum_p3')
    pl.suptitle('Rx270')
    pl.xlabel('atmogram std')
    pl.legend()
    #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
    pl.savefig(pf+'PairSum_wvrcleaned_std_date_wp3_'+bk_tag[:8]+'_rx_270_alpha'+str(alpha)+'.png')
    pl.close()

    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_psum_270, bins=np.linspace(0, 0.4, 400), histtype='step', linewidth=2, color='k', facecolor='r', fill=True)
    # pl.suptitle('Rx270 Psum std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_psum_270)))
    # pl.savefig(pf+'PairSum_std_date_'+bk_tag[:8]+'_rx_270.png')
    # pl.close()
    #
    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_pdiff_270, bins=np.linspace(0, 0.1, 100), histtype='step', linewidth=2, color='k', facecolor='blue', fill=True)
    # pl.suptitle('Rx270 Pdiff std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_pdiff_270)))
    # pl.savefig(pf+'PairDiff_std_date_'+bk_tag[:8]+'_rx_270.png')
    # pl.close()
    #
    # pl.figure(figsize=(12,10))
    # n, x, _ = pl.hist(std_psum_p3_270, bins=np.linspace(0, 0.1, 100), histtype='step', linewidth=2, color='k', facecolor='orange', fill=True)
    # pl.suptitle('Rx270 Psum p3 filtered std')
    # pl.xlabel('std')
    # pl.title('avg = '+str(np.nanmean(std_psum_p3_270)))
    # pl.savefig(pf+'PairSum_p3_std_date_'+bk_tag[:8]+'_rx_270.png')
    # pl.close()
    #




pf='../../../Postings/WVR_postings/20220310_BK_WVR_datacleaning/plots/'

wvr_scan='20200921_150134_scanAz_fast.txt'
bk_tag='20200921B09_dk293'


wvr_scan1='20200415_140134_scanAz_fast.txt'
bk_tag1='20200415D03_dk293'

wvr_scan2='20200418_190135_scanAz_fast.txt'
bk_tag2='20200418B01_dk203'

wvr_scan3='20200723_190134_scanAz_fast.txt'
bk_tag3='20200723B09_dk203'

wvr_scan4='20200921_150134_scanAz_fast.txt'
bk_tag4='20200921B09_dk293'

wvr_scan_list=[wvr_scan4, wvr_scan3, wvr_scan2]
bk_tag_list=[bk_tag4, bk_tag3, bk_tag2]

# wvr_scan_list=[wvr_scan2]
# bk_tag_list=[bk_tag2]

#
# make_std_histo(bk_tag2, wvr_scan2)
# sys.exit()

# for iscan in range(len(wvr_scan_list)):
#
#     wvr_scan = wvr_scan_list[iscan]
#     bk_tag=bk_tag_list[iscan]
#
#     raz.read_pwvatmo(wvr_scan, savefig_fn_posting=pf+bk_tag[:8]+'_pwvatmo.png', show=0)
#
# sys.exit()

for iscan in range(len(wvr_scan_list)):

    wvr_scan = wvr_scan_list[iscan]
    bk_tag=bk_tag_list[iscan]

    D_folder = 'bk_wvr_clean_data/'

    std_psum_clean_210=[]
    std_psum_210=[]
    std_pdiff_210=[]
    std_psum_p3_210=[]

    std_psum_clean_270=[]
    std_psum_270=[]
    std_pdiff_270=[]
    std_psum_p3_270=[]

    #raz.read_pwvatmo(wvr_scan, savefig_fn_posting=pf+bk_tag[:8]+'_pwvatmo.png')

    for fn in os.listdir(D_folder):
        if fn[29:46] == bk_tag:

            rx=int(fn[82:85])

            if (fn[-8:-4]=='True'):
                gs_filt=True
                if (fn[-23:-19]=='True'):
                    p3_filt=True
                else:
                    p3_filt=False
            else:
                gs_filt=False
                if (fn[-24:-20]=='True'):
                    p3_filt=True
                else:
                    p3_filt=False



            try:
                det=int(fn[91:93])
            except:
                det=int(fn[91:92])


            #make_maps(D_folder, fn, bk_tag, wvr_scan, rx, det, p3_filt, gs_filt)

            if p3_filt == True:
                try:
                    psum_clean, psum, pdiff, psum_p3 = make_maps_clean(D_folder, fn, bk_tag, wvr_scan, rx, det, gs_filt)
                    #make_ts(D_folder, fn, bk_tag, wvr_scan, rx, det, gs_filt)
                except Exception as e:
                    print(e)
    #             if rx == 210:
    #                 std_psum_clean_210.append(psum_clean)
    #                 std_psum_210.append(psum)
    #                 std_pdiff_210.append(pdiff)
    #                 std_psum_p3_210.append(psum_p3)
    #
    #             elif rx == 270:
    #                 std_psum_clean_270.append(psum_clean)
    #                 std_psum_270.append(psum)
    #                 std_pdiff_270.append(pdiff)
    #                 std_psum_p3_270.append(psum_p3)
    #
    #
    # std_dict = {'std_psum_clean_210':std_psum_clean_210,'std_psum_210': std_psum_210, 'std_pdiff_210': std_pdiff_210, 'std_psum_p3_210': std_pdiff_210,
    #             'std_psum_clean_270':std_psum_clean_270,'std_psum_270': std_psum_270, 'std_pdiff_270': std_pdiff_270, 'std_psum_p3_270': std_pdiff_270}
    #
    # D_folder = 'bk_wvr_clean_data/'
    # std_fn = D_folder+'BK-WVR_maps_std_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_alpha1.txt'
    #
    # f = open(std_fn,'wb')
    # pk.dump(std_dict, f)
    # f.close()
    #
    # make_std_histo(bk_tag, wvr_scan, alpha = 1)
