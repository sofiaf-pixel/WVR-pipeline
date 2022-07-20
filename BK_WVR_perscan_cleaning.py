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
import scipy.signal
from matplotlib.pyplot import cm

raz=pA.ReadAzscan()
x_am=am.AM_functions()


class struct(object):
    pass


def startup(pf, wvr_scan, bk_tag, rx, det, gs_filt=False):

    D_folder = 'bk_wvr_clean_data/'

    fn = 'BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtFalse_groundsub_'+str(gs_filt)+'.txt'
    fn_p3 = 'BK-WVR_data_cleaning_date_bk_'+bk_tag+'_wvr_'+wvr_scan[:-4]+'_rx_'+str(rx)+'_idet_'+str(det)+'_p3filtTrue_groundsub_'+str(gs_filt)+'.txt'
    #these files are created by Clean_BK_w_WVR.py

    f = open(D_folder+fn,'rb')
    D_dict = pk.load(f)
    f.close()

    a_det = D_dict['a_det']
    b_det = D_dict['b_det']
    D_sum = D_dict['D_sum']
    D_diff = D_dict['D_diff']
    D_polA = D_dict['D_polA']
    D_polB = D_dict['D_polB']
    D_wvr = D_dict['wvr_atmo']

    f_p3 = open(D_folder+fn_p3,'rb')
    D_dict_p3 = pk.load(f_p3)
    f_p3.close()

    a_det_p3 = D_dict_p3['a_det']
    b_det_p3 = D_dict_p3['b_det']
    D_sum_p3 = D_dict_p3['D_sum']
    D_diff_p3 = D_dict_p3['D_diff']
    D_polA_p3 = D_dict_p3['D_polA']
    D_polB_p3 = D_dict_p3['D_polB']

    wvr_atmo_Tplanck_tcut_azcut_resized = np.array(D_wvr.wvr_atmo_Tplanck_tcut_azcut_resized)
    D_sum_bk_tcut = np.array(D_sum.D_sum_bk_tcut)
    D_sum_bk_p3_tcut = np.array(D_sum_p3.D_sum_bk_tcut)
    D_diff_bk_tcut = np.array(D_diff.D_diff_bk_tcut)
    D_diff_bk_p3_tcut = np.array(D_diff_p3.D_diff_bk_tcut)


    fig,(ax1, ax2, ax3) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(D_sum_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_sum.bk_time_mask[0], D_sum.bk_time_mask[len(D_sum.bk_time_mask)-1], np.min(D_sum.x_sum), np.max(D_sum.x_sum)], origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title('BK - Psum\nstd = '+str(round(np.nanstd(D_sum_bk_tcut),3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    # ax1.set_xticklabels([])

    pos2=ax2.imshow(D_diff_bk_tcut, aspect='auto', interpolation='nearest', extent=[D_diff.bk_time_mask[0], D_diff.bk_time_mask[len(D_diff.bk_time_mask)-1], np.min(D_diff.x_diff), np.max(D_diff.x_diff)], origin='lower')
    ax2.set_title('BK - Pdiff\nstd = '+str(round(np.nanstd(D_diff_bk_tcut),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pos3=ax3.imshow(wvr_atmo_Tplanck_tcut_azcut_resized, aspect='auto', interpolation='nearest', extent=D_wvr.extent_wvr_cut, origin='lower')
    ax3.set_title('WVR - T_cmb\nstd = '+str(round(np.nanstd(wvr_atmo_Tplanck_tcut_azcut_resized),3)))
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('rx ='+str(rx)+' - idet = '+str(det)+' - groundsub ='+str(gs_filt))
    pl.savefig(pf+'zero-atmograms_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    return wvr_atmo_Tplanck_tcut_azcut_resized, D_sum_bk_tcut, D_sum_bk_p3_tcut, D_diff_bk_tcut, D_diff_bk_p3_tcut, D_sum, D_diff, D_wvr





def filter_p3az(pf, wvr_atmo, bk_atmo, bk_atmo_p3, extent_bk, gcp_pair, var='Psum'):

    scan_len = np.shape(bk_atmo)[1]
    az_len = np.shape(bk_atmo)[0]

    p3_scaling = []
    bk_p3_az_map = np.full(np.shape(bk_atmo), np.nan)
    wvr_p3_az_map = np.full(np.shape(bk_atmo), np.nan)
    bk_wvrclean = np.full(np.shape(bk_atmo), np.nan)

    for i in range(scan_len):

        in_bk = np.array(bk_atmo[:,i] - np.nanmean(bk_atmo[:,i])) #remove avg per scan
        in_wvr = np.array(wvr_atmo[:,i] - np.nanmean(wvr_atmo[:,i]))

        bk_par = np.polyfit(np.arange(az_len), in_bk, 3)
        p_bk = np.poly1d(bk_par)
        bk_p3_model = p_bk(np.arange(az_len))
        wvr_par = np.polyfit(np.arange(az_len)[~np.isnan(in_wvr)], in_wvr[~np.isnan(in_wvr)], 3)
        p_wvr = np.poly1d(wvr_par)
        wvr_p3_model = p_wvr(np.arange(az_len))

        pl.figure(figsize=(10,8))
        pl.plot(in_bk, c='r', label='bk_data')
        pl.plot(p_bk(np.arange(az_len)), c='r', alpha=0.6, ls='--', label='bk_p3fit')
        pl.plot(in_wvr, c='blue', alpha=0.6, label='wvr_data')
        pl.plot(p_wvr(np.arange(az_len)), c='blue', alpha=0.6, ls='--', label='wvr_p3fit')
        pl.xlabel('Az_i')
        pl.legend()
        pl.title(var+'-WVR\n nscan = '+str(i))
        pl.suptitle('GCP_pair = '+str(gcp_pair))
        pl.savefig(pf+var+'-WVR_TOD_az_nscan'+str(i)+'_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
        pl.close()

        bk_p3_az_map[:,i] = bk_p3_model
        wvr_p3_az_map[:,i] = wvr_p3_model

        alpha, offset = np.polyfit(wvr_p3_model, bk_p3_model, 1)
        p3_scaling.append(alpha)


        pl.figure(figsize=(10,8))
        pl.plot(in_bk, c='r', label='bk_data')
        pl.plot(p_bk(np.arange(az_len)), c='r', alpha=0.6, ls='--', label='bk_p3fit')
        pl.plot(alpha*in_wvr, c='blue', alpha=0.6, label='wvr_data- scaled:\nalpha='+str(round(alpha,4)))
        pl.plot(alpha*p_wvr(np.arange(az_len)), c='blue', alpha=0.6, ls='--', label='wvr_p3fit - scaled')
        pl.plot(in_bk - alpha*p_wvr(np.arange(az_len)), c='c', label='bk_data - wvr_pmodel')
        pl.xlabel('Az_i')
        pl.legend()
        pl.title(var+'-WVR\n nscan = '+str(i))
        pl.suptitle('GCP_pair = '+str(gcp_pair))
        pl.savefig(pf+var+'-WVR_TOD_Nscan_i'+str(i)+'_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
        pl.close()


        bk_wvrclean[:,i] = (bk_atmo[:,i] - alpha*wvr_p3_model)




    fig,(ax3, ax2, ax1) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+'-WVR_p3Az\nstd = '+str(round(np.nanstd(bk_wvrclean),3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    # ax1.set_xticklabels([])

    pos2=ax2.imshow(bk_atmo, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+'\nstd = '+str(round(np.nanstd(bk_atmo),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+'_p3\nstd = '+str(round(np.nanstd(bk_atmo_p3),3)))
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('Filter = P3_az\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    pl.savefig(pf+var+'-WVR_p3az_atmo_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    pl.figure(figsize=(12,10))
    pl.plot(p3_scaling, c='k')
    pl.axhline(y=np.nanmean(p3_scaling), c='r', ls='--',label='avg ='+str(round(np.nanmean(p3_scaling),4)))
    pl.legend()
    pl.title('p3 scaling factor (bk/wvr)')
    pl.suptitle('GCP_pair = '+str(gcp_pair))
    pl.savefig(pf+var+'-WVR_p3az_scaling_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()

    return bk_wvrclean, bk_p3_az_map, wvr_p3_az_map, p3_scaling




def filter_p3az_double(pf, wvr_atmo, bk_atmo, bk_atmo_p3, x, calib_az, extent_bk, gcp_pair, fit_deg, var='Psum'):


    bk_wvrclean_p3az, bk_p3_map_az, wvr_p3_map_az, p3_scaling_az = filter_p3az(pf, wvr_atmo, bk_atmo, bk_atmo_p3, extent_bk, gcp_pair, var=var)


    bk_atmo = bk_wvrclean_p3az

    scan_len = np.shape(bk_atmo)[1]
    az_len = np.shape(bk_atmo)[0]

    p3_scaling = []
    bk_p3_az_map = np.full(np.shape(bk_atmo), np.nan)
    wvr_p3_az_map = np.full(np.shape(bk_atmo), np.nan)
    bk_wvrclean = np.full(np.shape(bk_atmo), np.nan)

    for i in range(scan_len):

        in_bk = np.array(bk_atmo[:,i] - np.nanmean(bk_atmo[:,i])) #remove avg per scan
        in_wvr = np.array(wvr_atmo[:,i] - np.nanmean(wvr_atmo[:,i]))

        bk_par = np.polyfit(np.arange(az_len), in_bk, 3)
        p_bk = np.poly1d(bk_par)
        bk_p3_model = p_bk(np.arange(az_len))
        wvr_par = np.polyfit(np.arange(az_len)[~np.isnan(in_wvr)], in_wvr[~np.isnan(in_wvr)], 3)
        p_wvr = np.poly1d(wvr_par)
        wvr_p3_model = p_wvr(np.arange(az_len))

        pl.figure(figsize=(10,8))
        pl.plot(in_bk, c='r', label='bk_data')
        pl.plot(p_bk(np.arange(az_len)), c='r', alpha=0.6, ls='--', label='bk_p3fit')
        pl.plot(in_wvr, c='blue', alpha=0.6, label='wvr_data')
        pl.plot(p_wvr(np.arange(az_len)), c='blue', alpha=0.6, ls='--', label='wvr_p3fit')
        pl.xlabel('Az_i')
        pl.legend()
        pl.title(var+'-WVR\n nscan = '+str(i))
        pl.suptitle('GCP_pair = '+str(gcp_pair))
        pl.savefig(pf+var+'-WVR_TOD_time+az_nscan'+str(i)+'_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
        pl.close()

        bk_p3_az_map[:,i] = bk_p3_model
        wvr_p3_az_map[:,i] = wvr_p3_model

        alpha, offset = np.polyfit(wvr_p3_model, bk_p3_model, 1)
        p3_scaling.append(alpha)


        pl.figure(figsize=(10,8))
        pl.plot(in_bk, c='r', label='bk_data')
        pl.plot(p_bk(np.arange(az_len)), c='r', alpha=0.6, ls='--', label='bk_p3fit')
        pl.plot(alpha*in_wvr, c='blue', alpha=0.6, label='wvr_data- scaled:\nalpha='+str(round(alpha,4)))
        pl.plot(alpha*p_wvr(np.arange(az_len)), c='blue', alpha=0.6, ls='--', label='wvr_p3fit - scaled')
        pl.plot(in_bk - alpha*p_wvr(np.arange(az_len)), c='c', label='bk_data - wvr_pmodel')
        pl.xlabel('Az_i')
        pl.legend()
        pl.title(var+'-WVR\n nscan = '+str(i))
        pl.suptitle('GCP_pair = '+str(gcp_pair))
        pl.savefig(pf+var+'-WVR_TOD_time+az_Nscan_i'+str(i)+'_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
        pl.close()


        bk_wvrclean[:,i] = (bk_atmo[:,i] - alpha*wvr_p3_model)




    fig,(ax3, ax2, ax1) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+'-WVR_p3Az\nstd = '+str(round(np.nanstd(bk_wvrclean),3)))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    # ax1.set_xticklabels([])

    pos2=ax2.imshow(bk_atmo, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+'\nstd = '+str(round(np.nanstd(bk_atmo),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+'_p3\nstd = '+str(round(np.nanstd(bk_atmo_p3),3)))
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('Filter = P3_az\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    pl.savefig(pf+var+'-WVR_time+az_atmo_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()


    pl.figure(figsize=(12,10))
    pl.plot(p3_scaling, c='k')
    pl.axhline(y=np.nanmean(p3_scaling), c='r', ls='--',label='avg ='+str(round(np.nanmean(p3_scaling),4)))
    pl.legend()
    pl.title('p3 scaling factor (bk/wvr)')
    pl.suptitle('GCP_pair = '+str(gcp_pair))
    pl.savefig(pf+var+'-WVR_time+az_scaling_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'.png')
    pl.close()

    return bk_wvrclean, bk_p3_az_map, wvr_p3_az_map, p3_scaling








def filter_ptime(pf, wvr_atmo, bk_atmo, bk_atmo_p3, x, calib_az, extent_bk, gcp_pair, fit_deg, var, sf_string, sfig = 0):

    az_len = np.shape(bk_atmo)[0]
    scan_len = np.shape(bk_atmo)[1]

    az_wvr_max = np.max(calib_az)
    az_wvr_min = np.min(calib_az)

    az_wvr = np.linspace(az_wvr_min, az_wvr_max, np.shape(wvr_atmo)[0])

    az_labels = [round(azi, 1) for azi in az_wvr]

    bk_p_t_map = np.full(np.shape(bk_atmo), np.nan)
    wvr_p_t_map = np.full(np.shape(wvr_atmo), np.nan)
    bk_wvrclean = np.full(np.shape(bk_atmo), np.nan)

    p_scaling = []

    bk_dataset=[]
    wvr_dataset=[]
    rho_list=[]
    for i in range(az_len):
        try:
            in_bk = bk_atmo[i,:] - np.nanmean(bk_atmo[i,:])
            in_wvr = wvr_atmo[i,:] - np.nanmean(wvr_atmo[i,:])

            rho = np.corrcoef(in_bk, in_wvr)
            rho_list.append(rho[0,1])

            param_p_t_bk = np.polyfit(np.arange(len(in_bk)),in_bk, fit_deg)
            p_t_bk = np.poly1d(param_p_t_bk)
            p_t_bk_model = p_t_bk(np.arange(scan_len))
            bk_p_t_map[i,:] = p_t_bk_model

            param_p_t_wvr = np.polyfit(np.arange(len(in_wvr)),in_wvr, fit_deg)
            p_t_wvr = np.poly1d(param_p_t_wvr)
            p_t_wvr_model = p_t_wvr(np.arange(scan_len))
            wvr_p_t_map[i,:] = p_t_wvr_model

            alpha, offset = np.polyfit(p_t_wvr_model, p_t_bk_model, 1)
            lin_corr = np.poly1d([alpha, offset])
            p_scaling.append(alpha)

            bk_dataset.append(in_bk)
            wvr_dataset.append(in_wvr)

            pl.figure(figsize=(12,8))

            pl.plot(in_bk, c='r', alpha=0.6, label='bk_data')
            pl.plot(p_t_bk_model, c='r', ls='--', label='bk_data_pfit_t')

            pl.plot(alpha*in_wvr, c='blue', alpha=0.6, label='wvr_data - scaled\nalpha = '+str(round(alpha,4)))
            pl.plot(alpha*p_t_wvr_model, c='blue', ls='--', label='wvr_pmodel:\nwvr_pfit - alpha scaled')

            pl.plot(in_bk - alpha*p_t_wvr_model, c='c', label='bk_data - wvr_pmodel')

            bk_wvrclean[i,:] = in_bk - alpha*p_t_wvr_model

            pl.legend()
            pl.xlabel('Nscan')
            pl.title(var+'-WVR\n Az = '+str(az_labels[i]))
            pl.suptitle('GCP_pair = '+str(gcp_pair))
            if sfig==1:
                pl.savefig(pf+var+'-WVR_TOD_Az_i'+str(i)+'_'+sf_string+'.png')
            pl.close()



        except Exception as e:
            print(e)


    wvr_dataset = np.array(wvr_dataset)
    bk_dataset = np.array(bk_dataset)

    alpha_all, offset_all = np.polyfit(wvr_dataset.ravel(), bk_dataset.ravel(), 1)
    lin_corr = np.poly1d([alpha_all, offset_all])

    pl.figure()
    pl.scatter(wvr_dataset, bk_dataset, alpha=0.5)
    pl.plot(wvr_dataset.ravel(), lin_corr(wvr_dataset.ravel()), ls='--', c='k', label='alpha='+str(round(alpha_all,2)))
    pl.legend()
    pl.close()




    for i in range(scan_len):
        bk_wvrclean[:,i] = np.array((bk_wvrclean[:,i]) - np.nanmean(bk_wvrclean[:,i])) #remove avg per Az scan

    maxs=np.array([np.max(bk_wvrclean), np.max(bk_atmo_p3)])
    c_max = np.max(maxs)

    mins=np.array([np.min(bk_wvrclean), np.min(bk_atmo_p3)])
    c_min = np.min(mins)

    x1 = np.nanstd(bk_wvrclean)
    x2 = np.nanstd(bk_atmo)
    x3 = np.nanstd(bk_atmo_p3)

    fig,(ax2, ax1, ax3) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+' - WVR filtered [time] - std = '+f"{x1:.2e}")
    if var == 'Pdiff':
        pos1.set_clim(c_min, c_max)
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    ax1.set_xticklabels([])

    pos2=ax2.imshow(bk_atmo, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+' - no filter - std = '+f"{x2:.2e}")
    if var == 'Pdiff':
        pos2.set_clim(c_min, c_max)
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels([])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+' - p3 filtered - std = '+f"{x3:.2e}")
    if var == 'Pdiff':
        pos3.set_clim(c_min, c_max)
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    ax3.set_xticks(D_sum.xticks_mask[::10])
    ax3.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('Filter = P_time\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    if sfig==1:
        pl.savefig(pf+var+'-WVR_pt_atmo_'+sf_string+'.png')
    pl.close()



    pl.figure(figsize=(12,10))
    pl.plot(p_scaling, c='k')
    pl.axhline(y=np.nanmean(p_scaling), c='r', ls='--',label='avg ='+str(round(np.nanmean(p_scaling),4)))
    pl.legend()
    pl.title('P scaling factor (bk/wvr)')
    pl.suptitle('GCP_pair = '+str(gcp_pair))
    if sfig==1:
        pl.savefig(pf+var+'-WVR_pt_scaling_'+sf_string+'.png')
    pl.close()



    fig,(ax1, ax3, ax2) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+' - WVR filtered [time] - std = '+str(round(np.nanstd(bk_wvrclean),3)))
    if var == 'Pdiff':
        pos1.set_clim(c_min, c_max)
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    ax1.set_xticklabels([])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+' - p3 filtered [time] - std = '+str(round(np.nanstd(bk_atmo_p3),3)))
    if var == 'Pdiff':
        pos3.set_clim(c_min, c_max)
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels([])

    pos2=ax2.imshow(bk_wvrclean - bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+'_wvrfilt - '+var+'_p3filt - std = '+str(round(np.nanstd(bk_wvrclean - bk_atmo_p3),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    if var == 'Pdiff':
        pos2.set_clim(c_min, c_max)
    ax2.set_ylabel('Az[deg]')
    ax2.set_xticks(D_sum.xticks_mask[::10])
    ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('WVR Filter vs P3 Filter\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    if sfig==1:
        pl.savefig(pf+var+'-WVRfilter_vs_p3filter_'+sf_string+'.png')
    pl.close()


    pl.figure(figsize=(12,10))
    pl.plot(p_scaling, c='k')
    pl.axhline(y=np.nanmean(p_scaling), c='r', ls='--',label='avg ='+str(round(np.nanmean(p_scaling),4)))
    pl.legend()
    pl.title('P scaling factor (bk/wvr)')
    pl.suptitle('GCP_pair = '+str(gcp_pair))
    if sfig==1:
        pl.savefig(pf+var+'-WVR_pt_scaling_'+sf_string+'.png')
    pl.close()


    return bk_wvrclean, bk_p_t_map, wvr_p_t_map, p_scaling, rho_list






def makeplots_ptime(pf, bk_wvrclean, p_scaling,  bk_atmo, bk_atmo_p3, x, calib_az, extent_bk, gcp_pair, fit_deg, sf_string, var):

    az_len = np.shape(bk_atmo)[0]
    scan_len = np.shape(bk_atmo)[1]

    az_wvr_max = np.max(calib_az)
    az_wvr_min = np.min(calib_az)

    az_wvr = np.linspace(az_wvr_min, az_wvr_max, np.shape(wvr_atmo)[0])

    az_labels = [round(azi, 1) for azi in az_wvr]


    fig,(ax2, ax1, ax3) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+' - WVR filtered [time]\nstd = '+str(round(np.nanstd(bk_wvrclean),3)))
    pos1.set_clim(np.min(bk_wvrclean), np.max(bk_wvrclean))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    # ax1.set_xticklabels([])

    pos2=ax2.imshow(bk_atmo, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+' - no filter\nstd = '+str(round(np.nanstd(bk_atmo),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+' - p3 filtered\nstd = '+str(round(np.nanstd(bk_atmo_p3),3)))
    pos3.set_clim(np.min(bk_wvrclean), np.max(bk_wvrclean))
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('Filter = P_time\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    pl.savefig(pf+var+'-WVR_pt_atmo_'+sf_string+'.png')
    pl.close()



    fig,(ax1, ax3, ax2) = pl.subplots(3,1, figsize=(12,8))

    pos1=ax1.imshow(bk_wvrclean, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    #ax1.set_title('GCP idx a/b = '+str(a_det)+'/'+str(b_det)+'\nrx'+str(rx)+' Pair Sum \nstd = '+str(round(std_sum,3)))
    ax1.set_title(var+' - WVR filtered [time]\nstd = '+str(round(np.nanstd(bk_wvrclean),3)))
    pos1.set_clim(np.min(bk_wvrclean), np.max(bk_wvrclean))
    cbar1 = pl.colorbar(pos1, ax=ax1)
    cbar1.set_label('T[K]')
    ax1.set_ylabel('Az[deg]')
    # ax1.set_xticks(D_sum.xticks_mask[::10])
    # ax1.set_xticklabels([])

    pos3=ax3.imshow(bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax3.set_title(var+' - p3 filtered [time]\nstd = '+str(round(np.nanstd(bk_atmo_p3),3)))
    pos3.set_clim(np.min(bk_wvrclean), np.max(bk_wvrclean))
    cbar3 = pl.colorbar(pos3, ax=ax3)
    cbar3.set_label('T[K]')
    ax3.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pos2=ax2.imshow(bk_wvrclean - bk_atmo_p3, aspect='auto', interpolation='nearest', extent=extent_bk, origin='lower')
    ax2.set_title(var+'_wvrcleaned - '+var+'_p3filt\nstd = '+str(round(np.nanstd(bk_wvrclean - bk_atmo_p3),3)))
    cbar2 = pl.colorbar(pos2, ax=ax2)
    cbar2.set_label('T[K]')
    ax2.set_ylabel('Az[deg]')
    # ax2.set_xticks(D_sum.xticks_mask[::10])
    # ax2.set_xticklabels(D_sum.xlabels_mask[::10])

    pl.suptitle('WVR Filter vs P3 Filter\nrx ='+str(rx)+' - GCP_pair = '+str(gcp_pair)+' - groundsub ='+str(gs_filt))
    pl.savefig(pf+var+'-WVRfilter_vs_p3filter_'+sf_string+'.png')
    pl.close()


    pl.figure(figsize=(12,10))
    pl.plot(p_scaling, c='k')
    pl.axhline(y=np.nanmean(p_scaling), c='r', ls='--',label='avg ='+str(round(np.nanmean(p_scaling),4)))
    pl.legend()
    pl.title('P scaling factor (bk/wvr)')
    pl.suptitle('GCP_pair = '+str(gcp_pair))
    pl.savefig(pf+var+'-WVR_pt_scaling_'+sf_string+'.png')
    pl.close()



    pl.figure(figsize=(8,6))
    n, x, _ = pl.hist(p_scaling, bins=np.linspace(np.min(p_scaling), np.max(p_scaling), 100), histtype='step', alpha = 0.5, linewidth=2, color='c', facecolor='c', fill=True)
    #pl.legend()
    pl.suptitle('Scale Factor\n'+var)
    pl.xlabel('alpha')
    #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
    pl.savefig(pf+var+'-scale_factor_histo_'+sf_string+'.png')
    pl.close()





#main

pf='../../../Postings/WVR_postings/20220315_BK_WVR_cleaning_III/plots/'

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

# wvr_scan_list=[wvr_scan3, wvr_scan4, wvr_scan2]
# bk_tag_list=[bk_tag3, bk_tag4, bk_tag2]

wvr_scan_list=[wvr_scan2]
bk_tag_list=[bk_tag2]

#
# wvr_scan = wvr_scan2
# bk_tag=bk_tag2

#rx = 270
#det = 14
#det = 15
#det = 3
#gs_filt = False

make_file = False

for iscan in range(len(wvr_scan_list)):

    wvr_scan = wvr_scan_list[iscan]
    bk_tag=bk_tag_list[iscan]

    ets=bts.extract_ts(bk_tag, wvr_scan)

    # for rx in [210, 270]:
    for rx in [270]:
        # for gs_filt in [True, False]:
        for gs_filt in [False]:

            if make_file == True:


                x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol, el_fpu_center= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png', el_lim=10)
                offset_from_bkel = 55 - el_fpu_center
                x_b = x_pair[np.where(det_pol=='b')]
                x_a = x_pair[np.where(det_pol=='a')]
                y_b = y_pair[np.where(det_pol=='b')]
                y_a = y_pair[np.where(det_pol=='a')]

                rho_map = np.zeros(len(det_a_list))

                scale_sum = []
                scale_diff = []

                Psum_p3filt_std = []
                Pdiff_p3filt_std = []

                Psum_wvrfilt_std = []
                Pdiff_wvrfilt_std = []

                print('len det list = ', len(det_a_list))

                for det in range(len(det_a_list[150:200])):

                    fdeg = 12
                    outfile = 'clean_data/D_clean_psum_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det_a_list[det])+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt'

                    # if not os.path.exists(outfile):

                    a_det = det_a_list[det]
                    i_det_b=np.where(x_b==x_a[det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    gcp_pair = [a_det, b_det]

                    D_clean_sum = struct();
                    D_clean_diff = struct();

                    wvr_atmo, bk_atmo_sum, bk_atmo_sum_p3, bk_atmo_diff, bk_atmo_diff_p3, D_sum, D_diff, D_wvr = startup(pf, wvr_scan, bk_tag, rx, det, gs_filt)

                    x_d = D_diff.x_diff
                    extent_bk = [D_diff.bk_time_mask[0], D_diff.bk_time_mask[len(D_diff.bk_time_mask)-1], np.min(D_diff.x_diff), np.max(D_diff.x_diff)]
                    waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, savefig_fn_posting=pf+bk_tag[:8]+'_pwvatmo.png', show=0)

                    #
                    D_clean_sum.map_x = x
                    D_clean_sum.map_y = y
                    D_clean_sum.x_a = x_a[det]
                    D_clean_sum.y_a = y_a[det]
                    D_clean_sum.gcp_idx = (a_det, b_det)

                    #
                    #
                    # #p3 filtering in az
                    # D_clean_sum.bk_wvrclean_p3az, D_clean_sum.bk_p3_map_az, D_clean_sum.wvr_p3_map_az, D_clean_sum.p3_scaling_az = filter_p3az(pf, wvr_atmo, bk_atmo_sum, bk_atmo_sum_p3, extent_bk, gcp_pair, var='Psum')
                    # D_clean_diff.bk_wvrclean_p3az, D_clean_diff.bk_p3_map_az, D_clean_diff.wvr_p3_map_az, D_clean_diff.p3_scaling_az = filter_p3az(pf, wvr_atmo, bk_atmo_diff, bk_atmo_diff_p3, extent_bk,gcp_pair, var='Pdiff')
                    #
                    #


                    # for fdeg in [3,4,6,8]:

                    savefigdata = bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)

                    #p3 filtering in time
                    D_clean_sum.bk_wvrclean_pt, D_clean_sum.bk_p_map_t, D_clean_sum.wvr_p_map_t, D_clean_sum.p_scaling_t, D_clean_sum.rho = filter_ptime(pf, wvr_atmo, bk_atmo_sum, bk_atmo_sum_p3, x_d, calib_az, extent_bk, gcp_pair, fit_deg = fdeg, var='Psum', sf_string=savefigdata, sfig = 1)
                    D_clean_diff.bk_wvrclean_pt, D_clean_diff.bk_p_map_t, D_clean_diff.wvr_p_map_t, D_clean_diff.p_scaling_t, D_clean_diff.rho = filter_ptime(pf, wvr_atmo, bk_atmo_diff, bk_atmo_diff_p3, x_d, calib_az, extent_bk, gcp_pair, fit_deg = fdeg, var='Pdiff', sf_string=savefigdata, sfig = 1)

                    # scale_diff = np.abs(D_clean_diff.p_scaling_t)
                    # scale_sum = np.abs(D_clean_sum.p_scaling_t)

                    rho_map[det]=np.nanmean(D_clean_diff.rho)

                    pl.figure(figsize=(10,6))
                    pl.scatter(x,y)
                    pl.scatter(x_a[det],y_a[det], s=100, marker='*', c=np.nanmean(D_clean_diff.rho))
                    pl.title('FPU Map\nGCP_idx='+str(a_det)+', '+str(b_det))
                    pl.colorbar()
                    pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'_rho_diff.png')
                    pl.close()

                    #
                    # #saving output in pk files
                    f = open('clean_data/D_clean_psum_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','wb')
                    pk.dump(D_clean_sum, f)
                    f.close()

                    f = open('clean_data/D_clean_pdiff_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','wb')
                    pk.dump(D_clean_diff, f)
                    f.close()

                    #loading from pk files
                    # f = open('clean_data/D_clean_psum_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','rb')
                    # D_clean_sum = pk.load(f)
                    # f.close()
                    #
                    # f = open('clean_data/D_clean_pdiff_'+bk_tag[:8]+'_rx_'+str(rx)+'_idet_'+str(det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','rb')
                    # D_clean_diff = pk.load(f)
                    # f.close()

                    scale_diff.append(np.abs(D_clean_diff.p_scaling_t))
                    scale_sum.append(np.abs(D_clean_sum.p_scaling_t))

                    Psum_p3filt_std.append(np.nanstd(bk_atmo_sum_p3))
                    Pdiff_p3filt_std.append(np.nanstd(bk_atmo_diff_p3))

                    Psum_wvrfilt_std.append(np.nanstd(D_clean_sum.bk_wvrclean_pt))
                    Pdiff_wvrfilt_std.append(np.nanstd(D_clean_diff.bk_wvrclean_pt))


            else:

                print('Reading pickle file')

                x, y, det_a_list, det_b_list, x_pair, y_pair, det_pol, el_fpu_center= ets.det_fpu_location(rx, fn_save = pf+bk_tag+'det_fpu_location_'+str(rx)+'.png', el_lim=10)
                offset_from_bkel = 55 - el_fpu_center
                x_b = x_pair[np.where(det_pol=='b')]
                x_a = x_pair[np.where(det_pol=='a')]
                y_b = y_pair[np.where(det_pol=='b')]
                y_a = y_pair[np.where(det_pol=='a')]

                rho_map = np.zeros(len(det_a_list))

                scale_sum = []
                scale_diff = []

                Psum_p3filt_std = []
                Pdiff_p3filt_std = []

                Psum_wvrfilt_std = []
                Pdiff_wvrfilt_std = []

                print('len det list = ', len(det_a_list))

                for a_det in (det_a_list[:40]):

                    fdeg = 12
                    mask_det = np.where(det_a_list == a_det)[0]
                    det=mask_det[0]

                    outfile = 'clean_data/D_clean_psum_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(a_det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt'

                    # if not os.path.exists(outfile):


                    i_det_b=np.where(x_b==x_a[det])[0]
                    b_det = det_b_list[i_det_b[0]]

                    gcp_pair = [a_det, b_det]

                    savefigdata = bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)

                    #
                    # #saving output in pk files
                    f = open('clean_data/D_clean_psum_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(a_det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','rb')
                    D_clean_sum=pk.load(f)
                    f.close()

                    f = open('clean_data/D_clean_pdiff_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(a_det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','rb')
                    D_clean_diff=pk.load(f)
                    f.close()

                    # pl.figure(figsize=(10,6))
                    # pl.scatter(x,y)
                    # pl.scatter(x_a[det],y_a[det], s=100, marker='*', c=np.nanmean(D_clean_diff.rho))
                    # pl.title('FPU Map\nGCP_idx='+str(a_det)+', '+str(b_det))
                    # pl.colorbar()
                    # pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'_rho_diff.png')
                    # pl.show()


                    wvr_atmo, bk_atmo_sum, bk_atmo_sum_p3, bk_atmo_diff, bk_atmo_diff_p3, D_sum, D_diff, D_wvr = startup(pf, wvr_scan, bk_tag, rx, det, gs_filt)

                    x_d = D_diff.x_diff
                    extent_bk = [D_diff.bk_time_mask[0], D_diff.bk_time_mask[len(D_diff.bk_time_mask)-1], np.min(D_diff.x_diff), np.max(D_diff.x_diff)]
                    waz, az, calib_az, fs, idx, pwv_ts, D_pwv, tot = raz.read_pwvatmo(wvr_scan, savefig_fn_posting=pf+bk_tag[:8]+'_pwvatmo.png', show=0)

                    rho_map[det]=np.nanmean(D_clean_diff.rho)
                    # for fdeg in [3,4,6,8]:

                    scale_diff.append(np.abs(D_clean_diff.p_scaling_t))
                    scale_sum.append(np.abs(D_clean_sum.p_scaling_t))

                    Psum_p3filt_std.append(np.nanstd(bk_atmo_sum_p3))
                    Pdiff_p3filt_std.append(np.nanstd(bk_atmo_diff_p3))

                    Psum_wvrfilt_std.append(np.nanstd(D_clean_sum.bk_wvrclean_pt))
                    Pdiff_wvrfilt_std.append(np.nanstd(D_clean_diff.bk_wvrclean_pt))



                rho_map_abs = [np.abs(rho_i) for rho_i in rho_map]

                fig = pl.figure(figsize=(10,6))
                ax = fig.add_subplot(111)
                pl.scatter(x,y)
                pl.scatter(x_a, y_a, s=100, marker='*', c=rho_map_abs)
                ax.fill_between(np.linspace(np.min(x), np.max(x), 1000), offset_from_bkel - 1.5, offset_from_bkel + 1.5, alpha=0.2)
                pl.title('FPU Map\nGCP_idx='+str(a_det)+', '+str(b_det))
                clb=pl.colorbar()
                clb.ax.set_title('BK-WVR correlation coeff')
                pl.savefig(pf+bk_tag+'det_fpu_location_'+str(rx)+'_det'+str(det)+'_rho_diff_map.png')
                pl.show()




                histo = struct();

                histo.scale_diff = scale_diff
                histo.scale_sum = scale_sum

                histo.Psum_p3filt_std = Psum_p3filt_std
                histo.Pdiff_p3filt_std = Pdiff_p3filt_std

                histo.Psum_wvrfilt_std = Psum_wvrfilt_std
                histo.Pdiff_wvrfilt_std = Pdiff_wvrfilt_std


                f = open('clean_data/histograms_'+bk_tag[:8]+'_rx_'+str(rx)+'_adet_'+str(det)+'_groundsub_'+str(gs_filt)+'_tfitdeg_'+str(fdeg)+'.txt','wb')
                pk.dump(histo, f)
                f.close()

                scale_sum = np.array(scale_sum)
                scale_diff = np.array(scale_diff)

                print('scale_sum =', scale_sum)

                ndet = len(scale_sum)
                color_list = cm.rainbow(np.linspace(0, 1, ndet))

                pl.figure(figsize=(8,6))
                for i in range(ndet):
                    color = color_list[i]
                    n, xh, _ = pl.hist(scale_sum[i], bins=np.linspace(np.min(scale_sum[i]), np.max(scale_sum[i]), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor=color, fill=True)
                #pl.axvline(x=np.nanmedian(scale_sum.ravel()), c='k', linestyle='--', label='median='+str(round(np.nanmedian(scale_sum.ravel()),3)))
                pl.legend()
                pl.suptitle('Scale Factor\nRx = '+str(rx)+' - PSum')
                pl.xlabel('alpha')
                #pl.xlim(np.nanmedian(scale_sum.ravel())-1, np.nanmedian(scale_sum.ravel())+1)
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Psum-scale_factor_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'_colordet.png')
                pl.show()

                pl.figure(figsize=(8,6))
                for i in range (ndet):
                    color = color_list[i]
                    n, xh, _ = pl.hist(scale_diff[i], bins=np.linspace(np.min(scale_diff[i]), np.max(scale_diff[i]), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor=color, fill=True)
                #pl.axvline(x=np.nanmedian(scale_diff), c='k', linestyle='--', label='median='+str(round(np.nanmedian(scale_diff),3)))
                pl.legend()
                #pl.xlim([np.nanmedian(scale_diff.ravel())-0.05, np.nanmedian(scale_diff.ravel())+0.05])
                pl.suptitle('Scale Factor\nRx = '+str(rx)+' - PDiff')
                pl.xlabel('alpha')
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Pdiff-scale_factor_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'_colordet.png')
                pl.show()





                #histograms for std
                pl.figure(figsize=(8,6))
                n, xh, _ = pl.hist(Psum_p3filt_std, bins=np.linspace(np.min(Psum_p3filt_std), np.max(Psum_p3filt_std), 200), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='r', fill=True)
                pl.axvline(x=np.nanmedian(Psum_p3filt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Psum_p3filt_std),3)))
                #pl.xlim([np.nanmedian(Psum_p3filt_std)-1, np.nanmedian(Psum_p3filt_std)+1])
                pl.legend()
                pl.suptitle('std Psum - p3 filtered \nRx = '+str(rx))
                pl.xlabel('std')
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Psum-p3_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
                pl.show()


                pl.figure(figsize=(8,6))
                n, xh, _ = pl.hist(Pdiff_p3filt_std, bins=np.linspace(np.min(Pdiff_p3filt_std), np.max(Pdiff_p3filt_std), 200), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='c', fill=True)
                pl.axvline(x=np.nanmedian(Pdiff_p3filt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Pdiff_p3filt_std),3)))
                pl.legend()
                #pl.xlim([np.nanmedian(Pdiff_p3filt_std)-1, np.nanmedian(Pdiff_p3filt_std)+1])
                pl.suptitle('std Pdiff - p3 filtered \nRx = '+str(rx))
                pl.xlabel('std')
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Pdiff-p3_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
                pl.show()


                pl.figure(figsize=(8,6))
                n, xh, _ = pl.hist(Psum_wvrfilt_std, bins=np.linspace(np.min(Psum_wvrfilt_std), np.max(Psum_wvrfilt_std), 200), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='r', fill=True)
                pl.axvline(x=np.nanmedian(Psum_wvrfilt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Psum_wvrfilt_std),3)))
                pl.legend()
                #pl.xlim([np.nanmedian(Psum_wvrfilt_std)-1, np.nanmedian(Psum_wvrfilt_std)+1])
                pl.suptitle('std Psum - wvr filtered \nRx = '+str(rx))
                pl.xlabel('std')
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Psum-wvr_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
                pl.show()


                pl.figure(figsize=(8,6))
                n, xh, _ = pl.hist(Pdiff_wvrfilt_std, bins=np.linspace(np.min(Pdiff_wvrfilt_std), np.max(Pdiff_wvrfilt_std), 200), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='c', fill=True)
                pl.axvline(x=np.nanmedian(Pdiff_wvrfilt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Pdiff_wvrfilt_std),3)))
                pl.legend()
                #pl.xlim([np.nanmedian(Pdiff_wvrfilt_std)-1, np.nanmedian(Pdiff_wvrfilt_std)+1])
                pl.suptitle('std Pdiff - wvr filtered \nRx = '+str(rx))
                pl.xlabel('std')
                #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
                pl.savefig(pf+'Pdiff-wvr_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
                pl.show()


                #try:

                scale_sum_median = [np.nanmedian(ss_i) for ss_i in scale_sum]
                scale_diff_median = [np.nanmedian(sd_i) for sd_i in scale_diff]



            sys.exit()



            #
            # pl.figure(figsize=(10,6))
            # pl.scatter(x,y)
            # pl.scatter(x_pair, y_pair, s=100, c=scale_sum_median)
            # pl.title('Scale Factor\nRx = '+str(rx)+' - PSum')
            # pl.colorbar()
            # pl.savefig(pf+'Psum-scale_factor_FPU_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()
            #
            #
            # pl.figure(figsize=(10,6))
            # pl.scatter(x,y)
            # pl.scatter(x_pair,y_pair, s=100, c=scale_diff_median)
            # pl.title('Scale Factor\nRx = '+str(rx)+' - PDiff')
            # pl.colorbar()
            # pl.savefig(pf+'Psum-scale_factor_FPU_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()

            #
            # #plot for alpha
            # ndet = len(scale_sum)
            # color_list = cm.rainbow(np.linspace(0, 10, ndet))
            #
            # pl.figure(figsize=(8,6))
            # for i in range(ndet):
            #     color = color_list[i]
            #     n_bins, x_bins, _ = pl.hist(scale_sum[i], bins=np.linspace(np.min(scale_sum[i]), np.max(scale_sum[i]), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor=color, fill=True)
            # pl.axvline(x=np.nanmedian(scale_sum[i]), c='k', linestyle='--', label='median='+str(round(np.nanmedian(scale_sum[i]),3)))
            # pl.legend()
            # pl.suptitle('Scale Factor\nRx = '+str(rx)+' - PSum')
            # pl.xlabel('alpha')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Psum-scale_factor_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'_colordet.png')
            # pl.show()
            #
            #
            # pl.figure(figsize=(8,6))
            # for i in range (ndet):
            #     color = color_list[i]
            #     n_bins, x_bins, _ = pl.hist(scale_diff[i], bins=np.linspace(np.min(scale_diff[i]), np.max(scale_diff[i]), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor=color, fill=True)
            # pl.axvline(x=np.nanmedian(scale_diff[i]), c='k', linestyle='--', label='median='+str(round(np.nanmedian(scale_diff[i]),3)))
            # pl.legend()
            # pl.suptitle('Scale Factor\nRx = '+str(rx)+' - PDiff')
            # pl.xlabel('alpha')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Pdiff-scale_factor_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'_colordet.png')
            # pl.show()
            #
            #
            # #histograms for std
            # pl.figure(figsize=(8,6))
            # n, x, _ = pl.hist(Psum_p3filt_std, bins=np.linspace(np.min(Psum_p3filt_std), np.max(Psum_p3filt_std), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='r', fill=True)
            # pl.axvline(x=np.nanmedian(Psum_p3filt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Psum_p3filt_std),3)))
            # pl.legend()
            # pl.suptitle('std Psum - p3 filtered \nRx = '+str(rx))
            # pl.xlabel('std')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Psum-p3_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()
            #
            #
            # pl.figure(figsize=(8,6))
            # n, x, _ = pl.hist(Pdiff_p3filt_std, bins=np.linspace(np.min(Pdiff_p3filt_std), np.max(Pdiff_p3filt_std), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='c', fill=True)
            # pl.axvline(x=np.nanmedian(Pdiff_p3filt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Pdiff_p3filt_std),3)))
            # pl.legend()
            # pl.suptitle('std Pdiff - p3 filtered \nRx = '+str(rx))
            # pl.xlabel('std')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Pdiff-p3_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()
            #
            #
            #
            #
            #
            # pl.figure(figsize=(8,6))
            # n, x, _ = pl.hist(Psum_wvrfilt_std, bins=np.linspace(np.min(Psum_wvrfilt_std), np.max(Psum_wvrfilt_std), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='r', fill=True)
            # pl.axvline(x=np.nanmedian(Psum_wvrfilt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Psum_wvrfilt_std),3)))
            # pl.legend()
            # pl.suptitle('std Psum - wvr filtered \nRx = '+str(rx))
            # pl.xlabel('std')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Psum-wvr_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()
            #
            #
            # pl.figure(figsize=(8,6))
            # n, x, _ = pl.hist(Pdiff_wvrfilt_std, bins=np.linspace(np.min(Pdiff_wvrfilt_std), np.max(Pdiff_wvrfilt_std), 20), histtype='step', alpha = 0.5, linewidth=2, color='k', facecolor='c', fill=True)
            # pl.axvline(x=np.nanmedian(Pdiff_wvrfilt_std), c='k', linestyle='--', label='median='+str(round(np.nanmedian(Pdiff_wvrfilt_std),3)))
            # pl.legend()
            # pl.suptitle('std Pdiff - wvr filtered \nRx = '+str(rx))
            # pl.xlabel('std')
            # #pl.title('avg = '+str(np.nanmean(std_psum_clean_270)))
            # pl.savefig(pf+'Pdiff-wvr_std_histo_'+bk_tag[:8]+'_rx_'+str(rx)+'_groundsub_'+str(gs_filt)+'_fitdeg_'+str(fdeg)+'.png')
            # pl.show()
            # #
            # # except Exception as e:
            # #
            # #     print(e)
