import plots_for_paper as pfp
import pickle as pk
import os, sys
import numpy as np
import pylab as pl
import plotAzscan as pAz
import AM_functions as am
import matplotlib.gridspec as gridspec


x_paper=pfp.make_plots()
x_am=am.AM_functions()
x_az=pAz.ReadAzscan()


years_list=[] #to collect in on eplot points for a full year
months_list=[] #to collect in one plot points for a full month
days_list=[] #to collect in one plot points for a full day

months_range=range(1,13)
#months_range=[6,7,8]
day_range=range(1,31)
#day_range=[10,11,12]

year='2020'
years_list.append(year)
for month in months_range:
    if month < 10:
        month_str='0'+str(month)
    else:
        month_str=str(month)
    months_list.append(year+month_str)
    for day in day_range:
        if day < 10:
            day_str='0'+str(day)
        else:
            day_str=str(day)
        days_list.append(year+month_str+day_str)


if sys.argv[1]=='zenith_pwv':

    if sys.argv[2]=='w':
        #write
        year='2020'
        #month_list=['05', '06', '08', '09', '10', '11', '12']
        month_list=['08']
        failed_months=[]

        for month_string in month_list:
            #try:
            date_string=year+month_string
            print('Extracting Zenith PWV for '+date_string)
            x_paper.extract_zenith_PWV(date=date_string, out_folder='paper_plots/', posting_folder='None')
            # except:
            #     failed_months.append(month_string)
        print('Months that failed: ', failed_months)

    if sys.argv[2]=='r':
        target_date_list=['202001', '202002', '202003', '202004', '202007']
        #read
        x_paper.extract_monthly_PS(target_date_list)



#to produce scatter plots
elif sys.argv[1]=='scatter_plots':

    single_day_list=['20200822', '20200823', '20200824']
    single_months_list=['202008', '202009']

    date_list_list=[days_list, months_list, years_list]

    for date_list in date_list_list:
        pf_scatter='../../Postings/WVR_postings/20211019_Zenith_T_channels_correlation/plots/'
        x_paper.make_ch_scatterplots_fit(date_list, pk_fn='paper_plots/zenithT_LinFit_'+date_list[0]+'-'+date_list[len(date_list)-1]+'_pk.txt', posting_folder=pf_scatter, write_pk=1)


#to make skyDip PWV plots
elif sys.argv[1]=='skydip_pwv':
    if sys.argv[2]=='w':
        pf_skydip='../../Postings/WVR_postings/20210625_SkyDips_Template_Comparison/plots/'
        date=sys.argv[3]
        x_paper.extract_PWV_skydip(year=date[:4], month=date[4:6], day=date[6:8], posting_folder=pf_skydip)

    elif sys.argv[2]=='r':
        date=sys.argv[3]
        filename=date+'_140002_skyDip_fast.txt'
        temp_merra = 'MERRA_'+date[:4]+'-'+date[4:6]+'-'+date[6:8]+'_14.amc'
        temp_standard = 'SPole_annual_50.amc'
        el_list_merra, pwv_list_merra=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template=temp_merra, spline=2, show=0)
        el_list, pwv_list=x_am.plot_am_fit_2(filename, var='pwv', pwv_layer='tropo', template=temp_standard, spline=2, show=0)


        cut_start=20
        cut_end=70

        el_list_merra, pwv_list_merra = el_list_merra[cut_start:-cut_end], pwv_list_merra[cut_start:-cut_end]
        el_list, pwv_list = el_list[cut_start:-cut_end], pwv_list[cut_start:-cut_end]


        fig = pl.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = pl.subplot(gs[0])
        ax2 = pl.subplot(gs[1])

        ax1.plot(el_list_merra, pwv_list_merra, color='k', alpha=0.2, label='Daily Template - MERRA data\nPWV_mean[um]='+str(int(np.nanmean(pwv_list_merra))))
        ax1.scatter(el_list_merra, pwv_list_merra, color='k', alpha=0.2, s=40, marker='o')
        ax1.axhline(y=np.nanmean(pwv_list_merra), color='k', ls='--')
        ax1.plot(el_list, pwv_list, color='red', alpha=0.6, label='Yearly Template - SPole 50 percentile\nPWV_mean[um]='+str(int(np.nanmean(pwv_list))))
        ax1.scatter(el_list, pwv_list, color='red', s=40, marker='o')
        ax1.axhline(y=np.nanmean(pwv_list), color='red', ls='--')

        print('el_list_merra =', el_list_merra )
        print('el_list = ', el_list)



        # ax1.scatter(el_list_merra[len(el_list_merra)-1], pwv_list_merra[len(el_list_merra)-1], color='c', s=70, marker='*')
        # ax1.scatter(el_list[len(el_list)-1], pwv_list[len(el_list)-1], color='c', s=70, marker='*')

        #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-3., 'pwv_mean[um]={:.2f}'.format(np.mean(pwv_list))+'±{:.2f}'.format(np.std(pwv_list)), fontsize=16)
        #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-6., 'pwv_mean_up[um]={:.2f}'.format(np.mean(pwv_up))+'±{:.2f}'.format(np.std(pwv_up)), color='r', fontsize=16)
        #pl.text(np.max(El_list)-30.,  np.max(pwv_list)-9., 'pwv_mean_dn[um]={:.2f}'.format(np.mean(pwv_dn))+'±{:.2f}'.format(np.std(pwv_dn)), color='b', fontsize=16)
        #pl.axhline(y=np.mean(pwv_list), c=col,alpha=0.6, linestyle='--', label=('pwv_avg='+str(round(np.mean(pwv_list),2))))

        ax1.set_ylabel('pwv[um]', fontsize=14)
        ax1.legend(loc='upper right', fontsize=14)
        ax1.set_xticklabels('')
        ax1.tick_params(axis='both', labelsize=14)
        #ax1.set_title(filename[:-4], fontsize=16)#+'\nTemplate: '+template[:-4], fontsize=16)
        ax2.set_ylabel('pwv[um]', fontsize=14)
        ax2.plot(el_list_merra, (pwv_list_merra-pwv_list[:-1]), '.', color='k', alpha=0.6, label='SPole Yearly Template')
        ax2.set_title('PWV_dailytemp - PWV_yearlytemp', fontsize=16)
        ax2.set_xlabel('El[deg]', fontsize=14)
        ax2.tick_params(axis='both', labelsize=14)
        #pl.savefig(pf2+filename[:-4]+'_fitoutput_pwv_'+pwv_layer+'_wlineshape_corr.png')
        pl.show()



#to plot power spectra from Az scans
elif sys.argv[1]=='azscan_ps':
    pf_ps='../../Postings/WVR_postings/20211025_PowerSpectra/plots/'
    data_folder='../../wvr1_data_local/'

    amfit_data_folder='am_datfiles_Az/SPole_annual_50/'

    test_file0='20200418_140135_scanAz_fast.txt'
    test_file1='20200101_140135_scanAz_fast.txt'

    #x_paper.plot_WVR_ps(test_file=test_file1, posting_folder=pf_ps)
    failed_scans=[]
    no_skydip_files=[]

    for flist in os.listdir(amfit_data_folder):
        pickle_fn='am_datfiles_Az/SPole_annual_50/'+flist+'/'+flist+'_clean_mod3_method_import_model_fitoutput.txt'
        if os.path.exists(pickle_fn):
            try:
                x_paper.plot_WVR_ps(flist+'.txt', posting_folder=pf_ps)
            except:
                failed_scans.append(flist)
                print(flist+' failed.')
        else:
            no_skydip_files.append(flist)


    print('List of files that failed:', failed_scans)
    print('No SkyDip Files:', no_skydip_files)


# failed_scans_to_relaunch= ['20200418_130135_scanAz_fast']
#
# for file in failed_scans_to_relaunch:
#     pickle_fn='am_datfiles_Az/SPole_annual_50/'+flist+'/'+flist+'_clean_mod3_method_import_model_fitoutput.txt'
#     x_paper.plot_WVR_ps(flist+'.txt', posting_folder=pf_ps)


#to make PWV atmo
elif sys.argv[1]=='azscan_pwv':
    if sys.argv[2]=='w':
        month_list=['202001']
        for month_str in month_list:
            x_paper.extract_PWV_atmo(date=month_str, posting_folder='None')
    elif sys.argv[2]=='r':
        date=sys.argv[3]
        try:
            x_az.read_pwvatmo(date+'_140135_scanAz_fast.txt')
        except:
            x_az.read_pwvatmo(date+'_140134_scanAz_fast.txt')




elif sys.argv[1]=='slopes':
    if sys.argv[2]=='w':
        x_paper.extract_PWV_skydip_slopes(days_list)
    elif sys.argv[2]=='r':

        f = open('SkyDips_PWV_slopes_2020.txt','rb')
        slopes=pk.load(f)
        f.close()

        date=np.array(slopes['date'])
        annual_lowel=np.array(slopes['SPole_annual_LowEl'])
        annual_highel=np.array(slopes['SPole_annual_HighEl'])
        merra_lowel=np.array(slopes['MERRA_LowEl'])
        merra_highel=np.array(slopes['MERRA_HighEl'])
        pwv_avg=np.array(slopes['PWV_avg'])
        pwv_avg_merra=np.array(slopes['PWV_avg_MERRA'])


        xticks_pos=[]
        xticks_lab=[]
        count=0
        for flist in date:
            if flist[4:6] != month:
                month=flist[4:6]
                xticks_pos.append(count)
                xticks_lab.append(flist[4:8])
            count+=1

        fig, ax = pl.subplots(2,1, figsize=(12,10))

        ax[0].scatter(date, annual_lowel, c='b', label='El<45')
        ax[0].plot(date, annual_lowel, c='k', alpha=0.5)

        ax[0].scatter(date, annual_highel, c='r', label='El>55')
        ax[0].plot(date, annual_highel, c='k', alpha=0.5)
        ax[0].set_ylabel('d(PWV)/d(EL)[um/deg]')
        ax[0].legend()
        ax[0].set_title('SPole_annual_50')
        ax[0].set_xticks(ticks=[])

        ax[1].scatter(date, merra_lowel, c='b', label='El<45')
        ax[1].plot(date, merra_lowel, c='k', alpha=0.5)

        ax[1].scatter(date, merra_highel, c='r', label='El>55')
        ax[1].plot(date, merra_highel, c='k', alpha=0.5)
        ax[1].set_ylabel('d(PWV)/d(EL)[um/deg]')
        ax[1].legend()
        ax[1].set_title('MERRA_daily')

        pl.suptitle('SkyDips PWV Slopes\n2020')

        pl.xticks(ticks=xticks_pos, labels=xticks_lab)

        pl.savefig('paper_plots/SkyDips_PWV_Slopes_2020.png')

        pl.show()


        fig, ax = pl.subplots(2,1, figsize=(12,10))

        coef_low=np.polyfit(pwv_avg, annual_lowel, 1)
        fit_low=np.poly1d(coef_low)
        pwv_fit_low=fit_low(pwv_avg)

        coef_high=np.polyfit(pwv_avg, annual_highel, 1)
        fit_high=np.poly1d(coef_high)
        pwv_fit_high=fit_high(pwv_avg)

        ax[0].scatter(pwv_avg, annual_lowel, c='c', alpha='0.6', label='El<45')
        ax[0].plot(pwv_avg, pwv_fit_low, c='c', ls='--', label='LinFit-Slope='+str(round(coef_low[0],4)))
        ax[0].scatter(pwv_avg, annual_highel, c='orange', alpha='0.6', label='El>55')
        ax[0].plot(pwv_avg, pwv_fit_high, c='orange', ls='--', label='LinFit-Slope='+str(round(coef_high[0],4)))
        ax[0].legend()
        ax[0].set_title('SPole_annual_50')
        ax[0].set_xticks(ticks=[])
        ax[0].set_ylabel('d(PWV)/d(EL)[um/deg]')

        coef_merra_low=np.polyfit(pwv_avg_merra, merra_lowel, 1)
        fit_merra_low=np.poly1d(coef_merra_low)
        pwv_fit_merra_low=fit_merra_low(pwv_avg_merra)

        coef_merra_high=np.polyfit(pwv_avg_merra, merra_highel, 1)
        fit_merra_high=np.poly1d(coef_merra_high)
        pwv_fit_merra_high=fit_merra_high(pwv_avg_merra)

        ax[1].scatter(pwv_avg_merra, merra_lowel, c='c', alpha='0.6', label='El<45')
        ax[1].plot(pwv_avg_merra, pwv_fit_merra_low, c='c', ls='--', label='LinFit-Slope='+str(round(coef_merra_low[0],4)))
        ax[1].scatter(pwv_avg_merra, merra_highel, c='orange',  alpha='0.6', label='El>55')
        ax[1].plot(pwv_avg_merra, pwv_fit_merra_high, c='orange', ls='--', label='LinFit-Slope='+str(round(coef_merra_high[0],4)))
        ax[1].set_ylabel('d(PWV)/d(EL)[um/deg]')
        ax[1].set_xlabel('PWV[um]')
        ax[1].legend()
        ax[1].set_title('MERRA_daily')

        pl.suptitle('Slopes vs Skydip PWV avg\n2020')


        pl.savefig('paper_plots/SkyDips_PWV_Slopes_2020_vs_PWVavg.png')

        pl.show()




else:
    print('No option selected.')
