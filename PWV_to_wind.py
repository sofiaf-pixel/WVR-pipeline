import pylab as pl
import numpy as np
import Read_BICEP_ts as bts
import plotAzscan as pA
import numpy as np
import pylab as pl
import os, sys
import pickle as pk
import AM_functions as am
import math
import scipy.interpolate as interp
import scipy.optimize as sp
import math
from scipy.optimize import curve_fit,least_squares
import scipy.stats as stats


raz=pA.ReadAzscan()
ets=bts.extract_ts()
x_am=am.AM_functions()



class extract_wind(object):

    def __init__(self, wvr_scan, pf, show, unit=None, verb=True):

        '''
        '''
        self.wvr_scan=wvr_scan
        self.pf=pf
        self.show=show
        self.el=55. #degrees

    def sim_sky(self, size=(5,20), x_c=4., y_c=2., az_steps=361):

        sky_model = {'theta':[], 'x':[], 'y': [], 's':[]}

        s=np.random.normal(size=size)
        sky_model['s']=s
        fig = pl.figure(figsize=(12,4))
        ax = fig.add_subplot(111)

        ax.imshow(s, interpolation='nearest')
        ax.set_aspect(1)
        theta=np.linspace(0,2*np.pi,az_steps) #why 6. and not 6.28 (2pi)?
        sky_model['theta'] = theta
        x=2.*np.cos(theta)+x_c #x_c and y_c are the coordinates of the center of the wvr circle.
        y=2.*np.sin(theta)+y_c
        sky_model['x'] = x
        sky_model['y'] = y
        pl.plot(x,y,'w')
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        #pl.xlabel()
        pl.title('Simulated Random Sky', fontsize=18)
        pl.savefig(self.pf+'/Sky_Sim.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()
        return sky_model


    def scan_sky_oneturn(self, sky_model):
        theta=sky_model['theta']
        x=sky_model['x']
        y=sky_model['y']
        s=sky_model['s']
        signal=np.zeros(len(theta))
        for n in range (len(theta)):
            k=np.int(np.round(x[n]))
            j=np.int(np.round(y[n]))
            signal[n]=s[j,k]
            #signal is what the wvr sees in one (almost) full turn.
        pl.plot(signal)
        pl.xlabel('Az[deg]', fontsize=18)
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        pl.title('Az Signal - One Turn', fontsize=18)
        pl.savefig(self.pf+'/Sky_Sim_WVRmapped_1Azturn.png', fontsize=18)
        if self.show==1:
            pl.show()
        else:
            pl.close()

    def scan_sky_map(self, sky_model, v=20.*(1./50.), nscans=110):
        theta=sky_model['theta']
        x=sky_model['x']
        y=sky_model['y']
        s=sky_model['s']

        atmogram=np.zeros(shape=(len(theta),nscans))
        for m in range(nscans):
            for n in range (len(theta)):
                k=np.int(np.round(x[n]+m*v))%20
                j=np.int(np.round(y[n]))
                atmogram[n,m]=s[j,k]


        fig = pl.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        ax.imshow(atmogram, aspect='auto', interpolation='nearest', origin='lower')
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        pl.title('Simulated Sky mapped by the WVR', fontsize=18)
        pl.savefig(self.pf+'/Sky_Sim_WVRmapped_v_'+str(round(v,1))+'.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        return atmogram



    def load_real_atmo(self, wvr_scan):

        waz, az_wvr, az_corr, fs, idx, pwv, D_pwv = raz.read_pwvatmo(wvr_scan)

        print(np.shape(D_pwv)[0])

        #az_ticks=np.arange(np.shape(D_pwv)[0])
        az_lab_raw=az_wvr[fs.s[10]:fs.e[10]]
        az_lab_corr_raw=az_corr[fs.s[10]:fs.e[10]]

        az_lab=[int(a) for a in az_lab_raw]
        az_lab_corr=[int(ac) for ac in az_lab_corr_raw]

        print('fs.s=', fs.s)
        print('fs.e=', fs.e)
        print('az_wvr=', az_wvr)

        pl.plot(az_lab, az_lab_corr)
        pl.show()

        fig = pl.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        im=ax.imshow(D_pwv,aspect='auto', interpolation='nearest', origin='lower')
        pl.title(wvr_scan[:-4])
        pl.xticks(fontsize=18)
        ax.set_yticks(az_lab[::50])
        ax.set_yticklabels(az_lab[::50], fontsize=18)
        cbar = fig.colorbar(im, extend='both')
        cbar.set_label('PWV[um]', fontsize=18)
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_atmo_data_systematics_clean.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        fig = pl.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        im=ax.imshow(D_pwv,aspect='auto', interpolation='nearest', origin='lower')
        pl.title(wvr_scan[:-4], fontsize=18)
        pl.xticks(fontsize=18)
        pl.yticks(ticks=az_lab[::20], labels=az_lab_corr[::20], fontsize=18)
        cbar = fig.colorbar(im, extend='both')
        cbar.set_label('PWV[um]', fontsize=18)
        cbar.ax.tick_params(labelsize=18)
        #pl.ylim(-120, -70)
        pl.xlabel('Nscan', fontsize=18)
        pl.ylabel('Az[deg]', fontsize=18)
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_atmo_data_systematics_clean_point.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        return D_pwv, waz, pwv



    def remove_atmo_dipole(self, D_pwv, pk_fn):
    #Removing Atmospheric Dipole
        mod=1
        def modulation(x, a, phi, C):
            return C + (a * np.sin(np.radians(mod*x) + phi))

        def err_modulation(p_j, x, y):
            return modulation(x, p_j[0], p_j[1], p_j[2]) - y

        def err_global_modulation(p, x, y_matrix):
            err_list=[]
            for i in range(len(y_matrix[1,:])): #=nscans
                p_i= p[1], p[i+2], p[0] #A,phi,C
                y_i=y_matrix[:,i]
                y_i[~np.isfinite(y_i)]=np.nanmean(y_i)
                err_i=err_modulation(p_i, x, y_i)
                err_list.append(err_i)

            err_array=np.array(err_list)
            return np.concatenate(err_array)

        def remove_dipole(theta_az, D_matrix):
            model=np.zeros(np.shape(D_matrix))
            mod_removed_data=np.zeros(np.shape(D_matrix))

            amp=5.
            offs=100.
            phi=0.
            p=[offs, amp]

            for i in range(len(D_matrix[1,:])):
                p.append(phi)

            bounds_inf=np.full(len(p), -np.pi*2)
            bounds_inf[0]=0.
            bounds_inf[1]=0.
            bounds_sup=np.full(len(p), np.pi*2)
            bounds_sup[0]=np.inf
            bounds_sup[1]=np.inf

            res = sp.least_squares(err_global_modulation, p, bounds=(bounds_inf, bounds_sup), args=(theta_az, D_matrix))

            p_best=res.x
            cov_x=res.jac
            p_err=np.sqrt(np.diag(cov_x))
            amp_out=p_best[0]

            for j in range (len(D_matrix[1,:])):
                p_out=p_best[1], p_best[2+j], p_best[0]
                model[:,j] = modulation(theta_az, *p_out)
                mod_removed_data[:,j] = D_matrix[:,j] - model[:,j]

            return p_best, p_err, mod_removed_data, model



        x_az= np.arange(0,361)
        p_best, p_err, mod_removed_data, model = remove_dipole(x_az, D_pwv)

        fig, axs = pl.subplots(3, 1, figsize=(20,20))
        im0=axs[0].imshow(D_pwv, aspect='auto',interpolation='nearest', origin='lower')
        axs[0].set_title('PWV Atmogram \n(Instrumental dipole and quadrupole removed)', fontsize=18)
        axs[0].set_xticks()
        cbar = fig.colorbar(im0, extend='both', ax=axs[0])
        cbar.set_label('PWV[um]', fontsize=18)
        im1=axs[1].imshow(model, aspect='auto',interpolation='nearest', origin='lower')
        axs[1].set_title('Atmospheric Dipole Model', fontsize=18)
        axs[1].set_xticks()
        cbar = fig.colorbar(im1, extend='both', ax=axs[1])
        cbar.set_label('PWV[um]', fontsize=18)
        im2=axs[2].imshow(mod_removed_data, aspect='auto',interpolation='nearest', origin='lower')
        axs[2].set_title('PWV Atmogram - Atmospheric Dipole Model', fontsize=18)
        cbar = fig.colorbar(im2, extend='both', ax=axs[2])
        cbar.set_label('PWV[um]', fontsize=18)
        subplots_adjust(hspace=0.01)
        pl.colorbar()
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_atmo_data_systematics+atmo_dipole_clean.png')

        if self.show==1:
            pl.show()
        else:
            pl.close()

        dipole_removed_atmo = {'D_pwv':[], 'model':[], 'mod_removed_data':[]}
        dipole_removed_atmo['D_pwv']=D_pwv
        dipole_removed_atmo['model']=model
        dipole_removed_atmo['mod_removed_data']=mod_removed_data

        f = open(pk_fn,'wb')
        pk.dump(dipole_removed_atmo,f)
        f.close()

        return dipole_removed_atmo

    def remove_offs(self, double_mod_removed_data, wvr_scan):

        nscans=len(double_mod_removed_data[10,:])
        double_mod_removed_data_nooffs=np.full(np.shape(double_mod_removed_data), np.nan)

        for i in range(nscans):
            y=double_mod_removed_data[:,i]
            double_mod_removed_data_nooffs[:,i]=y-np.mean(y)

        fig=pl.figure(figsize=(10,8))
        im=pl.imshow(double_mod_removed_data_nooffs, aspect='auto', vmin=-1.5, vmax=1.5, interpolation='nearest', origin='lower')
        # pl.title('PWV Atmogram \nAtmospheric Dipole Model and Offset subtracted', fontsize=18)
        #pl.title('PWV Atmogram', fontsize=18)
        pl.title('PWV Fluctuations Map', fontsize=30)
        pl.xlabel('Nscan[t/30s]', fontsize=25)
        pl.ylabel('Az[deg]', fontsize=25)
        pl.xticks(fontsize=25)
        pl.yticks(fontsize=25)
        cbar = fig.colorbar(im, extend='both')
        cbar.set_label('PWV[um]', fontsize=25)
        cbar.ax.tick_params(labelsize=25)
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_atmo_data_systematics+atmo_dipole+offset_clean.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        return double_mod_removed_data_nooffs

    def find_ws_wd(self, atmogram, wvr_scan, delta_az=20):

        def gauss(x, *p):
            A, mu, sigma, offs = p
            return A*np.exp(-(x-mu)**2/(2.*sigma**2))+(offs)

        nscans=len(atmogram[10,:])

        lag=np.arange(2*nscans-1)-nscans
        npts=len(atmogram[:,1])-delta_az
        #npts=5 #just for test

        delay=np.zeros(npts)
        phase=np.arange(npts)

        c_set=['r', 'g', 'y', 'c', 'b']

        pl.figure(figsize=(15,8))
        for k in range(1,npts-1):
            rowcor1=np.correlate(atmogram[k,:],atmogram[k+delta_az,:], mode='full')/(atmogram[k,:].std()*atmogram[k+delta_az,:].std())
            #delay[k]=rowcor1.argmax()
            pl.plot(lag, rowcor1, label='Az'+str(k)+'-'+str(k+20))
            try:
                fit_mask=np.where((lag>=-6) & (lag<=6))
                p0=[100., 0., 5., 0.]
                coeff, var_matrix = curve_fit(gauss, lag[fit_mask], rowcor1[fit_mask], p0=p0, bounds=((0, -100, 0, 0), (1000, 100, 10, 100)))
                x_model=np.linspace(np.min(lag[fit_mask]), np.max(lag[fit_mask]), 100)
                g=gauss(x_model, *coeff)
                delay[k]=coeff[1]
                #delay[k]=nscans-rowcor1.argmax()
            except:
                delay[k]=nscans-rowcor1.argmax()
            #pl.plot(x_model, g, c=c_set[k], label='gaussian model:\nmu[deg] ='+str(round(coeff[1],1)))
            #pl.xlim(-200,200)

        #pl.legend()
        pl.xlabel('lag', fontsize=18)
        pl.ylabel('rowcor/nscans', fontsize=18)
        pl.suptitle('Correlation Function', fontsize=18)
        pl.title('delta_az = '+str(delta_az), fontsize=18)
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        pl.suptitle(wvr_scan[:-4], fontsize=18)
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_correlation_plots_deltaAz_'+str(delta_az)+'.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        #mask_good=np.where(((nscans-delay)<7) & ((nscans-delay)>-7))
        #mask_good=np.where((delay<10) & (delay>-10))

        def sin_func(x, A, phi, C):
            return A*np.sin((x*2.*np.pi/np.max(phase))+phi)+C

        p0=[1., 0., 0.]
        #coeff, var_matrix = curve_fit(sin_func, phase[mask_good], (nscans-delay)[mask_good], p0=p0, bounds=((0., -2.*np.pi, 0.), (np.inf, 2.*np.pi, np.inf)))
        coeff, var_matrix = curve_fit(sin_func, phase, delay, p0=p0, bounds=((0., -2.*np.pi, -np.inf), (np.inf, 2.*np.pi, np.inf)))
        model=sin_func(phase, *coeff)


        noise = delay - model
        sigma = np.std(noise)
        pl.figure()
        n, x, _ = pl.hist(noise, bins=np.linspace(-10, 10, 40), histtype=u'step', density=True)
        pl.suptitle('Noise 1')
        pl.title('sigma = '+ str(round(sigma,2)))
        #pl.show()
        pl.close()

        pl.figure(figsize=(15,8))
        pl.suptitle('Noise Cut', fontsize=18)
        pl.plot(phase, noise, label='noise')
        sigma2_cut = np.where((noise <= 2*sigma) & (noise >= (-2*sigma)))
        pl.axhline(y=1*sigma, ls='--', alpha=0.5, label='2sigma')
        noise2= noise[sigma2_cut]
        sigma2=np.std(noise2)
        pl.close()
        pl.figure()
        n, x, _ = pl.hist(noise2, bins=np.linspace(-10, 10, 40), histtype=u'step', density=True)
        pl.suptitle('Noise 2')
        pl.title('sigma = '+ str(round(sigma2,2)))
        pl.close()

        pl.plot(phase[sigma2_cut], noise[sigma2_cut], label='noise_cut1')
        pl.axhline(y=1*sigma2, ls='--', alpha=0.5, label='2sigma2')
        stage2_sigma2_cut = np.where((noise <= 2*sigma2) & (noise >= (-2*sigma2)))
        pl.plot(phase[stage2_sigma2_cut], noise[stage2_sigma2_cut], label='noise_cut2')
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        # noise3=noise[stage2_sigma2_cut]
        # sigma3=np.std(noise3)
        # pl.axhline(y=2*sigma3, ls='--', alpha=0.5, label='2sigma2')
        # new_sigma2_cut = np.where((noise <= 2*sigma2) & (noise >= (-2*sigma2)))
        # pl.plot(phase[new_sigma2_cut], noise[new_sigma2_cut], label='noise_cut2')

        pl.legend(fontsize=18)
        #pl.show()
        pl.close()

        data_fit_mask = stage2_sigma2_cut
        p0=[1., 0., 0.]
        coeff, var_matrix = curve_fit(sin_func, phase[data_fit_mask], delay[data_fit_mask], p0=p0, bounds=((0., -2.*np.pi, -np.inf), (np.inf, 2.*np.pi, np.inf)))
        print('offset=', coeff[2])
        print('avg=', np.mean(delay[data_fit_mask]))
        model=sin_func(phase, *coeff)
        tau = (coeff[0]*30) #1scan=30s
        h_clouds=1900#m
        delta_x_sky = h_clouds*np.sin(np.radians(delta_az))*(1./np.tan(math.radians(self.el)))
        ws = delta_x_sky/tau
        print('wind speed =', ws)
        wd = np.degrees(coeff[1])+(delta_az/2.) #midpoint between 2 correlated timestreams
        if wd<0:
            wd=360.+wd
        print('wind dir =', wd)


        pl.figure(figsize=(15,8))
        #pl.plot(phase[mask_good],(nscans-delay)[mask_good],'-bs', markersize=5, linewidth=0.7, c='blue')
        pl.plot(phase[data_fit_mask], delay[data_fit_mask],'-bs', markersize=5, linewidth=0.7, c='blue')
        pl.ylim(-10,10)
        pl.xlabel('Azimuthal Phase [deg]', fontsize=18)
        pl.ylabel('Delay [s/30]', fontsize=18)
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        pl.suptitle(wvr_scan[:-4], fontsize=18)
        pl.title('ws[m/s] ='+str(round(ws,2))+'\nwd[deg] = '+str(round(wd,2)), fontsize=18)

        label_real='Amp[s/30]='+str(round(coeff[0],2))+'\n phase[deg]='+str(round(np.degrees(coeff[1]),2))
        pl.plot(phase, model, c='r', label=label_real)
        pl.legend(fontsize=18)
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_delay_deltaAz_'+str(delta_az)+'.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()

        return ws, wd






    def atmo_to_sky(self, waz, pwv, atmogram, ws, wd, h_clouds=1900.):#sky_model, v=20.*(1./50.), nscans=110):

        scan_v=360./30.
        t=waz/scan_v

        w_y=ws*np.sin(np.radians(wd))
        w_x=ws*np.cos(np.radians(wd))

        x=np.zeros(len(waz))
        y=np.zeros(len(waz))

        y_size=round(h_clouds*(1./np.tan(math.radians(self.el)))+ w_y*np.max(t))
        x_size=round(h_clouds*(1./np.tan(math.radians(self.el)))+ w_x*np.max(t))

        sky=np.full((int(x_size)+1,int(y_size)+1), np.nan)

        for i in range (len(waz)):
            waz_i=waz[i]
            t_i=waz_i/scan_v
            (scannum_i, az_i) = divmod(waz_i,360)
            y[i]=h_clouds*(1./np.tan(math.radians(self.el)))*np.sin(math.radians(az_i))+ w_y*t_i
            x[i]=h_clouds*(1./np.tan(math.radians(self.el)))*np.cos(math.radians(az_i))+ w_x*t_i
            k=np.int(np.round(x[i]))
            j=np.int(np.round(y[i]))
            sky[k,j]=pwv[i]
            #sky[int(x[i]), int(y[i])]=atmogram[int(az_i),int(scannum_i)]
        #
        #
        pl.scatter(x, y, s=4, c=pwv)#, cmap='plasma')
        pl.colorbar()
        pl.suptitle('PWV atmogram projected')
        #pl.title(wvr_scan[:-4])
        pl.show()



        fig = pl.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        ax.imshow(sky, aspect='auto', interpolation='gaussian', origin='lower', cmap='viridis')

        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        pl.title('Simulated Sky mapped by the WVR', fontsize=18)
        # pl.savefig(self.pf+'/Sky_Sim_WVRmapped_v_'+str(round(v,1))+'.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()



        return sky


    def NOAA_wind_data(self, wvr_scan):
        #Measure of how much the wind varies in one hour
        import datetime

        pkfn='SP_windspeed_042020.pk'
        f=open(pkfn, "rb")
        wind_dict=pk.load(f)
        f.close()


        t=wind_dict['t']
        ws=wind_dict['ws']

        yr=int(wvr_scan[:4])
        mm= int(wvr_scan[4:6])
        dd=int(wvr_scan[6:8])

        H_start=int(wvr_scan[9:11])
        H_end=H_start+1


        pl.figure(figsize=(12,8))
        pl.scatter(t,ws, s=6, c='r')
        #pl.axhline(y=np.mean(ws), ls='--', c='k', alpha=0.5, label='ws_avg[m/s]='+str(round(np.mean(ws),2)))
        pl.plot(t,ws)
        pl.xlim([datetime.datetime(yr,mm,dd,H_start,0,0),datetime.datetime(yr,mm,dd,H_end,0,0)])
        pl.ylim(8,10)
        pl.xlabel('2020-04-dd H:M', fontsize=18)
        pl.ylabel('wind speed [m/s]', fontsize=18)
        pl.title('wind speed from NOAA', fontsize=18)
        pl.xticks(fontsize=18)
        pl.yticks(fontsize=18)
        #pl.legend()
        pl.savefig(self.pf+'/'+wvr_scan[:-4]+'_NOAA_wind.png')
        if self.show==1:
            pl.show()
        else:
            pl.close()
