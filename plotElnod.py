import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import numpy as np
import pickle
import pylab as pl
import scipy.optimize as sp
import scipy.interpolate as interp
from matplotlib.gridspec import GridSpec


class ReadElnod(object):
    '''

    '''



    def __init__(self, unit=None, verb=True):

        '''


        '''

        self.T_hot=363.15
        self.T_cold=283.15
        self.LO=91.65
        self.IF=2*self.LO

        self.ch=[1.05,1.9,3.2,6.1]
        self.freqs=[self.IF-self.ch[3],self.IF-self.ch[2],self.IF-self.ch[1],self.IF-self.ch[0],self.IF+self.ch[0],self.IF+self.ch[1],self.IF+self.ch[2],self.IF+self.ch[3]]


    def read_elnod_slow(self, filename, data_folder='', plot_data=1):
        with open(data_folder+filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
            FH = np.loadtxt(lines, delimiter=' ', skiprows=1, usecols = (1,19,20,21,22,23,24))
        if plot_data==1:
            pl.plot(FH[:,0], FH[:,1], label="T_src-Ch0")
            pl.plot(FH[:,0], FH[:,2], label="T_src-Ch1")
            pl.plot(FH[:,0], FH[:,3], label="T_src-Ch2")
            pl.plot(FH[:,0], FH[:,4], label="T_src-Ch3")
            pl.legend()
            pl.title(filename)
            pl.show()
        return FH

    def interp_fastdata(self, time, data_w_zeros):
        good_ind=np.where(data_w_zeros > 0.)
        data_wout_zeros=data_w_zeros[good_ind]
        x_axis=time[good_ind]
        f=interp.interp1d(x_axis,data_wout_zeros, kind='cubic', fill_value='extrapolate')
        data_interp=f(time)
        return data_interp

    def calibrate_fast_data(self, interp_data):
        def calibrate_one_ch(data_C, data_A, data_H, data_B):
            calib_ch_data=np.zeros(len(data_A))
            T_hot_arr=np.full(len(data_A), self.T_hot)
            T_cold_arr=np.full(len(data_A), self.T_cold)

            T_ref=(T_hot_arr+T_cold_arr)/2.
            V_ref=(data_H+data_C)/2.
            G=(data_H-data_C)/(T_hot_arr-T_cold_arr)
            V_sky=(data_A+data_B)/2.
            calib_ch_data=T_ref+((V_sky-V_ref)/G)
            return calib_ch_data, G

        data_calibrated=np.zeros((len(interp_data[:,0]),7))
        G_array=np.zeros((len(interp_data[:,0]),4))
        data_calibrated[:,0]=interp_data[:,0]
        data_calibrated[:,5]=interp_data[:,17]
        data_calibrated[:,6]=interp_data[:,18]
        for i in range (4):
            data_calibrated[:,1+i],G_array[:,i]=calibrate_one_ch(interp_data[:,4*i+1], interp_data[:,4*i+2], interp_data[:,4*i+3], interp_data[:,4*i+4])
        return data_calibrated, G_array




    def read_elnod_fast(self, filename, test_plot_int=1, test_plot_cal=0, data_folder=''):
    #this function read fast files and returns as an array the
    #interpolated values for [t, T0_A, T1_A, T2_A, T3_A, El, Az]
    #plus it writes them into a file *_interp.txt
        with open(data_folder+filename) as f:
            lines = (line for line in f if not (line.startswith('#') or line.startswith('TIME')))
            FH = np.loadtxt(lines, delimiter=' ', skiprows=1, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
            #The data order is: TIMEWVR CH0C CH0A CH0H CH0B CH1C CH1A CH1H CH1B CH2C CH2A CH2H CH2B CH3C CH3A CH3H CH3B EL AZ

        interp_data=np.zeros((len(FH[:,0]),19))
        interp_data[:,0]=FH[:,0]
        interp_data[:,17]=FH[:,17]
        interp_data[:,18]=FH[:,18]
        for i in range(len(interp_data[0,:])-3):
            interp_data[:,i+1]=self.interp_fastdata(FH[:,0], FH[:,1+i])

        if test_plot_int==1:
            #fig, axes = plt.subplots(2, 2)
            fig=plt.figure()
            gs = GridSpec(4, 2, width_ratios=[1, 1], height_ratios=[3, 1, 3, 1], hspace=0.5)
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
            ax5 = fig.add_subplot(gs[4])
            ax6 = fig.add_subplot(gs[5])
            ax7 = fig.add_subplot(gs[6])
            ax8 = fig.add_subplot(gs[7])

            T0=np.concatenate((interp_data[:,2], interp_data[:,4]), axis=None)
            el_fit=np.concatenate((interp_data[:,17], interp_data[:,17]), axis=None)
            T0_coef=np.polyfit(el_fit, T0, 7)
            T0_fit=np.poly1d(T0_coef)
            T0_plot=T0_fit(interp_data[:,17])

            ax1.plot(interp_data[:,17], interp_data[:,2], c='r', label="V_A")
            ax1.plot(interp_data[:,17], interp_data[:,4], c='c', label="V_B")
            ax1.plot(interp_data[:,17], T0_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax1.set_xlabel('El[deg]')
            ax1.set_ylabel('V_sky [ADC units]')
            ax1.legend(loc='upper right')
            ax1.set_title('Channel 0')

            ax3.plot(interp_data[:,17], interp_data[:,2]-T0_plot, c='r', label="V_A - PolyFit")
            ax3.plot(interp_data[:,17], interp_data[:,4]-T0_plot, c='c', label="V_B - PolyFit")
            ax3.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax3.legend(loc='upper right')

            T1=np.concatenate((interp_data[:,6], interp_data[:,8]), axis=None)
            T1_coef=np.polyfit(el_fit, T1, 7)
            T1_fit=np.poly1d(T1_coef)
            T1_plot=T1_fit(interp_data[:,17])

            ax2.plot(interp_data[:,17], interp_data[:,6], c='r', label="V_A")
            ax2.plot(interp_data[:,17], interp_data[:,8], c='c', label="V_B")
            ax2.plot(interp_data[:,17], T1_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax2.set_xlabel('El[deg]')
            ax2.set_ylabel('V_sky [ADC units]')
            ax2.legend(loc='upper right')
            ax2.set_title('Channel 1')

            ax4.plot(interp_data[:,17], interp_data[:,6]-T1_plot, c='r', label="V_A - PolyFit")
            ax4.plot(interp_data[:,17], interp_data[:,8]-T1_plot, c='c', label="V_B - PolyFit")
            ax4.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax4.legend(loc='upper right')

            T2=np.concatenate((interp_data[:,10], interp_data[:,12]), axis=None)
            T2_coef=np.polyfit(el_fit, T2, 7)
            T2_fit=np.poly1d(T2_coef)
            T2_plot=T2_fit(interp_data[:,17])

            ax5.plot(interp_data[:,17], interp_data[:,10], c='r', label="V_A")
            ax5.plot(interp_data[:,17], interp_data[:,12], c='c', label="V_B")
            ax5.plot(interp_data[:,17], T2_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax5.set_xlabel('El[deg]')
            ax5.set_ylabel('V_sky [ADC units]')
            ax5.legend(loc='upper right')
            ax5.set_title('Channel 2')

            ax7.plot(interp_data[:,17], interp_data[:,10]-T2_plot, c='r', label="V_A - PolyFit")
            ax7.plot(interp_data[:,17], interp_data[:,12]-T2_plot, c='c', label="V_B - PolyFit")
            ax7.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax7.legend(loc='upper right')

            T3=np.concatenate((interp_data[:,14], interp_data[:,16]), axis=None)
            T3_coef=np.polyfit(el_fit, T3, 7)
            T3_fit=np.poly1d(T3_coef)
            T3_plot=T3_fit(interp_data[:,17])

            ax6.plot(interp_data[:,17], interp_data[:,14], c='r', label="V_A")
            ax6.plot(interp_data[:,17], interp_data[:,16], c='c', label="V_B")
            ax6.plot(interp_data[:,17], T3_plot, linestyle='--', linewidth=0.5, c='k', label="Poly Fit")
            ax6.set_xlabel('El[deg]')
            ax6.set_ylabel('V_sky [ADC units]')
            ax6.legend(loc='upper right')
            ax6.set_title('Channel 3')

            ax8.plot(interp_data[:,17], interp_data[:,14]-T3_plot, c='r', label="V_A - PolyFit")
            ax8.plot(interp_data[:,17], interp_data[:,16]-T3_plot, c='c', label="V_B - PolyFit")
            ax8.axhline(y=0, linestyle='--', linewidth=0.5, c='k')
            ax8.legend(loc='upper right')

            pl.suptitle('Raw Temperature readings at mirror position A and B\n'+filename[-31:-4])
            #pl.savefig('Output/'+filename[:-4]+'_Backlash-corrected_data.png')
            #pl.show()
            pl.close()



        calib_data,G=self.calibrate_fast_data(interp_data)

        if test_plot_cal==1:
            fig = pl.figure()
            ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
            ax2 = pl.subplot2grid((3,1), (2,0))
            ax1.plot(calib_data[:,0], calib_data[:,1], c='k')
            ax1.set_title('Calibrated Temperature')
            ax2.plot(calib_data[:,0], G[:,0], c='r')
            ax2.set_title('Gain')
            # pl.plot(calib_data[:,0], calib_data[:,2], label="T_src-Ch1")
            # pl.plot(calib_data[:,0], calib_data[:,3], label="T_src-Ch2")
            # pl.plot(calib_data[:,0], calib_data[:,4], label="T_src-Ch3")
            pl.suptitle('Calibrated T0 data\n'+filename[:-4])
            pl.show()

        mat = np.matrix(interp_data)
        with open(data_folder+filename[:-4]+'_interp.txt','wb') as f:
            for line in mat:
                np.savetxt(f, line, fmt='%.5f')
        f.close()

        mat_c = np.matrix(calib_data)
        with open(data_folder+filename[:-4]+'_calib.txt','wb') as f:
            for line in mat_c:
                np.savetxt(f, line, fmt='%.5f')
        f.close()

        return calib_data

    def read_elnod(self, filename):
        if filename[-8:-4]=='slow':
            FH=self.read_elnod_slow(filename)
        elif filename[-8:-4]=='fast':
            FH=self.read_elnod_fast(filename)
        else:
            print("ERROR in the given filename.")
            FH=0
        return FH


    def filter_data(self, ymodel_unfilt, alpha):
        #returns data filtered including the averaging over a time constant
        X=np.array([ymodel_unfilt[0]])
        for j in range(len(ymodel_unfilt)):
             X=np.append(X, np.abs(alpha)*ymodel_unfilt[j] + (1.-np.abs(alpha))*X[j])
        return X[1:]


    def fit_elnod_vsfreq(self, filename, el_range=(0,90.), plot=1):
        low_el, high_el = el_range
        def line_gauss(x, c, a, sigma):
            return c + a * np.exp(-(x - 183.3)**2 / (2 * sigma**2))  #to understand why it doesn't work if i define it out of here
        function=line_gauss
        FH=self.read_elnod(filename)
        fit_par=np.zeros((3,len(FH[:,1])))
        T_data=[]
        for i in range (0, len(FH[:,1])):
            T=[FH[i][4],FH[i][3],FH[i][2],FH[i][1],FH[i][1],FH[i][2],FH[i][3],FH[i][4]]
            T_data.append([self.freqs, T])
            popt,pcov = sp.curve_fit(function, self.freqs, T, p0 = [1, 10, 10])
            if FH[i][5]>=90.:
                print("El=", FH[i][5])
                print("Ts=", T)
            if (plot==1) and (low_el<= FH[i][5] <=high_el):
                freqs_fit=np.linspace(170,190,2000)
                pl.plot(self.freqs, T, '.', label="El="+str(FH[i][5]))
                pl.xlabel('f[GHz]')
                pl.ylabel('T[K]')
                pl.plot(freqs_fit, line_gauss(freqs_fit, *popt), linewidth=1, label="El="+str(FH[i][5])+'_fit')
                pl.legend()
                pl.xlim(175,190)
            fit_par[:,i]=popt
        if plot==1:
            pl.show()

        return fit_par, T_data


    def fit_elnod_vsel(self, filename, layers=1, plot=1):
        def line_shape_1lay(x, T, beta, alpha):
            tau=np.exp(-beta/np.sin((x)*np.pi/180.))
            T_meas=T*(1.-tau)
            return self.filter_data(T_meas, alpha) #best bet on time constant alpha=0.6

        def line_shape_2lay(x, Th, Tl, betah, betal):#, alpha):
            tau_l=np.exp(-betal/np.sin((x)*np.pi/180.))
            tau_h=np.exp(-betah/np.sin((x)*np.pi/180.))
            T_meas=Th*(1.-tau_h)+Tl*tau_h*(1.-tau_l) #putting d=1
            return self.filter_data(T_meas, alpha=0.6)

        FH=self.read_elnod(filename)

        t=FH[:,0]
        T1=FH[:,1]
        T2=FH[:,2]
        T3=FH[:,3]
        T4=FH[:,4]
        El=FH[:,5]

        p_0=[250., 0.5, 0]
        p_1=[250., 0.5, 0]
        p_2=[250., 0.5, 0]
        p_3=[250., 0.5, 0]

        p_best0, cov_x0 = sp.curve_fit(line_shape_1lay, El, T1, p0=p_0)
        p_best1, cov_x1 = sp.curve_fit(line_shape_1lay, El, T2, p0=p_1)
        p_best2, cov_x2 = sp.curve_fit(line_shape_1lay, El, T3, p0=p_2)
        p_best3, cov_x3 = sp.curve_fit(line_shape_1lay, El, T4, p0=p_3)

        fit_par0=np.zeros((len(p_best0),2))
        fit_par1=np.zeros((len(p_best1),2))
        fit_par2=np.zeros((len(p_best2),2))
        fit_par3=np.zeros((len(p_best3),2))

        fit_par0[:,0]=p_best0
        fit_par0[:,1]=np.sqrt(np.diag(cov_x0))
        fit_par1[:,0]=p_best1
        fit_par1[:,1]=np.sqrt(np.diag(cov_x1))
        fit_par2[:,0]=p_best2
        fit_par2[:,1]=np.sqrt(np.diag(cov_x2))
        fit_par3[:,0]=p_best3
        fit_par3[:,1]=np.sqrt(np.diag(cov_x3))


        fig, axes = plt.subplots(2, 2)
        axes[0, 0].plot(T1, '.', ms=3, c='k', label="Ch_0")
        axes[0, 0].plot(line_shape_1lay(El, fit_par0[0,0], fit_par0[1,0], fit_par0[2,0]), linewidth=1, c='r', label="Fit")
        axes[0, 0].text(0.,  T1.max()-10., 'T ={:.2f}'.format(round(fit_par0[0,0],2))+'±{:.2f}'.format(round(fit_par0[0,1],2)), fontsize=10)
        axes[0, 0].text(0.,  T1.max()-20., 'beta ={:.4f}'.format(round(fit_par0[1,0],3))+'±{:.4f}'.format(round(fit_par0[1,1],4)), fontsize=10)
        axes[0, 0].text(0.,  T1.max()-30., 'alpha ={:.3f}'.format(round(fit_par0[2,0],3))+'±{:.3f}'.format(round(fit_par0[2,1],3)), fontsize=10)
        axes[0, 0].legend(loc='upper right', fontsize='large')
        axes[0, 0].set_ylabel('T[K]')
        axes[0, 0].set_xlabel('sample #')

        axes[0, 1].plot(T2, '.', ms=3, c='k',label="Ch_1")
        axes[0, 1].plot(line_shape_1lay(El, fit_par1[0,0], fit_par1[1,0], fit_par1[2,0]), linewidth=1, c='r', label="Fit")
        axes[0, 1].text(0.,  T2.max()-10., 'T ={:.2f}'.format(round(fit_par1[0,0],2))+'±{:.2f}'.format(round(fit_par1[0,1],2)), fontsize=10)
        axes[0, 1].text(0.,  T2.max()-20., 'beta ={:.4f}'.format(round(fit_par1[1,0],3))+'±{:.4f}'.format(round(fit_par1[1,1],3)), fontsize=10)
        axes[0, 1].text(0.,  T2.max()-30., 'alpha ={:.3f}'.format(round(fit_par1[2,0],3))+'±{:.3f}'.format(round(fit_par1[2,1],3)), fontsize=10)
        axes[0, 1].legend(loc='upper right', fontsize='large')
        axes[0, 1].set_ylabel('T[K]')
        axes[0, 1].set_xlabel('sample #')

        axes[1, 0].plot(T3, '.', ms=3, c='k',label="Ch_2")
        axes[1, 0].plot(line_shape_1lay(El, fit_par2[0,0], fit_par2[1,0], fit_par2[2,0]), linewidth=1, c='r', label="Fit")
        axes[1, 0].text(0.,  T3.max()-10., 'T ={:.2f}'.format(round(fit_par2[0,0],2))+'±{:.2f}'.format(round(fit_par2[0,1],2)), fontsize=10)
        axes[1, 0].text(0.,  T3.max()-20., 'beta ={:.4f}'.format(round(fit_par2[1,0],3))+'±{:.4f}'.format(round(fit_par2[1,1],3)), fontsize=10)
        axes[1, 0].text(0.,  T3.max()-30., 'alpha ={:.3f}'.format(round(fit_par2[2,0],3))+'±{:.3f}'.format(round(fit_par2[2,1],3)), fontsize=10)
        axes[1, 0].legend(loc='upper right', fontsize='large')
        axes[1, 0].set_ylabel('T[K]')
        axes[1, 0].set_xlabel('sample #')

        axes[1, 1].plot(T4, '.', ms=3, c='k',label="Ch_3")
        axes[1, 1].plot(line_shape_1lay(El, fit_par3[0,0], fit_par3[1,0], fit_par3[2,0]), linewidth=1, c='r', label="Fit")
        axes[1, 1].text(0.,  T4.max()-10., 'T ={:.2f}'.format(round(fit_par3[0,0],2))+'±{:.2f}'.format(round(fit_par3[0,1],2)), fontsize=10)
        axes[1, 1].text(0.,  T4.max()-20., 'beta ={:.4f}'.format(round(fit_par3[1,0],3))+'±{:.4f}'.format(round(fit_par3[1,1],3)), fontsize=10)
        axes[1, 1].text(0.,  T4.max()-30., 'alpha ={:.3f}'.format(round(fit_par3[2,0],3))+'±{:.3f}'.format(round(fit_par3[2,1],3)), fontsize=10)
        axes[1, 1].legend(loc='upper right', fontsize='large')
        axes[1, 1].set_ylabel('T[K]')
        axes[1, 1].set_xlabel('sample #')

        pl.suptitle(filename[:-4])#+'\nTatm free to vary with the channel')
        pl.show()

        if layers==1:
            model1=line_shape_1lay
            fit_par=np.zeros((6,2))

        elif layers==2:
             model1=line_shape_2lay
             #fit_par=np.zeros((11,2)) #if we also want to fit alpha
             fit_par=np.zeros((10,2))

        def quad_fit(p_array, xdata, ydata1, ydata2, ydata3, ydata4, fit_model, layers=layers):
            def err(p, x, y):
                if layers==1:
                    return model1(x, p[0],p[1],p[2]) - y
                if layers==2:
                    return model1(x, p[0],p[1],p[2],p[3]) - y

            def err_global(p, x, y1, y2, y3, y4):
                # p is now a_1, b, c_1, a_2, c_2, with b shared between the two
                if layers==1:
                    p1 = p[0], p[1], p[5]
                    p2 = p[0], p[2], p[5]
                    p3 = p[0], p[3], p[5]
                    p4 = p[0], p[4], p[5]
                if layers==2:
                    p1 = p[0], p[1], p[2], p[3]#, p[10]
                    p2 = p[0], p[1], p[4], p[5]#, p[10]
                    p3 = p[0], p[1], p[6], p[7]#, p[10]
                    p4 = p[0], p[1], p[8], p[9]#, p[10]
                    # p1 = p[0], p[0], p[2], p[2]#, p[10]
                    # p2 = p[0], p[0], p[4], p[4]#, p[10]
                    # p3 = p[0], p[0], p[6], p[6]#, p[10]
                    # p4 = p[0], p[0], p[8], p[8]#, p[10]


                err1 = err(p1, x, y1)
                err2 = err(p2, x, y2)
                err3 = err(p3, x, y3)
                err4 = err(p4, x, y4)

                return np.concatenate((err1, err2, err3, err4))

            p_best, cov_x, infodict, mesg, ier = sp.leastsq(err_global, p_array, args=(xdata, ydata1, ydata2, ydata3, ydata4), full_output=1)
            print("p_best=",p_best)
            print("cov_x=",cov_x)
            #cov_x=np.zeros((6,6))
            fit_par[:,0]=p_best
            fit_par[:,1]=np.sqrt(np.diag(cov_x))
            print("p_best=", p_best)
            print("p_err=", np.sqrt(np.diag(cov_x)))
            if ier==1:
                print("The Fit converged.")
            else:
                print("The Fit did not converge.")

            return fit_par

        if layers==1:
            p_array=[250., 0.5, 0.5, 0.5, 0.5, 0.]
        if layers==2:
            p_array=[250., 20., 1.2, 33., 0.5, 20., 0.2, 1.1, 0.1, 1.]

        fit_par = quad_fit(p_array, El, T1, T2, T3, T4, model1)

        if (plot==1 and layers==1):
            fig, axes = plt.subplots(2, 2)
            axes[0, 0].plot(T1, '.', ms=3, c='k', label="Ch_0")
            axes[0, 0].plot(model1(El, fit_par[0,0], fit_par[1,0], fit_par[5,0]), linewidth=1, c='r', label="Fit")
            axes[0, 0].text(0.,  T1.max()-10., 'T ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=10)
            axes[0, 0].text(0.,  T1.max()-20., 'beta ={:.4f}'.format(round(fit_par[1,0],4))+'±{:.4f}'.format(round(fit_par[1,1],4)), fontsize=10)
            axes[0, 0].text(0.,  T1.max()-30., 'alpha ={:.3f}'.format(round(fit_par[5,0],3))+'±{:.3f}'.format(round(fit_par[5,1],3)), fontsize=10)
            axes[0, 0].legend(loc='upper right', fontsize='large')
            axes[0, 0].set_ylabel('T[K]')
            axes[0, 0].set_xlabel('sample #')

            axes[0, 1].plot(T2, '.', ms=3, c='k',label="Ch_1")
            axes[0, 1].plot(model1(El, fit_par[0,0], fit_par[2,0], fit_par[5,0]), linewidth=1, c='r', label="Fit")
            axes[0, 1].text(0.,  T2.max()-10., 'T ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=10)
            axes[0, 1].text(0.,  T2.max()-20., 'beta ={:.4f}'.format(round(fit_par[2,0],4))+'±{:.4f}'.format(round(fit_par[2,1],4)), fontsize=10)
            axes[0, 1].text(0.,  T2.max()-30., 'alpha ={:.3f}'.format(round(fit_par[5,0],3))+'±{:.3f}'.format(round(fit_par[5,1],3)), fontsize=10)
            axes[0, 1].legend(loc='upper right', fontsize='large')
            axes[0, 1].set_ylabel('T[K]')
            axes[0, 1].set_xlabel('sample #')

            axes[1, 0].plot(T3, '.', ms=3, c='k',label="Ch_2")
            axes[1, 0].plot(model1(El, fit_par[0,0], fit_par[3,0], fit_par[5,0]), linewidth=1, c='r', label="Fit")
            axes[1, 0].text(0.,  T3.max()-10., 'T ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=10)
            axes[1, 0].text(0.,  T3.max()-20., 'beta ={:.4f}'.format(round(fit_par[3,0],4))+'±{:.4f}'.format(round(fit_par[3,1],4)), fontsize=10)
            axes[1, 0].text(0.,  T3.max()-30., 'alpha ={:.3f}'.format(round(fit_par[5,0],3))+'±{:.3f}'.format(round(fit_par[5,1],3)), fontsize=10)
            axes[1, 0].legend(loc='upper right', fontsize='large')
            axes[1, 0].set_ylabel('T[K]')
            axes[1, 0].set_xlabel('sample #')

            axes[1, 1].plot(T4, '.', ms=3, c='k',label="Ch_3")
            axes[1, 1].plot(model1(El, fit_par[0,0], fit_par[4,0], fit_par[5,0]), linewidth=1, c='r', label="Fit")
            axes[1, 1].text(0.,  T4.max()-10., 'T ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=10)
            axes[1, 1].text(0.,  T4.max()-20., 'beta ={:.4f}'.format(round(fit_par[4,0],4))+'±{:.4f}'.format(round(fit_par[4,1],4)), fontsize=10)
            axes[1, 1].text(0.,  T4.max()-30., 'alpha ={:.3f}'.format(round(fit_par[5,0],3))+'±{:.3f}'.format(round(fit_par[5,1],3)), fontsize=10)
            axes[1, 1].legend(loc='upper right', fontsize='large')
            axes[1, 1].set_ylabel('T[K]')
            axes[1, 1].set_xlabel('sample #')

            pl.suptitle(filename[:-4])#+'\n1 Tatm for all the channels')

            pl.show()

        if (plot==1 and layers==2):
            fig, axes = plt.subplots(2, 2)
            axes[0, 0].plot(T1, '.', label="Ch_0")
            axes[0, 0].plot(model1(El, fit_par[0,0], fit_par[1,0], fit_par[2,0], fit_par[3,0]), linewidth=1, label="Fit")
            axes[0, 0].text(0.,  T1.max()-10, 'T_h ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=6)
            axes[0, 0].text(0.,  T1.max()-20, 'T_l ={:.2f}'.format(round(fit_par[1,0],2))+'±{:.2f}'.format(round(fit_par[1,1],2)), fontsize=6)
            axes[0, 0].text(0.,  T1.max()-30, 'beta_h ={:.3f}'.format(round(fit_par[2,0],3))+'±{:.3f}'.format(round(fit_par[2,1],3)), fontsize=6)
            axes[0, 0].text(0.,  T1.max()-40, 'beta_l ={:.3f}'.format(round(fit_par[3,0],3))+'±{:.3f}'.format(round(fit_par[3,1],3)), fontsize=6)
            # axes[0, 0].text(0.,  T1.max()-50, 'alpha ={:.3f}'.format(round(fit_par[10,0],3))+'±{:.3f}'.format(round(fit_par[10,1],3)), fontsize=6)
            axes[0, 0].legend(loc='upper right')

            axes[0, 1].plot(T2, '.', label="Ch_1")
            axes[0, 1].plot(model1(El, fit_par[0,0], fit_par[1,0], fit_par[4,0], fit_par[5,0]), linewidth=1, label="Fit")
            axes[0, 1].text(0.,  T2.max()-10, 'T_h ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=6)
            axes[0, 1].text(0.,  T2.max()-20, 'T_l ={:.2f}'.format(round(fit_par[1,0],2))+'±{:.2f}'.format(round(fit_par[1,1],2)), fontsize=6)
            axes[0, 1].text(0.,  T2.max()-30, 'beta_h ={:.3f}'.format(round(fit_par[4,0],3))+'±{:.3f}'.format(round(fit_par[4,1],3)), fontsize=6)
            axes[0, 1].text(0.,  T2.max()-40, 'beta_l ={:.3f}'.format(round(fit_par[5,0],3))+'±{:.3f}'.format(round(fit_par[5,1],3)), fontsize=6)
            axes[0, 1].legend(loc='upper right')

            axes[1, 0].plot(T3, '.', label="Ch_2")
            axes[1, 0].plot(model1(El, fit_par[0,0], fit_par[1,0], fit_par[6,0], fit_par[7,0]), linewidth=1, label="Fit")
            axes[1, 0].text(0.,  T3.max()-10, 'T_h ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=6)
            axes[1, 0].text(0.,  T3.max()-20, 'T_l ={:.2f}'.format(round(fit_par[1,0],2))+'±{:.2f}'.format(round(fit_par[1,1],2)), fontsize=6)
            axes[1, 0].text(0.,  T3.max()-30, 'beta_h ={:.3f}'.format(round(fit_par[6,0],3))+'±{:.3f}'.format(round(fit_par[6,1],3)), fontsize=6)
            axes[1, 0].text(0.,  T3.max()-40, 'beta_l ={:.3f}'.format(round(fit_par[7,0],3))+'±{:.3f}'.format(round(fit_par[7,1],3)), fontsize=6)
            axes[1, 0].legend(loc='upper right')

            axes[1, 1].plot(T3, '.', label="Ch_3")
            axes[1, 1].plot(model1(El, fit_par[0,0], fit_par[1,0], fit_par[8,0], fit_par[9,0]), linewidth=1, label="Fit")
            axes[1, 1].text(0.,  T3.max()-10, 'T_h ={:.2f}'.format(round(fit_par[0,0],2))+'±{:.2f}'.format(round(fit_par[0,1],2)), fontsize=6)
            axes[1, 1].text(0.,  T3.max()-20, 'T_l ={:.2f}'.format(round(fit_par[1,0],2))+'±{:.2f}'.format(round(fit_par[1,1],2)), fontsize=6)
            axes[1, 1].text(0.,  T3.max()-30, 'beta_h ={:.3f}'.format(round(fit_par[8,0],3))+'±{:.3f}'.format(round(fit_par[8,1],3)), fontsize=6)
            axes[1, 1].text(0.,  T3.max()-40, 'beta_l ={:.3f}'.format(round(fit_par[9,0],3))+'±{:.3f}'.format(round(fit_par[9,1],3)), fontsize=6)
            axes[1, 1].legend(loc='upper right')

            pl.show()

        return FH, fit_par

    def extract_pwv(self, FH, beta):
        f=self.freqs
        El=np.array(FH[:,5])
        print(El)
        tau=np.zeros((len(beta),len(El)))
        for i in range(len(beta)):
            tau[i,:]=np.exp(-beta[i]/np.sin((El)*np.pi/180.))
            pl.plot(El,tau[i,:],label='Ch_'+str(i))
        pl.legend()
        pl.xlabel('El[deg]')
        pl.ylabel('Tau')
        pl.show()

        for e in range (len(El)):
            tau_doubleband=np.array([tau[3,e],tau[2,e],tau[1,e],tau[0,e],tau[0,e],tau[1,e],tau[2,e],tau[3,e]])
            pl.plot(f,tau_doubleband,label='El='+str(El[e]))
        pl.legend()
        pl.xlabel('freq[GHz]')
        pl.ylabel('Tau')
        pl.show()



    def plotElnod(self, filename, el_range=(0,90)):
        freqs_fit=np.linspace(170,190,2000)
        fit_par, FH = self.fit_elnod(filename)
        for i in range (0, len(fit_par[0,:])):
            T=[FH[i][3],FH[i][2],FH[i][1],FH[i][0],FH[i][0],FH[i][1],FH[i][2],FH[i][3]]


        pl.show()
