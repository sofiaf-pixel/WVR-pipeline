import os, sys
from pylab import *
from scipy.optimize import curve_fit,least_squares
import scipy.integrate as integrate
from scipy.integrate import simps
import scipy.interpolate as sint
import numpy as np
import pickle
import pylab as pl
from matplotlib import dates
import datetime
import scipy.optimize as sp
from pathlib import Path
import plotElnod as pE
import plotAzscan as pA
from time import perf_counter
import glob
import re as r

re=pE.ReadElnod()
raz=pA.ReadAzscan()




    def plot_full_folder(self, folder_name, picklefn='', template= 'SPole_annual_50.amc', spline=0):
        full_folder_pwv = {'filename': [], 'day':[], 'time':[], 'pwv':[], 'pwv_err':[], 'Nscale':[], 'Nscale_err':[], 'fit_time':[]}
        full_folder_pwv_red = {'filename': [], 'day':[], 'time':[], 'pwv':[], 'pwv_err':[], 'Nscale':[], 'Nscale_err':[], 'fit_time':[]}
        print('start')
        path='wvr1_data/'+folder_name+'/'
        print(path)
        if not os.path.exists(path):
            print('ERROR. Path does not exist')
        for folderlist in os.listdir(path):
            if folderlist[-6:]=='skyDip':
                for filelist in os.listdir(f'{path+folderlist}'):
                    if filelist[-8:]=='fast.txt':
                        print('creating am datfile')
                        print(path+folderlist+'/'+filelist)
                        fit_output = {'filename': [], 'El': [], 'Nscale': [], 'Nscale_err': [], 'pwv':[], 'pwv_err':[]} #defined just not to delete all the parts that call that
                        #self.create_am_datfile(filelist, path_to_data=folder_name+'/'+folderlist+'/') #it's included in fit_w_am
                        #print(fit_output)
            #            try:
                        #if not os.path.exists('am_datfiles/'+template[:-4]+'/'+filelist[:-4]???+'.ams'):
                        #t1=perf_counter()
                        #fit_output, zenith=self.fit_w_am(filelist, path_to_data=folder_name+'/'+folderlist+'/', template=template, spline=0)
                        #t2=perf_counter()
                        #print('fit time=', t2-t1)
                        #full_folder_pwv['fit_time'].append(t2-t1)

                        t1_red=perf_counter()
                        fit_output_red=self.fit_w_am(filelist, path_to_data=folder_name+'/'+folderlist+'/', template=template, spline=spline)
                        t2_red=perf_counter()
                        print('fit time red=', t2_red-t1_red)
                        full_folder_pwv_red['fit_time'].append(t2_red-t1_red)

                        pwv=np.mean(fit_output['pwv'])
                        pwv_err=np.std(fit_output['pwv'])
                        full_folder_pwv['filename'].append(filelist)
                        full_folder_pwv['day'].append(filelist[:8])
                        full_folder_pwv['time'].append(filelist[9:15])
                        full_folder_pwv['pwv'].append(pwv)
                        full_folder_pwv['pwv_err'].append(pwv_err)

                        pwv_red=np.mean(fit_output_red['pwv'])
                        pwv_err_red=np.std(fit_output_red['pwv'])
                        full_folder_pwv_red['filename'].append(filelist)
                        full_folder_pwv_red['day'].append(filelist[:8])
                        full_folder_pwv_red['time'].append(filelist[9:15])
                        full_folder_pwv_red['pwv'].append(pwv_red)
                        full_folder_pwv_red['pwv_err'].append(pwv_err_red)

                        print('filename:', filelist)
                        print('day:', filelist[:8])
                        print('time:', filelist[9:15])
                        print('pwv:', pwv)
                        print('pwv_err=', pwv_err)


                        #pl.scatter(fit_output['El'], fit_output['pwv'], marker='.', c='k', alpha=0.6)
                        #pl.axhline(y=pwv, color='k', linestyle='--', label='pwv='+str(round(pwv,2))+'+-'+str(round(pwv_err,2)))
                        pl.scatter(fit_output_red['El'], fit_output_red['pwv'], marker='*', c='r')
                        pl.axhline(y=pwv_red, color='r', linestyle='--', label='pwv_fit='+str(round(pwv_red,2))+'+-'+str(round(pwv_err_red,2)))
                        pl.suptitle('PWV - One Elnod')
                        pl.title(filelist[:-4])
                        pl.xlabel('El')
                        pl.ylabel('pwv[um]')

                        pl.legend()
                        pl.savefig('fit_output_plots/'+filelist+'_pwvscatter_corr.png')
                        pl.close()


            #            except:
            #                print('ERROR for file '+filelist)

        print(full_folder_pwv_red)
        f=open(picklefn, "wb")
        pickle.dump(full_folder_pwv_red, f)
        f.close()
        return full_folder_pwv_red





def plot_Apr_pwv(self, fn):
    f = open("pwv_april.txt","rb")
    dict = pickle.load(f)

    filename=dict["filename"]
    pwv=dict["pwv"]
    pwv_err=dict["pwv_err"]

    #f_new=['']

    # for i in range (len(f)):
    #     string=str(f[i])
    #     f_new.append(string[4:8]+'_'+string[8:10])
    #     print(f_new[i])
    #
    # for i in range (len(f_new)):
    #     f_new[i]=f_new[i+1]
    #
    # f_new=fnew[0:(len(f_new-1))]
    #
    # for i in range (len(f_new)):
    #     string=f_new[i]
    #     print(string)
    #     date=string[0:2]+'-'+string[2:4]
    #     print(date)

    for i in range(len(filename)):
        filename[i]='2020-'+filename[i]

    converted_dates = list(map(datetime.datetime.strptime, datelist, len(datelist)*['%Y-%m-%d']))
    x_axis = converted_dates

    fig, ax = pl.subplots()
    pl.errorbar(x_axis, pwv, yerr=pwv_err, fmt='.')
    pl.text(0.95, 0.90, 'pwv_median[um]='+str(round(np.median(pwv),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
    pl.text(0.95, 0.85, 'pwv_std[um]='+str(round(np.std(pwv),2)), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
    pl.xticks(rotation=45)
    pl.title('PWV variations as measured with the WVR')
    pl.xlabel('date')
    pl.ylabel('pwv[um]')

    pl.show() #make it bigger if text not visible
