import AM_functions as am
import pickle
import numpy as np
import pylab as pl

x_am=am.AM_functions()

wvr_scan='20200527_020002_skyDip_fast.txt'
#
# data=x_am.fit_w_am(wvr_scan)
#
# print(data)
#


template= 'SPole_annual_50.amc'
path='am_datfiles/'+template[:-4]+'/'+wvr_scan[:-4]
pickle_fn=path+'/'+wvr_scan[:-4]+'_fitoutput_corr.txt'


f = open(pickle_fn,'rb')
data=pickle.load(f)
f.close()

el=np.array(data['El'])
PWV_z=np.array(data['pwv_tropo'])
PWV_los=np.array(data['pwv_los_total'])


pl.plot(el, PWV_z, label='pwv zenith')
pl.plot(el, PWV_los, label='pwv line of sight')
pl.legend()
pl.title(wvr_scan[:-9])
pl.suptitle('Skydip PWV')
pl.show()
