import urllib
from pylab import *
from datetime import datetime
import pylab as pl
import pickle as pk


# download and read hourly files
#year = arange(1975,2020)
wind_dict={'t':[], 'ws':[], 'wd':[]}
year=[2020]
ws=[]
wd=[]
t=[]
for y in year:
    fname = 'https://gml.noaa.gov/aftp/data/meteorology/in-situ/spo/met_spo_insitu_1_obop_hour_%s.txt'%y
    print(y,fname)
    try:
        r = urllib.request.urlopen(fname)
    except:
        print('file missing')
    lines = r.read().decode('utf-8')
    for k,line in enumerate(lines.split('\n')[:-1]):
        print(k,line)
        sline=line.split()
        print(sline[6])
        wd.append(float(sline[5]))
        ws.append(float(sline[6]))# maybe also this has to be 7
        t.append(datetime.strptime('%s %s %s %s'%(sline[1],sline[2],sline[3],sline[4]),'%Y %m %d %H'))

wind_dict['t']=t
wind_dict['ws']=ws
wind_dict['wd']=wd

pkfn='SP_windspeed_042020_hour.pk'
f=open(pkfn, "wb")
pk.dump(wind_dict, f)
f.close()


wind_dict={'t':[], 'ws':[], 'wd':[]}
#year = arange(1979,2020)
year=[2020]
ws=[]
wd=[]
t=[]
for y in year:
    for mo in arange(4,5):
    #for mo in arange(1,13):
        fname = 'https://gml.noaa.gov/aftp/data/meteorology/in-situ/spo/%s/met_spo_insitu_1_obop_minute_%s_%02d.txt'%(y,y,mo)
        print(y,mo,fname)
        try:
            r = urllib.request.urlopen(fname)
        except:
            print('file missing')
        lines = r.read().decode('utf-8')
        for k,line in enumerate(lines.split('\n')[:-1]):
            sline=line.split()
            print(sline[7])
            wd.append(float(sline[5]))
            #ws.append(float(sline[6]))
            ws.append(float(sline[7]))
            t.append(datetime.strptime('%s %s %s %s %s'%(sline[1],sline[2],sline[3],sline[4],sline[5]),'%Y %m %d %H %M'))

pl.plot(t,ws)
pl.show()

wind_dict['t']=t
wind_dict['ws']=ws
wind_dict['wd']=wd

pkfn='SP_windspeed_042020.pk'
f=open(pkfn, "wb")
pk.dump(wind_dict, f)
f.close()
