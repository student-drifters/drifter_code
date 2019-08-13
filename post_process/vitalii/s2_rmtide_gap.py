"""
s2_rmtide_gap.py

remove tide and find gaps in data

compare drifter dataset velocity with fvcom hincast
following the drifter trajectory.
fvcom data pydap access from monthly files.

To do: 
1. spacial interpolation of fvcom velocity data
currently, the nearest neighbor
@author: Vitalii Sheremet, FATE Project
"""

# -*- coding: utf-8 -*-
import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
from scipy import interpolate
import sys
from SeaHorseTide import *
import shutil

SOURCEDIR='driftfvcom_data3/'

k=0
FNs=['ID_100390731.npz']

FN=FNs[k]
#ID_19965381.npz
FN1=SOURCEDIR+FN
print k, FN1
Z=np.load(FN1) 
tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
udh=Z['udh'];vdh=Z['vdh'];
umoh=Z['umoh'];vmoh=Z['vmoh'];
tgap=Z['tgap'];flag=Z['flag'];
udm=Z['udm'];vdm=Z['vdm'];
udti=Z['udti'];vdti=Z['vdti'];
umom=Z['umom'];vmom=Z['vmom'];
umoti=Z['umoti'];vmoti=Z['vmoti'];
Z.close()
        
flag=umoh*0.+1. # 1 inside GOM3 NaN outside

flag_gap=tdh*0.+1.
i=np.argwhere(tgap>12./24.) # exceeds 12h
flag_gap[i]=np.NaN

flag=flag*flag_gap

v=np.sqrt(udh*udh+vdh*vdh)    
i=np.argwhere(v>2.7) # 2 m/s
flag_v=tdh*0.+1.
flag_v[i]=np.NaN

latdhf=latdh*flag
i=np.argwhere(~np.isnan(latdhf)).flatten()
latref=latdhf[i].mean()
A=1./np.cos(latref*np.pi/180.) # axes aspect ratio

plt.close('all')
plt.figure()
plt.plot(londh,latdh,'r.-')
plt.plot(londh*flag,latdh*flag,'m.-')
flag=flag*flag_v
plt.plot(londh*flag,latdh*flag,'b.-')
plt.axes().set_aspect(A)
plt.grid('on')
plt.title(FN)
plt.show()

td0=datetime.fromordinal(tdh[0].astype('int'))
t0=datetime(td0.year,1,1).toordinal()-1# beginning of year as day 1   
plt.figure()
plt.plot(tdh-t0,udh,'b.-',tdh-t0,umoh,'r.-',tdh-t0,udm,'k-',tdh-t0,umom,'g-')
plt.title('U '+FN)
plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
plt.xlabel('Yearday of '+str(td0.year)+',UTC')
plt.grid('on')
plt.show()

plt.figure()
plt.plot(tdh-t0,vdh,'b.-',tdh-t0,vmoh,'r.-',tdh-t0,vdm,'k-',tdh-t0,vmom,'g-')
plt.title('V '+FN)
plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
plt.xlabel('Yearday of '+str(td0.year)+',UTC')
plt.grid('on')
plt.show()
    
plt.figure()
sepx=(udh-umoh)*24.*60*60/1000. # km/day
sepy=(vdh-vmoh)*24.*60*60/1000. # km/day
sepmx=(udm-umom)*24.*60*60/1000. # km/day
sepmy=(vdm-vmom)*24.*60*60/1000. # km/day

N=len(sepx)
M=24
sepx1d=sepx[0:N/M*M].reshape((N/M,M)).mean(1)
sepy1d=sepy[0:N/M*M].reshape((N/M,M)).mean(1)
sep1d=np.sqrt(sepx1d*sepx1d+sepy1d*sepy1d)

plt.plot(tdh-t0,sepx,'g.',tdh-t0,sepy,'b.')
plt.plot(tdh-t0,sepmx,'g-',tdh-t0,sepmy,'b-')
plt.plot((tdh-t0)[0:N/M*M:M],sep1d,'r-o')
plt.title('Separation km/day '+FN)
plt.legend(['X delta','Y delta','X rmtide','Y rmtide','Separ 1day'])
plt.xlabel('Yearday of '+str(td0.year)+',UTC')
plt.ylabel('km')
plt.grid('on')
plt.show()


