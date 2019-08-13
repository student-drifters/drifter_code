"""
s4_plt.py

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
DESTINDIR='driftfvcom_data4/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])

k=0
try:
    P = np.genfromtxt('kpic.txt',dtype=None,names=['INDEX'],delimiter=',')
    kpic=P['INDEX']
    k=kpic
except IOError:
    k=0


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

plt.close('all')
plt.figure(1)
plt.plot(londh,latdh,'r.-')
plt.plot(londh*flag,latdh*flag,'m.-')
flag=flag*flag_v
plt.plot(londh*flag,latdh*flag,'b.-')
plt.title(FN1)
plt.show()

i=np.argwhere(flag == 1).flatten()
if len(i)>0:
    FN2=DESTINDIR+FN
    print FN2
    #        np.savez(FN2,tdh=tdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh,umoh=umoh,vmoh=vmoh)
    shutil.copyfile(FN1,FN2)


"""
plt.figure()
plt.plot(tdhgap,udh,'b.-',tdhgap,umoh,'r.-',tdhgap,udm,'k-',tdhgap,umom,'g-')
plt.title('U '+FN0)
plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
plt.show()

plt.figure()
plt.plot(tdhgap,vdh,'b.-',tdhgap,vmoh,'r.-',tdhgap,vdm,'k-',tdhgap,vmom,'g-')
plt.title('V '+FN0)
plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
plt.show()
"""

k=k+1
f=open('kpic.txt','w')
f.write(str(k))
f.close()
