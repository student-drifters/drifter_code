"""
s5_windslip.py

analyze wind slip of drifters
remove tide and find gaps in data

@author: Vitalii Sheremet, FATE Project
"""

# -*- coding: utf-8 -*-
import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
#from scipy import interpolate
import sys
from SeaHorseTide import *
import shutil
import matplotlib.mlab as mlab
import matplotlib.cm as cm


###########################################################

SOURCEDIR='driftfvcom_data4/'
DESTINDIR='driftfvcom_data4/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])

lath=np.array([])
lonh=np.array([])
th=np.array([])
flagh=np.array([])
udr=np.array([])
vdr=np.array([])
umo=np.array([])
vmo=np.array([])
uws=np.array([])
vws=np.array([])

#for k in range(10):
for k in range(len(FNs)):
    FN=FNs[k]
    #ID_19965381.npz
    FN1=SOURCEDIR+FN
#    print k, FN1
    Z=np.load(FN1) 
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    umoh=Z['umoh'];vmoh=Z['vmoh'];
    tgap=Z['tgap'];flag=Z['flag'];
    udm=Z['udm'];vdm=Z['vdm'];
    udti=Z['udti'];vdti=Z['vdti'];
    umom=Z['umom'];vmom=Z['vmom'];
    umoti=Z['umoti'];vmoti=Z['vmoti'];
#    uvwsh=Z['uwsh'];vwsh=Z['vwsh'];
    uwsm=Z['uwsm'];vwsm=Z['vwsm'];
    Z.close()
    
    lath=np.append(lath,latdh)
    lonh=np.append(lonh,londh)
    th=np.append(th,tdh)
    flagh=np.append(flagh,flag)

    udr=np.append(udr,udm*flag)            
    vdr=np.append(vdr,vdm*flag)
    umo=np.append(umo,umom*flag)            
    vmo=np.append(vmo,vmom*flag)
    uws=np.append(uws,uwsm)
    vws=np.append(vws,vwsm)
                
    if len(np.argwhere(umom==0.0).flatten())>21:
        print 'zeros',k,FN1            
  
i=np.argwhere(np.isnan(udr-umo)==False).flatten()
lath=lath[i]
lonh=lonh[i]
th=th[i]
udr=udr[i]
vdr=vdr[i]
umo=umo[i]
vmo=vmo[i]
uws=uws[i]
vws=vws[i]

np.save('lonh.npy',lonh)
np.save('lath.npy',lath)
np.save('th.npy',th)

np.save('umo.npy',umo)
np.save('vmo.npy',vmo)
np.save('udr.npy',udr)
np.save('vdr.npy',vdr)
np.save('uws.npy',uws)
np.save('vws.npy',vws)


lonh=np.load('lonh.npy')
lath=np.load('lath.npy')
th  =np.load('th.npy')

umo=np.load('umo.npy')
vmo=np.load('vmo.npy')
udr=np.load('udr.npy')
vdr=np.load('vdr.npy')
uws=np.load('uws.npy')
vws=np.load('vws.npy')

veldiff=(umo-udr)+(vmo-vdr)*1j
wind=uws+vws*1j


plt.figure()
plt.plot(udr,umo,'b.')
plt.xlabel('u drifter')
plt.ylabel('u model')
plt.show()

plt.figure()
plt.plot(vdr,vmo,'b.')
plt.xlabel('v drifter')
plt.ylabel('v model')
plt.show()


plt.figure()
plt.plot(veldiff.real,veldiff.imag,'b.')
plt.xlabel('u diff')
plt.ylabel('v diff')
plt.show()



plt.figure()
plt.plot(np.abs(wind),np.abs(veldiff),'b.')
plt.xlabel('wind stress')
plt.ylabel('vel diff')
plt.show()

