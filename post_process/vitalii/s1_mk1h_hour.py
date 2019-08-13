"""
s1_mk1h_hour.py
compare drifter dataset velocity with fvcom hincast
following the drifter trajectory

step 1
fvcom data from local hourly files

To do: 
1. spacial interpolation of fvcom velocity data
currently, the nearest neighbor
fixed to polygonal baricentric coordinate interpolation

before running
reset kdhour.txt file

@author: Vitalii Sheremet, FATE Project
"""

import numpy as np
#import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
#from scipy import interpolate
import sys
from SeaHorseTide import *


def nearlonlat(lon,lat,lonp,latp):
    """
i=nearlonlat(lon,lat,lonp,latp)
find the closest node in the array (lon,lat) to a point (lonp,latp)
input:
lon,lat - np.arrays of the grid nodes, spherical coordinates, degrees
lonp,latp - point on a sphere
output:
i - index of the closest node
min_dist - the distance to the closest node, degrees
For coordinates on a plane use function nearxy

Vitalii Sheremet, FATE Project
"""
    cp=np.cos(latp*np.pi/180.)
# approximation for small distance
    dx=(lon-lonp)*cp
    dy=lat-latp
    dist2=dx*dx+dy*dy
# dist1=np.abs(dx)+np.abs(dy)
    i=np.argmin(dist2)
#    min_dist=np.sqrt(dist2[i])
    return i 


def polygonal_barycentric_coordinates(xp,yp,xv,yv):
    """
Calculate generalized barycentric coordinates within an N-sided polygon.

    w=polygonal_barycentric_coordinates(xp,yp,xv,yv)
    
    xp,yp - a point within an N-sided polygon
    xv,yv - vertices of the N-sided polygon, length N
    w     - polygonal baricentric coordinates, length N,
            normalized w.sum()=1
   
Used for function interpolation:
    fp=(fv*w).sum()
    where fv - function values at vertices,
    fp the interpolated function at the point (xp,yp)
    
Vitalii Sheremet, FATE Project    
    """
    N=len(xv)   
    j=np.arange(N)
    ja=(j+1)%N # next vertex in the sequence 
    jb=(j-1)%N # previous vertex in the sequence
# area of the chord triangle j-1,j,j+1
    Ajab=np.cross(np.array([xv[ja]-xv[j],yv[ja]-yv[j]]).T,np.array([xv[jb]-xv[j],yv[jb]-yv[j]]).T) 
# area of triangle p,j,j+1
    Aj=np.cross(np.array([xv[j]-xp,yv[j]-yp]).T,np.array([xv[ja]-xp,yv[ja]-yp]).T)  

# In FVCOM A is O(1.e7 m2) .prod() may result in inf
# to avoid this scale A
    Aj=Aj/max(abs(Aj))
    Ajab=Ajab/max(abs(Ajab))
    
    w=xv*0.
    j2=np.arange(N-2)
    
    for j in range(N):
# (j2+j+1)%N - list of triangles except the two adjacent to the edge pj
# For hexagon N=6 j2=0,1,2,3; if j=3  (j2+j+1)%N=4,5,0,1
        w[j]=Ajab[j]*Aj[(j2+j+1)%N].prod()
# timing [s] per step:  1.1976 1.478
# timing [s] per step:  1.2048 1.4508 
        
    
#    w=np.array([Ajab[j]*Aj[(j2+j+1)%N].prod() for j in range(N)])
# timing [s] per step:  1.2192 1.4572
# list comprehension does not affect speed

# normalize w so that sum(w)=1       
    w=w/w.sum() 
       
    return w,Aj


def VelInterp_lonlat(lonp,latp,Grid,u,v):
    """
Velocity interpolating function

    ui,vi,kv=VelInterp_lonlat(lonp,latp,Grid,u,v)
    
    lonp,latp - arrays of points where the interpolated velocity is desired
    Grid - parameters of the triangular grid
    u,v - velocity field defined at the triangle baricenters
    
    urci - interpolated u/cos(lat)
    ui   - intepolated u
    vi   - interpolated v
    The Lame coefficient cos(lat) of the spherical coordinate system
    is needed for RungeKutta4_lonlat: dlon = u/cos(lat)*tau, dlat = vi*tau
    kv - the vertex corresponding to the dual mesh polygon
    
    """
    
# find the nearest vertex    
    kv=nearlonlat(Grid['lon'],Grid['lat'],lonp,latp)
#    print kv
# list of triangles surrounding the vertex kv    
    kfv=Grid['kfv'][0:Grid['nfv'][kv],kv]
#    print kfv
# coordinates of the (dual mesh) polygon vertices: the centers of triangle faces
    lonv=Grid['lonc'][kfv];latv=Grid['latc'][kfv] 
    w,Aj=polygonal_barycentric_coordinates(lonp,latp,lonv,latv)
# baricentric coordinates are invariant wrt coordinate transformation (xy - lonlat), check! 

# Check whether any Aj are negative, which would mean that a point is outside the polygon.
# Otherwise, the polygonal interpolation will not be continous.
# This check is not needed if the triangular mesh and its dual polygonal mesh
# are Delaunay - Voronoi. 

# normalize subareas by the total area 
# because the area sign depends on the mesh orientation.    
    Aj=Aj/Aj.sum()
    if np.argwhere(Aj<0).flatten().size>0:
# if point is outside the polygon try neighboring polygons
#        print kv,kfv,Aj
        for kv1 in Grid['kvv'][0:Grid['nvv'][kv],kv]:
            kfv1=Grid['kfv'][0:Grid['nfv'][kv1],kv1]
            lonv1=Grid['lonc'][kfv1];latv1=Grid['latc'][kfv1] 
            w1,Aj1=polygonal_barycentric_coordinates(lonp,latp,lonv1,latv1)
            Aj1=Aj1/Aj1.sum()
            if np.argwhere(Aj1<0).flatten().size==0:
                w=w1;kfv=kfv1;kv=kv1;Aj=Aj1
#                print kv,kfv,Aj

# Now there should be no negative w
# unless the point is outside the triangular mesh
    if np.argwhere(w<0).flatten().size>0:
#        print kv,kfv,w
        
# set w=0 -> velocity=0 for points outside 
#        w=w*0.
        kv=-1        
# set w=NaN -> velocity=NaN for points outside 
        w=w*np.NaN
        kv=-1        


# interpolation within polygon, w - normalized weights: w.sum()=1.    
# use precalculated Lame coefficients for the spherical coordinates
# coslatc[kfv] at the polygon vertices
# essentially interpolate u/cos(latitude)
# this is needed for RungeKutta_lonlat: dlon = u/cos(lat)*tau, dlat = vi*tau

# In this version the resulting interpolated field is continuous, C0.
#    cv=Grid['coslatc'][kfv]    
#    urci=(u[kfv]/cv*w).sum()
#    vi=(v[kfv]*w).sum()

# cos factor is not needed here
    ui=(u[kfv]*w).sum()
    vi=(v[kfv]*w).sum()
        
    return ui,vi,kv
    
def RataDie(yr,mo=1,da=1):
    """

RD = RataDie(yr,mo=1,da=1,hr=0,mi=0,se=0)
RD = RataDie(yr,mo=1,da=1)

returns the serial day number in the (proleptic) Gregorian calendar
or elapsed time in days since 0001-01-00.

Vitalii Sheremet, SeaHorse Project, 2008-2013.
"""
#
#    yr+=(mo-1)//12;mo=(mo-1)%12+1; # this extends mo values beyond the formal range 1-12
    RD=367*yr-(7*(yr+((mo+9)//12))//4)-(3*(((yr+(mo-9)//7)//100)+1)//4)+(275*mo//9)+da-396;
    return RD    
    
    

###############################################################################################
# FVCOM GOM3 triangular grid
"""
from pydap.client import open_url
from netCDF4 import Dataset
URL='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
#http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?
#a1u[0:1:3][0:1:90414],a2u[0:1:3][0:1:90414],art1[0:1:48450],art2[0:1:48450],
#aw0[0:1:2][0:1:90414],awx[0:1:2][0:1:90414],awy[0:1:2][0:1:90414],cc_hvc[0:1:90414],
#h[0:1:48450],lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],
#nbe[0:1:2][0:1:90414],nbsn[0:1:10][0:1:48450],nbve[0:1:8][0:1:48450],
#nn_hvc[0:1:48450],nprocs,ntsn[0:1:48450],ntve[0:1:48450],nv[0:1:2][0:1:90414],
#partition[0:1:90414],siglay[0:1:44][0:1:48450],siglev[0:1:45][0:1:48450],
#x[0:1:48450],xc[0:1:90414],y[0:1:48450],yc[0:1:90414],z0b[0:1:90414],
#Itime[0:1:171882],Itime2[0:1:171882],Times[0:1:171882],file_date[0:1:171882],
#iint[0:1:171882],kh[0:1:171882][0:1:45][0:1:48450],
#km[0:1:171882][0:1:45][0:1:48450],kq[0:1:171882][0:1:45][0:1:48450],
#l[0:1:171882][0:1:45][0:1:48450],net_heat_flux[0:1:171882][0:1:48450],
#omega[0:1:171882][0:1:45][0:1:48450],q2[0:1:171882][0:1:45][0:1:48450],
#q2l[0:1:171882][0:1:45][0:1:48450],salinity[0:1:171882][0:1:44][0:1:48450],
#short_wave[0:1:171882][0:1:48450],temp[0:1:171882][0:1:44][0:1:48450],
#time[0:1:171882],u[0:1:171882][0:1:44][0:1:90414],ua[0:1:171882][0:1:90414],
#uwind_stress[0:1:171882][0:1:90414],v[0:1:171882][0:1:44][0:1:90414],
#va[0:1:171882][0:1:90414],vwind_stress[0:1:171882][0:1:90414],
#ww[0:1:171882][0:1:44][0:1:90414],zeta[0:1:171882][0:1:48450]

#ds=open_url(URL)                 # pydap version 
ds = Dataset(URL,'r').variables   # netCDF4 version

#xxx=ds['xxx']; np.save('gom3.xxx.npy',np.array(xxx))
a1u=ds['a1u']; np.save('gom3.a1u.npy',np.array(a1u))
a2u=ds['a2u']; np.save('gom3.a2u.npy',np.array(a2u))
art1=ds['art1']; np.save('gom3.art1.npy',np.array(art1))
art2=ds['art2']; np.save('gom3.art2.npy',np.array(art2))
aw0=ds['aw0']; np.save('gom3.aw0.npy',np.array(aw0))
awx=ds['awx']; np.save('gom3.awx.npy',np.array(awx))
awy=ds['awy']; np.save('gom3.awy.npy',np.array(awy))
cc_hvc=ds['cc_hvc']; np.save('gom3.cc_hvc.npy',np.array(cc_hvc))
    
h=ds['h']; np.save('gom3.h.npy',np.array(h))

lat=ds['lat']; np.save('gom3.lat.npy',np.array(lat))
lon=ds['lon']; np.save('gom3.lon.npy',np.array(lon))
latc=ds['latc']; np.save('gom3.latc.npy',np.array(latc))
lonc=ds['lonc']; np.save('gom3.lonc.npy',np.array(lonc))

nbe=ds['nbe']; np.save('gom3.nbe.npy',np.array(nbe))
nbsn=ds['nbsn']; np.save('gom3.nbsn.npy',np.array(nbsn))
nbve=ds['nbve']; np.save('gom3.nbve.npy',np.array(nbve))
nn_hvc=ds['nn_hvc']; np.save('gom3.nn_hvc.npy',np.array(nn_hvc))
nprocs=ds['nprocs']; np.save('gom3.nprocs.npy',np.array(nprocs))
ntsn=ds['ntsn']; np.save('gom3.ntsn.npy',np.array(ntsn))
ntve=ds['ntve']; np.save('gom3.ntve.npy',np.array(ntve))
nv=ds['nv']; np.save('gom3.nv.npy',np.array(nv))
partition=ds['partition']; np.save('gom3.partition.npy',np.array(partition))
siglay=ds['siglay']; np.save('gom3.siglay.npy',np.array(siglay))
siglev=ds['siglev']; np.save('gom3.siglev.npy',np.array(siglev))

x=ds['x']; np.save('gom3.x.npy',np.array(x))
xc=ds['xc']; np.save('gom3.xc.npy',np.array(xc))
y=ds['y']; np.save('gom3.y.npy',np.array(y))
yc=ds['yc']; np.save('gom3.yc.npy',np.array(yc))
"""    
    
x=np.load('gom3.x.npy')
y=np.load('gom3.y.npy')
xc=np.load('gom3.xc.npy')
yc=np.load('gom3.yc.npy')

lon=np.load('gom3.lon.npy')
lat=np.load('gom3.lat.npy')
lonc=np.load('gom3.lonc.npy')
latc=np.load('gom3.latc.npy')

# precalculate Lame coefficients for the spherical coordinates
coslat=np.cos(lat*np.pi/180.)
coslatc=np.cos(latc*np.pi/180.)


#nv: Array of 32 bit Integers [three = 0..2][nele = 0..90414] 
#long_name: nodes surrounding element
#standard_name: face_node_connectivity
#start_index: 1
nv=np.load('gom3.nv.npy')
nv-=1 # convert from FORTRAN to python 0-based indexing
#kvf=nv

#nbe: Array of 32 bit Integers [three = 0..2][nele = 0..90414] 
# long_name: elements surrounding each element
nbe=np.load('gom3.nbe.npy')
nbe-=1 # convert from FORTRAN to python 0-based indexing
#kff=nbe

#nbsn: Array of 32 bit Integers [maxnode = 0..10][node = 0..48450]
#long_name: nodes surrounding each node
 # list of nodes surrounding a given node, 1st and last entries identical to make a closed loop
nbsn=np.load('gom3.nbsn.npy')
nbsn-=1 # convert from FORTRAN to python 0-based indexing
#kvv=nbsn

#ntsn: Array of 32 bit Integers [node = 0..48450]
#long_name: #nodes surrounding each node
 # the number of nodes surrounding a given node + 1, because 1st and last entries identical to make a closed loop
ntsn=np.load('gom3.ntsn.npy')
#nvv=ntsn

#nbve: Array of 32 bit Integers [maxelem = 0..8][node = 0..48450] 
#long_name: elems surrounding each node
# list of elements surrounding a given node, 1st and last entries identical to make a closed loop
nbve=np.load('gom3.nbve.npy')
nbve-=1 # convert from FORTRAN to python 0-based indexing
#kfv=nbve

#ntve: Array of 32 bit Integers [node = 0..48450] 
#long_name: #elems surrounding each node
# the number of elements surrounding a given node + 1, because 1st and last entries identical to make a closed loop
ntve=np.load('gom3.ntve.npy')
#nfv=ntve

Grid={'x':x,'y':y,'xc':xc,'yc':yc,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,'coslat':coslat,'coslatc':coslatc,'kvf':nv,'kff':nbe,'kvv':nbsn,'nvv':ntsn,'kfv':nbve,'nfv':ntve}


######################################

FList = np.genfromtxt('drift_data/FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])
kdr=0
try:
    P = np.genfromtxt('kdrhour.txt',dtype=None,names=['INDEX'],delimiter=',')
    kdr=P['INDEX']
#    kdr = np.load('kdr.npy')
    print 'continue with drifter '+str(kdr)
except IOError:
    kdr=0
    print 'start with drifter '+str(kdr)

while kdr in range(len(FNs)):
#while kdr in range(2):
     
    FN='drift_data/'+FNs[kdr]
#    FN='drift_data/ID_96101.csv'
    print kdr, FN
# load file 'ID_19954471.csv' from cleaned ascii dataset
    """
    # ID - drifter ID
    # TimeRD - RataDie, serial day number, days since 0001-00-00 00:00:00, GMT
    # TIME_GMT - timestamp, yyyy-mm-dd, GMT
    # YRDAY0_GMT - yearday, zero based indexing, days since Jan-01 00:00:00, GMT
    # LON_DD - longitude, degrees east
    # LAT_DD - latitude, degrees north
    # TEMP - temperature, degrees Celsius
    # DEPTH_I - depth of instrument drogue, meters
    #Vars: ID, TimeRD, TIME_GMT, YRDAY0_GMT, LON_DD, LAT_DD, TEMP, DEPTH_I
    19954471,728387.5687,1995-04-04,93.5687,-67.39,44.167,9.99,40.0
    19954471,728387.9056,1995-04-04,93.9056,-67.429,44.213,9.99,40.0
    19954471,728387.975,1995-04-04,93.975,-67.428,44.205,9.99,40.0
    19954471,728388.0451,1995-04-05,94.0451,-67.431,44.202,9.99,40.0
    19954471,728388.2576,1995-04-05,94.2576,-67.44,44.224,9.99,40.0
    """
    D = np.genfromtxt(FN,dtype=None,names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH_I'],delimiter=',')
        
    td=np.array(D['TimeRD'])
    latd=np.array(D['LAT_DD'])
    lond=np.array(D['LON_DD'])
    depd=np.median(np.array(D['DEPTH_I']))
    
    #start and end times close to a whole hour
    t1=np.ceil(np.min(td)*24.)/24.
    t2=np.floor(np.max(td)*24.)/24.

# available range of FVCOM GOM3 30yr hindcast
# record must be long enogh for interpolation to work   
    if (t1>RataDie(1978,1,1)) and (t2<RataDie(2011,1,1)) and len(td)>3:
    
        # interpolate drifter trajectory to hourly 
        #tdh=np.arange(t1,t2,1./24.)
        tdh=np.arange(t1,t2,1./24.)
        #latdh=np.interp(tdh,td,latd)
        #londh=np.interp(tdh,td,lond)

        # optionally cubic interpolation: scipy interpolate interp1d
        #fip1=interpolate.interp1d(td,latd,kind='cubic',bounds_error=False) #returns function not array
        #latdh=fip1(tdh)
        #fip2=interpolate.interp1d(td,lond,kind='cubic',bounds_error=False) #returns function not array
        #londh=fip2(tdh)         
        # the above algorithm gets very slow for long arrays
        # and fails for many long tracks because of lack of memory
        # it must have used matrix inversion to calculate spline over all data points.
        
        # I use cubic Hermite polynomial spline, SeaHorseTide.sh_interp3.
        # It is a local algorithm only using neighbor data
        latdh=sh_interp3(tdh,td,latd)
        londh=sh_interp3(tdh,td,lond)
        
        
        
        """
        plt.figure()
        plt.plot(D['LON_DD'],D['LAT_DD'],'b.-')
        plt.plot(londh,latdh,'r.')
        plt.show()
        """
        
        ND=len(tdh)
        print 'points in track', ND 
        kvdh=np.zeros(ND,dtype=int)
        u0=np.zeros(ND,dtype=float)   
        v0=np.zeros(ND,dtype=float)   
        uws=np.zeros(ND,dtype=float)   
        vws=np.zeros(ND,dtype=float)   

        for i in range(ND):
            lonp=londh[i];latp=latdh[i]

            t=tdh[i]
                   
            tt=np.round(t*24.)/24.
            ti=datetime.fromordinal(int(tt))
            YEAR=str(ti.year)
            MO=str(ti.month).zfill(2)
            DA=str(ti.day).zfill(2)
            hr=(tt-int(tt))*24
            HR=str(int(np.round(hr))).zfill(2)            
            TS=YEAR+MO+DA+HR+'0000'

            tchk=RataDie(int(YEAR),int(MO),int(DA))+int(HR)/24.
            if (t-tchk)*24. > 0.5 :
                print 'warning tchk ', kdr,TS
            PATH0='/home/vsheremet/FATE/'            
            #PATH0='F:/'
            FNU=PATH0+'GOM3_DATA/GOM3_'+YEAR+'/u0/'+TS+'_u0.npy'
            FNV=PATH0+'GOM3_DATA/GOM3_'+YEAR+'/v0/'+TS+'_v0.npy'
            FNUWS=PATH0+'GOM3_DATA/GOM3_'+YEAR+'/uwind_stress/'+TS+'_uwind_stress.npy'
            FNVWS=PATH0+'GOM3_DATA/GOM3_'+YEAR+'/vwind_stress/'+TS+'_vwind_stress.npy'

#            fu=np.load(FNU)
#            fv=np.load(FNV)
            fu=np.load(FNU).flatten()
            fv=np.load(FNV).flatten()
            fuws=np.load(FNUWS).flatten()
            fvws=np.load(FNVWS).flatten()

               
            u0[i],v0[i],kvdh[i]=VelInterp_lonlat(lonp,latp,Grid,fu,fv)   
            uws[i],vws[i],kvdh[i]=VelInterp_lonlat(lonp,latp,Grid,fuws,fvws)   
# note that kedh is in fact kv, vertex corresponding to a polygon of the dual mesh 
       
        udh=np.zeros(ND,dtype=float)   
        vdh=np.zeros(ND,dtype=float)
        Coef=111111./86400. # deg/day -> m/s
        #print Coef
        for i in range(1,ND-1):
            udh[i]=(londh[i+1]-londh[i-1])/(tdh[i+1]-tdh[i-1])*Coef*np.cos(latdh[i]*np.pi/180.)
            vdh[i]=(latdh[i+1]-latdh[i-1])/(tdh[i+1]-tdh[i-1])*Coef
        
        if np.isnan(u0).all():
            pass # the whole track is outside gom3
        else:
            FN1=FNs[kdr]
            FN2=FN1[0:-4]
            FN2='driftfvcom_data1/'+FN2+'.npz'
            print FN2
            np.savez(FN2,tdh=tdh,kvdh=kvdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh,umoh=u0,vmoh=v0,uwsh=uws,vwsh=vws)
        
    kdr=kdr+1 # next drifter
#    np.save('kdr.npy',kdr)
    f=open('kdrhour.txt','w')
    f.write(str(kdr))
    f.close()    
