# -*- coding: utf-8 -*-
"""
Local functions for all orientation codes
A. Doran and G. Laske
"""

#############################
# Import all packages for entire program
#############################

from numpy import *
from geographiclib.geodesic import *
import subprocess
import numpy as np
import scipy.signal as sig
import obspy 
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from scipy.stats import circmean as cmean
from scipy.stats import circstd as cstd
from scipy.stats import hmean as hm
import os

#############################
# Functions for driver program
#############################


# First pass at detrending
def detr(T):
    T.detrend()
    T.detrend('linear')
    return T
    
# Organize channels by coordinate system
def org(T,nameconv):
    # Need to add to this in case other type (eg LHZ) is used
    if nameconv==1:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??N")
        bh2=T.select(channel="??E")    
    if nameconv==2:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??1")
        bh2=T.select(channel="??2")
    if nameconv==3:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??2")
        bh2=T.select(channel="??1")
    
    s=bh1+bh2+bhz

    # decimate just to make things easier
    # for now assume all channels are same sampling rate
    while s[0].stats.sampling_rate >= 10.0 and s[0].stats.sampling_rate / 10.0 == int(s[0].stats.sampling_rate / 10.0):
        s.decimate(10)

    while s[0].stats.sampling_rate / 5.0 == int(s[0].stats.sampling_rate / 5.0):
        s.decimate(5)

    return s
  
# used by catclean
def close(x1,x2):
    if abs(x1-x2)<0.8:
        return True
    else:
        return False

# look for repeat events
def catclean(stime,lat,lon,mag):
    rep=array([],dtype=int)
    for k in arange((len(stime))):
        for kk in arange((len(stime))):
            if stime[kk]-stime[k]<60*15 and close(lat[kk],lat[k]) and close(lon[kk],lon[k]) and mag[kk]<mag[k]-0.3:
                rep=append(rep,kk)
    return rep
    

#############################
# PHASE VELOCITY CALCULATIONS
#############################

def getf(freq,A):
    for i in arange((len(A))):
        if A[i,0]==freq:
            v=A[i,1]
    return v            

###
# find nearest value in an array
def nv(x,v):
    # x is array
    # v is value
    idx = (abs(x-v)).argmin()
    return x[idx]

# Overall function to get path-averaged group velocity
def pathvels(lat1,lon1,lat2,lon2,map10,map15,map20,map25,map30,map35,map40):

    Rearth=6371.25;
    circE=2*np.pi*Rearth;
    
    # Get distance and azimuth
    p=Geodesic.WGS84.Inverse(lat1,lon1,lat2,lon2)
    
    minor=float(p['s12']) / 1000.00
    major=circE-minor
            
    l1=Geodesic.WGS84.Line(lat1,lon1,p['azi1'])         
    l2=Geodesic.WGS84.Line(lat1,lon1,p['azi1']-180)     
    
    deg=111.17
    D1=zeros([1,2])
    for i in arange((361)):
        b=l1.Position(deg*1000*i)
        p2=array([b['lat2'],b['lon2']])
    
        if i==0:
            D1[0,0]=p2[0]
            D1[0,1]=p2[1]
        else:
            D1=vstack((D1,p2))
        
        bb=Geodesic.WGS84.Inverse(lat2,lon2,b['lat2'],b['lon2'])
        if bb['s12'] <= deg*1000.0 :
            break
    
    
    D2=zeros([1,2])
    for i in arange((361)):
        b=l2.Position(deg*1000*i)
        p2=array([b['lat2'],b['lon2']])
    
        if i==0:
            D2[0,0]=p2[0]
            D2[0,1]=p2[1]
        else:
            D2=vstack((D2,p2))
        
        bb=Geodesic.WGS84.Inverse(lat2,lon2,b['lat2'],b['lon2'])
        if bb['s12'] <= deg*1000.0 :
            break
    
    """
    We now have lat and lon points along the major and minor great circles.
    We calcaulte the group velocity of at each point, and then
    find the average velocities. 
    """
    
    for k in arange((len(D1))):
        if D1[k,1]<0:
            D1[k,1]+=360
    for k in arange((len(D2))):
        if D2[k,1]<0:
            D2[k,1]+=360

    def Ray(D):
        U1=zeros([len(D),7])
        for k in arange((len(D))):
            # do latitude first
            
            ## get correct precision
            ## designed to match results of Ma et al codes
            if abs(D[k,1]) < 10:
                D[k,1]=round(D[k,1],5)
            if abs(D[k,1]) >= 10 and abs(D[k,1]) < 100:
                D[k,1]=round(D[k,1],4)
            if abs(D[k,1]) > 100:
                D[k,1]=round(D[k,1],3)
            if abs(D[k,0]) < 10:
                D[k,0]=round(D[k,0],5)
            if abs(D[k,0]) >= 10:
                D[k,0]=round(D[k,0],4)
            #    
            # find right index
            q=where( map10[:,1] == nv(map10[:,1], (D[k,0])  ))[0]
            qq=where( map10[q,0]==nv(map10[q,0],  (D[k,1])  ))[0]
            idx=q[qq]
            
            # update path
            U1[k,0]= map10[ idx, 2]
            U1[k,1]= map15[ idx, 2]
            U1[k,2]= map20[ idx, 2]
            U1[k,3]= map25[ idx, 2]
            U1[k,4]= map30[ idx, 2]
            U1[k,5]= map35[ idx, 2]
            U1[k,6]= map40[ idx, 2]
        mhz=array([10,15,20,25,30,35,40])
        return array((mhz, hm(U1,axis=0))).T

    return Ray(D1),Ray(D2)

#####################
## ANGLE CALCULATION FUNCTIONS
#####################

# keep eqs above certain cc limit
# also keep which earthquakes were kept
# Different version of C1
def C1_2(phi,cc,clim):
    PHI=array([]); C=array([]); ix=array([])
    for i in arange((len(phi))):
        if cc[i]>clim:
            PHI=append(PHI,phi[i])
            C=append(C,cc[i])
            ix=append(ix,i)
    return PHI, C, ix

# C2 culling - keep values within 95% circ conf of circ mean
# also keep which earthquakes were kept

# Get unique events used in final calculation
def uniqueevents(phis,ccs,n,R1cc,R2cc):
    L=len(phis)/2
    ii=zeros((len(n)))
    for i in arange((len(n))):
        if n[i]<L:
            ii[i]=where(R1cc==ccs[int(n[i])])[0][0]
        else:
            ii[i]=where(R2cc==ccs[int(n[i])])[0][0]

    return unique(ii)


# Plotting function
def centerat(phi,m=0):
    phinew=copy(phi)
    if len(shape(phi))==1:
        for i in arange((len(phi))):
            if phi[i]>=m+180:
                phinew[i]-=360
            if phi[i]<=m-180:
                phinew[i]+=360
    else:
        for k in arange((shape(phi)[1])):
            for i in arange((shape(phi)[0])):
                if phi[i,k]>=m+180:
                    phinew[i,k]-=360
                if phi[i,k]<=m-180:
                    phinew[i,k]+=360
    return phinew


# Function to flatten result arrays
def flatten(X):
    return reshape(X,[X.shape[0]*X.shape[1],1])

# median absolute deviation
def mad(x):
    return median(abs(x-median(x)))

# remove outliars
def outlier1(Tphi,ix,lim=5.0):
    devs=abs(Tphi-median(Tphi))/mad(Tphi)
    ixs=where(devs<lim)
    return Tphi[ixs],ix[ixs]

# bootstrap mean
def boot1(phi,bootnum):
    m=zeros((bootnum)); L=len(phi)
    for i in arange((bootnum)):
        a=np.random.choice(phi,size=L,replace=True)
        m[i]=cmean(a,high=360)
    return m

# reorganize results
def resort(phi,col2):
    phi2=centerat(phi,m=cmean(phi,high=360))
    t=zeros((len(phi2),2))
    t[:,0]=phi2; t[:,1]=col2
    t = t[t[:,0].argsort()]
    return t[:,0], t[:,1]

# final Doran-Laske calculation
def fcalc1(phi,cc,lim,R1cc,R2cc):
    # keep cc over limit
    Tphi,Tcc,ii=C1_2(phi,cc,lim)
    if len(Tphi)==0:
        return 0,180,array([0]),0
    if mad(Tphi)==0:
        return mean(Tphi),90,array([1]),len(Tphi)
    # remove outliers using MAD
    Tphi,ii=resort(Tphi,ii)
    Tphi2,ii2=outlier1(Tphi,ii)
    # bootstrap results for statistic
    m=boot1(Tphi2,5000)
    
    
    return cmean(m,high=360),2*1.96*cstd(m,high=360),uniqueevents(phi,cc,ii2,R1cc,R2cc),len(Tphi2) 

# create histogram of results
def fhist1(phi,cc,lim):
    Tphi,Tcc,ii=C1_2(phi,cc,lim)
    Tphi,ii=resort(Tphi,ii)
    Tphi2,ii2=outlier1(Tphi,ii)
#    plt.figure()
    plt.hist(Tphi2,bins=25)
    plt.xlabel('Orientation Estimate')
    plt.ylabel('Counts')
    return 

    
    
    
#############################
# A few other random necessary ones
#############################

# Define rotation function
#   -rotates horizontal components CW from N by alpha (in degrees)
def rot2d(N,E,alpha):
    a=deg2rad(alpha)  # convert from degrees to radians
    r=cos(a)*N - sin(a)*E
    t=sin(a)*N + cos(a)*E
    return(r,t)

# root mean square
def rms(x):
    return sqrt(mean(abs(x)**2))

def find_nearest(array,value):
    return (abs(array-value)).argmin()

def checklen(st,hrs):
    # checks to see if there is enough downloaded data to run program
    L=len(st)
    for i in arange((L)):
        if (UTCDateTime(st[i].stats.endtime)-UTCDateTime(st[i].stats.starttime))+100 < hrs:
            return True        
    if var(st[0].data)<1 or var(st[1].data)<1 or var(st[2].data)<1:
        return True
    return False

# save resutls
def saved(R1cc,R2cc,R1phi,R2phi,loc='temp.dir'):
    if not os.path.exists(loc):
        os.mkdir(loc)
    savetxt(loc+'/'+'R1cc',R1cc)
    savetxt(loc+'/'+'R2cc',R2cc)
    savetxt(loc+'/'+'R1phi',R1phi)
    savetxt(loc+'/'+'R2phi',R2phi)
    return

"""
Functions from file formerly called Phases
Mostly deal with computing correlations
and rotating data
"""



# Resize arrays to all identical shapes
def resiz(x1,x2,x3):
    a1=len(x1); a2=len(x2); a3=len(x3)
    L=min(array([a1,a2,a3]))
    return x1[0:L], x2[0:L], x3[0:L]
    
# preprocess segments of data
# taper, zerophase filter, detrend
def sw1proc(T,LPF,HPF,corn=4):
    T.taper(type='hann',max_percentage=0.05)
    T.filter("lowpass",freq=LPF,corners=corn,zerophase=True)
    T.filter("highpass",freq=HPF,corners=corn,zerophase=True)
    T.detrend()
    
    return T

# DORAN-LASKE calculation for one freq, one orbit of surface wave
def SW1(TT,Rf,LPF,HPF,daz1,A,nameconv,winlen=10.0,ptype=0):
    # event info
    daz=daz1[0]/1000.0 # convert to KM
    baz=daz1[1] # angle from station to event        

    Rvel=getf(Rf,A) # Group velocity at Rf
    R1window=(1.0/(Rf/1000.0))*winlen

    # Process
    T=sw1proc(TT.copy(),LPF,HPF)

    # Window info
    arv=1.0/Rvel * daz
    r1=arv-R1window/2.0
    r2=arv+R1window/2.0

    dt=T[0].stats.starttime
    P=T.slice(starttime=dt+r1,endtime=dt+r2)
    
    rdat=P[0].data
    rdat2=P[1].data
    vdat=imag(sig.hilbert(P[2].data))
    
    # Ensure all data vectors are same length
    rdat,rdat2,vdat=resiz(rdat,rdat2,vdat)
    
    # rotate through and find max cc
    degs=360*4
    ang=arange((degs))
    cc=zeros((degs)); cc2=zeros((degs)); cc3=zeros((degs))
    for k in ang:
        n,e=rot2d(rdat,rdat2,k/4.0)
        covmat=corrcoef(n,vdat)
        cc[k]=covmat[0,1]
        cstar=cov(vdat,n)/cov(vdat)
        cc2[k]=cstar[0,1]
        cstar=cov(vdat,e)/cov(vdat)
        cc3[k]=cstar[0,1]

    # Keep angle determined by cstar, but use rating from corrcoef
    #   Formulas in Stachnik paper
    ANG=cc2.argmax(); #CC[j]=cc[ANG]
    # correct for angles above 360
    or_ang= (baz- (360 - ANG/4.0) )

    # ADJUST FOR NAMING CONVENTION
    if nameconv==3: or_ang+=180

    if or_ang<0: or_ang+=360
    if or_ang>=360: or_ang-=360
    # Can plot xc

    # plotting:
    # ptype=0, no plot
    # ptype=1, Rayleigh plot
    # ptype=2, Love plot

    if ptype==1:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        ZZ=dat.epoch2num(T)
        Z=dat.num2date(ZZ)
        n,e=rot2d(rdat,rdat2,ANG/4.0)
#        plt.figure()
        plt.plot(Z,vdat,label='Vetical')
#        savetxt('/Users/adoran/Desktop/T.txt',T)
#        savetxt('/Users/adoran/Desktop/vdat.txt',vdat)
#        savetxt('/Users/adoran/Desktop/n.txt',n)
        plt.hold("on")
        plt.plot(Z,n,label='BH1')    
        plt.legend(loc=4)
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.title('D-L Results (%i mHz)' %(Rf))
    elif ptype==2:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        ZZ=dat.epoch2num(T)
        Z=dat.num2date(ZZ)
        n,e=rot2d(rdat,rdat2,ANG/4.0)
        plt.figure()
        plt.subplot(121)
        plt.plot(Z,vdat,label='Vetical')
        plt.hold("on")
        plt.plot(Z,n,label='BH1')    
        plt.legend(loc=4)
        plt.xlabel('Time')
        plt.suptitle('D-L Results (%i mHz)' %(Rf))
        plt.subplot(122)
        plt.plot(Z,e,label='BH2')
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.legend(loc=4)
    elif ptype==3:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        n,e=rot2d(rdat,rdat2,ANG/4.0)
        plt.figure()
        plt.plot(T,vdat,label='Vetical')
        
#
    return or_ang, cc[ANG]





