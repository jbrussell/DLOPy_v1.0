# -*- coding: utf-8 -*-
"""
PRIMARY ORIENTATION PROGRAM
ADRIAN. K. DORAN
GABI LASKE
VERSION 1.0
RELEASED APRIL 2017
"""

#########################
# Import necessary packages, functions, and parameter files
from Orient_PF import *
from locfuns import *
#########################

# Define variables using param file inputs1
client1=Client(cat_client)      # Catalog Client
client2=Client(wf_client)     # Waveform Client
t1=UTCDateTime(time1)           # start date
t2=UTCDateTime(time2)           # end date
cat=client1.get_events(starttime=t1,endtime=t2,minmagnitude=minmag)#,maxmagnitude=maxmag)

# get station info
if localdata==1 or localcoords==1:
    sta_lat=inputlat
    sta_lon=inputlon
else:
    # Get station and event data
    inv=client1.get_stations(network=NET,station=STA)
    sta_lat=inv[0][0].latitude
    sta_lon=inv[0][0].longitude
    
# Organize station and event data
# organize data from catalog
L=len(cat.events)
lat=zeros(L); lon=zeros(L); mag=zeros(L); stime=zeros(L); dep=zeros(L); deg=zeros(L); baz=zeros(L)
for i in arange(L):
    lat[i]=cat.events[i].origins[0].latitude                            # latitude
    lon[i]=cat.events[i].origins[0].longitude                            # longitude
    dep[i]=cat.events[i].origins[0].depth                            # depth
    stime[i]=UTCDateTime(cat.events[i].origins[0].time)                # event start time
    mag[i]=cat.events[i].magnitudes[0].mag                            # magnitude
    #daz1=obspy.core.util.gps2DistAzimuth(sta_lat,sta_lon,lat[i],lon[i])   # depricated version
    daz1=obspy.geodetics.gps2dist_azimuth(sta_lat,sta_lon,lat[i],lon[i])   # distance b/t station and event
    #deg[i]=obspy.core.util.kilometer2degrees(daz1[0]/1000)                # depricated version
    deg[i]=obspy.geodetics.kilometer2degrees(daz1[0]/1000)               
    baz[i]=daz1[1]                                                     # angle from station to event
# get index of repeat events, save for later
reps=unique(catclean(stime,lat,lon,mag))

# to save catalog:
if savecat==1:
    ts=array([],dtype=object)
    h1=array(['Time','Lon','Lat','Dep(km)','Mag'],dtype=object)
    for i in arange(L):
        ts=append(ts,UTCDateTime(stime[i]))
    catprint2=array((ts,lon,lat,dep/1000,mag),dtype=object).T
    catprint=vstack((h1,catprint2))
    savetxt(catname,catprint,fmt="%s")
    
#   INITIALIZE INITIALIZE
# Initialize surface wave arrays
numsurfcalcs=7
R1phi=zeros([L,numsurfcalcs]); R1cc=zeros([L,numsurfcalcs])
R2phi=zeros([L,numsurfcalcs]); R2cc=zeros([L,numsurfcalcs]);
# Initialize Stachnik arrays
R4phi=zeros((L)); R4cc=zeros((L))

hrs=4*60*60     # Length of data to download

# load group velocity maps
map10=loadtxt('R.gv.10.txt'); map15=loadtxt('R.gv.15.txt')
map20=loadtxt('R.gv.20.txt'); map25=loadtxt('R.gv.25.txt')
map30=loadtxt('R.gv.30.txt'); map35=loadtxt('R.gv.35.txt')
map40=loadtxt('R.gv.40.txt')

#      LOOP OVER ALL EVENTS
for j in arange((L)):

#     GET WAVEFORMS, including protections
    try:
        if localdata==0:
            # download data from client
            s=client2.get_waveforms(NET,STA,LOC,CHA,UTCDateTime(stime[j]),UTCDateTime(stime[j]+hrs))
        else:
            # access local data
            s=LF(UTCDateTime(stime[j]))
            
        # merge waveforms (sometimes downloaded in several segments)
        s.merge()
        
        # don't want any masked data or data with nans
        for q in arange((len(s))):
            if ma.count_masked(s[q].data)>0:
                continue
            
        # remove mean and trend
        s.detrend()
        s.detrend('linear')
        if len(s)<3:
            continue
        # [0] bh1, [1] bh2, [2] bhz
        st=org(s.copy(),nameconv) 
            # organizes data by coordinate system
            # also downsamples to <10 Hz
    
        # check data length, data quality        
        if checklen(st,hrs):
            continue

    except:
        continue
        
    # get some additional parameters
    #daz1=obspy.core.util.gps2DistAzimuth(sta_lat,sta_lon,lat[j],lon[j])
    daz1=obspy.geodetics.gps2dist_azimuth(sta_lat,sta_lon,lat[j],lon[j])
    daz2=copy(daz1)
    Rearth=6371.25*1000; circE=2*np.pi*Rearth;
    daz2[0]=circE-daz2[0]; daz2[1]=daz2[1]+180  # major & minor arc calculation
    if daz2[1]>=360: daz2[1]-=360
    
    # SURFACE WAVE CALCULATIONS
    
    # conditions
    # minimum distance, maximum distance, and maximum depth
    if deg[j]<mindeg_sw or deg[j]>maxdeg_sw or dep[j]>=maxdep_sw*1000:
        continue
    # clean catalog of repeats (repeat conditions set in catclean)
    if j in reps:
        continue
    
    # get path-averaged group velocities
    Ray1,Ray2=pathvels(sta_lat,sta_lon,lat[j],lon[j],map10,map15,map20,map25,map30,map35,map40)  
#        
    # FOR EACH FREQUENCY AND ORBIT, calculate arrival angle

##    # freq 1 (40 mHz)
    Rf=40.0; HPF=0.035; LPF=0.045
    
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=20.0,ptype=0)
    R1phi[j,0]=ANG; R1cc[j,0]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=24.0,ptype=0)
    R2phi[j,0]=ANG; R2cc[j,0]=cc

##    # freq 2 (35 mHz)
    Rf=35.0; HPF=0.030; LPF=0.040
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=17.0,ptype=0)
    R1phi[j,1]=ANG; R1cc[j,1]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=20.0,ptype=0)
    R2phi[j,1]=ANG; R2cc[j,1]=cc

#
###    # freq 3 (30 mHz)
    Rf=30.0; HPF=0.025; LPF=0.035
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=14.0,ptype=0)
    R1phi[j,2]=ANG; R1cc[j,2]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=16.0,ptype=0)
    R2phi[j,2]=ANG; R2cc[j,2]=cc


# # # freq 4 (25 mHz)
    Rf=25.0; HPF=0.020; LPF=0.030
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=12.0,ptype=0)
    R1phi[j,3]=ANG; R1cc[j,3]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=13.0,ptype=0)
    R2phi[j,3]=ANG; R2cc[j,3]=cc

###    # freq 5 (20 mHz)
    Rf=20.0; HPF=0.015; LPF=0.025
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=10.0,ptype=0)
    R1phi[j,4]=ANG; R1cc[j,4]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=10.0,ptype=0)
    R2phi[j,4]=ANG; R2cc[j,4]=cc

###    # freq 6 (15 mHz)
    Rf=15.0; HPF=0.020; LPF=0.010
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=10.0,ptype=0)
    R1phi[j,5]=ANG; R1cc[j,5]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=10.0,ptype=0)
    R2phi[j,5]=ANG; R2cc[j,5]=cc

###    # freq 7 (10 mHz)
    Rf=10.0; HPF=0.005; LPF=0.015
#
    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=7.0,ptype=0)
    R1phi[j,6]=ANG; R1cc[j,6]=cc

    ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=7.0,ptype=0)
    R2phi[j,6]=ANG; R2cc[j,6]=cc

# save up that Data
    if constsave==1:
        saved(R1cc,R2cc,R1phi,R2phi,loc=str(saveloc))          
        
    # WHAT TO OUTPUT AT THE END OF EACH ITERATION
    if verb==1:
        # Just output number
        print("%s: %i / %i" %(STA,j+1,L))
    elif verb==2:
        print("%s: %i / %i" %(STA,j+1,L))
        print("R1-30 cc: %.2f   R1-30 phi: %.2f" %(R1cc[j,2],R1phi[j,2]))
    

# PLOT ALL RESULTS
if finplot==1:
    plt.figure()
    plt.subplot(1,1,1)
    plt.title('Surf Waves')
    plt.plot(R1cc,R1phi,'x',R2cc,R2phi,'x')
    plt.ylim([0,360]); plt.xlim([0,1])

# SAVE DATA
saved(R1cc,R2cc,R1phi,R2phi,loc=str(saveloc))


