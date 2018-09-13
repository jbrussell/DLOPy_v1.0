# -*- coding: utf-8 -*-
"""
Parameter file for Orient
A. Doran and G. Laske
"""
#################
# REQUIRED INPUTS
#################
from numpy import *

# STATION INFO
# ------------
# NETWORK CODE
NET="XM"      
# STATION NAME  
STA="OSN1B"    
# CHANNELS TO ORIENT 
CHA="BH?"     
# CHANNEL LOCATION  
LOC="*"         

####
##    IF USING LOCAL DATA
####

localdata=0
if localdata==1:
    from readlocal import *
    LF=getmsd1
####
##    
####

# If want to input own station lat and lon 
localcoords=0
if localcoords==1:
    inputlat=32.533617
    inputlon=-120.49978    

####
##
####

# DATE INFO
# Time frame in which to calculation orientaiton
#       Must be of form YYYY-MM-DD HH:MM:SS

time1="1998-05-01 00:00:00"     # Start date
time2="1998-06-01 00:00:00"     # End date

# CLIENT INFO
#       Where to download catalog and waveform data
cat_client="IRIS"       # catalog data
wf_client="IRIS"        # waveform data
# Other options according to OBSPY man page:
# ‘BGR’, ‘EMSC’, ‘ETH’, ‘GEONET’, ‘GFZ’, ‘INGV’, ‘IPGP’, 
# ‘IRIS’, ‘ISC’, ‘KOERI’, ‘LMU’, ‘NCEDC’, ‘NIEP’, ‘NOA’, 
# ‘ODC’, ‘ORFEUS’, ‘RESIF’, ‘SCEDC’, ‘USGS’, ‘USP’

# COORDINATE SYSTEM
nameconv=2
# Options for channel naming parameter
# 1 - HZ, HN, HE
# 2 - HZ + Left-handed system: H2 is 90 degs CW of H1 
# 3 - HZ + Right-handed system: H2 is 90 degs CCW of H1

# EQ SPECIFICATIONS
minmag=5.5          # Minimum magniutde EQ used
mindeg_sw=5.0      # Minimum event degree distance for surface waves
maxdeg_sw=175.0    # Maximum event degree distance for surface waves
maxdep_sw=150.0     # Maximum event depth (km) for surface waves
    
# How much information to output while the program is running
verb=2
# 0 - outputs nothing
# 1 - outputs current event being analyzed out of total (eg 4 / 79)
# 2 - outputs event number and one result (R1-30)

# Save results along the way?
constsave=1
# 0 - no
# 1 - yes
saveloc=copy(STA)
# Different save location? must be string
# saveloc=

# plot results at end?
finplot=0
# 0 - no
# 1 - yes

# save event catalog?
savecat=1
# 0 - no
# 1 - yes
catname=str(STA+'.cat.txt')


