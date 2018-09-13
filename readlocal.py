# -*- coding: utf-8 -*-
"""
File with functions to read local data
A. Doran and G. Laske
"""

# import necessary functions
import obspy
from obspy import read as r
from obspy.core import UTCDateTime
from numpy import *
from scipy import signal


def getmsd1(t1):
    dir1='/Volumes/AKD1/AKD-data/'
    t1=UTCDateTime(t1)

    d=r(dir1+'AKD.CH?.%i.%03d.*.msd'%(t1.year,t1.julday))

	# include next day of data if close to edge
    if t1.hour>=20:
        d2=r(dir1+'AKD.CH?.%i.%03d.*.msd'%(t1.year,t1.julday+1))
        d3=d+a2
        d=d3.merge()

	# rename channels
    C1=d.select(channel='CH1')
    C1[0].stats.channel='BH1'
    
    C2=d.select(channel='CH2')
    C2[0].stats.channel='BH2'
    
    C3=d.select(channel='CH3')
    C3[0].stats.channel='BHZ'
    
    st1=C1+C2+C3

    # set sampling rate
    sps=20
    # define four hours of data
    t2=60*60*4

    return st1.slice(starttime=t1,endtime=t1+t2)

    


