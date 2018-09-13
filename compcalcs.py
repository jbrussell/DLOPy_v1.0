# -*- coding: utf-8 -*-
"""
Final Orientation Calculation File
A. Doran and G. Laske
"""

#####################
## ANGLE CALCULATION PARAMETERS
#####################

# Get necessary functions
from locfuns import *

# location of result files
loc1='OSN1B/'

LIM=0.8         # CC limit for Surface wave calculations

#
### Specify phases to use
R1use=1
R1_40=1; R1_35=1; R1_30=1; R1_25=1; R1_20=1; R1_15=1; R1_10=1
R2use=1
R2_40=1; R2_35=1; R2_30=1; R2_25=1; R2_20=1; R2_15=1; R2_10=1
#
#

## Load files
R1phi=loadtxt(loc1+'R1phi')
R1cc=loadtxt(loc1+'R1cc')

R2phi=loadtxt(loc1+'R2phi')
R2cc=loadtxt(loc1+'R2cc')

######################
### FINAL ANGLE CALCULATIONS
######################
#
# Initialize arrays
L=len(R1phi)
phis=array([])
ccs=array([])
finval=array([]); finerr=array([])
N=array([]); 

N=full((L,L),-1.0)
LN=zeros((L))
phases=array([]);

startL=0
endL=0

A=array([L])
if endL!=0:
    A=array([endL])
    

# If not all calculations are desired, adjust accordingly
sha=shape(R1phi)
if R1use==0: R1cc=zeros(sha)
if R1use==1 and R1_40==0: R1cc[:,0]=zeros((sha[0]))
if R1use==1 and R1_35==0: R1cc[:,1]=zeros((sha[0]))
if R1use==1 and R1_30==0: R1cc[:,2]=zeros((sha[0]))
if R1use==1 and R1_25==0: R1cc[:,3]=zeros((sha[0]))
if R1use==1 and R1_20==0: R1cc[:,4]=zeros((sha[0]))
if R1use==1 and R1_15==0: R1cc[:,5]=zeros((sha[0]))
if R1use==1 and R1_10==0: R1cc[:,6]=zeros((sha[0]))

if R2use==0: R2cc=zeros(sha)    
if R2use==1 and R2_40==0: R2cc[:,0]=zeros((sha[0]))
if R2use==1 and R2_35==0: R2cc[:,1]=zeros((sha[0]))
if R2use==1 and R2_30==0: R2cc[:,2]=zeros((sha[0]))
if R2use==1 and R2_25==0: R2cc[:,3]=zeros((sha[0]))
if R2use==1 and R2_20==0: R2cc[:,4]=zeros((sha[0]))
if R2use==1 and R2_15==0: R2cc[:,5]=zeros((sha[0]))
if R2use==1 and R2_10==0: R2cc[:,6]=zeros((sha[0]))

for i in A:
    # create one massive list with necessary angles and cc values
    phis=concatenate((flatten(R1phi[startL:i,:]),flatten(R2phi[startL:i,:])))
    ccs=concatenate((flatten(R1cc[startL:i,:]),flatten(R2cc[startL:i,:])))
    
    # Doran-Laske calculation
    val,err,n,ph=fcalc1(phis,ccs,LIM,R1cc,R2cc)
    finval=append(finval,val)
    finerr=append(finerr,err)
    phases=append(phases,ph)
    for k in arange((len(n))):
        N[k,i-1]=n[k]
    LN[i-1]=len(n)
    
# output results to termianl
print "Station %s" %(loc1)
print "D-L mean, error, data included, unique events: %.2f, %.2f, %i, %i" %(finval[-1],finerr[-1],phases[-1],max(LN))
print "D-L CC level: %f" %(LIM)

###
###
#####
###
###

# create figure

CEN=finval[-1]
YLIM1=[-10+CEN,10+CEN]

plt.figure()
plt.subplot(1,1,1)
plt.title('DLOPy results',fontsize=16)
plt.plot(R1cc,centerat(R1phi,m=CEN),'x',R2cc,centerat(R2phi,m=CEN),'x')
plt.ylabel('BH1 Orientation \n Angle ($^\circ$)',fontsize=16)
plt.ylim([CEN-180,CEN+180]); plt.xlim([0,1])
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)

# save figure
plt.savefig(loc1+'cluster.eps',fmt='eps')

    

