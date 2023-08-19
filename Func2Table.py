import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

savefile = False # Will save output files if True
trajName = 'test'
filename = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+trajName+'.txt'; # Output file name
currentwdr = os.getcwd()
if os.path.exists('./'+datetime.now().strftime('%Y-%m-%d'))==False:
    os.mkdir(datetime.now().strftime('%Y-%m-%d'));

numpoints = ; #Number of time steps
tinc = 1; #Increment time length, in integer multiples of 100 us
t_tot = numpoints * tinc * 0.0001 #Total frequency sweep time in units of seconds 
tvec = np.arange(numpoints)
defaultfreq0='080.0000000'; #frequency when table is not running
defaultphase0='00000'; #phase when table is not running
defaultamp0='0696'; #amplitude when table is not running (N/1023) default 0696
defaultfreq1='080.0000000'; #frequency when table is not running
defaultphase1='00000'; #phase when table is not running
defaultamp1='0820'; #amplitude when table is not running (N/1023) default 0820

def freqfunc0(t):
    f=[80.0 for i in t]
    #f=["{0:.7f}".format(x) for x in func] #Put freq function in format(__) in MHz
    return f

def phasefunc0(t):
    ph=np.zeros(len(t)).astype(int).astype('str')
    return ph

def ampfunc0(t):
    a=['0696' for i in np.arange(len(t))] #default 0696
    return a

a=0.3;
stepamp=a*1.05391;
t1=a*1250;
#t2=(numpoints-4*t1)/2;
t2=3784;

def freqfunc1(t):
    f=[stepamp*(1-np.cos(2*np.pi/(2*t1)*i))*(np.heaviside(i,1)-np.heaviside(i-t1,1))
    +2*stepamp*(np.heaviside(i-t1,1)-np.heaviside(i-(t1+t2),1))
    +stepamp*(1+np.cos(2*np.pi/(2*t1)*(i-(t1+t2))))*(np.heaviside(i-(t1+t2),1)-np.heaviside(i-(2*t1+t2),1))
    +80.0 for i in t]
    
    #f=[-0.3*(1-np.cos(2*np.pi/(numpoints)*i))+0.3*(1-np.cos(2*np.pi/(numpoints)*i))*np.heaviside(i-numpoints/2,1)-2*0.3*np.heaviside(i-numpoints/2,1)+80.0 for i in t]
    
    #f=[-0.805585*(1-np.cos(2*np.pi/(numpoints/2)*i))*(2*np.heaviside(i-numpoints/2,1)-1)+80.0 for i in t]
    
    #f=[-0.064092*np.sin(np.sin(np.pi/(numpoints)*i)**2)**2+80.0 for i in t]
    
    #f=[0.9581*(1-np.cos(2*np.pi/(numpoints)*i))+80.0 for i in t]
    
    #f=[-0.037594*np.sin(np.pi/(numpoints)*i)**2+80.0 for i in t]
    
    #f=[-0.029526*np.sin(np.pi/(numpoints)*i)+80.0 for i in t]
    
    #f=[80.0-0.15 for i in t]
    
    #f=["{0:.7f}".format(x) for x in func] #Put freq function in format(__) in MHz
    
    return f

def phasefunc1(t):
    ph=np.zeros(len(t)).astype(int).astype('str')
    return ph

def ampfunc1(t):
    a=['0820' for i in np.arange(len(t))] #default 0820
    return a

fnum0=freqfunc0(tvec); #freq float
freq0=["{0:.7f}".format(x) for x in fnum0]; #freq string
phase0=phasefunc0(tvec);
amp0=ampfunc0(tvec);

fnum1=freqfunc1(tvec); #freq float
freq1=["{0:.7f}".format(x) for x in fnum1]; #freq string
phase1=phasefunc1(tvec);
amp1=ampfunc1(tvec);

#format strings to get fixed number of characters per line
for i in np.arange(numpoints):
    while len(freq0[i][:-8])<3:
        freq0[i]='0'+freq0[i];
    while len(freq1[i][:-8])<3:
        freq1[i]='0'+freq1[i];
    while len(phase0[i])<5:
        phase0[i]=phase0[i]+'0';
    while len(phase1[i])<5:
        phase1[i]=phase1[i]+'0';


#create properly formatted table of strings of all values        
table=['']*((numpoints+1)*2);
table[0]='0 '+defaultfreq0+' '+defaultphase0+' '+defaultamp0+' '+'0255';
table[1]='1 '+defaultfreq1+' '+defaultphase1+' '+defaultamp1+' '+'0255';
for i in (np.arange(numpoints)+1):
    if i<numpoints:
        table[i*2]= '0 '+freq0[i-1]+' '+phase0[i-1]+' '+amp0[i-1]+' '+format(tinc,"04")
        table[2*i+1]= '1 '+freq1[i-1]+' '+phase1[i-1]+' '+amp1[i-1]+' '+format(tinc,"04")
    else:
        table[i*2]= '0 '+freq0[i-1]+' '+phase0[i-1]+' '+amp0[i-1]+' '+'0000'
        table[2*i+1]= '1 '+freq1[i-1]+' '+phase1[i-1]+' '+amp1[i-1]+' '+'0000'

#save to text file
if savefile==True:        
    txtfile=open(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+filename,'w+')
    for i in np.arange(len(table)):
        txtfile.write(table[i]+'\n')
    txtfile.close()

#compute distance travelled
lattice_spacing=0.000000532;
diff=[0]*len(fnum0);
for i in np.arange(len(fnum1)):
    diff[i]=2*1000000*(fnum1[i]-fnum0[i]) #factor of 2 because of quadruple pass AOM configuration
dist=np.trapz(diff,tvec*tinc*0.1/1000)*lattice_spacing;
print('Distance travelled (m) = '+str(dist))

#compute acceleration, max acceleration
accel=np.gradient(diff,tinc*0.0001)*lattice_spacing;
max_accel=np.max(accel)
print('Max acceleration (m/s^2) = '+str(max_accel))

#plot freq
fig,ax0=plt.subplots()
fig.suptitle(filename,fontweight='bold')
fig.subplots_adjust(top=0.8)
ax0.set_title('Distance (m) = '+str(dist)+'\n'+'Max accel (m/s^2) = '+str(max_accel))
ax0.plot(tvec*tinc*0.1,fnum0,marker='.');
ax0.plot(tvec*tinc*0.1,fnum1,marker='.');

ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Frequency (MHz)')
ax0.legend(['freq0','freq1'])
#plt.ylim(np.min(fnum1),80)

if savefile==True:
    fig.savefig(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+filename+'.png')
