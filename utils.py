import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

wavelength = 1064e-9 # Wavelength of the laser in units of meters
defaultfreq0 = '080.0000000' # frequency when table is not running
defaultphase0 = '00000' # phase when table is not running
defaultamp0 = '0696' # amplitude when table is not running (N/1023) default 0696
defaultfreq1 = '080.0000000' # frequency when table is not running
defaultphase1 = '00000' # phase when table is not running
defaultamp1 = '0820' # amplitude when table is not running (N/1023) default 0820

def distance_traveled(freqdiff, wavelength=wavelength, interval=1e-4):
    # freqdiff is in units of MHz
    # wavelength is in units of meters
    # distance is in units of meters
    return np.sum(abs(freqdiff)) * wavelength * 1e6 * interval

def freqdiff_flattop(t_inc, t_tot, dist, round_trip, wait_time, t1=375):
    if type(t1) == float:
        T1 = int(t1 * t_tot / (1e-4*t_inc))
    elif type(t1) == int:
        T1 = t1

    t = np.arange(t_tot / (1e-4*t_inc), dtype=int)
    f = (1 - np.cos(np.pi * t / T1)) * (t < T1) \
        + 2 * ((t >= T1) & (t < len(t) - T1)) \
        + (1 + np.cos(np.pi * (t - (len(t) - T1)) / T1)) * (t >= len(t) - T1)
    
    amp = dist / distance_traveled(f, interval=t_inc*1e-4)
    if round_trip:
        f = np.concatenate((f, np.zeros(int(wait_time / (t_inc*1e-4))), -f[::-1]))
        t = np.arange(len(f), dtype=int)
    return t, np.around(f * amp, decimals=7)

def save_table(freq0, freq1, amp0, amp1, phase0, phase1, t_inc, trajName, fig):
    if os.path.exists('./' + datetime.now().strftime('%Y-%m-%d')) == False:
        os.mkdir(datetime.now().strftime('%Y-%m-%d'))
    
    currentwdr = os.getcwd()
    fileName = datetime.now().strftime('%Y-%m-%d_%H-%M-%S_') + trajName + '.txt' # Output file name

    with open(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+fileName,'w+') as textFile:
        textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + ' 0255\n')
        textFile.write('1 ' + defaultfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + ' 0255\n')

        for i in range(len(freq0)):
            textFile.write('0 ' + freq0[i] + ' ' + phase0[i] + ' ' + amp0[i] + ' ' + format(t_inc,"04") + '\n')
            textFile.write('1 ' + freq1[i] + ' ' + phase1[i] + ' ' + amp1[i] + ' ' + format(t_inc,"04") + '\n')
        
        textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + ' 0000\n')
        textFile.write('1 ' + defaultfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + ' 0000')

    print('Table saved as ' + fileName)

    fig.savefig(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+fileName+'.png')
    print('Figure saved as ' + fileName + '.png')