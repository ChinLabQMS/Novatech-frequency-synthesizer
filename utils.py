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

# Wrapper for round trip trajectory
def round_trip(func):
    def inner(wait_time, t_inc, *args, **kwargs):
        t, f = func(t_inc, *args, **kwargs)
        f = np.concatenate((f, np.zeros(int(wait_time / (t_inc*1e-4))), -f[::-1]))
        t = np.arange(len(f), dtype=int)
        return t, f
    return inner

def freqdiff_flattop(t_inc, t_tot, t_ramp, dist):
    if type(t_ramp) == float:
        T1 = int(t_ramp * t_tot / (1e-4*t_inc))
    elif type(t_ramp) == int:
        T1 = t_ramp

    t = np.arange(t_tot / (1e-4*t_inc), dtype=int)
    f = (1 - np.cos(np.pi * t / T1)) * (t < T1) \
        + 2 * ((t >= T1) & (t < len(t) - T1)) \
        + (1 + np.cos(np.pi * (t - (len(t) - T1)) / T1)) * (t >= len(t) - T1)
    
    amp = dist / (np.sum(abs(f)) * wavelength * t_inc * 1e2)
    return t, np.around(amp * f, decimals=7)

def make_figure(freqdiff, t_inc, trajName):
    dist_actual = np.cumsum(freqdiff) * wavelength * t_inc *1e2
    accel = np.diff(freqdiff) * 1e6 / (t_inc*0.0001) * wavelength
    vel = freqdiff * 1e6 * wavelength
    description = f'Max distance (m) = {dist_actual.max():.4f}\nMax acceleration (m/s^2) = {accel.max():.2f}\nMax velocity (m/s) = {vel.max():.3f}'

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(trajName, fontweight='bold')
    fig.subplots_adjust(top=0.85)
    ax.set_title(description)
    ax.plot(np.arange(len(freqdiff)) * t_inc * 0.1, freqdiff)
    ax.axhline(0, c='k', linewidth=0.5)

    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Frequency diff (MHz)')
    ax.legend(['freq1 - freq0'])
    plt.show()

    return fig

def save_table(freqdiff, t_inc, trajName, fig):
    freq1 = ["0%#10.7f" % (x + float(defaultfreq1)) for x in (freqdiff)]

    if os.path.exists('./' + datetime.now().strftime('%Y-%m-%d')) == False:
        os.mkdir(datetime.now().strftime('%Y-%m-%d'))
    
    currentwdr = os.getcwd()
    fileName = datetime.now().strftime('%Y-%m-%d_%H-%M-%S_') + trajName + '.txt' # Output file name

    counter = 1
    prevfreq1 = defaultfreq1
    with open(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+fileName,'w+') as textFile:
        textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + ' 0255\n')
        textFile.write('1 ' + defaultfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + ' 0255\n')

        for i in range(len(freq1)):
            if freq1[i] != prevfreq1:
                prevfreq1, counter = freq1[i], 1
            else:
                counter += 1
                continue
            textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + ' ' + format(t_inc * counter,"04") + '\n')
            textFile.write('1 ' + prevfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + ' ' + format(t_inc * counter,"04") + '\n')
        textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + format(t_inc * counter,"04") + '\n')
        textFile.write('1 ' + prevfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + format(t_inc * counter,"04") + '\n')
        textFile.write('0 ' + defaultfreq0 + ' ' + defaultphase0 + ' ' + defaultamp0 + ' 0000\n')
        textFile.write('1 ' + defaultfreq1 + ' ' + defaultphase1 + ' ' + defaultamp1 + ' 0000')

    print('Table saved as ' + fileName)

    fig.savefig(currentwdr+'\\'+datetime.now().strftime('%Y-%m-%d')+'\\'+fileName+'.png')
    print('Figure saved as ' + fileName + '.png')