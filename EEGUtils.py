import struct
import numpy as np
from matplotlib import pyplot as plt
from bisect import bisect_left
import os
import math
from scipy.signal import freqz
from scipy.signal import butter, lfilter
from scipy.signal import decimate
#import peakutils


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
   
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
    
def takeClosest(myList, myNumber):
  pos = bisect_left(myList, myNumber)
  if pos == 0:
     return myList[0]
  if pos == len(myList):
     return myList[-1]
     
  before = myList[pos - 1]
  after = myList[pos]
  rpos = pos

  if after - myNumber < myNumber - before:
     return after
  else:
     return before


def readEEG(filename):
           
	with open(filename, 'rb') as f:  
		eegdata = f.read()[16384:]
	
	nrecs=int(len(eegdata)/1044)
	eegstarts = [20+i*1044 for i in range(nrecs)]	
	tsstarts = [i*1044 for i in range(nrecs)]
	
	eeg = struct.unpack(str((512*len(eegdata)/1044))+'h', "".join([eegdata[start:start+1024] for start in eegstarts]))
	ts = struct.unpack(str((len(eegdata)/1044))+'Q', "".join([eegdata[start:start+8] for start in tsstarts]))

	return eeg, ts

def returnStability(eegdata):
    decimated = decimate(eegdata, 1000, ftype='fir')
    return decimated

def returnThetaIndex(eegdata, starts, stops):
	#fs = 512/.257552
	fs = 2003
        time_step = 1./fs

	theta_index = []
	gamma_index = []

        pstot = np.zeros([fs,len(starts)])

        j = 0

	for start,stop in zip(starts,stops):		
		#freqs = np.fft.fftfreq(stop-start, time_step)
                freqs = np.fft.fftfreq(fs)*fs
        
		ps = abs(np.fft.fft(eegdata[start:stop],fs))**2
		theta = sum([ps[i] for i,freq in enumerate(freqs) if freqs[i] >= 6 and freqs[i] <= 10])
		gamma = sum([ps[i] for i, freq in enumerate(freqs) if freqs[i] >= 25 and freqs[i] <= 160])
		delta = sum([ps[i] for i,freq in enumerate(freqs) if freqs[i] >= 1 and freqs[i] <= 3])
		theta_index.append(theta/delta)
		gamma_index.append(gamma/delta)
                pstot[:,j] = ps[:]
                j+=1

	return theta_index, gamma_index, pstot, freqs


#def returnPeakIndices(eegvals, m, t):
#     ind_peaks = peakutils.indexes(eegvals, thres=t, min_dist=m)
#     ind_troughs = peakutils.indexes(eegvals*-1, thres=t, min_dist=m)
#     return ind_peaks, ind_troughs
     
def returnPhases(peaks, troughs, flags):
  phases = []
  for flag in flags:

    nearest_peak = takeClosest(peaks, flag)
    nearest_trough = takeClosest(troughs, flag)
    #print("peak {}, trough {}, event {}".format(nearest_peak, nearest_trough, event))
    offset = .5*(flag - nearest_peak)/(nearest_trough-nearest_peak)
    #print(offset)
    if (nearest_trough > nearest_peak):
      phase = .5+offset
    else:
      phase = .5-offset
    phases.append(phase)
    
  return phases
