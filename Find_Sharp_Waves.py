#compare theta band for walking/not walking

import EEGUtils as eeg
from matplotlib import pyplot as plt
import matplotlib.lines as lines
import numpy as np
import sys
from scipy.signal import hilbert, chirp
import peakutils

def returnPeakIndices(eegvals, m, t):
     ind_peaks = peakutils.indexes(eegvals, thres=t, min_dist=m)
     ind_troughs = peakutils.indexes(eegvals*-1, thres=t, min_dist=m)
     return ind_peaks, ind_troughs

def return_grouping_indicator(ts, index, thresh):
  indicator = np.zeros(len(ts))
  indicator[index] = 1
  print(np.shape(index))
  grouped = group_consecutives(index)
  for g in grouped:
#    print(ts[g[-1]]-ts[g[0]])
    if (ts[g[-1]]-ts[g[0]] < thresh):
      indicator[g] = 0
  return indicator

def group_consecutives(vals, step=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result


fs = 512/.257552

eegfile = sys.argv[1]

if len(sys.argv) == 3:
   plotit = sys.argv[2]
else:
   plotit = 'noplot'

eegdata, eegtimestamps = eeg.readEEG(eegfile)
#ts = np.linspace(eegtimestamps[0], eegtimestamps[-1], len(eegdata))
#print(1000000*len(eegdata)/(ts[-1]-ts[0]))

eegdata = eeg.butter_bandpass_filter(eegdata, 2, 500, fs, order=4)

npzfile = np.load("./RawData/NOTMOVING.npz")
starts=npzfile['arr_0'].astype(int)
stops=npzfile['arr_1'].astype(int)

#find closest eegtimestamps to tstart, tstop
EEGstart = [eeg.takeClosest(eegtimestamps, start) for start in starts]
EEGstop = [eeg.takeClosest(eegtimestamps, stop) for stop in stops] 

startidx = [eegtimestamps.index(start)*512 for start in EEGstart]
stopidx =  [eegtimestamps.index(stop)*512 for stop in EEGstop]

epoch = [list(eegdata[start:stop]) for start, stop in zip(startidx,stopidx)]
#plt.plot(epoch[0])
#plt.show()

big_e = sum(epoch,[])
ts = np.linspace(0, 1000*len(big_e)/(1988), len(big_e))
#plt.plot(ts,big_e)
#plt.show()
theta = eeg.butter_bandpass_filter(big_e, 2,15, fs, order=4)
filtered = eeg.butter_bandpass_filter(big_e, 120, 240, fs, order=4)
analytic_signal = hilbert(filtered)
amp_env = np.abs(analytic_signal)
instantaneous_phase = np.unwrap(np.angle(analytic_signal))
instantaneous_frequency = (np.diff(instantaneous_phase) / (2.0*np.pi) * fs)
std_amp_env = np.std(amp_env)
   
upper_rip_thresh = np.mean(amp_env) + (np.std(amp_env) * 2)
lower_rip_thresh = np.mean(amp_env) + (np.std(amp_env) * 1)

#plt.plot(ts,amp_env)
#plt.plot(filtered, 'k')
#plt.plot([0, len(amp_env)], [upper_rip_thresh, upper_rip_thresh], 'c')
#plt.plot([0, ts[-1]], [lower_rip_thresh, lower_rip_thresh], 'r')
#plt.show()

above_index = np.where(amp_env > upper_rip_thresh)
below_index = np.where(amp_env > lower_rip_thresh)

mirror_above = np.ones(len(amp_env))
mirror_above[above_index[0]] = 0
mirror_index = np.where(mirror_above == 1)

mirror_above = np.abs(return_grouping_indicator(ts, mirror_index[0],50) -1)

#plt.plot(ts,mirror_above*upper_rip_thresh, 'k')
#plt.show()
above_index = np.where(mirror_above == 1)
ind_above = return_grouping_indicator(ts, above_index[0],35)

#plt.plot(ts, ind_above*upper_rip_thresh, 'k')
#plt.show()

mirror_below = np.ones(len(amp_env))
mirror_below[below_index[0]] = 0
mirror_index = np.where(mirror_below == 1)

mirror_below = np.abs(return_grouping_indicator(ts, mirror_index[0],50) -1)
below_index = np.where(mirror_below == 1)
ind_below = return_grouping_indicator(ts, below_index[0], 35)


ripbottom = np.zeros(len(ts), dtype=np.uint16)

grouped_below = group_consecutives(np.where(ind_below==1)[0])

for i,g in enumerate(grouped_below):
  ripbottom[g] = i
   
rips = ind_above * ripbottom
#plt.plot(ts,rips,'r')
#plt.show()


vals = np.unique(rips)

centered_rips = np.empty([len(vals), 4000])
centered_rips[:] = np.nan


for i,val in enumerate(vals):
  print("i: {}, val: {}".format(i,val))
  start = grouped_below[int(val)][0]
  stop = grouped_below[int(val)][-1]
#  peak = np.max(
  rip = filtered[start:stop]
  ripts = np.arange(len(rip))
  peaks, troughs = returnPeakIndices(rip ,10, .25)
  print(peaks)
  print(troughs)
  center_trough = troughs[np.round(len(troughs)/2)]
  print(center_trough)
  c_to_end = rip[center_trough:]
  c_to_begin = rip[:center_trough]

  centered_rips[i,2000:2000+len(c_to_end)] = c_to_end
  centered_rips[i,2000-len(c_to_begin):2000] = c_to_begin

 
  print("rip: {},  duration -- start: {}, stop: {}".format(int(val), start, stop))
  #plt.plot(ripts, rip)
  #plt.plot(ripts[troughs], rip[troughs], 'r.')
   
# plt.plot(theta[start-100:stop+100], 'r')
  #plt.show()
  #plt.plot(centered_rips[i,:])
  #plt.show()
#if plotit == 'plot':
#  plt.plot(filtered)
#  plt.plot(big_e)
#  plt.plot(rips*np.max(big_e), 'k') 
#  plt.plot(rips*upper_rip_thresh,'k')
#  plt.show()   
 
#outfile = eegfile.replace('.Ncs', '_RIPS')
#np.savez(outfile, outfiltered, outrips, outupper, outlower, outts)

rip_mean = np.nanmean(centered_rips,0)

plt.plot(rip_mean)
plt.show()
