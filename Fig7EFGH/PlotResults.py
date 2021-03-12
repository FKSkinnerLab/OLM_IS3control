### Plot Results
from __future__ import division
import numpy
import efel
import time
import random
import scipy
from scipy import signal
from scipy import stats
import matplotlib
from matplotlib import pyplot
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import Axes3D

# Cell_Name = 'Gormadoc'
Cell_Name = 'Isembard'
Case = Cell_Name + '_E_COM_I_COM'

# HC Treshold Measurement Values
tstop = 10000 # seconds
dt = 0.1
font_size = 13
numrerands = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'numrerands.npy')

Examples = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelParams.npy')
ExampleStrings = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelStrings.npy')
ExampleString = ExampleStrings[0]

timevec = numpy.arange(1000,int(tstop),dt,dtype=numpy.float)
modfreqs = numpy.array([8,8],dtype=numpy.float)
addSUPs = numpy.array([0,1,1,1])
ModShifts = numpy.array([0,-15.625,0,15.625]) # i.e. [Baseline,EC-Modulated,Balanced Modulation,CA3-Modulated]

volts = numpy.zeros((numrerands,len(addSUPs),len(timevec)), dtype=numpy.float)
sptim = numpy.zeros((numrerands,len(addSUPs),len(timevec)), dtype=numpy.float)
for rep in range(1,numrerands+1):
	for l in range(len(addSUPs)):
		volts[rep-1][l] = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_AddSupp_' + str(addSUPs[l]) + '_ModShift_' + str(ModShifts[l]) + '_vvec.npy')
		sptim[rep-1][l] = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_AddSupp_' + str(addSUPs[l]) + '_ModShift_' + str(ModShifts[l]) + '_vvecbinary.npy')

baseIS3SynNum = 30
RateOfSynIncrease = 7 # synapses/cycle
SupStartTime = 2000 + (15.625) # ms
Period = 125 # ms
cycle = 0
SpikesPerCycleVec = numpy.zeros((len(timevec),), dtype=numpy.float)
for t in range(0,len(timevec)):
	if ((timevec[t] >= SupStartTime) & (Period >= 125)):
		cycle = cycle + 1
		Period = 0
	elif ((timevec[t] >= SupStartTime) & (Period < 125)):
		Period = Period + dt
	
	SpikesPerCycleVec[t] = SpikesPerCycleVec[t] + baseIS3SynNum + RateOfSynIncrease*cycle

fig, (ax0) = pyplot.subplots(nrows=1)
ax0.plot(timevec/1000,SpikesPerCycleVec,color='k')
ylims = ax0.get_ylim()
ax0.set_xlim(numpy.amin(timevec/1000),numpy.amax(timevec/1000))
ax0.set_ylim(ylims[0],ylims[1])
ax0.set_xlabel('Time [s]')
ax0.set_ylabel('IS3 Spikes Per Cycle')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_SpikesPerCycle.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_SpikesPerCycle.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

fmin = 0
fmax = 30
for rep in range(0,numrerands):
	fig, (ax0, ax1, ax2, ax3) = pyplot.subplots(nrows=4)
	Fs = 1/(dt/1000)
	f, t, Sxx = scipy.signal.spectrogram(volts[rep][0], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax0.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=5.5)
	ax0.set_ylabel('Frequency [Hz]')
	cbar0 = fig.colorbar(im, ax=ax0)
	cbar0.set_label('PSD')
	ax0.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(volts[rep][1], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax1.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=5.5)
	ax1.set_ylabel('Frequency [Hz]')
	cbar1 = fig.colorbar(im, ax=ax1)
	cbar1.set_label('PSD')
	ax1.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(volts[rep][2], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax2.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=5.5)
	ax2.set_ylabel('Frequency [Hz]')
	cbar2 = fig.colorbar(im, ax=ax2)
	cbar2.set_label('PSD')
	ax2.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(volts[rep][3], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax3.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=5.5)
	ax3.set_ylabel('Frequency [Hz]')
	cbar3 = fig.colorbar(im, ax=ax3)
	cbar3.set_label('PSD')
	ax3.axis(ymin=fmin, ymax=fmax)
	ax3.set_xlabel('Time [sec]')
	fig.tight_layout()
	# pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltSpec.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltSpec.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	fig, (ax0, ax1, ax2, ax3) = pyplot.subplots(nrows=4)
	Fs = 1/(dt/1000)
	f, t, Sxx = scipy.signal.spectrogram(sptim[rep][0], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax0.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=0.000004)
	ax0.set_ylabel('Frequency [Hz]')
	cbar0 = fig.colorbar(im, ax=ax0, format='%.0e')
	cbar0.set_label('PSD')
	ax0.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(sptim[rep][1], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax1.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=0.000004)
	ax1.set_ylabel('Frequency [Hz]')
	cbar1 = fig.colorbar(im, ax=ax1, format='%.0e')
	cbar1.set_label('PSD')
	ax1.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(sptim[rep][2], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax2.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=0.000004)
	ax2.set_ylabel('Frequency [Hz]')
	cbar2 = fig.colorbar(im, ax=ax2, format='%.0e')
	cbar2.set_label('PSD')
	ax2.axis(ymin=fmin, ymax=fmax)
	
	f, t, Sxx = scipy.signal.spectrogram(sptim[rep][3], fs=Fs, window=('gaussian', 25000), nperseg=500, noverlap=None, nfft=int(Fs*5), scaling='density', mode='psd')
	im = ax3.pcolormesh(t+1, f[(f>fmin) & (f<fmax)], Sxx[(f>fmin) & (f<fmax)], cmap='inferno', vmin=0, vmax=0.000004)
	ax3.set_ylabel('Frequency [Hz]')
	cbar3 = fig.colorbar(im, ax=ax3, format='%.0e')
	cbar3.set_label('PSD')
	ax3.axis(ymin=fmin, ymax=fmax)
	ax3.set_xlabel('Time [sec]')
	fig.tight_layout()
	# pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_SpikeSpec.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_SpikeSpec.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	# Instantaneous Frequency Plots
	spiketimes0 = (timevec[sptim[rep][0]==1])/1000
	spiketimes1 = (timevec[sptim[rep][1]==1])/1000
	spiketimes2 = (timevec[sptim[rep][2]==1])/1000
	spiketimes3 = (timevec[sptim[rep][3]==1])/1000
	bin_width = 1 # seconds
	interval = 0.25 # seconds
	Bins = numpy.arange(1+bin_width/2,10-bin_width/2,interval,dtype=numpy.float)
	inst_freqs0 = numpy.zeros((len(Bins),), dtype=numpy.float)
	inst_freqs1 = numpy.zeros((len(Bins),), dtype=numpy.float)
	inst_freqs2 = numpy.zeros((len(Bins),), dtype=numpy.float)
	inst_freqs3 = numpy.zeros((len(Bins),), dtype=numpy.float)
	count = 0
	for b in Bins:
		inst_freqs0[count] = len(spiketimes0[(spiketimes0>(b-bin_width/2)) & (spiketimes0<(b+bin_width/2))])/(bin_width)
		inst_freqs1[count] = len(spiketimes1[(spiketimes1>(b-bin_width/2)) & (spiketimes1<(b+bin_width/2))])/(bin_width)
		inst_freqs2[count] = len(spiketimes2[(spiketimes2>(b-bin_width/2)) & (spiketimes2<(b+bin_width/2))])/(bin_width)
		inst_freqs3[count] = len(spiketimes3[(spiketimes3>(b-bin_width/2)) & (spiketimes3<(b+bin_width/2))])/(bin_width)
		count = count + 1
	
	fig, (ax0) = pyplot.subplots(nrows=1)
	ax0.plot(Bins,inst_freqs0,color='k',label='Baseline')
	ax0.plot(Bins,inst_freqs1,color='r',label='EC-Mod')
	ax0.plot(Bins,inst_freqs2,color='g',label='Even-Mod')
	ax0.plot(Bins,inst_freqs3,color='b',label='CA3-Mod')
	ylims = ax0.get_ylim()
	ax0.plot(numpy.array([2,2]),numpy.array([ylims[0],ylims[1]]),linestyle='--',color='k')
	# ax0.legend(loc='upper left')
	ax0.set_xlim(Bins[0],Bins[-1:])
	ax0.set_ylim(ylims[0],ylims[1])
	ax0.set_xlabel('Time [sec]')
	ax0.set_ylabel('Spike Rate [Hz]')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_InstaFreq.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_InstaFreq.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	fig, (ax0, ax1, ax2, ax3) = pyplot.subplots(nrows=4)
	ax0.plot(timevec/1000,volts[rep][0],color='k',label='Baseline')
	ax0.plot(numpy.array([2,2]),numpy.array([numpy.amin(volts[rep][0]),numpy.amax(volts[rep][0])]),linestyle='--',color='k')
	ax0.set_xlim(1,10)
	ax0.legend(loc='upper left')
	# ax0.set_ylabel('Voltage [mV]')
	ax1.plot(timevec/1000,volts[rep][1],color='r',label='EC-Mod')
	ax1.plot(numpy.array([2,2]),numpy.array([numpy.amin(volts[rep][1]),numpy.amax(volts[rep][1])]),linestyle='--',color='k')
	ax1.set_xlim(1,10)
	ax1.legend(loc='upper left')
	# ax1.set_ylabel('Voltage [mV]')
	ax2.plot(timevec/1000,volts[rep][2],color='g',label='Even-Mod')
	ax2.plot(numpy.array([2,2]),numpy.array([numpy.amin(volts[rep][2]),numpy.amax(volts[rep][2])]),linestyle='--',color='k')
	ax2.set_xlim(1,10)
	ax2.legend(loc='upper left')
	# ax2.set_ylabel('Voltage [mV]')
	ax3.plot(timevec/1000,volts[rep][3],color='b',label='CA3-Mod')
	ax3.plot(numpy.array([2,2]),numpy.array([numpy.amin(volts[rep][3]),numpy.amax(volts[rep][3])]),linestyle='--',color='k')
	ax3.set_xlim(1,10)
	ax3.legend(loc='upper left')
	# ax3.set_ylabel('Voltage [mV]')
	ax3.set_xlabel('Time [sec]')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_Voltage.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_Voltage.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	fig, (ax0, ax1, ax2, ax3) = pyplot.subplots(nrows=4)
	ax0.plot(timevec/1000,volts[rep][0],color='k',label='Baseline')
	ax0.set_xlim(9,10)
	ax0.legend(loc='upper left')
	# ax0.set_ylabel('Voltage [mV]')
	ax1.plot(timevec/1000,volts[rep][1],color='r',label='EC-Mod')
	ax1.set_xlim(9,10)
	ax1.legend(loc='upper left')
	# ax1.set_ylabel('Voltage [mV]')
	ax2.plot(timevec/1000,volts[rep][2],color='g',label='Even-Mod')
	ax2.set_xlim(9,10)
	ax2.legend(loc='upper left')
	# ax2.set_ylabel('Voltage [mV]')
	ax3.plot(timevec/1000,volts[rep][3],color='b',label='CA3-Mod')
	ax3.set_xlim(9,10)
	ax3.legend(loc='upper left')
	# ax3.set_ylabel('Voltage [mV]')
	ax1.set_xlabel('Time [sec]')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltageZoomed.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltageZoomed.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	f0, Pxx_den0 = signal.welch(sptim[rep][0], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f1, Pxx_den1 = signal.welch(sptim[rep][1], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f2, Pxx_den2 = signal.welch(sptim[rep][2], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f3, Pxx_den3 = signal.welch(sptim[rep][3], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	fig, (ax0) = pyplot.subplots(nrows=1)
	ax0.plot(f0,Pxx_den0,color='k',label='Baseline')
	ax0.plot(f1,Pxx_den1,color='r',label='EC-Mod')
	ax0.plot(f2,Pxx_den2,color='g',label='Even-Mod')
	ax0.plot(f3,Pxx_den3,color='b',label='CA3-Mod')
	ylims = ax0.get_ylim()
	ax0.plot(numpy.array([8,8]),numpy.array([ylims[0],ylims[1]]),linestyle='--',color='k')
	ax0.set_xlim(0,30)
	ax0.set_ylim(ylims[0],ylims[1])
	ax0.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax0.legend(loc='upper left')
	ax0.set_xlabel('Frequency [Hz]')
	ax0.set_ylabel('Power')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTSpikeSpec.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTSpikeSpec.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	f0, Pxx_den0 = signal.welch(volts[rep][0], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f1, Pxx_den1 = signal.welch(volts[rep][1], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f2, Pxx_den2 = signal.welch(volts[rep][2], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	f3, Pxx_den3 = signal.welch(volts[rep][3], fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	fig, (ax0) = pyplot.subplots(nrows=1)
	ax0.plot(f0,Pxx_den0,color='k',label='Baseline')
	ax0.plot(f1,Pxx_den1,color='r',label='EC-Mod')
	ax0.plot(f2,Pxx_den2,color='g',label='Even-Mod')
	ax0.plot(f3,Pxx_den3,color='b',label='CA3-Mod')
	ylims = ax0.get_ylim()
	ax0.plot(numpy.array([8,8]),numpy.array([ylims[0],ylims[1]]),linestyle='--',color='k')
	ax0.set_xlim(0,30)
	ax0.set_ylim(ylims[0],ylims[1])
	ax0.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	# ax0.legend(loc='upper left')
	ax0.set_xlabel('Frequency [Hz]')
	ax0.set_ylabel('Power')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTVolt.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTVolt.png')
	ax0.set_xlim(7,9)
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTVoltZoomed.pdf')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_FFTVoltZoomed.png')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
	
	tvec0 = timevec # ms
	vvec0 = volts[rep][0]
	vvec1 = volts[rep][1]
	vvec2 = volts[rep][2]
	vvec3 = volts[rep][3]
	tvec0split = numpy.split(tvec0,72)
	volt0split = numpy.split(vvec0,72)
	volt1split = numpy.split(vvec1,72)
	volt2split = numpy.split(vvec2,72)
	volt3split = numpy.split(vvec3,72)
	volt0mean = numpy.mean(volt0split,axis=0)
	volt1mean = numpy.mean(volt1split,axis=0)
	volt2mean = numpy.mean(volt2split,axis=0)
	volt3mean = numpy.mean(volt3split,axis=0)
	volt0std = numpy.std(volt0split,axis=0)
	volt1std = numpy.std(volt1split,axis=0)
	volt2std = numpy.std(volt2split,axis=0)
	volt3std = numpy.std(volt3split,axis=0)
	f, axarr = matplotlib.pyplot.subplots(4)
	axarr[0].fill_between(tvec0split[0],volt0mean-volt0std,volt0mean+volt0std,facecolor='k', alpha=0.5)
	axarr[0].plot(tvec0split[0],volt0mean,color='k')
	axarr[1].fill_between(tvec0split[0],volt1mean-volt1std,volt1mean+volt1std,facecolor='red', alpha=0.5)
	axarr[1].plot(tvec0split[0],volt1mean,color='red')
	axarr[2].fill_between(tvec0split[0],volt2mean-volt2std,volt2mean+volt2std,facecolor='green', alpha=0.5)
	axarr[2].plot(tvec0split[0],volt2mean,color='green')
	axarr[3].fill_between(tvec0split[0],volt3mean-volt3std,volt3mean+volt3std,facecolor='blue', alpha=0.5)
	axarr[3].plot(tvec0split[0],volt3mean,color='blue')
	axarr[0].tick_params(labelsize=font_size)
	axarr[1].tick_params(labelsize=font_size)
	axarr[2].tick_params(labelsize=font_size)
	axarr[3].tick_params(labelsize=font_size)
	axarr[3].set_xlabel('Time (ms)',fontsize=font_size)
	# axarr[1].set_ylabel('Voltage (mV)',fontsize=font_size)
	axarr[0].set_xlim(1000,1125)
	axarr[1].set_xlim(1000,1125)
	axarr[2].set_xlim(1000,1125)
	axarr[3].set_xlim(1000,1125)
	axarr[0].set_xticks([1000, 1031.25, 1062.5, 1093.75, 1125])
	axarr[1].set_xticks([1000, 1031.25, 1062.5, 1093.75, 1125])
	axarr[2].set_xticks([1000, 1031.25, 1062.5, 1093.75, 1125])
	axarr[3].set_xticks([1000, 1031.25, 1062.5, 1093.75, 1125])
	axarr[0].set_xticklabels(['', '', '', '', ''],fontsize=font_size)
	axarr[1].set_xticklabels(['', '', '', '', ''],fontsize=font_size)
	axarr[2].set_xticklabels(['', '', '', '', ''],fontsize=font_size)
	axarr[3].set_xticklabels(['0', '31.25', '62.5', '93.75', '125'],fontsize=font_size)
	
	axarr2 = axarr[0].twiny()
	axarr2.set_xlim(axarr[0].get_xlim())
	axarr2.set_xticks([1000, 1031.25, 1062.5, 1093.75, 1125])
	axarr2.set_xticklabels([r"$0^\circ$", r"$90^\circ$", r"$180^\circ$", r"$270^\circ$", r"$360^\circ$"],fontsize=font_size)
	# axarr2.axis["right"].major_ticklabels.set_visible(False)
	# axarr2.axis["top"].major_ticklabels.set_visible(True)
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltAvg.pdf', bbox_inches='tight')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Example_' + str(rep) + '_VoltAvg.png', bbox_inches='tight')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()
