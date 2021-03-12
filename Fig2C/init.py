from __future__ import division
import numpy
import neuron
from neuron import h, gui
import matplotlib.pyplot as plt
import efel
import glob
import IPython, os
import scipy
from scipy import signal
from scipy import stats
from currents_visualization import plotCurrentscape

font_size = 15

### Instantiate Model ###
h.load_file("init_model.hoc")
h.load_file("SynParamSearch.hoc")

cell = h.cell
h.ic_hold.amp = 0
find7HzHolding = False
if cell == 1:
	cellname = 'Gormadoc'
	currhold = 0.049 # 7.25 Hz
elif cell == 2:
	cellname = 'Isembard'
	currhold = 0.045 # 7.25 Hz

modelnum = h.modelnum
if modelnum == 1:
	modelname = '1st'
elif modelnum == 2:
	modelname = '2nd'
elif modelnum == 3:
	modelname = '3rd'
elif modelnum == 4:
	modelname = '4th'
elif modelnum == 5:
	modelname = '5th'
elif modelnum == 6:
	modelname = 'FT'

### Load Experimental ###
Exp_Dir = cellname + '_exp_data'

files_list0 = glob.glob1(Exp_Dir, "*120pA_dnqx.dat")
files_list1 = glob.glob1(Exp_Dir, "*30pA_dnqx.dat")
files_list2 = glob.glob1(Exp_Dir, "*60pA_dnqx.dat")
files_list3 = glob.glob1(Exp_Dir, "*90pA_dnqx.dat")

data0 = numpy.loadtxt(os.path.join(Exp_Dir,files_list0[0]))
data1 = numpy.loadtxt(os.path.join(Exp_Dir,files_list1[0]))
data2 = numpy.loadtxt(os.path.join(Exp_Dir,files_list2[0]))
data3 = numpy.loadtxt(os.path.join(Exp_Dir,files_list3[0]))
data0 = data0[1:80000,:]
data1 = data1[1:80000,:]
data2 = data2[1:80000,:]
data3 = data3[1:80000,:]

### Parameters ###
Nsyns = numpy.array([0,6,12,18,24,30,36,42,48,54,60])
freqs = numpy.array([1,5,8,10,20])

### Setup holding currents ###

h.cvode_active(1) # i.e. to make it run faster for this portion
h.cvode.atol(1e-05)
h.tstop = 5000

spikerate_prev = 0
if find7HzHolding:
	for stepcurr in range(40,101,1):
		h.ic_hold.amp = stepcurr*0.001 # holding current
		h.f(0,0)
		apctimes = numpy.array(h.apctimes)
		apctimes = apctimes[apctimes > 1000] # Only spike times > 1s
		spikerate = len(apctimes)/4
		print('Spike Rate = ' + str(spikerate) + ' Hz at ' + str(stepcurr*0.001) + ' nA')
		if (spikerate > 7) & (spikerate < 7.6):
			print('Holding Current = ' + str(stepcurr*0.001) + ' nA')
			print('Baseline Spike Rate = ' + str(spikerate) + ' Hz')
			break
		elif spikerate > 7.6:
			h.ic_hold.amp = (stepcurr-1)*0.001 # holding current
			print('Spike Rate Passed Acceptable Range')
			print('Holding Current = ' + str((stepcurr-1)*0.001) + ' nA')
			print('Baseline Spike Rate = ' + str(spikerate_prev) + ' Hz')
			break
		spikerate_prev = spikerate
else:
	h.ic_hold.amp = currhold

h.cvode_active(0) # i.e turn cvode off for simulations
h.tstop = 10000

## For Gormadoc the holding current was 0.049 nA (gives a ~7.25 Hz spike rate)
## For Isembard the holding current was 0.045 nA (gives a ~7.25 Hz spike rate)

h.ic_step.amp = 0 # i.e. only use this if injecting current steps

### Setup Recordings and Distance Vector ###
t_vec = h.Vector()
t_vec.record(h._ref_t)

v_vec = h.Vector()
v_vec.record(h.soma[0](0.5)._ref_v)

v_vec_d = h.Vector()
v_vec_d.record(h.dend[0](0)._ref_v)

Ih = h.Vector()
Ih.record(h.dend[0](0)._ref_ih_Ih)

INa = h.Vector()
INa.record(h.dend[0](0)._ref_ina_Nadend)

IKa = h.Vector()
IKa.record(h.dend[0](0)._ref_ik_Ika)

IKdrf = h.Vector()
IKdrf.record(h.dend[0](0)._ref_ik_Ikdrf)

IKdrs = h.Vector()
IKdrs.record(h.dend[0](0)._ref_ik_Ikdrs)

Im = h.Vector()
Im.record(h.dend[0](0)._ref_ik_IM)

IKCa = h.Vector()
IKCa.record(h.dend[0](0)._ref_ik_kca)

ICaL = h.Vector()
ICaL.record(h.dend[0](0)._ref_ica_cal)

ICaT = h.Vector()
ICaT.record(h.dend[0](0)._ref_ica_cat)

Il = h.Vector()
Il.record(h.dend[0](0)._ref_i_passsd)

labels = ('gKA','gKdrf','gKdrs','gKCa','gM','gL','gNa','gCaT','gCaL','gH')

Areas = numpy.zeros((len(freqs),len(Nsyns)))
Power = numpy.zeros((len(freqs),len(Nsyns)))
AreasSubtracted = numpy.zeros((len(freqs),len(Nsyns)))
PowerSubtracted = numpy.zeros((len(freqs),len(Nsyns)))

rcount = 0
for rate in freqs:
	scount = 0
	fi, axarr = plt.subplots(11, sharex = True, sharey = True)
	fi3, axarr3 = plt.subplots(11, sharex = True, sharey = True)
	for synnum in Nsyns:
		h.f(synnum,rate) # Sets synaptic inputs
		if scount == 0 or scount == 10:
			print('Plotting Currentscape')
			TStart = 90001
			vvec0 = numpy.array(v_vec_d)
			tvec0 = numpy.array(t_vec)
			pIKa = numpy.array(IKa)
			pIKdrf = numpy.array(IKdrf)
			pIKdrs = numpy.array(IKdrs)
			pIKCa = numpy.array(IKCa)
			pIm = numpy.array(Im)
			pIl = numpy.array(Il)
			pINa = numpy.array(INa)
			pICaT = numpy.array(ICaT)
			pICaL = numpy.array(ICaL)
			pIh = numpy.array(Ih)
			
			r0 = [pIKa[TStart:],pIKdrf[TStart:],pIKdrs[TStart:],pIKCa[TStart:],pIm[TStart:],pIl[TStart:],pINa[TStart:],pICaT[TStart:],pICaL[TStart:],pIh[TStart:]]
			
			fig0 = plotCurrentscape(vvec0[TStart:], r0)
			fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_Frequency_' + str(rate) + 'Hz_SynapseNum_' + str(synnum) + '.pdf',dpi=500)
			fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_Frequency_' + str(rate) + 'Hz_SynapseNum_' + str(synnum) + '.png',dpi=500)
			fig0.clf()
			plt.close(fig0)
		
		vvec = numpy.array(v_vec)
		tvec = numpy.array(t_vec)
		numpy.save('Output/' + cellname + modelname + '_NumberOfSyns_' + str(synnum) + '_Frequency_' + str(rate) + '_Voltage.npy',vvec)
		numpy.save('Output/' + cellname + modelname + '_NumberOfSyns_' + str(synnum) + '_Frequency_' + str(rate) + '_Time.npy',tvec)
		apctimes = numpy.array(h.apctimes)
		spikebinvec = numpy.zeros(len(tvec))
		for x in apctimes: spikebinvec[numpy.where(tvec == x)] = 1
		
		dt = h.dt
		f, Pxx_den = signal.welch(spikebinvec, 1/(dt/1000), nperseg=25000)
		if synnum == 0:
			Pxx_den0 = Pxx_den
			f0 = f
		Areas[rcount][scount] = numpy.trapz(Pxx_den[(f>4) & (f<12)],x=f[(f>4) & (f<12)])
		Power[rcount][scount] = numpy.max(Pxx_den[(f>rate-1) & (f<rate+1)])
		# Power[rcount][scount] = numpy.max(Pxx_den[f==rate])
		if synnum != 0:
			Pxx_den_Subtracted = Pxx_den-Pxx_den0
			AreasSubtracted[rcount][scount] = numpy.trapz(Pxx_den[(f>4) & (f<12)],x=f[(f>4) & (f<12)])-numpy.trapz(Pxx_den0[(f0>4) & (f0<12)],x=f0[(f0>4) & (f0<12)])
			PowerSubtracted[rcount][scount] = numpy.max(Pxx_den_Subtracted[(f>rate-1) & (f<rate+1)])
			# Power[rcount][scount] = numpy.max(Pxx_den[f==rate])
		
		fi2, axarr2 = plt.subplots(1)
		axarr2.loglog(f, Pxx_den,'k')
		axarr2.set_xlim(0,100)
		axarr2.set_xlabel('frequency (Hz)')
		axarr2.set_ylabel(r'$PSD (Spikes^2 / Hz)$')
		axarr2.spines['right'].set_visible(False)
		axarr2.spines['top'].set_visible(False)
		fi2.tight_layout()
		fi2.savefig('PLOTFFTfiles/' + cellname + modelname + '_PSD_' + str(rate) + 'Hz_Rate_' + str(synnum) + '_NumSynapses.pdf', bbox_inches='tight')
		fi2.savefig('PLOTFFTfiles/' + cellname + modelname + '_PSD_' + str(rate) + 'Hz_Rate_' + str(synnum) + '_NumSynapses.png', bbox_inches='tight')
		# pyplot.gcf().clear()
		# plt.cla()
		fi2.clf()
		plt.close(fi2)
		
		axarr[scount].plot(tvec, vvec, color='k')
		axarr[scount].set_xlim(9000,10000)
		axarr[scount].set_ylim(-85,30)
		axarr[scount].spines['right'].set_visible(False)
		axarr[scount].spines['top'].set_visible(False)
		if scount != 10:
			axarr[scount].spines['bottom'].set_visible(False)
			for tic in axarr[scount].xaxis.get_major_ticks():
				tic.tick1On = tic.tick2On = False
		
		tvec1 = tvec[10001:]
		vvec1 = vvec[10001:]
		tvecsplit = numpy.split(tvec1,9000/((1/rate)*1000))
		vvecsplit = numpy.split(vvec1,9000/((1/rate)*1000))
		vvecmean = numpy.mean(vvecsplit,axis=0)
		vvecstd = numpy.std(vvecsplit,axis=0)
		axarr3[scount].fill_between(tvecsplit[0],vvecmean-vvecstd,vvecmean+vvecstd,facecolor='black', alpha=0.5)
		axarr3[scount].plot(tvecsplit[0],vvecmean,color='black')
		axarr3[scount].set_xlim(1000,1000+(1/rate)*1000)
		if scount != 0:
			axarr3[scount].spines['top'].set_visible(False)
		if scount == 0:
			axarr4 = axarr3[scount].twiny()
			axarr4.spines['bottom'].set_visible(False)
			axarr4.set_xlim(axarr3[scount].get_xlim())
			axarr4.set_xticks([1000, 1000+((1/rate)*1000)*1/4, 1000+((1/rate)*1000)*1/2, 1000+((1/rate)*1000)*3/4, 1000+(1/rate)*1000])
			axarr4.set_xticklabels([r"$0^\circ$", r"$90^\circ$", r"$180^\circ$", r"$270^\circ$", r"$360^\circ$"],fontsize=font_size)
		if scount != 10:
			axarr3[scount].spines['bottom'].set_visible(False)
			for tic in axarr3[scount].xaxis.get_major_ticks():
				tic.tick1On = tic.tick2On = False
		if scount == 10:
			# axarr3[scount].tick_params(labelsize=font_size)
			axarr3[scount].set_xlabel('Time (ms)',fontsize=font_size)
			axarr3[scount].set_xticks([1000, 1000+((1/rate)*1000)*1/4, 1000+((1/rate)*1000)*1/2, 1000+((1/rate)*1000)*3/4, 1000+(1/rate)*1000])
			axarr3[scount].set_xticklabels(['0', str(((1/rate)*1000)*1/4), str(((1/rate)*1000)*1/2), str(((1/rate)*1000)*3/4), str((1/rate)*1000)],fontsize=font_size)
		
		scount = scount + 1
	
	axarr[10].set_xlabel('Time (ms)',fontsize = font_size)
	fi.savefig('PLOTfiles/' + cellname + modelname + '_' + str(rate) + 'Hz_Traces.pdf', bbox_inches='tight')
	fi.savefig('PLOTfiles/' + cellname + modelname + '_' + str(rate) + 'Hz_Traces.png', bbox_inches='tight')
	# plt.gcf().clear()
	# plt.cla()
	fi.clf()
	plt.close(fi)
	
	fi3.savefig('PLOTfiles/' + cellname + modelname + '_' + str(rate) + 'Hz_AveragedCycles.pdf', bbox_inches='tight')
	fi3.savefig('PLOTfiles/' + cellname + modelname + '_' + str(rate) + 'Hz_AveragedCycles.png', bbox_inches='tight')
	# plt.gcf().clear()
	# plt.cla()
	fi3.clf()
	plt.close(fi3)
	
	rcount = rcount + 1

plt.plot(Nsyns,Areas[0], color = 'b', label = '1 Hz')
plt.plot(Nsyns,Areas[1], color = 'r', label = '5 Hz')
plt.plot(Nsyns,Areas[2], color = 'g', label = '8 Hz')
plt.plot(Nsyns,Areas[3], color = 'c', label = '10 Hz')
plt.plot(Nsyns,Areas[4], color = 'm', label = '20 Hz')
plt.legend()
plt.xlabel('Number of IS3 Synapses',fontsize = font_size)
plt.ylabel('Power Area (4-12 Hz)',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAreas.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAreas.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(Nsyns,Power[0], color = 'b', label = '1 Hz')
plt.plot(Nsyns,Power[1], color = 'r', label = '5 Hz')
plt.plot(Nsyns,Power[2], color = 'g', label = '8 Hz')
plt.plot(Nsyns,Power[3], color = 'c', label = '10 Hz')
plt.plot(Nsyns,Power[4], color = 'm', label = '20 Hz')
plt.legend()
plt.xlabel('Number of IS3 Synapses',fontsize = font_size)
plt.ylabel('Power at Input Frequency',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAtInputFrequencies.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAtInputFrequencies.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(Nsyns,AreasSubtracted[0], color = 'b', label = '1 Hz')
plt.plot(Nsyns,AreasSubtracted[1], color = 'r', label = '5 Hz')
plt.plot(Nsyns,AreasSubtracted[2], color = 'g', label = '8 Hz')
plt.plot(Nsyns,AreasSubtracted[3], color = 'c', label = '10 Hz')
plt.plot(Nsyns,AreasSubtracted[4], color = 'm', label = '20 Hz')
plt.legend()
plt.xlabel('Number of IS3 Synapses',fontsize = font_size)
plt.ylabel('Power Area (4-12 Hz)',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAreasSubtracted.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAreasSubtracted.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(Nsyns,PowerSubtracted[0], color = 'b', label = '1 Hz')
plt.plot(Nsyns,PowerSubtracted[1], color = 'r', label = '5 Hz')
plt.plot(Nsyns,PowerSubtracted[2], color = 'g', label = '8 Hz')
plt.plot(Nsyns,PowerSubtracted[3], color = 'c', label = '10 Hz')
plt.plot(Nsyns,PowerSubtracted[4], color = 'm', label = '20 Hz')
plt.legend()
plt.xlabel('Number of IS3 Synapses',fontsize = font_size)
plt.ylabel('Power at Input Frequency',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAtInputFrequenciesSubtracted.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + '_PowerAtInputFrequenciesSubtracted.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()
