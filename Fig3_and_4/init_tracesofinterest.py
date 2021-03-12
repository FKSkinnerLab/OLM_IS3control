from __future__ import division
import numpy
import neuron
from neuron import h, gui
import matplotlib
import matplotlib.pyplot as plt
import efel
import glob
import IPython, os
import scipy
from scipy import signal
from scipy import stats
from currents_visualization import plotCurrentscape

font_size = 20

font = {'family' : 'normal',
		'weight' : 'normal',
		'size'   : font_size}

matplotlib.rc('font', **font)


### Instantiate Model ###
h.load_file("init_model.hoc")
h.load_file("SynParamSearch.hoc")

currsize = 'rheo'
cell = h.cell
h.ic_hold.amp = 0
findHolding = False # Finds the current needed (i.e. by increasing the current incrementally) to get the spike rate listed by TargetFreq
if cell == 1:
	cellname = 'Gormadoc'
	if currsize == '15Hz': currhold = 0.073 # 15 Hz
	elif currsize == '7.25Hz': currhold = 0.049 # 7.25 Hz
	elif currsize == '4Hz': currhold = 0.037 # ~4 Hz
	elif currsize == 'rheo': currhold = 0.030 # ~1-2 Hz (almost rheobase)
	else :
		currhold = 0.071
		findHolding = True # Finds the current needed (i.e. by increasing the current incrementally) to get the spike rate listed by TargetFreq
		if findHolding: TargetFreq = 15
elif cell == 2:
	cellname = 'Isembard'
	if currsize == '15Hz': currhold = 0.0715 # 15 Hz
	elif currsize == '7.25Hz': currhold = 0.045 # 7.25 Hz
	elif currsize == '4Hz': currhold = 0.031 # ~4 Hz
	elif currsize == 'rheo': currhold = 0.022 # ~1-2 Hz (almost rheobase)
	else :
		currhold = 0.0715
		findHolding = True # Finds the current needed (i.e. by increasing the current incrementally) to get the spike rate listed by TargetFreq
		if findHolding: TargetFreq = 15
else:
	print('No cell selected')

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

### Parameters ###
Nsyns = numpy.array([30])
freqs = numpy.array([1])
# starttimes = numpy.array([1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130])
# starttimes = numpy.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140])
starttimes = numpy.array([5,12])
base = 14

### Setup holding currents ###

h.cvode_active(1) # i.e. to make it run faster for this portion
h.cvode.atol(1e-05)
h.tstop = 5000

spikerate_prev = 0
if findHolding:
	for stepcurr in numpy.linspace(currhold*1000,currhold*1000+101,101):
		h.ic_hold.amp = stepcurr*0.001 # holding current
		h.f(0,0,h.tstop)
		apctimes = numpy.array(h.apctimes)
		apctimes = apctimes[apctimes > 1000] # Only spike times > 1s
		spikerate = len(apctimes)/4
		print('Spike Rate = ' + str(spikerate) + ' Hz at ' + str(stepcurr*0.001) + ' nA')
		if (spikerate > TargetFreq-0.2) & (spikerate < TargetFreq+0.2):
			print('Holding Current = ' + str(stepcurr*0.001) + ' nA')
			print('Baseline Spike Rate = ' + str(spikerate) + ' Hz')
			break
		elif spikerate > TargetFreq+0.2:
			h.ic_hold.amp = (stepcurr-1)*0.001 # holding current
			print('Spike Rate Passed Acceptable Range')
			print('Holding Current = ' + str((stepcurr-1)*0.001) + ' nA')
			print('Baseline Spike Rate = ' + str(spikerate_prev) + ' Hz')
			break
		spikerate_prev = spikerate
else:
	h.ic_hold.amp = currhold

h.cvode_active(0) # i.e turn cvode off for simulations
h.tstop = 5000
h.ic_hold.dur = h.tstop

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

fi, axarr = plt.subplots()

# First find first spike time after 1000 ms
h.f(Nsyns[0],freqs[0],h.tstop)
apctimes = numpy.array(h.apctimes)
spikes = apctimes[apctimes>1000]
spike0 = spikes[0]
ISIBase = spikes[1]-spikes[0]
pertimeres = ISIBase/base
starttimes = pertimeres*starttimes

# Plot base
TStart = int(spike0*10) - int(ISIBase*15)
TEnd = int(spike0*10) + int(ISIBase*25)
vvec = numpy.array(v_vec)
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

r0 = [pIKa[TStart:TEnd],pIKdrf[TStart:TEnd],pIKdrs[TStart:TEnd],pIKCa[TStart:TEnd],pIm[TStart:TEnd],pIl[TStart:TEnd],pINa[TStart:TEnd],pICaT[TStart:TEnd],pICaL[TStart:TEnd],pIh[TStart:TEnd]]

vvec_init = vvec

cols = ['orange','crimson']
startcount = 0
for starttime in starttimes:
	
	h.f(Nsyns[0],freqs[0],spike0+starttime) # Sets synaptic inputs
	
	TStart = int(spike0*10) - int(ISIBase*15)
	TEnd = int(spike0*10) + int(ISIBase*25)
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
	
	r0 = [pIKa[TStart:TEnd],pIKdrf[TStart:TEnd],pIKdrs[TStart:TEnd],pIKCa[TStart:TEnd],pIm[TStart:TEnd],pIl[TStart:TEnd],pINa[TStart:TEnd],pICaT[TStart:TEnd],pICaL[TStart:TEnd],pIh[TStart:TEnd]]
	
	vvec = numpy.array(v_vec)
	tvec = numpy.array(t_vec)
	apctimes = numpy.array(h.apctimes)
	
	axarr.plot(tvec, vvec, color=cols[startcount])
	# axarr.plot(numpy.array([spike0+starttime,spike0+starttime]), numpy.array([-85,30]), color=cols[startcount],linestyle=':')
	
	startcount = startcount + 1

axarr.plot(tvec, vvec_init, color='grey',linestyle='dashed')
axarr.set_xlim(spike0*1.0 - ISIBase*0.1,spike0*1.0 + ISIBase*1.5)
# axarr.set_ylim(-85,30)
# axarr.set_xlabel('Time (ms)',fontsize = font_size)
# axarr.spines['right'].set_visible(False)
# axarr.spines['top'].set_visible(False)
# if startcount != 14:
# 	axarr.spines['bottom'].set_visible(False)
# 	for tic in axarr.xaxis.get_major_ticks():
# 		tic.tick1On = tic.tick2On = False
fi.savefig('PLOTfiles/' + cellname + modelname + currsize + '_TracesOfInterest.pdf', bbox_inches='tight')
fi.savefig('PLOTfiles/' + cellname + modelname + currsize + '_TracesOfInterest.png', bbox_inches='tight')
# plt.gcf().clear()
# plt.cla()
fi.clf()
plt.close(fi)
