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
starttimes = numpy.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])
PhaseShift = numpy.zeros(len(starttimes))
IhShift = numpy.zeros(len(starttimes))
ICaTShift = numpy.zeros(len(starttimes))
IKAShift = numpy.zeros(len(starttimes))
IKdrfShift = numpy.zeros(len(starttimes))
IKdrsShift = numpy.zeros(len(starttimes))
IKCaShift = numpy.zeros(len(starttimes))
IMShift = numpy.zeros(len(starttimes))
ILShift = numpy.zeros(len(starttimes))
INaShift = numpy.zeros(len(starttimes))
ICaLShift = numpy.zeros(len(starttimes))

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

fi, axarr = plt.subplots(15, sharex = True, sharey = True)

# First find first spike time after 1000 ms
h.f(Nsyns[0],freqs[0],h.tstop)
apctimes = numpy.array(h.apctimes)
spikes = apctimes[apctimes>1000]
spike0 = spikes[0]
ISIBase = spikes[1]-spikes[0]
pertimeres = ISIBase/starttimes[-1:]
starttimes = pertimeres*starttimes

# Plot base
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

fig0 = plotCurrentscape(vvec0[TStart:TEnd], r0)
# fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_StartTime_' + str(spike0+starttime) + currsize + '.pdf',dpi=500)
fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_StartTime_' + currsize + '_' + str(0) + '.png',dpi=500)
fig0.clf()
plt.close(fig0)

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
	
	fig0 = plotCurrentscape(vvec0[TStart:TEnd], r0)
	# fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_StartTime_' + str(spike0+starttime) + currsize + '.pdf',dpi=500)
	fig0.savefig('PLOTfiles/' + cellname + '_Currentscape_' + modelname + '_StartTime_' + currsize + '_' + str(spike0+starttime) + '.png',dpi=500)
	fig0.clf()
	plt.close(fig0)
	
	vvec = numpy.array(v_vec)
	tvec = numpy.array(t_vec)
	apctimes = numpy.array(h.apctimes)
	
	# Compute Phase Shift Caused by Perturbation
	apctimes1 = apctimes[apctimes<=(spike0+starttime)]
	apctimes2 = apctimes[apctimes>(spike0+starttime)]
	ISIVal0 = apctimes1[-1:]-apctimes1[-2:-1] # Find last ISI before perturbation
	ISIVal1 = apctimes2[0]-apctimes1[-1:] # Find ISI caused by perturbation
	PhaseShift[startcount] = (ISIVal1 - ISIVal0)/ISIVal0 # Normalized to 1
	# NOTE: Positive phase shift values means a shortened ISI and is called a phase advance
	
	# Compute Change in Ih Caused by Perturbation
	pIhMAX0 = numpy.amax(abs(pIh[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max Ih current between 2nd last spike before perturbation and perturbation
	pIhMAX1 = numpy.amax(abs(pIh[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max Ih current between perturbation and 2nd spike following
	IhShift[startcount] = (pIhMAX1 - pIhMAX0)/pIhMAX0 # Normalized to 1
	
	# Compute Change in ICaT Caused by Perturbation
	pICaTMAX0 = numpy.amax(abs(pICaT[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pICaTMAX1 = numpy.amax(abs(pICaT[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	ICaTShift[startcount] = (pICaTMAX1 - pICaTMAX0)/pICaTMAX0 # Normalized to 1
	# Compute Change in ICaT Caused by Perturbation
	
	pICaLMAX0 = numpy.amax(abs(pICaL[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pICaLMAX1 = numpy.amax(abs(pICaL[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	ICaLShift[startcount] = (pICaLMAX1 - pICaLMAX0)/pICaLMAX0 # Normalized to 1
	
	pIKCaMAX0 = numpy.amax(abs(pIKCa[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pIKCaMAX1 = numpy.amax(abs(pIKCa[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	IKCaShift[startcount] = (pIKCaMAX1 - pIKCaMAX0)/pIKCaMAX0 # Normalized to 1
	
	pIKAMAX0 = numpy.amax(abs(pIKa[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pIKAMAX1 = numpy.amax(abs(pIKa[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	IKAShift[startcount] = (pIKAMAX1 - pIKAMAX0)/pIKAMAX0 # Normalized to 1
	
	pIKdrfMAX0 = numpy.amax(abs(pIKdrf[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pIKdrfMAX1 = numpy.amax(abs(pIKdrf[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	IKdrfShift[startcount] = (pIKdrfMAX1 - pIKdrfMAX0)/pIKdrfMAX0 # Normalized to 1
	
	pIKdrsMAX0 = numpy.amax(abs(pIKdrs[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pIKdrsMAX1 = numpy.amax(abs(pIKdrs[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	IKdrsShift[startcount] = (pIKdrsMAX1 - pIKdrsMAX0)/pIKdrsMAX0 # Normalized to 1
	
	pIMMAX0 = numpy.amax(abs(pIm[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pIMMAX1 = numpy.amax(abs(pIm[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	IMShift[startcount] = (pIMMAX1 - pIMMAX0)/pIMMAX0 # Normalized to 1
	
	pILMAX0 = numpy.amax(abs(pIl[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pILMAX1 = numpy.amax(abs(pIl[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	ILShift[startcount] = (pILMAX1 - pILMAX0)/pILMAX0 # Normalized to 1
	
	pINaMAX0 = numpy.amax(abs(pINa[(tvec>apctimes1[-2:-1]) & (tvec<(spike0+starttime))])) # Find max ICaT current between 2nd last spike before perturbation and perturbation
	pINaMAX1 = numpy.amax(abs(pINa[(tvec>(spike0+starttime)) & (tvec<apctimes2[1])])) # Find max ICaT current between perturbation and 2nd spike following
	INaShift[startcount] = (pINaMAX1 - pINaMAX0)/pINaMAX0 # Normalized to 1
	
	axarr[startcount].plot(tvec, vvec, color='k')
	axarr[startcount].plot(numpy.array([spike0+starttime,spike0+starttime]), numpy.array([-85,30]), color='k',linestyle=':')
	axarr[startcount].set_xlim(spike0*1.0 - ISIBase*1.5,spike0*1.0 + ISIBase*2.5)
	axarr[startcount].set_ylim(-85,30)
	axarr[startcount].spines['right'].set_visible(False)
	axarr[startcount].spines['top'].set_visible(False)
	if startcount != 14:
		axarr[startcount].spines['bottom'].set_visible(False)
		for tic in axarr[startcount].xaxis.get_major_ticks():
			tic.tick1On = tic.tick2On = False
	
	startcount = startcount + 1

axarr[14].set_xlabel('Time (ms)',fontsize = font_size)
fi.savefig('PLOTfiles/' + cellname + modelname + currsize + '_Traces.pdf', bbox_inches='tight')
fi.savefig('PLOTfiles/' + cellname + modelname + currsize + '_Traces.png', bbox_inches='tight')
# plt.gcf().clear()
# plt.cla()
fi.clf()
plt.close(fi)


# Plot PRCs
PerturbationPhase = (starttimes/ISIVal0)*100 # Normalize to the most recently indexed ISI value preceding a perturbation
plt.plot(PerturbationPhase,PhaseShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Phase Shift',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_PhaseShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_PhaseShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IhShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max Ih Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IhShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IhShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,ICaTShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max ICaT Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ICaTShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ICaTShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,ICaLShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max ICaL Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ICaLShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ICaLShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IKCaShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max IKCa Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IKCaShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IKCaShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IKAShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max IA Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IAShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IAShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IKdrfShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max Ikdrf Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IkdrfShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IkdrfShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IKdrsShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max Ikdrs Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IkdrsShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IkdrsShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,IMShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max IM Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IMShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_IMShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,ILShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max IL Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ILShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_ILShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

plt.plot(PerturbationPhase,INaShift*100, color = 'k')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max INa Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_INaShift.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_INaShift.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

# Choos similar colors to currentscapes
plt.plot(PerturbationPhase,IKAShift*100, color = 'tab:blue', label = 'IKa')
plt.plot(PerturbationPhase,IKdrfShift*100, color = 'tab:orange', label = 'IKdrf')
plt.plot(PerturbationPhase,IKdrsShift*100, color = 'tab:green', label = 'IKdrs')
plt.plot(PerturbationPhase,IKCaShift*100, color = 'tab:red', label = 'IKCa')
plt.plot(PerturbationPhase,IMShift*100, color = 'tab:purple', label = 'IM')
plt.plot(PerturbationPhase,ILShift*100, color = 'tab:brown', label = 'IL')
plt.plot(PerturbationPhase,INaShift*100, color = 'tab:pink', label = 'INa')
plt.plot(PerturbationPhase,ICaTShift*100, color = 'tab:gray', label = 'ICaT')
plt.plot(PerturbationPhase,ICaLShift*100, color = 'tab:olive', label = 'ICaL')
plt.plot(PerturbationPhase,IhShift*100, color = 'tab:cyan', label = 'IH')
plt.legend(loc = 'upper left')
plt.xlabel('Perturbation % Phase',fontsize = font_size)
plt.ylabel('% Max Current Change',fontsize = font_size)
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_AllIShifts.pdf', bbox_inches='tight')
plt.savefig('PLOTfiles/' + cellname + modelname + currsize + '_AllIShifts.png', bbox_inches='tight')
plt.gcf().clear()
plt.cla()
plt.clf()
plt.close()

numpy.save('Output/' + cellname + modelname + currsize + '_' + 'PhaseShift.npy',PhaseShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IKAShift.npy',IKAShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IKdrfShift.npy',IKdrfShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IKdrsShift.npy',IKdrsShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IKCaShift.npy',IKCaShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IMShift.npy',IMShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'ILShift.npy',ILShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'INaShift.npy',INaShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'ICaTShift.npy',ICaTShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'ICaLShift.npy',ICaLShift)
numpy.save('Output/' + cellname + modelname + currsize + '_' + 'IhShift.npy',IhShift)
