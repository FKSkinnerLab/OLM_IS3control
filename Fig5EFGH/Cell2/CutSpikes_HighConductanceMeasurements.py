### Test Script for a file found in the SDprox2 results from initial Parallel Simulations
import efel
import numpy
import time
import random
from scipy import signal
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

t = time.time()

data = numpy.zeros((100001,), dtype=numpy.float64)
tstop = h.tstop
dt = h.dt

h.randomize_syns(1,1)

def getMeasures(inhsyn,excsyn,inhspikes,excspikes,examplenum,repnum,modfreq,modfreqnum,currentstep,yind):	
	print 'Exc Num = ' + str(excsyn) + ', Inh Num = ' + str(inhsyn) + ', Exc Rate = ' + str(excspikes/10) + ', Inh Rate = ' + str(inhspikes/10)
	# rep = int(repnum)
	# example = int(examplenum)
	# modfreq = int(modfreq)
	# modfreqnum = int(modfreqnum)
	# currentstep = currhold + repnum*0.0025 - 0.0025
	# t_vec = h.Vector()
	# t_vec.record(h._ref_t)
	# v_vec = h.Vector()
	# v_vec.record(h.soma[0](0.5)._ref_v)
	
	saverep = repnum
	savemod = modfreqnum
	savemodfreq = modfreq
	savecurstep = currentstep
	# h.randomize_syns(1,1) # Randomizes synapse location with a different seed on each repetition
	h.f(0,0,0,0,0,1,1,0,0,0,30,0,modfreq,currentstep,int(yind)) # Runs Simulation
	h.run()
	# data1 = numpy.array(t_vec)
	# data = numpy.array(v_vec)
	# timevec = data1[10001:len(data1)]
	data = numpy.array(h.recV[int(yind)],dtype=numpy.float)
	timevec = numpy.arange(1000,int(tstop),dt,dtype=numpy.float)
	voltage = data[10001:len(data)] # Cut out first second of simulation since there are transient effects still present
	trace = {}
	trace['T'] = timevec
	trace['V'] = voltage
	trace['stim_start'] = [1000]
	trace['stim_end'] = [10000]
	traces = [trace]
	
	traces_results = efel.getFeatureValues(traces,['AP_begin_indices','AP_end_indices'])
	traces_mean_results = efel.getMeanFeatureValues(traces,['AP_amplitude','ISI_CV','Spikecount','mean_frequency'])
	SpikeRates_i3 = traces_mean_results[0]['mean_frequency']
	
	#### Create Trace Where Spikes Are Removed ####
	AP_begin = traces_results[0]['AP_begin_indices']
	AP_end = traces_results[0]['AP_end_indices']
	spikecut_voltage = []
    
	## Create spiketime 0/1 vector using acptimes and compute PSD using pwelch function
	apctimes = numpy.array(h.apctimes[int(yind)])
	apctimes = apctimes[apctimes>1000]
	NumSpikes_i2 = len(apctimes)
	if NumSpikes_i2 > 0:
		SpikeRates_i2 = NumSpikes_i2/((apctimes[-1:][0]-1000)/1000)
	else:
		SpikeRates_i2 = 0
	
	HC_SpikeTimes = numpy.zeros((len(voltage),), dtype=numpy.float)
	for i in range(0,len(apctimes)): HC_SpikeTimes[int(apctimes[i]/dt)-10001] = 1 # apctimes[i]
	f1, Pxx_den1 = signal.welch(HC_SpikeTimes, fs=1/(dt/1000), scaling='density', nperseg=20000) # fs/nperseg = frequency resolution
	e0Hz = Pxx_den1[f1==0]
	e0Hz = e0Hz[0]
	e05Hz = Pxx_den1[f1==0.5]
	e05Hz = e05Hz[0]
	e1Hz = Pxx_den1[f1==1]
	e1Hz = e1Hz[0]
	e2Hz = Pxx_den1[f1==2]
	e2Hz = e2Hz[0]
	e3Hz = Pxx_den1[f1==3]
	e3Hz = e3Hz[0]
	e4Hz = Pxx_den1[f1==4]
	e4Hz = e4Hz[0]
	e5Hz = Pxx_den1[f1==5]
	e5Hz = e5Hz[0]
	e8Hz = Pxx_den1[f1==8]
	e8Hz = e8Hz[0]
	e9Hz = Pxx_den1[f1==9]
	e9Hz = e9Hz[0]
	e10Hz = Pxx_den1[f1==10]
	e10Hz = e10Hz[0]
	e12Hz = Pxx_den1[f1==12]
	e12Hz = e12Hz[0]
	e15Hz = Pxx_den1[f1==15]
	e15Hz = e15Hz[0]
	e16Hz = Pxx_den1[f1==16]
	e16Hz = e16Hz[0]
	e20Hz = Pxx_den1[f1==20]
	e20Hz = e20Hz[0]
	e25Hz = Pxx_den1[f1==25]
	e25Hz = e25Hz[0]
	e30Hz = Pxx_den1[f1==30]
	e30Hz = e30Hz[0]
    # print e0Hz
    # print e05Hz
    # print e1Hz
    # print e2Hz
    # print e3Hz
    # print e4Hz
    # print e5Hz
    # print e8Hz
    # print e9Hz
    # print e10Hz
    # print e12Hz
    # print e15Hz
    # print e16Hz
    # print e20Hz
    # print e25Hz
    # print e30Hz
	
	if AP_begin is not None:
		for i in range(0,len(AP_begin)):
			# Cut out action potentials + 100 index points preceding each spike
			if i == 0:
				spikecut_tempvoltage = voltage[0:AP_begin[i]-80]
				spikecut_voltage.append(spikecut_tempvoltage)
			elif i == len(AP_begin):
				spikecut_tempvoltage = [voltage[AP_end[i-1]:AP_begin[i]]-80, voltage[AP_end[i]:len(voltage)]]
				spikecut_voltage.append(spikecut_tempvoltage)
			else:
				spikecut_tempvoltage = voltage[AP_end[i-1]:AP_begin[i]-80]
				spikecut_voltage.append(spikecut_tempvoltage)
			# Find lengths of appended arrays and rebuild voltage trace array
			x = []
			for i in range(0,len(spikecut_voltage)):
				newlength = len(spikecut_voltage[i])
				x.append(newlength)
			totallength = numpy.sum(x)
			spv = numpy.zeros((totallength,), dtype=numpy.int)
			count = 0
			for i in range(0,len(spikecut_voltage)):
				for j in range(0,len(spikecut_voltage[i])):
					spv[count] = spikecut_voltage[i][j]
					count = count + 1
	else:
		spv = voltage
	
	spv = spv[spv < -50] # Remove all voltage instinces greater than -50 mV
	spt = numpy.arange(0,len(spv),1)*0.1 # Build new tvec for trace with spike cut
	
	### Generate Measurements ###
	if len(spv) == 0:
		StdVolt_i = 0
		MeanVolt_i = -50 # i.e. set to highest possible average if all data points get cut
	else:
		StdVolt_i = numpy.std(spv)
		MeanVolt_i = numpy.mean(spv)
	NumSpikes_i = traces_mean_results[0]['Spikecount']
	SpikeRates_i = NumSpikes_i/9
	if traces_mean_results[0]['AP_amplitude'] is not None:
		MeanAPamp_i = traces_mean_results[0]['AP_amplitude']
	else:
		MeanAPamp_i = 0
	if traces_mean_results[0]['ISI_CV'] is not None:
		ISICV_i = traces_mean_results[0]['ISI_CV']
	else:
		ISICV_i = 0
		print 'FAILED AT I = ' + str(currentstep) + ' nA & Mod = ' + str(modfreq) + ' Hz'
		reload(efel) # added due to previous error following which all models had ISICV values of 0 regardless of spikes
	
	if modfreq == 0:
		fig = pyplot.figure()
		pyplot.plot(timevec,voltage)
		pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_Mod_' + str(h.synfreqINH) + 'Hz_NEURONSpikeRateEstimate_' + str(SpikeRates_i2) + 'Hz_I_' + str(currentstep) + 'nA_yind_' + str(yind) + '.png')
		pyplot.gcf().clear()
		pyplot.cla()
		pyplot.clf()
		pyplot.close()
	
	AvgPot_Thresh = -70.588
	StdPot_Thresh = 2.2
	ISICV_Thresh = 0.8
	AMP_DB_Thresh = 40
	Rate_maxthresh = 25 # Greater than ~3 Hz and less than ~25 Hz during resting (Varga et al, 2012 - though up to 50 Hz during runs; Katona et al, 2014)
	Rate_minthresh = 3 # Greater than ~3 Hz and less than ~25 Hz during resting (Varga et al, 2012 - though up to 50 Hz during runs; Katona et al, 2014)
	HC_Metric = 0
	HC_Metric = HC_Metric + (MeanVolt_i >= AvgPot_Thresh) + (StdVolt_i >= StdPot_Thresh) + (ISICV_i >= ISICV_Thresh) + ((SpikeRates_i <= Rate_maxthresh) & (SpikeRates_i >= Rate_minthresh)) - 5*((MeanAPamp_i <= AMP_DB_Thresh) & (NumSpikes_i > 0))
	print 'IRepNum = ' + str(repnum) + ', IRepNum_OG = ' + str(saverep) + ', I = ' + str(currentstep) + ' nA, I_OG = ' + str(savecurstep) + ' nA, I_NEURON = ' + str(h.ic_hold.amp) + ' nA, ModNum = ' + str(modfreqnum) + ', ModNum_OG = ' + str(savemod) + ', Mod = ' + str(modfreq) + ' Hz, Mod_OG = ' + str(savemodfreq) + ' Hz, Mod_NEURON = ' + str(h.synfreqINH) + ' Hz, Spike Rate = ' + str(SpikeRates_i) + ', Spike Rate 2 = ' + str(SpikeRates_i2) + ', Spike Rate 3 = ' + str(SpikeRates_i3)
	outputresults = [int(repnum), HC_Metric, modfreq, SpikeRates_i, int(modfreqnum), e0Hz,e05Hz,e1Hz,e2Hz,e3Hz,e4Hz,e5Hz,e8Hz,e9Hz,e10Hz,e12Hz,e15Hz,e16Hz,e20Hz,e25Hz,e30Hz,currentstep,SpikeRates_i2,ISICV_i,SpikeRates_i3]
	return outputresults

