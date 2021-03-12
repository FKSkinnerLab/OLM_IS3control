from __future__ import division
import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import efel
import glob
import IPython, os
import scipy
from scipy import signal
from scipy import stats
from currents_visualization import plotCurrentscape
font = {'family' : 'normal',
		'weight' : 'normal',
		'size'   : 20}

matplotlib.rc('font', **font)

font_size = 15

modelname = '1st'
cellnames = ['Gormadoc', 'Isembard']
Rates = ['rheo', '4Hz', '7.25Hz', '15Hz']
PerturbationPhase = (numpy.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])/14)*100

figPhase, axPhase = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIKA, axIKA = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIKdrf, axIKdrf = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIKdrs, axIKdrs = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIKCa, axIKCa = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIM, axIM = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIL, axIL = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figINa, axINa = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figICaT, axICaT = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figICaL, axICaL = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))
figIh, axIh = plt.subplots(1, 4, sharey=True, sharex=True,figsize=(15,6))

for i in range(0,len(axPhase)):
	axPhase[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIKA[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIKdrf[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIKdrs[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIKCa[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIM[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIL[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axINa[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axICaT[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axICaL[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)
	axIh[i].plot([0,100],[0,0], color = 'dimgray', ls = 'dotted', alpha = 0.5)

for cellname in cellnames:
	if cellname == 'Gormadoc':
		wls = 'solid'
		labelc = 'Cell 1'
	elif cellname == 'Isembard':
		wls = 'dashed'
		labelc = 'Cell 2'
	for j, currsize in enumerate(Rates):
		PhaseShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'PhaseShift.npy')
		IKAShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IKAShift.npy')
		IKdrfShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IKdrfShift.npy')
		IKdrsShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IKdrsShift.npy')
		IKCaShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IKCaShift.npy')
		IMShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IMShift.npy')
		ILShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'ILShift.npy')
		INaShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'INaShift.npy')
		ICaTShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'ICaTShift.npy')
		ICaLShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'ICaLShift.npy')
		IhShift = numpy.load('Output/' + cellname + modelname + currsize + '_' + 'IhShift.npy')
		
		axPhase[j].plot(PerturbationPhase,PhaseShift*100, color = 'k', ls = wls, label=labelc)
		axIKA[j].plot(PerturbationPhase,IKAShift*100, color = 'tab:blue', ls = wls, label=labelc)
		axIKdrf[j].plot(PerturbationPhase,IKdrfShift*100, color = 'tab:orange', ls = wls, label=labelc)
		axIKdrs[j].plot(PerturbationPhase,IKdrsShift*100, color = 'tab:green', ls = wls, label=labelc)
		axIKCa[j].plot(PerturbationPhase,IKCaShift*100, color = 'tab:red', ls = wls, label=labelc)
		axIM[j].plot(PerturbationPhase,IMShift*100, color = 'tab:purple', ls = wls, label=labelc)
		axIL[j].plot(PerturbationPhase,ILShift*100, color = 'tab:brown', ls = wls, label=labelc)
		axINa[j].plot(PerturbationPhase,INaShift*100, color = 'tab:pink', ls = wls, label=labelc)
		axICaT[j].plot(PerturbationPhase,ICaTShift*100, color = 'tab:gray', ls = wls, label=labelc)
		axICaL[j].plot(PerturbationPhase,ICaLShift*100, color = 'tab:olive', ls = wls, label=labelc)
		axIh[j].plot(PerturbationPhase,IhShift*100, color = 'tab:cyan', ls = wls, label=labelc)

axPhase[0].set_xlim(0,100)
axPhase[0].set_title('Rheobase')
axPhase[1].set_title('4 Hz')
axPhase[2].set_title('7.25 Hz')
axPhase[3].set_title('15 Hz')
axPhase[0].legend(loc='upper left')
figPhase.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_PhaseShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figPhase.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_PhaseShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figPhase)

axIKA[0].set_xlim(0,100)
# axIKA[0].set_title('Rheobase')
# axIKA[1].set_title('4 Hz')
# axIKA[2].set_title('7.25 Hz')
# axIKA[3].set_title('15 Hz')
axIKA[0].legend(loc='lower left')
figIKA.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKAShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIKA.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKAShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIKA)

axIKdrf[0].set_xlim(0,100)
# axIKdrf[0].set_title('Rheobase')
# axIKdrf[1].set_title('4 Hz')
# axIKdrf[2].set_title('7.25 Hz')
# axIKdrf[3].set_title('15 Hz')
axIKdrf[0].legend(loc='lower left')
figIKdrf.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKdrfShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIKdrf.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKdrfShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIKdrf)

axIKdrs[0].set_xlim(0,100)
# axIKdrs[0].set_title('Rheobase')
# axIKdrs[1].set_title('4 Hz')
# axIKdrs[2].set_title('7.25 Hz')
# axIKdrs[3].set_title('15 Hz')
axIKdrs[0].legend(loc='lower left')
figIKdrs.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKdrsShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIKdrs.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKdrsShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIKdrs)

axIKCa[0].set_xlim(0,100)
# axIKCa[0].set_title('Rheobase')
# axIKCa[1].set_title('4 Hz')
# axIKCa[2].set_title('7.25 Hz')
# axIKCa[3].set_title('15 Hz')
axIKCa[0].legend(loc='lower left')
figIKCa.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKCaShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIKCa.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IKCaShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIKCa)

axIM[0].set_xlim(0,100)
# axIM[0].set_title('Rheobase')
# axIM[1].set_title('4 Hz')
# axIM[2].set_title('7.25 Hz')
# axIM[3].set_title('15 Hz')
axIM[0].legend(loc='lower left')
figIM.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IMShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIM.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IMShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIM)

axIL[0].set_xlim(0,100)
# axIL[0].set_title('Rheobase')
# axIL[1].set_title('4 Hz')
# axIL[2].set_title('7.25 Hz')
# axIL[3].set_title('15 Hz')
axIL[0].legend(loc='lower left')
figIL.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ILShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIL.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ILShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIL)

axINa[0].set_xlim(0,100)
# axINa[0].set_title('Rheobase')
# axINa[1].set_title('4 Hz')
# axINa[2].set_title('7.25 Hz')
# axINa[3].set_title('15 Hz')
axINa[0].legend(loc='lower left')
figINa.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_INaShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figINa.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_INaShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figINa)

axICaT[0].set_xlim(0,100)
# axICaT[0].set_title('Rheobase')
# axICaT[1].set_title('4 Hz')
# axICaT[2].set_title('7.25 Hz')
# axICaT[3].set_title('15 Hz')
axICaT[0].legend(loc='lower left')
figICaT.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ICaTShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figICaT.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ICaTShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figICaT)

axICaL[0].set_xlim(0,100)
# axICaL[0].set_title('Rheobase')
# axICaL[1].set_title('4 Hz')
# axICaL[2].set_title('7.25 Hz')
# axICaL[3].set_title('15 Hz')
axICaL[0].legend(loc='lower left')
figICaL.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ICaLShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figICaL.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_ICaLShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figICaL)

axIh[0].set_xlim(0,100)
# axIh[0].set_title('Rheobase')
# axIh[1].set_title('4 Hz')
# axIh[2].set_title('7.25 Hz')
# axIh[3].set_title('15 Hz')
axIh[0].legend(loc='upper left')
figIh.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IhShiftAll.pdf', bbox_inches='tight', dpi = 300, transparent = True)
figIh.savefig('PLOTfiles_analysis/' + cellname + modelname + currsize + '_IhShiftAll.png', bbox_inches='tight', dpi = 300, transparent = True)
plt.close(figIh)
