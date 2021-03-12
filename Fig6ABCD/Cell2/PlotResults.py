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
font_size = 13
numrerands = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'numrerands.npy')

Examples = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelParams.npy')
ExampleStrings = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelStrings.npy')
ExampleString = ExampleStrings[0]

modfreqs = numpy.array([0,0.5,1,2,3,4,5,8,9,10,12,15,16,20,25,30])

numinh = Examples[0][0]
numexc = Examples[1][0]
inhspikes = Examples[2][0]
excspikes = Examples[3][0]

HCT = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'IVLMat.npy')
MODFREQ = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'MODFREQMat.npy')
RATES = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'RATESMat.npy')
PSDXX = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + 'PSDXXMat.npy')

baserates = numpy.zeros((numrerands,), dtype=numpy.float)
basehcts = numpy.zeros((numrerands,), dtype=numpy.float)
for y in range(0,numrerands):
	baserates[y] = RATES[y][0]
	basehcts[y] = HCT[y][0]

minrate = numpy.amin(baserates)
maxrate = numpy.amax(baserates)
percenthct = (sum(basehcts==4)/numrerands)*100

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l] # Find power ratio for each modulation frequency
	if BaseIVLval == 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='k')
	elif BaseIVLval < 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='r')
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(numrerands-1,-1,-1):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l] # Find power ratio for each modulation frequency
	
	ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

ax0 = pyplot.subplot2grid((10, 10), (0, 0), colspan=9)
ax1 = pyplot.subplot2grid((10, 10), (1, 9), rowspan=9)
ax2 = pyplot.subplot2grid((10, 10), (1, 0), rowspan=9, colspan=9)
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))
cmap = pl.cm.plasma
norm = matplotlib.colors.Normalize(vmin=numpy.amin(baserates), vmax=numpy.amax(baserates))
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label(r'$f_B$')

PeakFreqs = numpy.zeros((numrerands,), dtype=numpy.float)

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l] # Find power ratio for each modulation frequency
	
	ax2.plot(modfreqnums, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l] # Find power ratio for each modulation frequency
	
	PeakFreqs[y] = modfreqnums[numpy.argmax(PowerRatios[1:None])]
	num = numpy.random.rand(1)*.6-0.3
	ax2.scatter(modfreqnums[numpy.argmax(PowerRatios[1:None])]+num,maxpratio,color=colors[y], alpha=0.5)

ax0.hist(PeakFreqs,bins=numpy.arange(18)+1.5)
ax0.set_xlim(0,16)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.get_xaxis().set_ticks([])
ax0.yaxis.tick_right()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlim(0,16)
ax2.set_xticks(modfreqnums)
ax2.set_xticklabels(modfreqstrs, rotation = 45)
ax2.set_xlabel(r'$f_i$')
ax2.set_ylim(0,maxpratio)
ax2.set_ylabel('Power')
pyplot.axis('tight')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Power3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

avgpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecB = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]
		powvecB[y] = PSDXX[y][0][l]
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	avgpowB[l] = numpy.median(powvecB)
	stdpowB[l] = numpy.std(powvecB)
	upper_quartileB[l] = numpy.percentile(powvecB, 75)
	lower_quartileB[l] = numpy.percentile(powvecB, 25)
	minpowB[l] = numpy.amin(powvecB)
	maxpowB[l] = numpy.amax(powvecB)
	ttest = stats.ttest_rel(powvec, powvecB) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,maxpow[l]+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,maxpow[l]+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,maxpow[l]+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums,avgpow[1:None],color = 'tab:orange')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileB[v], upper_quartileB[v],color = 'k', linewidth=3)
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowB[v], maxpowB[v],color = 'k', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowB[1:None],color = 'k',label='Base', alpha=0.8)

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecB = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecNIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvecB[y] = PSDXX[y][0][l] # BASELINE POWER
		powvec[y] = PSDXX[y][l][l]
		if HCT[y][0] == 4:
			powvecIVL[y] = PSDXX[y][l][l]
		elif HCT[y][0] < 4:
			powvecNIVL[y] = PSDXX[y][l][l]
	
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	avgpowB[l] = numpy.median(powvecB)
	stdpowB[l] = numpy.std(powvecB)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	upper_quartileB[l] = numpy.percentile(powvecB, 75)
	lower_quartileB[l] = numpy.percentile(powvecB, 25)
	minpowB[l] = numpy.amin(powvecB)
	maxpowB[l] = numpy.amax(powvecB)
	if (numpy.sum(powvecIVL!=0) > 0):
		avgpowIVL[l] = numpy.median(powvecIVL[powvecIVL!=0])
		stdpowIVL[l] = numpy.std(powvecIVL[powvecIVL!=0])
		upper_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 75)
		lower_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 25)
		minpowIVL[l] = numpy.amin(powvecIVL[powvecIVL!=0])
		maxpowIVL[l] = numpy.amax(powvecIVL[powvecIVL!=0])
	if (numpy.sum(powvecNIVL!=0) > 0):
		avgpowNIVL[l] = numpy.median(powvecNIVL[powvecNIVL!=0])
		stdpowNIVL[l] = numpy.std(powvecNIVL[powvecNIVL!=0])
		upper_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 75)
		lower_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 25)
		minpowNIVL[l] = numpy.amin(powvecNIVL[powvecNIVL!=0])
		maxpowNIVL[l] = numpy.amax(powvecNIVL[powvecNIVL!=0])

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums+num,avgpow[1:None],color = 'tab:orange',label='All', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileB[v], upper_quartileB[v],color = 'k', linewidth=3)
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowB[v], maxpowB[v],color = 'k', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowB[1:None],color = 'k',label='Base', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileIVL[v], upper_quartileIVL[v],color = 'tab:red', linewidth=3)
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowIVL[v], maxpowIVL[v],color = 'tab:red', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowIVL[1:None],color = 'tab:red',label='IVL', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileNIVL[v], upper_quartileNIVL[v],color = 'tab:blue', linewidth=3)
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowNIVL[v], maxpowNIVL[v],color = 'tab:blue', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowNIVL[1:None],color = 'tab:blue',label='NIVL', alpha=0.8)

# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowB[1:None],stdpowB[1:None],color = 'k',ecolor='k',label='Base', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvecLow = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecHig = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		if baserates[y] < numpy.median(baserates):
			powvecLow[y] = PSDXX[y][l][l]
		elif baserates[y] >= numpy.median(baserates):
			powvecHig[y] = PSDXX[y][l][l]
	
	if (numpy.sum(powvecLow!=0) > 0):
		avgpowLow[l] = numpy.median(powvecLow[powvecLow!=0])
		stdpowLow[l] = numpy.std(powvecLow[powvecLow!=0])
		upper_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 75)
		lower_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 25)
		minpowLow[l] = numpy.amin(powvecLow[powvecLow!=0])
		maxpowLow[l] = numpy.amax(powvecLow[powvecLow!=0])
	if (numpy.sum(powvecHig!=0) > 0):
		avgpowHig[l] = numpy.median(powvecHig[powvecHig!=0])
		stdpowHig[l] = numpy.std(powvecHig[powvecHig!=0])
		upper_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 75)
		lower_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 25)
		minpowHig[l] = numpy.amin(powvecHig[powvecHig!=0])
		maxpowHig[l] = numpy.amax(powvecHig[powvecHig!=0])
	
	ttest = stats.ttest_ind(powvecLow, powvecHig) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileLow)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileLow[v], upper_quartileLow[v],color = 'tab:green', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowLow[v], maxpowLow[v],color = 'tab:green', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowLow[1:None],color = 'tab:green',label='Low ' + r'$f_B$', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileHig[v], upper_quartileHig[v],color = 'tab:purple', linewidth=3)
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowHig[v], maxpowHig[v],color = 'tab:purple', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowHig[1:None],color = 'tab:purple',label='High ' + r'$f_B$', alpha=0.8)

# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowB[1:None],stdpowB[1:None],color = 'k',ecolor='k',label='Base', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz; Med = ' + str(numpy.around(numpy.median(baserates),decimals=1)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerAvg3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


# Power Ratios
fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][l][0] # Find power ratio for each modulation frequency
	if BaseIVLval == 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='k')
	elif BaseIVLval < 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='r')
	maxpratio = numpy.amax([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(numrerands-1,-1,-1):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][l][0] # Find power ratio for each modulation frequency
	
	ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

ax0 = pyplot.subplot2grid((10, 10), (0, 0), colspan=9)
ax1 = pyplot.subplot2grid((10, 10), (1, 9), rowspan=9)
ax2 = pyplot.subplot2grid((10, 10), (1, 0), rowspan=9, colspan=9)
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))
cmap = pl.cm.plasma
norm = matplotlib.colors.Normalize(vmin=numpy.amin(baserates), vmax=numpy.amax(baserates))
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label(r'$f_B$')

PeakFreqs = numpy.zeros((numrerands,), dtype=numpy.float)

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][l][0] # Find power ratio for each modulation frequency
	
	ax2.plot(modfreqnums, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][l][0] # Find power ratio for each modulation frequency
	
	PeakFreqs[y] = modfreqnums[numpy.argmax(PowerRatios[1:None])]
	num = numpy.random.rand(1)*.6-0.3
	ax2.scatter(modfreqnums[numpy.argmax(PowerRatios[1:None])]+num,maxpratio,color=colors[y], alpha=0.5)

ax0.hist(PeakFreqs,bins=numpy.arange(18)+1.5)
ax0.set_xlim(0,16)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.get_xaxis().set_ticks([])
ax0.yaxis.tick_right()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlim(0,16)
ax2.set_xticks(modfreqnums)
ax2.set_xticklabels(modfreqstrs, rotation = 45)
ax2.set_xlabel(r'$f_i$')
ax2.set_ylim(0,maxpratio)
ax2.set_ylabel('Power Ratio')
pyplot.axis('tight')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatios3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

avgpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecB = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
		powvecB[y] = PSDXX[y][0][l]/PSDXX[y][0][0]
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	avgpowB[l] = numpy.median(powvecB)
	stdpowB[l] = numpy.std(powvecB)
	upper_quartileB[l] = numpy.percentile(powvecB, 75)
	lower_quartileB[l] = numpy.percentile(powvecB, 25)
	minpowB[l] = numpy.amin(powvecB)
	maxpowB[l] = numpy.amax(powvecB)
	ttest = stats.ttest_rel(powvec, powvecB) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,maxpow[l]+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,maxpow[l]+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,maxpow[l]+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums,avgpow[1:None],color = 'tab:orange')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileB[v], upper_quartileB[v],color = 'k', linewidth=3)
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowB[v], maxpowB[v],color = 'k', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowB[1:None],color = 'k',label='Base', alpha=0.8)

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecB = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecNIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvecB[y] = PSDXX[y][0][l]/PSDXX[y][0][0] # BASELINE POWER RATIO
		powvec[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
		if HCT[y][0] == 4:
			powvecIVL[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
		elif HCT[y][0] < 4:
			powvecNIVL[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
	
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	avgpowB[l] = numpy.median(powvecB)
	stdpowB[l] = numpy.std(powvecB)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	upper_quartileB[l] = numpy.percentile(powvecB, 75)
	lower_quartileB[l] = numpy.percentile(powvecB, 25)
	minpowB[l] = numpy.amin(powvecB)
	maxpowB[l] = numpy.amax(powvecB)
	if (numpy.sum(powvecIVL!=0) > 0):
		avgpowIVL[l] = numpy.median(powvecIVL[powvecIVL!=0])
		stdpowIVL[l] = numpy.std(powvecIVL[powvecIVL!=0])
		upper_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 75)
		lower_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 25)
		minpowIVL[l] = numpy.amin(powvecIVL[powvecIVL!=0])
		maxpowIVL[l] = numpy.amax(powvecIVL[powvecIVL!=0])
	if (numpy.sum(powvecNIVL!=0) > 0):
		avgpowNIVL[l] = numpy.median(powvecNIVL[powvecNIVL!=0])
		stdpowNIVL[l] = numpy.std(powvecNIVL[powvecNIVL!=0])
		upper_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 75)
		lower_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 25)
		minpowNIVL[l] = numpy.amin(powvecNIVL[powvecNIVL!=0])
		maxpowNIVL[l] = numpy.amax(powvecNIVL[powvecNIVL!=0])

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums+num,avgpow[1:None],color = 'tab:orange',label='All', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileB[v], upper_quartileB[v],color = 'k', linewidth=3)
for v in range(1,len(lower_quartileB)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowB[v], maxpowB[v],color = 'k', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowB[1:None],color = 'k',label='Base', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileIVL[v], upper_quartileIVL[v],color = 'tab:red', linewidth=3)
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowIVL[v], maxpowIVL[v],color = 'tab:red', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowIVL[1:None],color = 'tab:red',label='IVL', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileNIVL[v], upper_quartileNIVL[v],color = 'tab:blue', linewidth=3)
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowNIVL[v], maxpowNIVL[v],color = 'tab:blue', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowNIVL[1:None],color = 'tab:blue',label='NIVL', alpha=0.8)

# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowB[1:None],stdpowB[1:None],color = 'k',ecolor='k',label='Base', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvecLow = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecHig = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		if baserates[y] < numpy.median(baserates):
			powvecLow[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
		elif baserates[y] >= numpy.median(baserates):
			powvecHig[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
	
	if (numpy.sum(powvecLow!=0) > 0):
		avgpowLow[l] = numpy.median(powvecLow[powvecLow!=0])
		stdpowLow[l] = numpy.std(powvecLow[powvecLow!=0])
		upper_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 75)
		lower_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 25)
		minpowLow[l] = numpy.amin(powvecLow[powvecLow!=0])
		maxpowLow[l] = numpy.amax(powvecLow[powvecLow!=0])
	if (numpy.sum(powvecHig!=0) > 0):
		avgpowHig[l] = numpy.median(powvecHig[powvecHig!=0])
		stdpowHig[l] = numpy.std(powvecHig[powvecHig!=0])
		upper_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 75)
		lower_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 25)
		minpowHig[l] = numpy.amin(powvecHig[powvecHig!=0])
		maxpowHig[l] = numpy.amax(powvecHig[powvecHig!=0])
	
	ttest = stats.ttest_ind(powvecLow, powvecHig) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileLow)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileLow[v], upper_quartileLow[v],color = 'tab:green', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowLow[v], maxpowLow[v],color = 'tab:green', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowLow[1:None],color = 'tab:green',label='Low ' + r'$f_B$', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileHig[v], upper_quartileHig[v],color = 'tab:purple', linewidth=3)
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowHig[v], maxpowHig[v],color = 'tab:purple', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowHig[1:None],color = 'tab:purple',label='High ' + r'$f_B$', alpha=0.8)

# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowB[1:None],stdpowB[1:None],color = 'k',ecolor='k',label='Base', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz; Med = ' + str(numpy.around(numpy.median(baserates),decimals=1)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosAvg3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


# Baseline Ratios
fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')

modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][0][l] # Find power ratio for each modulation frequency
	if BaseIVLval == 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='k')
	elif BaseIVLval < 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='r')
	maxpratio = numpy.amax([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Baseline Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(numrerands-1,-1,-1):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Baseline Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

ax0 = pyplot.subplot2grid((10, 10), (0, 0), colspan=9)
ax1 = pyplot.subplot2grid((10, 10), (1, 9), rowspan=9)
ax2 = pyplot.subplot2grid((10, 10), (1, 0), rowspan=9, colspan=9)
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))
cmap = pl.cm.plasma
norm = matplotlib.colors.Normalize(vmin=numpy.amin(baserates), vmax=numpy.amax(baserates))
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label(r'$f_B$')

PeakFreqs = numpy.zeros((numrerands,), dtype=numpy.float)
PeakFreqs2 = numpy.zeros((numrerands,), dtype=numpy.float)
PeakFreqRates = numpy.zeros((numrerands,), dtype=numpy.float)
ax2.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	ax2.plot(modfreqnums, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	PeakFreqs[y] = modfreqnums[numpy.argmax(PowerRatios[1:None])]
	modfreqs2 = modfreqs[1:]
	PeakFreqs2[y] = modfreqs2[numpy.argmax(PowerRatios[1:None])]
	PeakFreqRates[y] = RATES[y][numpy.argmax(PowerRatios[1:None])]
	num = numpy.random.rand(1)*.6-0.3
	ax2.scatter(modfreqnums[numpy.argmax(PowerRatios[1:None])]+num,maxpratio,color=colors[y], alpha=0.5)

ax0.hist(PeakFreqs,bins=numpy.arange(18)+1.5)
ax0.set_xlim(0,16)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.get_xaxis().set_ticks([])
ax0.yaxis.tick_right()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlim(0,16)
ax2.set_xticks(modfreqnums)
ax2.set_xticklabels(modfreqstrs, rotation = 45)
ax2.set_xlabel(r'$f_i$')
ax2.set_ylim(0,maxpratio)
ax2.set_ylabel('Baseline Ratio')
ax2.autoscale(enable=True, axis='x', tight=True)
pyplot.axis('tight')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatio3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

pyplot.fill_between(numpy.array([0,37]), numpy.array([0,37]), 37, color='b',alpha=0.2)
pyplot.fill_between(numpy.array([0,37]), 0, numpy.array([0,37]), color='r',alpha=0.2)
pyplot.scatter(baserates,PeakFreqs2,color='k')
pyplot.plot(numpy.array([0,37]),numpy.array([0,37]),ls='--',color='k')
pyplot.xlim(0,37)
pyplot.ylim(0,37)
pyplot.xlabel(r'$f_B$ (Hz)')
pyplot.ylabel(r'$f_R$ (Hz)')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Mod_fbVSfr.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Mod_fbVSfr.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

pyplot.fill_between(numpy.array([0,37]), numpy.array([0,37]), 37, color='b',alpha=0.2)
pyplot.fill_between(numpy.array([0,37]), 0, numpy.array([0,37]), color='r',alpha=0.2)
pyplot.scatter(baserates,PeakFreqRates,color='k')
pyplot.plot(numpy.array([0,37]),numpy.array([0,37]),ls='--',color='k')
pyplot.xlim(0,37)
pyplot.ylim(0,37)
pyplot.xlabel(r'$f_B$ (Hz)')
pyplot.ylabel(r'$f_R$ Spike Rate (Hz)')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Mod_fbVSfrSpikeRate.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3Mod_fbVSfrSpikeRate.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

meanFb = numpy.mean(baserates)
meanFr = numpy.mean(PeakFreqs2)
meanFr_rate = numpy.mean(PeakFreqRates)
SDFb = numpy.std(baserates)
SDFr = numpy.std(PeakFreqs2)
SDFr_rate = numpy.std(PeakFreqRates)
print('fb = '+str(meanFb)+u"\u00B1"+str(SDFb)+' Hz')
print('fb range = '+str(numpy.min(baserates))+'-'+str(numpy.max(baserates))+' Hz')
print('fr = '+str(meanFr)+u"\u00B1"+str(SDFr)+' Hz')
print('Rate at fr = '+str(meanFr_rate)+u"\u00B1"+str(SDFr_rate)+' Hz')

fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	ttest = stats.ttest_1samp(powvec, 1)
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,maxpow[l]+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,maxpow[l]+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,maxpow[l]+0.02,'*')

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums,avgpow[1:None],color = 'tab:orange')
# pyplot.errorbar(modfreqnums,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',elinewidth=3)
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecNIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
		if HCT[y][0] == 4:
			powvecIVL[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
		elif HCT[y][0] < 4:
			powvecNIVL[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
	
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	if (numpy.sum(powvecIVL!=0) > 0):
		avgpowIVL[l] = numpy.median(powvecIVL[powvecIVL!=0])
		stdpowIVL[l] = numpy.std(powvecIVL[powvecIVL!=0])
		upper_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 75)
		lower_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 25)
		minpowIVL[l] = numpy.amin(powvecIVL[powvecIVL!=0])
		maxpowIVL[l] = numpy.amax(powvecIVL[powvecIVL!=0])
	if (numpy.sum(powvecNIVL!=0) > 0):
		avgpowNIVL[l] = numpy.median(powvecNIVL[powvecNIVL!=0])
		stdpowNIVL[l] = numpy.std(powvecNIVL[powvecNIVL!=0])
		upper_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 75)
		lower_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 25)
		minpowNIVL[l] = numpy.amin(powvecNIVL[powvecNIVL!=0])
		maxpowNIVL[l] = numpy.amax(powvecNIVL[powvecNIVL!=0])

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums+num,avgpow[1:None],color = 'tab:orange',label='All', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileIVL[v], upper_quartileIVL[v],color = 'tab:red', linewidth=3)
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowIVL[v], maxpowIVL[v],color = 'tab:red', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowIVL[1:None],color = 'tab:red',label='IVL', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileNIVL[v], upper_quartileNIVL[v],color = 'tab:blue', linewidth=3)
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowNIVL[v], maxpowNIVL[v],color = 'tab:blue', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowNIVL[1:None],color = 'tab:blue',label='NIVL', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvecLow = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecHig = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		if baserates[y] < numpy.median(baserates):
			powvecLow[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
		elif baserates[y] >= numpy.median(baserates):
			powvecHig[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
	
	if (numpy.sum(powvecLow!=0) > 0):
		avgpowLow[l] = numpy.median(powvecLow[powvecLow!=0])
		stdpowLow[l] = numpy.std(powvecLow[powvecLow!=0])
		upper_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 75)
		lower_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 25)
		minpowLow[l] = numpy.amin(powvecLow[powvecLow!=0])
		maxpowLow[l] = numpy.amax(powvecLow[powvecLow!=0])
	if (numpy.sum(powvecHig!=0) > 0):
		avgpowHig[l] = numpy.median(powvecHig[powvecHig!=0])
		stdpowHig[l] = numpy.std(powvecHig[powvecHig!=0])
		upper_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 75)
		lower_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 25)
		minpowHig[l] = numpy.amin(powvecHig[powvecHig!=0])
		maxpowHig[l] = numpy.amax(powvecHig[powvecHig!=0])
	
	ttest = stats.ttest_ind(powvecLow, powvecHig) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileLow)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileLow[v], upper_quartileLow[v],color = 'tab:green', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowLow[v], maxpowLow[v],color = 'tab:green', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowLow[1:None],color = 'tab:green',label='Low ' + r'$f_B$', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileHig[v], upper_quartileHig[v],color = 'tab:purple', linewidth=3)
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowHig[v], maxpowHig[v],color = 'tab:purple', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowHig[1:None],color = 'tab:purple',label='High ' + r'$f_B$', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz; Med = ' + str(numpy.around(numpy.median(baserates),decimals=1)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3ModToBaseRatioAvg3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


# Fourth Measurement: Power Ratio / Baseline Power Ratio
fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0]) # Find power ratio for each modulation frequency
	if BaseIVLval == 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='k')
	elif BaseIVLval < 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='r')
	maxpratio = numpy.amax([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power Ratio / Base Power Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(numrerands-1,-1,-1):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0]) # Find power ratio for each modulation frequency
	
	ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Power Ratio / Base Power Ratio')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

ax0 = pyplot.subplot2grid((10, 10), (0, 0), colspan=9)
ax1 = pyplot.subplot2grid((10, 10), (1, 9), rowspan=9)
ax2 = pyplot.subplot2grid((10, 10), (1, 0), rowspan=9, colspan=9)
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))
cmap = pl.cm.plasma
norm = matplotlib.colors.Normalize(vmin=numpy.amin(baserates), vmax=numpy.amax(baserates))
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label(r'$f_B$')

PeakFreqs = numpy.zeros((numrerands,), dtype=numpy.float)
ax2.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0]) # Find power ratio for each modulation frequency
	
	ax2.plot(modfreqnums, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0]) # Find power ratio for each modulation frequency
	
	PeakFreqs[y] = modfreqnums[numpy.argmax(PowerRatios[1:None])]
	num = numpy.random.rand(1)*.6-0.3
	ax2.scatter(modfreqnums[numpy.argmax(PowerRatios[1:None])]+num,maxpratio,color=colors[y], alpha=0.5)

ax0.hist(PeakFreqs,bins=numpy.arange(18)+1.5)
ax0.set_xlim(0,16)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.get_xaxis().set_ticks([])
ax0.yaxis.tick_right()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlim(0,16)
ax2.set_xticks(modfreqnums)
ax2.set_xticklabels(modfreqstrs, rotation = 45)
ax2.set_xlabel(r'$f_i$')
ax2.set_ylim(0,maxpratio)
ax2.set_ylabel('Power Ratio / Base Power Ratio')
ax2.autoscale(enable=True, axis='x', tight=True)
pyplot.axis('tight')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatio3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	ttest = stats.ttest_1samp(powvec, 1)
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,maxpow[l]+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,maxpow[l]+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,maxpow[l]+0.02,'*')

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums,avgpow[1:None],color = 'tab:orange')
# pyplot.errorbar(modfreqnums,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',elinewidth=3)
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio / Base Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowB = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecNIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
		if HCT[y][0] == 4:
			powvecIVL[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
		elif HCT[y][0] < 4:
			powvecNIVL[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
	
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	if (numpy.sum(powvecIVL!=0) > 0):
		avgpowIVL[l] = numpy.median(powvecIVL[powvecIVL!=0])
		stdpowIVL[l] = numpy.std(powvecIVL[powvecIVL!=0])
		upper_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 75)
		lower_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 25)
		minpowIVL[l] = numpy.amin(powvecIVL[powvecIVL!=0])
		maxpowIVL[l] = numpy.amax(powvecIVL[powvecIVL!=0])
	if (numpy.sum(powvecNIVL!=0) > 0):
		avgpowNIVL[l] = numpy.median(powvecNIVL[powvecNIVL!=0])
		stdpowNIVL[l] = numpy.std(powvecNIVL[powvecNIVL!=0])
		upper_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 75)
		lower_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 25)
		minpowNIVL[l] = numpy.amin(powvecNIVL[powvecNIVL!=0])
		maxpowNIVL[l] = numpy.amax(powvecNIVL[powvecNIVL!=0])

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums+num,avgpow[1:None],color = 'tab:orange',label='All', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileIVL[v], upper_quartileIVL[v],color = 'tab:red', linewidth=3)
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowIVL[v], maxpowIVL[v],color = 'tab:red', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowIVL[1:None],color = 'tab:red',label='IVL', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileNIVL[v], upper_quartileNIVL[v],color = 'tab:blue', linewidth=3)
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowNIVL[v], maxpowNIVL[v],color = 'tab:blue', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowNIVL[1:None],color = 'tab:blue',label='NIVL', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio / Base Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvecLow = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecHig = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		if baserates[y] < numpy.median(baserates):
			powvecLow[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
		elif baserates[y] >= numpy.median(baserates):
			powvecHig[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
	
	if (numpy.sum(powvecLow!=0) > 0):
		avgpowLow[l] = numpy.median(powvecLow[powvecLow!=0])
		stdpowLow[l] = numpy.std(powvecLow[powvecLow!=0])
		upper_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 75)
		lower_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 25)
		minpowLow[l] = numpy.amin(powvecLow[powvecLow!=0])
		maxpowLow[l] = numpy.amax(powvecLow[powvecLow!=0])
	if (numpy.sum(powvecHig!=0) > 0):
		avgpowHig[l] = numpy.median(powvecHig[powvecHig!=0])
		stdpowHig[l] = numpy.std(powvecHig[powvecHig!=0])
		upper_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 75)
		lower_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 25)
		minpowHig[l] = numpy.amin(powvecHig[powvecHig!=0])
		maxpowHig[l] = numpy.amax(powvecHig[powvecHig!=0])
	
	ttest = stats.ttest_ind(powvecLow, powvecHig) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileLow)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileLow[v], upper_quartileLow[v],color = 'tab:green', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowLow[v], maxpowLow[v],color = 'tab:green', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowLow[1:None],color = 'tab:green',label='Low ' + r'$f_B$', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileHig[v], upper_quartileHig[v],color = 'tab:purple', linewidth=3)
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowHig[v], maxpowHig[v],color = 'tab:purple', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowHig[1:None],color = 'tab:purple',label='High ' + r'$f_B$', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([1,1]),color='k',linestyle='--')
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Power Ratio / Base Power Ratio')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz; Med = ' + str(numpy.around(numpy.median(baserates),decimals=1)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3PowerRatiosOverBaseRatioAvg3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

# Fifth Measurement: Baseline Differences
fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')

modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]-PSDXX[y][0][l] # Find power ratio for each modulation frequency
	if BaseIVLval == 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='k')
	elif BaseIVLval < 4:
		ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color='r')
	maxpratio = numpy.amax([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Baseline Difference')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0
for y in range(numrerands-1,-1,-1):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]/PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	ax.plot(modfreqnums, BaseRate, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
ax.view_init(elev=20, azim=280)
ax.set_xlim(0,16)
ax.set_zlim(0,maxpratio)
ax.set_xticks(modfreqnums)
ax.set_xticklabels(modfreqstrs, rotation = 45)
ax.set_xlabel("\n" + "\n" + r'$f_i$')
ax.set_ylabel(r'$f_B$')
ax.set_zlabel('Baseline Difference')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

ax0 = pyplot.subplot2grid((10, 10), (0, 0), colspan=9)
ax1 = pyplot.subplot2grid((10, 10), (1, 9), rowspan=9)
ax2 = pyplot.subplot2grid((10, 10), (1, 0), rowspan=9, colspan=9)
PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
modfreqnums = numpy.arange(1,len(modfreqs))
modfreqstrs = numpy.array(modfreqs[1:None],dtype='|S3')
maxpratio = 0

n = numrerands
colors = pl.cm.plasma((baserates-minrate)/(maxrate-minrate))
cmap = pl.cm.plasma
norm = matplotlib.colors.Normalize(vmin=numpy.amin(baserates), vmax=numpy.amax(baserates))
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
cb1.set_label(r'$f_B$')

PeakFreqs = numpy.zeros((numrerands,), dtype=numpy.float)
ax2.plot(numpy.array([0,16]), numpy.array([0,0]),color='k',linestyle='--')

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	BaseIVLval = HCT[y][0]
	BaseRate = numpy.repeat(RATES[y][0],len(modfreqnums)) # repeat value for each mod frequency to create vector
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]-PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	ax2.plot(modfreqnums, PowerRatios[1:None],color=colors[y])
	maxpratio = numpy.max([maxpratio,numpy.amax(PowerRatios[1:None])])

for y in range(0,numrerands):
	PowerRatios = numpy.zeros((len(modfreqs),), dtype=numpy.float)
	for l in range(1,len(modfreqs)):
		PowerRatios[l] = PSDXX[y][l][l]-PSDXX[y][0][l] # Find power ratio for each modulation frequency
	
	PeakFreqs[y] = modfreqnums[numpy.argmax(PowerRatios[1:None])]
	num = numpy.random.rand(1)*.6-0.3
	ax2.scatter(modfreqnums[numpy.argmax(PowerRatios[1:None])]+num,maxpratio,color=colors[y], alpha=0.5)

ax0.hist(PeakFreqs,bins=numpy.arange(18)+1.5)
ax0.set_xlim(0,16)
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.get_xaxis().set_ticks([])
ax0.yaxis.tick_right()

ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlim(0,16)
ax2.set_xticks(modfreqnums)
ax2.set_xticklabels(modfreqstrs, rotation = 45)
ax2.set_xlabel(r'$f_i$')
ax2.set_ylim(0,maxpratio)
ax2.set_ylabel('Baseline Difference')
ax2.autoscale(enable=True, axis='x', tight=True)
pyplot.axis('tight')

pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifference3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	ttest = stats.ttest_1samp(powvec, 1)
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,maxpow[l]+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,maxpow[l]+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,maxpow[l]+0.02,'*')

pyplot.plot(numpy.array([0,16]), numpy.array([0,0]),color='k',linestyle='--')
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1],minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums,avgpow[1:None],color = 'tab:orange')
# pyplot.errorbar(modfreqnums,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',elinewidth=3)
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Difference')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartile = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowNIVL = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvec = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecNIVL = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		powvec[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
		if HCT[y][0] == 4:
			powvecIVL[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
		elif HCT[y][0] < 4:
			powvecNIVL[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
	
	avgpow[l] = numpy.median(powvec)
	stdpow[l] = numpy.std(powvec)
	upper_quartile[l] = numpy.percentile(powvec, 75)
	lower_quartile[l] = numpy.percentile(powvec, 25)
	minpow[l] = numpy.amin(powvec)
	maxpow[l] = numpy.amax(powvec)
	if (numpy.sum(powvecIVL!=0) > 0):
		avgpowIVL[l] = numpy.median(powvecIVL[powvecIVL!=0])
		stdpowIVL[l] = numpy.std(powvecIVL[powvecIVL!=0])
		upper_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 75)
		lower_quartileIVL[l] = numpy.percentile(powvecIVL[powvecIVL!=0], 25)
		minpowIVL[l] = numpy.amin(powvecIVL[powvecIVL!=0])
		maxpowIVL[l] = numpy.amax(powvecIVL[powvecIVL!=0])
	if (numpy.sum(powvecNIVL!=0) > 0):
		avgpowNIVL[l] = numpy.median(powvecNIVL[powvecNIVL!=0])
		stdpowNIVL[l] = numpy.std(powvecNIVL[powvecNIVL!=0])
		upper_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 75)
		lower_quartileNIVL[l] = numpy.percentile(powvecNIVL[powvecNIVL!=0], 25)
		minpowNIVL[l] = numpy.amin(powvecNIVL[powvecNIVL!=0])
		maxpowNIVL[l] = numpy.amax(powvecNIVL[powvecNIVL!=0])

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartile[v], upper_quartile[v],color = 'tab:orange', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpow[v], maxpow[v],color = 'tab:orange', linewidth=1)
pyplot.plot(modfreqnums+num,avgpow[1:None],color = 'tab:orange',label='All', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileIVL[v], upper_quartileIVL[v],color = 'tab:red', linewidth=3)
for v in range(1,len(lower_quartileIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowIVL[v], maxpowIVL[v],color = 'tab:red', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowIVL[1:None],color = 'tab:red',label='IVL', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileNIVL[v], upper_quartileNIVL[v],color = 'tab:blue', linewidth=3)
for v in range(1,len(lower_quartileNIVL)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowNIVL[v], maxpowNIVL[v],color = 'tab:blue', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowNIVL[1:None],color = 'tab:blue',label='NIVL', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([0,0]),color='k',linestyle='--')
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpow[1:None],stdpow[1:None],color = 'tab:orange',ecolor='tab:orange',label='All', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowIVL[1:None],stdpowIVL[1:None],color = 'tab:red',ecolor='tab:red',label='IVL', alpha=0.8)
# pyplot.errorbar(modfreqnums+numpy.random.rand(1)*.4-0.2,avgpowNIVL[1:None],stdpowNIVL[1:None],color = 'tab:blue',ecolor='tab:blue',label='NIVL', alpha=0.8)
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Difference')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg2.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg2.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


fig = pyplot.figure()
avgpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
avgpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
stdpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

upper_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowLow = numpy.zeros((len(modfreqs),), dtype=numpy.float)
upper_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
lower_quartileHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
minpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)
maxpowHig = numpy.zeros((len(modfreqs),), dtype=numpy.float)

pvals = numpy.zeros((len(modfreqs),), dtype=numpy.float)

for l in range(1,len(modfreqs)):
	powvecLow = numpy.zeros((numrerands,), dtype=numpy.float)
	powvecHig = numpy.zeros((numrerands,), dtype=numpy.float)
	for y in range(0,numrerands):
		if baserates[y] < numpy.median(baserates):
			powvecLow[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
		elif baserates[y] >= numpy.median(baserates):
			powvecHig[y] = PSDXX[y][l][l]-PSDXX[y][0][l]
	
	if (numpy.sum(powvecLow!=0) > 0):
		avgpowLow[l] = numpy.median(powvecLow[powvecLow!=0])
		stdpowLow[l] = numpy.std(powvecLow[powvecLow!=0])
		upper_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 75)
		lower_quartileLow[l] = numpy.percentile(powvecLow[powvecLow!=0], 25)
		minpowLow[l] = numpy.amin(powvecLow[powvecLow!=0])
		maxpowLow[l] = numpy.amax(powvecLow[powvecLow!=0])
	if (numpy.sum(powvecHig!=0) > 0):
		avgpowHig[l] = numpy.median(powvecHig[powvecHig!=0])
		stdpowHig[l] = numpy.std(powvecHig[powvecHig!=0])
		upper_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 75)
		lower_quartileHig[l] = numpy.percentile(powvecHig[powvecHig!=0], 25)
		minpowHig[l] = numpy.amin(powvecHig[powvecHig!=0])
		maxpowHig[l] = numpy.amax(powvecHig[powvecHig!=0])
	
	ttest = stats.ttest_ind(powvecLow, powvecHig) # Compare to mean of baseline power ratios
	pvals[l] = ttest[1]
	# if (ttest[1] < 0.001):
	#     pyplot.text(l-0.3,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'***')
	# elif (ttest[1] < 0.01):
	#     pyplot.text(l-0.2,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'**')
	# elif (ttest[1] < 0.05):
	#     pyplot.text(l-0.1,numpy.amax([maxpowLow[l],maxpowHig[l]])+0.02,'*')

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileLow)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileLow[v], upper_quartileLow[v],color = 'tab:green', linewidth=3)
for v in range(1,len(lower_quartile)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowLow[v], maxpowLow[v],color = 'tab:green', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowLow[1:None],color = 'tab:green',label='Low ' + r'$f_B$', alpha=0.8)

num = numpy.random.rand(1)*.4-0.2
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,lower_quartileHig[v], upper_quartileHig[v],color = 'tab:purple', linewidth=3)
for v in range(1,len(lower_quartileHig)):
	pyplot.vlines(modfreqnums[v-1]+num,minpowHig[v], maxpowHig[v],color = 'tab:purple', linewidth=1)
pyplot.plot(modfreqnums+num,avgpowHig[1:None],color = 'tab:purple',label='High ' + r'$f_B$', alpha=0.8)

pyplot.plot(numpy.array([0,16]), numpy.array([0,0]),color='k',linestyle='--')
pyplot.legend(loc='upper left')
pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Baseline Difference')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz; Med = ' + str(numpy.around(numpy.median(baserates),decimals=1)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg3.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaselineDifferenceAvg3.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

# Effects of fb
Rvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
for l in range(1,len(modfreqs)):
	rates = numpy.zeros((numrerands,), dtype=numpy.float)
	psdxx = numpy.zeros((numrerands,), dtype=numpy.float)
	fig, axarr = pyplot.subplots(nrows=1, ncols=1)
	for y in range(0,numrerands):
		rates[y] = RATES[y][0]
		psdxx[y] = PSDXX[y][l][l]
	
	[R, pval] = stats.pearsonr(rates,psdxx)
	Rvalues[l] = R
	pvalues[l] = pval
	axarr.scatter(rates,psdxx,color='tab:orange')
	axarr.set_yscale('log')
	axarr.set_xlim(minrate-1,maxrate+1)
	axarr.set_ylim(numpy.amin(psdxx)-numpy.amin(psdxx)/10,numpy.amax(psdxx)+numpy.amax(psdxx)/10)
	axarr.set_xlabel(r'$f_B$')
	axarr.set_ylabel(str(modfreqstrs[l-1]) + ' Hz PSD')
	axarr.set_title('R = ' + str(R) + '; p-value = ' + str(pval))
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPower.pdf', bbox_inches='tight')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPower.png', bbox_inches='tight')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()

fig = pyplot.figure()
pyplot.bar(modfreqnums,Rvalues[1:None])
for l in range(1,len(Rvalues)):
	if pvalues[l] < 0.001:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '***', ha='center', va='bottom')
	elif pvalues[l] < 0.01:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '**', ha='center', va='bottom')
	elif pvalues[l] < 0.05:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '*', ha='center', va='bottom')

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Pearson Correlation Coefficient')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_Power.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_Power.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()

Rvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
for l in range(1,len(modfreqs)):
	rates = numpy.zeros((numrerands,), dtype=numpy.float)
	psdxx = numpy.zeros((numrerands,), dtype=numpy.float)
	fig, axarr = pyplot.subplots(nrows=1, ncols=1)
	for y in range(0,numrerands):
		rates[y] = RATES[y][0]
		psdxx[y] = PSDXX[y][l][l]/PSDXX[y][l][0]
	
	[R, pval] = stats.pearsonr(rates,psdxx)
	Rvalues[l] = R
	pvalues[l] = pval
	axarr.scatter(rates,psdxx,color='tab:orange')
	# axarr.set_yscale('log')
	axarr.set_xlim(minrate-1,maxrate+1)
	axarr.set_ylim(numpy.amin(psdxx)-numpy.amin(psdxx)/10,numpy.amax(psdxx)+numpy.amax(psdxx)/10)
	axarr.set_xlabel(r'$f_B$')
	axarr.set_ylabel(str(modfreqstrs[l-1]) + ' Hz Power Ratio')
	axarr.set_title('R = ' + str(R) + '; p-value = ' + str(pval))
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPowerRatio.pdf', bbox_inches='tight')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPowerRatio.png', bbox_inches='tight')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()

fig = pyplot.figure()
pyplot.bar(modfreqnums,Rvalues[1:None])
for l in range(1,len(Rvalues)):
	if pvalues[l] < 0.001:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '***', ha='center', va='bottom')
	elif pvalues[l] < 0.01:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '**', ha='center', va='bottom')
	elif pvalues[l] < 0.05:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '*', ha='center', va='bottom')

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Pearson Correlation Coefficient')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_PowerRatio.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_PowerRatio.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


Rvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
for l in range(1,len(modfreqs)):
	rates = numpy.zeros((numrerands,), dtype=numpy.float)
	psdxx = numpy.zeros((numrerands,), dtype=numpy.float)
	fig, axarr = pyplot.subplots(nrows=1, ncols=1)
	for y in range(0,numrerands):
		rates[y] = RATES[y][0]
		psdxx[y] = PSDXX[y][l][l]/PSDXX[y][0][l]
	
	[R, pval] = stats.pearsonr(rates,psdxx)
	Rvalues[l] = R
	pvalues[l] = pval
	axarr.scatter(rates,psdxx,color='tab:orange')
	# axarr.set_yscale('log')
	axarr.set_xlim(minrate-1,maxrate+1)
	axarr.set_ylim(numpy.amin(psdxx)-numpy.amin(psdxx)/10,numpy.amax(psdxx)+numpy.amax(psdxx)/10)
	axarr.set_xlabel(r'$f_B$')
	axarr.set_ylabel(str(modfreqstrs[l-1]) + ' Hz Baseline Ratio')
	axarr.set_title('R = ' + str(R) + '; p-value = ' + str(pval))
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzModToBaseRatio.pdf', bbox_inches='tight')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzModToBaseRatio.png', bbox_inches='tight')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()

fig = pyplot.figure()
pyplot.bar(modfreqnums,Rvalues[1:None])
for l in range(1,len(Rvalues)):
	if pvalues[l] < 0.001:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '***', ha='center', va='bottom')
	elif pvalues[l] < 0.01:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '**', ha='center', va='bottom')
	elif pvalues[l] < 0.05:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '*', ha='center', va='bottom')

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Pearson Correlation Coefficient')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_ModToBaseRatio.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_ModToBaseRatio.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()


Rvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
pvalues = numpy.zeros((len(modfreqs),), dtype=numpy.float)
for l in range(1,len(modfreqs)):
	rates = numpy.zeros((numrerands,), dtype=numpy.float)
	psdxx = numpy.zeros((numrerands,), dtype=numpy.float)
	fig, axarr = pyplot.subplots(nrows=1, ncols=1)
	for y in range(0,numrerands):
		rates[y] = RATES[y][0]
		psdxx[y] = (PSDXX[y][l][l]/PSDXX[y][l][0])/(PSDXX[y][0][l]/PSDXX[y][0][0])
	
	[R, pval] = stats.pearsonr(rates,psdxx)
	Rvalues[l] = R
	pvalues[l] = pval
	axarr.scatter(rates,psdxx,color='tab:orange')
	# axarr.set_yscale('log')
	axarr.set_xlim(minrate-1,maxrate+1)
	axarr.set_ylim(numpy.amin(psdxx)-numpy.amin(psdxx)/10,numpy.amax(psdxx)+numpy.amax(psdxx)/10)
	axarr.set_xlabel(r'$f_B$')
	axarr.set_ylabel(str(modfreqstrs[l-1]) + ' Hz Power Ratio / Base Power Ratio')
	axarr.set_title('R = ' + str(R) + '; p-value = ' + str(pval))
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPowerRatioOverBasePowerRatio.pdf', bbox_inches='tight')
	pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulation_' + str(modfreqstrs[l-1]) + 'HzPowerRatioOverBasePowerRatio.png', bbox_inches='tight')
	pyplot.gcf().clear()
	pyplot.cla()
	pyplot.clf()
	pyplot.close()

fig = pyplot.figure()
pyplot.bar(modfreqnums,Rvalues[1:None])
for l in range(1,len(Rvalues)):
	if pvalues[l] < 0.001:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '***', ha='center', va='bottom')
	elif pvalues[l] < 0.01:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '**', ha='center', va='bottom')
	elif pvalues[l] < 0.05:
		pyplot.text(modfreqnums[l-1], Rvalues[l]+0.005*numpy.sign(Rvalues[l]), '*', ha='center', va='bottom')

pyplot.xlim(0,16)
pyplot.xticks(modfreqnums,modfreqstrs)
pyplot.xlabel(r'$f_i$')
pyplot.ylabel('Pearson Correlation Coefficient')
pyplot.title('N = ' + str(int(numrerands)) + '; Percent IVL = ' + str(int(percenthct)) + '%; Base Rates = ' + str(int(minrate)) + ' to ' + str(int(maxrate)) + ' Hz')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_PowerRatioOverBasePowerRatio.pdf')
pyplot.savefig('Plots_' + Cell_Name + '/' + Case + '_IS3BaseRateVModulationPearson_PowerRatioOverBasePowerRatio.png')
pyplot.gcf().clear()
pyplot.cla()
pyplot.clf()
pyplot.close()
