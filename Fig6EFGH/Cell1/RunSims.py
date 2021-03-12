### Test Script for a file found in the SDprox2 results from initial Parallel Simulations
from __future__ import division
import numpy
import os
# os.system("python -m pip install efel --user")
# os.system("python -m pip install --upgrade setuptools --user")
# os.system("python -m pip install --upgrade matplotlib --user")

Cell_Name = h.cellname
Case = Cell_Name + '_E_COM_I_COM'

execfile("CutSpikes_HighConductanceMeasurements.py")

# HC Treshold Measurement Values
tstop = 10000 # seconds
font_size = 13
numrerands = 50
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'numrerands.npy',numrerands)

# Examples = zip(LNI_LNE_LIS_LES,HNI_LNE_LIS_LES,LNI_HNE_LIS_LES,HNI_HNE_LIS_LES,LNI_LNE_HIS_LES,HNI_LNE_HIS_LES,LNI_HNE_HIS_LES,HNI_HNE_HIS_LES,LNI_LNE_LIS_HES,HNI_LNE_LIS_HES,LNI_HNE_LIS_HES,HNI_HNE_LIS_HES,LNI_LNE_HIS_HES,HNI_LNE_HIS_HES,LNI_HNE_HIS_HES,HNI_HNE_HIS_HES)
# ExampleStrings = ['LNI_LNE_LIS_LES','HNI_LNE_LIS_LES','LNI_HNE_LIS_LES','HNI_HNE_LIS_LES','LNI_LNE_HIS_LES','HNI_LNE_HIS_LES','LNI_HNE_HIS_LES','HNI_HNE_HIS_LES','LNI_LNE_LIS_HES','HNI_LNE_LIS_HES','LNI_HNE_LIS_HES','HNI_HNE_LIS_HES','LNI_LNE_HIS_HES','HNI_LNE_HIS_HES','LNI_HNE_HIS_HES','HNI_HNE_HIS_HES']
Examples = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelParams.npy')
ExampleStrings = numpy.load('NPYFiles_' + Cell_Name + '/' + Case + '_ExampleHCModelStrings.npy')
ExampleString = ExampleStrings[0]

modfreqs = numpy.array([0,0.5,1,2,3,4,5,8,9,10,12,15,16,20,25,30])

numinh = Examples[0][0]
numexc = Examples[1][0]
inhspikes = Examples[2][0]
excspikes = Examples[3][0]

# Set up parallel context
pc = h.ParallelContext()
pc.runworker()

print 'pc.nhost() = ' + str(pc.nhost())

HCT = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.int)
MODFREQ = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
RATES = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
PSDXX = numpy.zeros((numrerands,len(modfreqs),len(modfreqs)), dtype=numpy.float)

for y in range(1,numrerands+1):
	for modfreq in range(0,len(modfreqs)):
			pc.submit(getMeasures,numinh,numexc,inhspikes,excspikes,y*6,y,modfreqs[modfreq],modfreq)

while pc.working():
	results = pc.pyret()
	HCT[results[0]-1][results[4]] = results[1]
	MODFREQ[results[0]-1][results[4]] = results[2]
	RATES[results[0]-1][results[4]] = results[3]
	
	PSDXX[results[0]-1][results[4]][0] = results[5]
	PSDXX[results[0]-1][results[4]][1] = results[6]
	PSDXX[results[0]-1][results[4]][2] = results[7]
	PSDXX[results[0]-1][results[4]][3] = results[8]
	PSDXX[results[0]-1][results[4]][4] = results[9]
	PSDXX[results[0]-1][results[4]][5] = results[10]
	PSDXX[results[0]-1][results[4]][6] = results[11]
	PSDXX[results[0]-1][results[4]][7] = results[12]
	PSDXX[results[0]-1][results[4]][8] = results[13]
	PSDXX[results[0]-1][results[4]][9] = results[14]
	PSDXX[results[0]-1][results[4]][10] = results[15]
	PSDXX[results[0]-1][results[4]][11] = results[16]
	PSDXX[results[0]-1][results[4]][12] = results[17]
	PSDXX[results[0]-1][results[4]][13] = results[18]
	PSDXX[results[0]-1][results[4]][14] = results[19]
	PSDXX[results[0]-1][results[4]][15] = results[20] 

numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'IVLMat.npy',HCT)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'MODFREQMat.npy',MODFREQ)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'RATESMat.npy',RATES)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'PSDXXMat.npy',PSDXX)
