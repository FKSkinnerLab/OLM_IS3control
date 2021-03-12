### Test Script for a file found in the SDprox2 results from initial Parallel Simulations
from __future__ import division
import numpy
import os
# os.system("python -m pip install efel --user")
# os.system("python -m pip install --upgrade setuptools --user")
# os.system("python -m pip install --upgrade matplotlib --user")

Cell_Name = h.cellname
if h.cell == 1:
    # currhold = 0.049 # 7.25 Hz
    # currhold = 0.0289 # Below Threshold
    currhold = 0.030 # ~1-2 Hz (almost rheobase)
    # currhold = 0.037 # ~4 Hz
elif h.cell == 2:
    # currhold = 0.045 # 7.25 Hz
    # currhold = 0.0214 # Below Threshold
    currhold = 0.022 # ~1-2 Hz (almost rheobase)
    # currhold = 0.031 # ~4 Hz

Case = Cell_Name + '_E_COM_I_COM'

execfile("CutSpikes_HighConductanceMeasurements.py")

# HC Treshold Measurement Values
tstop = 10000 # seconds
font_size = 13
numrerands = 50
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'numrerands.npy',numrerands)

modfreqs = numpy.array([0,0.5,1,2,3,4,5,8,9,10,12,15,16,20,25,30])
currsteps = numpy.arange(currhold,currhold+0.0025*numrerands,0.0025)

numinh = 0
numexc = 0
inhspikes = 0
excspikes = 0

# Set up parallel context
HCT = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.int)
MODFREQ = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
RATES = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
RATES2 = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
RATES3 = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
CURRENTS = numpy.zeros((numrerands,len(modfreqs)), dtype=numpy.float)
PSDXX = numpy.zeros((numrerands,len(modfreqs),len(modfreqs)), dtype=numpy.float)

count1 = numrerands*len(modfreqs)
count2 = 0
ParamMat = numpy.zeros((count1,2), dtype=numpy.int)
for y2 in range(0,numrerands):
	for modfreq2 in range(0,len(modfreqs)):
		ParamMat[count2][0] = y2
		ParamMat[count2][1] = modfreq2
		count2 = count2 + 1

numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'ParamMat.npy',ParamMat)

pc = h.ParallelContext()
pc.runworker()

for y in range(0,count2):
	pc.submit(getMeasures,numinh,numexc,inhspikes,excspikes,ParamMat[y][0]*6,ParamMat[y][0],modfreqs[ParamMat[y][1]],ParamMat[y][1],currsteps[ParamMat[y][0]],y)
while pc.working():
	results = pc.pyret()
	HCT[results[0]][results[4]] = results[1]
	MODFREQ[results[0]][results[4]] = results[2]
	RATES[results[0]][results[4]] = results[3]
	CURRENTS[results[0]][results[4]] = results[21]
	RATES2[results[0]][results[4]] = results[22]
	RATES3[results[0]][results[4]] = results[24]
	
	PSDXX[results[0]][results[4]][0] = results[5]
	PSDXX[results[0]][results[4]][1] = results[6]
	PSDXX[results[0]][results[4]][2] = results[7]
	PSDXX[results[0]][results[4]][3] = results[8]
	PSDXX[results[0]][results[4]][4] = results[9]
	PSDXX[results[0]][results[4]][5] = results[10]
	PSDXX[results[0]][results[4]][6] = results[11]
	PSDXX[results[0]][results[4]][7] = results[12]
	PSDXX[results[0]][results[4]][8] = results[13]
	PSDXX[results[0]][results[4]][9] = results[14]
	PSDXX[results[0]][results[4]][10] = results[15]
	PSDXX[results[0]][results[4]][11] = results[16]
	PSDXX[results[0]][results[4]][12] = results[17]
	PSDXX[results[0]][results[4]][13] = results[18]
	PSDXX[results[0]][results[4]][14] = results[19]
	PSDXX[results[0]][results[4]][15] = results[20]

pc.done()

numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'IVLMat.npy',HCT)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'MODFREQMat.npy',MODFREQ)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'RATESMat.npy',RATES)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'RATES2Mat.npy',RATES2)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'RATES3Mat.npy',RATES3)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'PSDXXMat.npy',PSDXX)
numpy.save('NPYFiles_' + Cell_Name + '/' + Case + 'CURRENTSMat.npy',CURRENTS)
