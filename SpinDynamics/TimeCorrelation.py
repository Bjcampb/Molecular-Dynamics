import numpy as np

###############################################################################
# Load and format data
###############################################################################
'''Loaded data is presented as one column of many rows. If the total number of
   spins in the simulation is 10 and the number of time steps of the simulation
   is 100 then the file contains 1000 rows. The first 10 rows will be all the
   spins components for the first time step, the second 10 rows will be all the
   spin components for the second time step, etc.'''

#Load data
pos1 = np.loadtxt('C:\\Users\\bjcam\\Google Drive\\Molecular Dynamics\\'
                  'Project 5\\Data\\L10T01trail1.txt', delimiter=',')

pos2 = np.loadtxt('C:\\Users\\bjcam\\Google Drive\\Molecular Dynamics\\'
                  'Project 5\\Data\\L10T01trail2.txt', delimiter=',')

pos3 = np.loadtxt('C:\\Users\\bjcam\\Google Drive\\Molecular Dynamics\\'
                  'Project 5\\Data\\L10T01trail3.txt', delimiter=',')

pos4 = np.loadtxt('C:\\Users\\bjcam\\Google Drive\\Molecular Dynamics\\'
                  'Project 5\\Data\\L10T01trail4.txt', delimiter=',')

pos5 = np.loadtxt('C:\\Users\\bjcam\\Google Drive\\Molecular Dynamics\\'
                  'Project 5\\Data\\L10T01trail5.txt', delimiter=',')




NumberAtoms = 10 # Total number of spins
TotalSteps = 5000 # Total number of timesteps in simulation
halfSteps = int(TotalSteps / 2)

# Reshape loaded data to look like matrix (spin, timestep)
pos1 = np.transpose(np.reshape(pos1,(TotalSteps, NumberAtoms)))
pos2 = np.transpose(np.reshape(pos2,(TotalSteps, NumberAtoms)))
pos3 = np.transpose(np.reshape(pos3,(TotalSteps, NumberAtoms)))
pos4 = np.transpose(np.reshape(pos4,(TotalSteps, NumberAtoms)))
pos5 = np.transpose(np.reshape(pos5,(TotalSteps, NumberAtoms)))




###############################################################################
#      Define different functions for time correlation calculations
###############################################################################
'''Arrays for the following fucntions must be presented in the reshaped format
   given above.'''


def ensembleAverage(array1, array2):
    average = 0
    
    neigh1 = np.roll(array2,-1)
    neigh2 = np.roll(array2, 1)
    pair1 = array1 * neigh1
    pair2 = array1 * neigh2
    mult = (np.sum(pair1) + np.sum(pair2)) / np.shape(array1)[0]
    var = (np.mean(array1) * np.mean(array2))
    average = mult - var
    return average


def timeCorrelation(array, corrStep):
    ''' Conducts time correlation value averaged over entire dataset based on 
        time separation between averages.
        - array: reshaped array from loaded data
        - corrStep: the number of steps separating each average'''
    corrValue = 0 # initialize correlation value for output
    tmax = np.arange(int(np.shape(array)[1]/2))
    avgVal = int(np.shape(array)[1]/2)
    for i in tmax:
        j = i + corrStep
        val = ensembleAverage(array[:,i], array[:,j])
        corrValue += val
    corrValue /= avgVal
    return corrValue


def correlationFunction(array, tmax):
    CorrelationValues = np.zeros(tmax)
    steps = np.arange(tmax)
    for i in steps:
        CorrelationValues[i] = timeCorrelation(array, i+1)
        if i % 100 == 0:
            print("Step %d" % i)
    return CorrelationValues


###############################################################################
#      Run functions
###############################################################################

#Find Correlation Function
CorrelationValues1 = correlationFunction(pos1, halfSteps)
CorrelationValues2 = correlationFunction(pos2, halfSteps)
CorrelationValues3 = correlationFunction(pos3, halfSteps)
CorrelationValues4 = correlationFunction(pos4, halfSteps)
CorrelationValues5 = correlationFunction(pos5, halfSteps)

np.savetxt('CorrL10T01Trail1.txt', CorrelationValues1,delimiter=',')
np.savetxt('CorrL10T01Trail2.txt', CorrelationValues2,delimiter=',')
np.savetxt('CorrL10T01Trail3.txt', CorrelationValues3,delimiter=',')
np.savetxt('CorrL10T01Trail4.txt', CorrelationValues4,delimiter=',')
np.savetxt('CorrL10T01Trail5.txt', CorrelationValues5,delimiter=',')

#np.savetxt('CorrelationValues.txt', CorrelationValues, delimiter=',')   # Output Correlation
                   
## Perform FFT                     
#Y = np.fft.fft(CorrelationValues)
#np.savetxt('FFT.txt', Y, delimiter=',')   # Output FFT
#
##Find Absolute Value of FFT
#Absolute = np.absolute(Y)
#np.savetxt('MagFFT.txt', Absolute, delimiter=',')   # Output Mag of FFT

    
    
    
    
    
    
    