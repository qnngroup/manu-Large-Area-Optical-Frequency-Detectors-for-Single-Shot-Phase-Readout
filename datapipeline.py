#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import sys
import os
from scipy import optimize
from scipy import constants
import scipy as scp
from scipy.signal import butter,filtfilt
from numba import jit, njit, prange
import shelve as shv

import json
import scipy
from scipy import signal
from scipy import stats 
import h5py

eCharge= 1.602e-19





def loadData(path):
    # Load data from a file

    # Implement a flexible oscilloscope data loader by using the preambel!

    dataCH1 = (pd.read_hdf(path+'/dataCH1.h5', key= 'CH1').to_numpy(dtype="float32"))
    dataCH2 = (pd.read_hdf(path+'/dataCH2.h5', key= 'CH2').to_numpy(dtype="float32"))
    dataCH3 = pd.read_hdf(path+'/dataCH3.h5', key= 'CH3').to_numpy(dtype="int8")
    print('Data loaded')
    return dataCH1, dataCH2, dataCH3

def preprocessData(dataCH1, dataCH2, dataCH3, params):
    # Preprocess data
    # The output is the box integrated data of CH1 and CH2
    #
    dataN = dataCH1.shape[0]
    dataM = dataCH1.shape[1]

    print('Preprocessing data')
    # Preprocess data
    deltaW=np.gradient(dataCH3,axis = 1)
    lowestPeak=(np.amin(deltaW))

    if lowestPeak > -20:
       print('Lowest peak is too low. Check the data')
       return

    
    bstart =-3
    bStop = 140

    boxSum = np.zeros((dataN,dataM))
    boxSum2 = np.zeros((dataN,dataM))

    for l in range(0,dataN-1):
        peaks, _ = signal.find_peaks(-(deltaW[l,:]),distance=140,height=20)
        # Skipping first peak and last peak
        peaks = peaks[1:-2]
        for i,j in enumerate(peaks):    
            baseline= np.mean(dataCH1[l,j-(200-bStop):j])
            boxSum[l,i]=np.sum((dataCH1[l,j+bstart:j+bStop]-baseline)/1e9)*1e-7
            # 1e7 V/A TIA gain
            baseline= np.mean(dataCH2[l,j-(200-bStop):j])
            boxSum2[l,i]=np.sum((dataCH2[l,j+bstart:j+bStop]-baseline)/0.5e9)*1e-7

    

    return boxSum/eCharge, boxSum2/eCharge

def getTriggerSortedData(data, trigger,boxlen,startIndex, params):
    # Get the data sorted by the trigger
    # data is the data vector - expected to be a 1D vector and in units of V
    # boxlen is the length of the box in samples
    # trigger is the trigger position - expected to be a 1D vector
    deltaW=np.gradient(trigger)
    peaks, _ = signal.find_peaks(-(deltaW),distance=180,height=20)
    peaks = peaks[1:-2]

    dataTriggsorted = np.zeros((len(peaks),boxlen+np.abs(startIndex)))

    for i,j in enumerate(peaks):
        dataTriggsorted[i,:] = data[j+startIndex:j+boxlen]

    return dataTriggsorted


def getPulseShape(data,trigger,length,SubsetIDX):
    # Get the pulse shape of a single pulse
    # data is the data vector - expected to be a 1D vector and in units of V
    # trigger is the trigger position - expected to be a 1D vector
 

    deltaW=np.gradient(trigger)
    peaks, _ = signal.find_peaks(-(deltaW),distance=20,height=20)
    peaks = peaks[1:-2]
   
    for i,j in enumerate(peaks[SubsetIDX:-3]):

        pulseShape = pulseShape + (data[j:j+length])/1e9
        # 1e7 V/A TIA gain


    return pulseShape/len(peaks[SubsetIDX:-3])


def __main__():
    #print('This is not an executable script')

    dataPath = os.path.abspath(os.path.join(os.path.expanduser('~'),'Nextcloud/PHzElectronics/Data/cepDetector/04112022/'))
    measID = '4CMCT0'
    print(dataPath)

    measPath = os.path.abspath(os.path.join(dataPath,measID))

    dataCH1, dataCH2, dataCH3 = loadData(measPath)
    print('Data loaded ')
    
    n,m = dataCH1.shape
    lenbox = 200
    print(n,m)
    trigSorted = np.zeros((n,lenbox+60))
    listID = [0,10,20,25,30]
    for i in listID:
        trigSorted[i,:] = (-(np.mean(getTriggerSortedData(dataCH1[i,:], dataCH3[i,:], lenbox,startIndex=-60, params=0)[40000:-1,:],axis=0)))/1e9
        plt.plot(trigSorted[i,:]-np.min(trigSorted[i,:]))
    plt.xlabel('Time Index in (100ns)')
    plt.ylabel('Current (A)')
    plt.title('TIA current')
    plt.savefig(measPath+'/TIAcurrent.png')
    plt.show()

    trigSorted2 = np.zeros((n,lenbox+60))

    for i in listID:
        trigSorted2[i,:] = ((np.mean(getTriggerSortedData(dataCH2[i,:], dataCH3[i,:], lenbox,startIndex=-60, params=0)[40000:-1,:],axis=0)))/1e9
        plt.plot(trigSorted2[i,:]-np.min(trigSorted2[i,:]))
    plt.xlabel('Time Index in (100ns)')
    plt.ylabel('Current (A)')
    plt.title('Pyro current')
    plt.savefig(measPath+'/Pyrocurrent.png')
    plt.show()
    

    # params = []
    # boxSum, boxSum2 = preprocessData(dataCH1, dataCH2, dataCH3, params)

    # print('Data loaded and preprocessed')

    # plt.plot(boxSum[25,:])
    # plt.show()
    return

if __name__ == '__main__':
    __main__()  

# %%
