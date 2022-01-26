import numpy as np

def testing_simple(a):
  return a**3


import pandas as pd
from .curveprocessing import rubberband


def mini_bff(wavelengths, reflectance, polymer_type, sampleName, plot = False):
    
    ####step 1, remove baseline using rubber band method, added to frgtools.curveprocessing
	
    wl = np.array(wavelengths)
    raw_absorbance = np.log(1/np.array(reflectance))
    baseline = rubberband(wl,raw_absorbance)

    corrected_abs = raw_absorbance - baseline
    
    ##currently uses these fixed peaks for the gaussian fitting
    ##may change that later to automatically find the local peaks of the curve (especially for POE)

    peak1 = wl[15] #1730
    peak2 = wl[31] #1762
    peak3 = wl[51] #1802
    peak4 = wl[102] #1904
    peak5 = wl[124] #1948
    
    if polymer_type == 'EVA':
        peaks_list = [peak1,peak2,peak3,peak4,peak5]
    elif polymer_type == 'POE':
        peaks_list = [peak1,peak2,peak3,peak4]
    
    if plot == True:
        fig, ax = plt.subplots(1,2, figsize = (8,3))

        ax[0].plot(wl, raw_absorbance, label = 'Raw')
        ax[0].plot(wl, baseline, label = 'Baseline')
        ax[0].legend()
        ax[0].set_xlabel('Wavelengths (nm)')
        ax[0].set_ylabel('Absorbance (AU)')

        ax[1].plot(wl, corrected_abs, label = 'Corrected')
        ax[1].legend()
#         ax[1].set_ylim(top = 0.35)
        ax[1].set_xlabel('Wavelengths (nm)')
        ax[1].set_ylabel('Baseline-Removed Absorbance (AU)')
        fig.suptitle(sampleName, y = 1.05)
        
        for i in peaks_list:
            ypeak = corrected_abs[np.where(wl == i)][0]
            ax[1].plot(np.ones((2,))*i,[0,ypeak], linestyle = '--', color = 'g')

        plt.tight_layout()
        plt.show()
        
    return corrected_abs, baseline
