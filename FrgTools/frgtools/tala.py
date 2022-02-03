import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


def water_calibration(ratio, architecture, polymer):
	#glassglass
	#         y_406 = 0.02x + 0.016
	#         y_806 = 0.015x + 0.03
	#         y_TF4 = 0.009x + 0.013
	#         y_TF8 = -0.002x + 0.024
	#glassbs
	#         y_406 = 0.023x + 0.069
	#         y_806 = 0.031x + 0.216
	#         y_TF4 = 0.039x + 0.026
	#         y_TF8 = 0.071x + 0.096

	glassglassdict = {
	'406': [0.02,0.016],
	'806': [0.015,0.03],
	'TF4': [0.009,0.013],
	'TF8': [-0.002,0.024]
	}


	glassbsdict = {
	'406': [0.024,0.058],
	'806': [0.03,0.223],
	'TF4': [0.072,0.011],
	'TF8': [0.104,0.122]
	}
    
	if architecture == 'GPOLYG':
		x = (ratio - glassglassdict[polymer][1])/glassglassdict[polymer][0]
	elif architecture == 'GPOLYBS':
		x = (ratio - glassbsdict[polymer][1])/glassbsdict[polymer][0]
		
	return x

from tqdm import tqdm
def bscorrect_linescans(df):
	#currently using mini bff
	correctedAbs_list = []
	baselineAbs_list = []
	gaussian_fit_list = []

	for _,row in tqdm(df.iterrows()):
		for ii in row['r']:
			corrected, baseline = mini_bff(row['wl'], ii, row['polymer'], row['label'], plot = False)
			correctedAbs_list.append(corrected)
			baselineAbs_list.append(baseline)
	
	s = []
	for i in range(len(df)):
		s.append(len(df['position'][0])*i)
	print(s)
	
	
	correctedAbs_dict = {}
	baselineAbs_dict = {}
	# gaussianfit_dict = {}
	wardRatio_dict = {}

	if df['polymer'][0] == 'EVA':
		pol = '406'
	else:
		pol = 'TF4'

	for i in range(len(s)):
		x = i+1

		wardRatio_dict[i] = []

		if i+1 >= len(s):
			correctedAbs_dict[i] = correctedAbs_list[s[i]:]
			baselineAbs_dict[i] = baselineAbs_list[s[i]:]
			#         gaussianfit_dict[i] = gaussian_fit_list[s[i]:]
		else:
			correctedAbs_dict[i] = correctedAbs_list[s[i]:s[x]]
			baselineAbs_dict[i] = baselineAbs_list[s[i]:s[x]]
			#         gaussianfit_dict[i] = gaussian_fit_list[s[i]:s[x]]

		for pos in range(len(df['r'][0])):
			#         wardRatio_dict[i].append(cali_406(gaussianfit_dict[i][pos][102]/gaussianfit_dict[i][pos][15]))
			wardRatio_dict[i].append(water_calibration(
								correctedAbs_dict[i][pos][102]/correctedAbs_dict[i][pos][15],
								  df['architecture'][0],
								    pol
								    )
								)
	return correctedAbs_dict, baselineAbs_dict, wardRatio_dict


def getexpectedwater(t,rh,polymer_type):
	inv_t = 1/(t + 273.15)
	k = 8.617333e-5 #eV

	solubilitydict = {
	'406': [43.19284622, 0.28018506],
	'806': [7.64233379E04, 4.99533970E-01],
	'TF4': [0.87724685, 0.20267589],
	'TF8': [0.99962712, 0.20842967]
	}

	c0 = solubilitydict[polymer_type][0]
	dH = solubilitydict[polymer_type][1]

	c_h2o = c0*np.exp(-dH*inv_t/k) * (rh/100) #g/cm3
	return c_h2o*1000 #mg/cm3

def new_test():
	return "Hey! New Testing! Please work! If works, then make sure to commit changes and then restart VS Code, not just re-import :)"