import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .curveprocessing import rubberband


def mini_bff(wavelengths, reflectance, polymer_type, sampleName, plot = False):
    #creates a mini friend for tala
    	####step 1, remove baseline using rubber band method, added to frgtools.curveprocessing
	
	wl = np.array(wavelengths)
	raw_absorbance = np.log(1/np.array(reflectance))
	baseline = rubberband(wl,raw_absorbance)

	corrected_abs = raw_absorbance - baseline
    
	##currently uses these fixed peaks for the gaussian fitting
	##may change that later to automatically find the local peaks of the curve (especially for POE)

	peak1 = wl[np.where(wl == 1730)[0][0]] #1730
	peak2 = wl[np.where(wl == 1762)[0][0]] #1762
	peak3 = wl[np.where(wl == 1802)[0][0]] #1802
	peak4 = wl[np.where(wl == 1904)[0][0]] #1904
	peak5 = wl[np.where(wl == 1948)[0][0]] #1948

	if polymer_type == 'EVA':
		peaks_list = [peak1,peak2,peak3,peak4,peak5]
	elif polymer_type == 'POE':
		peaks_list = [peak1,peak2,peak3,peak4]
	# print(peaks_list)

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


def water_calibration(ratio, architecture, polymer, alt_cal = False):
	#this function takes the absorbance ratio between 1900/1730 as y and returns x, the equivalent water concentration 

	#glassglass
	#		y_406 = 0.014x + 0.021
	#		y_806 = 0.011x + 0.034
	#         y_TF4 = 0.009x + 0.013
	#         y_TF8 = -0.002x + 0.024
	#glassbs
	#         y_406 = 0.023x + 0.069
	#         y_806 = 0.031x + 0.216
	#         y_TF4 = 0.039x + 0.026
	#         y_TF8 = 0.071x + 0.096

	if alt_cal == True:
		p1 = 24.75
		p2 = 0.3461
		# x = (ratio - p2)/p1
		h2o_meas = ratio*p1 + p2
		return h2o_meas

	else:
		glassglassdict = {

		#GG dict as of 4/28/22
		'406': [0.02284,0.00317],
		'806': [0.01302,0.02943],
		'TF4': [0.00298,0.0056],
		'TF8': [0.0003,0.01496]
		}
		
		glassbsdict = {
		#y_406 = 0.01741x + 0.03712
		# y_806 = 0.03444x + 0.07498
		# y_TF4 = 0.00287x + 0.02262
		# y_TF8 = 0.01873x + 0.09765
		#GBS dict as of 4/28/22
		'406': [0.01741,0.03712],
		'806': [0.03444,0.07498],
		'TF4': [0.00287,0.02262],
		'TF8': [0.01873,0.0976]
		}

		if polymer in ['EVA', 'POE']:
			print('Incorrect input, specify type of EVA or POE')
			return
		
		if architecture == 'GPOLYG':
			x = (ratio - glassglassdict[polymer][1])/glassglassdict[polymer][0]
			return x
		elif architecture == 'GPOLYBS':
			x = (ratio - glassbsdict[polymer][1])/glassbsdict[polymer][0]
			return x
		else:
			return 'fuck you'
			

from tqdm import tqdm
def bscorrect_linescans(df,alt_cal = False, plot = False):
	#currently using mini bff
	correctedAbs_dict = {}
	baselineAbs_dict = {}
	wardRatio_dict = {}
	wardWater_dict = {}

	for _,row in tqdm(df.iterrows()):
		wl = row['wl']
		idx_ch2 = np.where(wl == 1730)[0][0]
		idx_h2o = np.where(wl == 1902)[0][0]
		name = row['label']
		correctedAbs_dict[name] = []
		baselineAbs_dict[name] = []
		wardRatio_dict[name] = []
		wardWater_dict[name] = []
		for e,ii in enumerate(row['r']):
			corrected, baseline = mini_bff(row['wl'], ii, row['polymer'], name, plot)
			correctedAbs_dict[name].append(corrected)
			baselineAbs_dict[name].append(baseline)
			wardRatio_dict[name].append(correctedAbs_dict[name][e][idx_h2o]/correctedAbs_dict[name][e][idx_ch2])
			wardWater_dict[name].append(water_calibration(
								correctedAbs_dict[name][e][idx_h2o]/correctedAbs_dict[name][e][idx_ch2],
								row['architecture'],
								row['polymer_type'],
								alt_cal
								)
								)
	return correctedAbs_dict, baselineAbs_dict, wardRatio_dict, wardWater_dict
	
#correctedAbs_dict,baselineAbs_dict,wardRatio_dict = tala.bscorrect_linescans(df)
# df['corrected_abs_avg'] = df.apply(lambda x: np.mean(correctedAbs_dict[x.name],axis = 1),axis = 1)
# df['baseline_abs_avg'] = df.apply(lambda x: np.mean(baselineAbs_dict[x.name],axis = 1), axis = 1)
# # df['gaussian_fit_avg'] = df.apply(lambda x: np.mean(gaussianfit_dict[x.name],axis = 1), axis = 1)
# df['WaRD Water'] = df.apply(lambda x: wardRatio_dict[x.name], axis = 1)


def getexpectedwater(t,rh,polymer_type):
	inv_t = 1/(t + 273.15)
	k = 8.617333e-5 #eV

	solubilitydict = {
	# '406': [10.01530793, 0.22649281],
	# '806': [70.55601517,  0.27894601],
	# 'TF4': [3.17877056, 0.20111016],
	# 'TF8': [115.04621941,   0.30469101]
	# }
	#new calibration below, as of 4/24/22
	'406': [10.01530793,0.22649281],
	'806': [70.55601517,0.27894601],
	'TF4': [0.48172078, 0.17849866],
	'TF8': [0.60779957, 0.18669524]
	}
	c0 = solubilitydict[polymer_type][0]
	dH = solubilitydict[polymer_type][1]

	c_h2o = c0*np.exp(-dH*inv_t/k) * (rh/100) #g/cm3
	return c_h2o*1000 #mg/cm3

def new_test():
	return "Hey! New Testing! Please work! If works, then make sure to commit changes and then restart VS Code, not just re-import :)"



def fit_threepointbifacial(wavelengths, reflectance, architecture, plot = False):
	wl_eva = 1730
	wl_h2o = 1902
	wl_ref = 1872

	if np.mean(reflectance) > 1:
		reflectance = reflectance / 100
	ab = -np.log(reflectance)	#convert reflectance values to absorbance

	allWavelengthsPresent = True
	missingWavelength = None
	for each in [wl_eva, wl_h2o, wl_ref]:
		if each not in wavelengths:
			allWavelengthsPresent = False
			missingWavelength = each
			break

	if not allWavelengthsPresent:
		print('Wavelength Error: Necessary wavelength {0} missing from dataset - cannot fit.'.format(missingWavelength))
		return

	evaIdx = np.where(wavelengths == wl_eva)[0][0]
	h2oIdx = np.where(wavelengths == wl_h2o)[0][0]
	refIdx = np.where(wavelengths == wl_ref)[0][0]
	
	ratio = np.divide(ab[h2oIdx]-ab[refIdx], ab[evaIdx]-ab[refIdx])

	h2o_meas = water_calibration(
					ratio,
					architecture,
					'406'
					)
	# h2o[h2o < 0] = 0	
	## Avg Reflectance Fitting
	# avgRef = np.mean(ref, axis = 2)
	return h2o_meas