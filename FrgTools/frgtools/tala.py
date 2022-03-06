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
	#         y_406 = 0.02x + 0.016
	#         y_806 = 0.015x + 0.03
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
	wl = df['wl'][0]

	idx_ch2 = np.where(wl == 1730)[0][0]
	idx_h2o = np.where(wl == 1902)[0][0]

	for _,row in tqdm(df.iterrows()):
		correctedAbs_dict[row['label']] = []
		baselineAbs_dict[row['label']] = []
		wardRatio_dict[row['label']] = []
		wardWater_dict[row['label']] = []
		for e,ii in enumerate(row['r']):
			corrected, baseline = mini_bff(row['wl'], ii, row['polymer'], row['label'], plot)
			correctedAbs_dict[row['label']].append(corrected)
			baselineAbs_dict[row['label']].append(baseline)
		wardRatio_dict[e].append(correctedAbs_dict[row['label']][e][idx_h2o]/correctedAbs_dict[row['label']][e][idx_ch2])
		wardWater_dict[e].append(water_calibration(
							correctedAbs_dict[row['label']][e][idx_h2o]/correctedAbs_dict[row['label']][e][idx_ch2],
							row['architecture'],
							row['polymer'],
							alt_cal
							)
							)
	return correctedAbs_dict, baselineAbs_dict, wardRatio_dict, wardWater_dict
	

	# 	for pos in range(len(df['r'][0])):
	# 		wardRatio_dict[i].append(correctedAbs_dict[i][pos][idx_h2o]/correctedAbs_dict[i][pos][idx_ch2])
	# 		wardWater_dict[i].append(water_calibration(
	# 							correctedAbs_dict[i][pos][idx_h2o]/correctedAbs_dict[i][pos][idx_ch2],
	# 							df['architecture'][0],
	# 							pol,
	# 							alt_cal
	# 							)
	# 							)
	# return correctedAbs_dict, baselineAbs_dict, wardRatio_dict, wardWater_dict

#correctedAbs_dict,baselineAbs_dict,wardRatio_dict = tala.bscorrect_linescans(df)
# df['corrected_abs_avg'] = df.apply(lambda x: np.mean(correctedAbs_dict[x.name],axis = 1),axis = 1)
# df['baseline_abs_avg'] = df.apply(lambda x: np.mean(baselineAbs_dict[x.name],axis = 1), axis = 1)
# # df['gaussian_fit_avg'] = df.apply(lambda x: np.mean(gaussianfit_dict[x.name],axis = 1), axis = 1)
# df['WaRD Water'] = df.apply(lambda x: wardRatio_dict[x.name], axis = 1)


def getexpectedwater(t,rh,polymer_type):
	inv_t = 1/(t + 273.15)
	k = 8.617333e-5 #eV

	solubilitydict = {
	'406': [0.013, 0.022],
	'806': [1.174, 0.148],
	'TF4': [0.518, 0.167],
	'TF8': [2916.762, 0.419]
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