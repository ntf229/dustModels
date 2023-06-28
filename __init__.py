# python version of the dust models as incorporated in FSPS 
# for all two-conponent dust models, the young stars are attenuated by a power-law

import numpy as np

def applyAttenuation(spectrum, A, fracNoDust=0):
	# applies attenuation A to spectrum
	return fracNoDust*spectrum + (1-fracNoDust)*spectrum*10**(-0.4*A)

def convertToMags(spectrum):
	# converts spectrum from Jansky to Pogson magnitudes
	return 22.5 - 2.5*np.log10((spectrum/3631) * 1e9)

def calzetti(wavelengths, spectrum, Av, fracNoDust):
	"""
	- old and young stars treated together
	- wavelengths in Angstroms
	- spectrum in Jansky
	- Av (attenuation in the V-band) in magnitudes
	- fracNoDust is the fraction of stars that don't get attenuated
	"""
	waveMask1 = wavelengths/1e4 <= 0.63
	waveMask2 = (wavelengths/1e4 >= 0.63) & (wavelengths/1e4 <= 2.20)
	Rv = 4.05
	k1 = 2.659*(-2.156 + (1.509/(wavelengths[waveMask1]/1e4)) - (
		 0.198/(wavelengths[waveMask1]/1e4)**2) + (
		 0.011/(wavelengths[waveMask1]/1e4)**3 )) + Rv
	k2 = 2.659*(-1.857 + (1.040/(wavelengths[waveMask2]/1e4))) + Rv
	k = np.append(k1, k2)
	A = k*Av/Rv # attenuation curve
	dustySpec = applyAttenuation(spectrum, A, fracNoDust)
	dustMags = convertToMags(dustySpec)
	noDustMags = convertToMags(spectrum)
	attenuationMags = dustMags - noDustMags
	return dustySpec, attenuationMags

def cardelli(wavelengths, youngSpec, oldSpec, Av, mwr, uvb, AvYoung, 
			 dustIndexYoung, fracNoDust, fracNoDustYoung):
	""" 
	- MW extinction curve R = Av / E(B-V)
	- old and young stars treated separate
	- wavelengths in Angstroms
	- youngSpec and oldSpec in Jansky
	- Av (attenuation in the V-band) in magnitudes, applies to young and old stars
	- mwr is the ratio of total to selective absorption which characterizes the MW extinction curve
	- uvb is the 2175A UV-bump strength
	- AvYoung is the Av applied to young stars
	- dustIndexYoung is the power law index for the young stars
	- fracNoDust is the fraction of yound and old stars that don't get attenuated
	- fracNoDustYoung is the fraction of young stars that don't get attenuated
	"""
	youngSpec = youngPowerLaw(wavelengths, youngSpec, AvYoung, dustIndexYoung, fracNoDustYoung)
	# combine old stars with attenuated young stars
	compositeSpec = youngSpec + oldSpec
	x = 1e4/wavelengths # x is in decreasing order
	y = x-1.82
	a = np.zeros(len(x))
	b = np.zeros(len(x))
	F_a = np.zeros(len(x))
	F_b = np.zeros(len(x))
	A = np.zeros(len(x)) # attenuation curve
	# wavelength indices
	mwdindex = np.zeros(6, dtype=int)
	mwdindex[0] = np.where(x == x[x >= 0.1][-1])[0][0]
	mwdindex[1] = np.where(x == x[x >= 1.1][-1])[0][0]
	mwdindex[2] = np.where(x == x[x >= 3.3][-1])[0][0]
	mwdindex[3] = np.where(x == x[x >= 5.9][-1])[0][0]
	mwdindex[4] = np.where(x == x[x >= 8.][-1])[0][0]
	if np.amax(x) >= 12:
		mwdindex[5] = np.where(x == x[x >= 12.][-1])[0][0] 
	else:
		mwdindex[5] = np.where(x == np.amax(x))[0][0]
	maskIR = slice(mwdindex[1], mwdindex[0])
	maskOpticalNIR = slice(mwdindex[2], mwdindex[1])
	maskNUV = slice(mwdindex[3], mwdindex[2])
	maskMUV = slice(mwdindex[4], mwdindex[3])
	maskFUV = slice(mwdindex[5], mwdindex[4])
	if np.amax(x) >= 12:
		maskExtendedUV = x >= 12. 
	# IR: 0.1 um^-1 <= x <= 1.1 um^-1 (paper only goes down to 0.3 um^-1, fsps goes to 0.1 um^-1)
	a[maskIR] = 0.574*x[maskIR]**(1.61)
	b[maskIR] = -0.527*x[maskIR]**(1.61)
	A[maskIR] = a[maskIR] + b[maskIR]/mwr
	# optical+near-IR: 1.1 um^-1 <= x <= 3.3 um^-1
	a[maskOpticalNIR] = 1 + 0.17699*y[maskOpticalNIR] - 0.50447*y[maskOpticalNIR]**2 - 0.02427*y[maskOpticalNIR]**3 + 0.72085*y[maskOpticalNIR]**4 + 0.01979*y[maskOpticalNIR]**5 - 0.77530*y[maskOpticalNIR]**6 + 0.32999*y[maskOpticalNIR]**7
	b[maskOpticalNIR] = 1.41338*y[maskOpticalNIR] + 2.28305*y[maskOpticalNIR]**2 + 1.07233*y[maskOpticalNIR]**3 - 5.38434*y[maskOpticalNIR]**4 - 0.62251*y[maskOpticalNIR]**5 + 5.30260*y[maskOpticalNIR]**6 - 2.09002*y[maskOpticalNIR]**7
	A[maskOpticalNIR] = a[maskOpticalNIR] + b[maskOpticalNIR]/mwr
	# near-UV: 3.3 um^-1 <= x <= 5.9 um^-1 (with fsps hack)
	a[maskNUV] = 1.752 - 0.316*x[maskNUV] - 0.104/((x[maskNUV] - 4.67)**2 + 0.341)*uvb
	b[maskNUV] = -3.090 + 1.825*x[maskNUV] + 1.206/((x[maskNUV] - 4.62)**2 + 0.263)*uvb
	aTemp = 1.752 - 0.316*x - 0.104/((x - 4.67)**2 + 0.341)*uvb
	bTemp = -3.090 + 1.825*x + 1.206/((x - 4.62)**2 + 0.263)*uvb
	hack = (x[mwdindex[2]]/x)**6 * (A[mwdindex[2]] - (aTemp[mwdindex[2]] + bTemp[mwdindex[2]]/mwr))
	A[maskNUV] = a[maskNUV] + b[maskNUV]/mwr + hack[maskNUV]
	# mid-UV: 5.9 um^-1 <= x <= 8. um^-1
	F_a[maskMUV] = -0.04473*(x[maskMUV] - 5.9)**2 - 0.009779*(x[maskMUV] - 5.9)**3
	F_b[maskMUV] = 0.2130*(x[maskMUV] - 5.9)**2 + 0.1207*(x[maskMUV] - 5.9)**3
	a[maskMUV] = 1.752 - 0.316*x[maskMUV] - 0.104/((x[maskMUV] - 4.67)**2 + 0.341)*uvb
	b[maskMUV] = -3.090 + 1.825*x[maskMUV] + 1.206/((x[maskMUV] - 4.62)**2 + 0.263)*uvb
	A[maskMUV] = (a[maskMUV]+F_a[maskMUV]) + (b[maskMUV]+F_b[maskMUV])/mwr
	# far-UV: 8. um^-1 <= x <= 12. um^-1 (paper only goes to 10 um^-1, fsps goes to 12 um^-1)
	a[maskFUV] = -1.073 - 0.628*(x[maskFUV] - 8) + 0.137*(x[maskFUV] - 8)**2 - 0.070*(x[maskFUV]-8)**3
	b[maskFUV] = 13.670 + 4.257*(x[maskFUV] - 8) - 0.420*(x[maskFUV] - 8)**2 + 0.374*(x[maskFUV] - 8)**3
	A[maskFUV] = a[maskFUV] + b[maskFUV]/mwr
	# extended UV: x >= 12. (fsps)
	if np.amax(x) >= 12:
		A[maskExtendedUV] = A[maskFUV][0]
	# total scaling
	A *= Av
	dustySpec = applyAttenuation(compositeSpec, A, fracNoDust)
	dustMags = convertToMags(dustySpec)
	noDustMags = convertToMags(compositeSpec)
	attenuationMags = dustMags - noDustMags
	return dustySpec, attenuationMags

def powerLaw(wavelengths, youngSpec, oldSpec, Av, dustIndex, AvYoung, 
			 dustIndexYoung, fracNoDust, fracNoDustYoung):
	""" 
	- old and young stars treated separate
	- wavelengths in Angstroms
	- youngSpec and oldSpec in Jansky
	- Av (attenuation in the V-band) in magnitudes, applies to young and old stars
	- dustIndex is the power law index for young and old stars
	- AvYoung is the Av applied to young stars
	- dustIndexYoung is the power law index for the young stars
	- fracNoDust is the fraction of yound and old stars that don't get attenuated
	- fracNoDustYoung is the fraction of young stars that don't get attenuated
	"""
	youngSpec = youngPowerLaw(wavelengths, youngSpec, AvYoung, dustIndexYoung, fracNoDustYoung)
	# combine old stars with attenuated young stars
	compositeSpec = youngSpec + oldSpec
	A = Av*(wavelengths/5500.)**dustIndex
	dustySpec = applyAttenuation(compositeSpec, A, fracNoDust)
	dustMags = convertToMags(dustySpec)
	noDustMags = convertToMags(compositeSpec)
	attenuationMags = dustMags - noDustMags
	return dustySpec, attenuationMags

def kriekAndConroy(wavelengths, youngSpec, oldSpec, Av, dustIndex, 
				   AvYoung, dustIndexYoung, fracNoDust, fracNoDustYoung):
	""" 
	- old and young stars treated separate
	- wavelengths in Angstroms
	- youngSpec and oldSpec in Jansky
	- Av (attenuation in the V-band) in magnitudes, applies to young and old stars
	- dustIndex is the power law index for young and old stars
	- AvYoung is the Av applied to young stars
	- dustIndexYoung is the power law index for the young stars
	- fracNoDust is the fraction of yound and old stars that don't get attenuated
	- fracNoDustYoung is the fraction of young stars that don't get attenuated
	"""
	youngSpec = youngPowerLaw(wavelengths, youngSpec, AvYoung, dustIndexYoung, fracNoDustYoung)
	compositeSpec = youngSpec + oldSpec
	eb = 0.85 - 1.9 * dustIndex
	Rv = 4.05
	waveMask1 = wavelengths/1e4 <= 0.63
	waveMask2 = (wavelengths/1e4 >= 0.63) & (wavelengths/1e4 <= 2.20)
	k1 = 2.659*(-2.156 + (1.509/(wavelengths[waveMask1]/1e4)) - (
		 0.198/(wavelengths[waveMask1]/1e4)**2) + (
		 0.011/(wavelengths[waveMask1]/1e4)**3 )) + Rv
	k2 = 2.659*(-1.857 + (1.040/(wavelengths[waveMask2]/1e4))) + Rv
	k = np.append(k1, k2)
	drude = eb*(wavelengths*350)**2 / ( (wavelengths**2-2175**2)**2 + (wavelengths*350)**2 )
	A = (Av/Rv)*(k+drude)*(wavelengths/5500)**dustIndex # attenuation curve
	# apply attenuation
	dustySpec = applyAttenuation(compositeSpec, A, fracNoDust)
	dustMags = convertToMags(dustySpec)
	noDustMags = convertToMags(compositeSpec)
	attenuationMags = dustMags - noDustMags
	return dustySpec, attenuationMags

def youngPowerLaw(wavelengths, spectrum, Av, dustIndex, fracNoDust):
	"""
	- this function is called from within the two-component models
	- it applies a power law exctinction curve to the young stars
	"""
	A = Av*(wavelengths/5500.)**dustIndex
	return fracNoDust*spectrum + (1-fracNoDust)*spectrum*10**(-0.4*A)

