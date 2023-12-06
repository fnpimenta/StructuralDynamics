import numpy as np
import streamlit as st
from scipy import signal
from scipy.signal import find_peaks

def LinearFitCoefficients(x,y):
	'''
	Computes the parameters for a linear fit and corresponding uncertainties.
	
	Parameters
	----------
	x : array_like 
		Independent variable time series.
	y : array_like 
		Dependent variable time series.
	
	Returns
	-------
	m : float
		Best fit slope.
	dm : float
		Best fit slope uncertainty.
	b : float
		Best fit offset for x=0
	db : float
		Best fit offset uncertainty
	'''
	n = len(x)
	
	Sx = np.sum(x)
	Sy = np.sum(y)
	Sxx = np.sum(x**2)
	Syy = np.sum(y**2)
	Sxy = np.sum(x*y)
	m = (n*Sxy-Sx*Sy)/(n*Sxx-Sx**2)
	b = Sy/n - m/n*Sx
	
	Se2 = 1/(n*(n-2)) * (n*Syy - Sy**2 - m**2*(n*Sxx-Sx**2))
	Sm2 = n*Se2/(n*Sxx-Sx**2)
	Sb2 = 1/n * Sm2 * Sxx
	
	dm = np.sqrt(Sm2)
	db = np.sqrt(Sb2)
	
	return m,dm,b,db

def offset_estimator(y):
	'''
	Estimates the time series offset.
	Only valid for well defined free decay simulations and not arbitrary time series. 
	
	The mean value is estimated from the expected zero crossings at 1/4 and 3/4 between consecutive peaks.
	Parameters
	----------
	y : array_like 
		Time series to be zeroed.
	
	Returns
	-------
	offset : float
		Estimated time series offset.

	'''
	
	all_peaks, _ = find_peaks(y, prominence=0)

	t_peaks = [int(x) for x in (3*all_peaks[:-1]/4 + all_peaks[1:]/4)]
	[t_peaks.append(int(x)) for x in (all_peaks[:-1]/4 + 3*all_peaks[1:]/4)]
	
	offset = y[t_peaks].mean()
	
	return offset


def peaks_estimator(t,y,peaks_type=0,add_1P=1,prominence=0):
	'''
	Estimates the time series peaks amplitude and location.
	
	Parameters
	----------
	t : array_like 
		Time array.
	y : array_like 
		Time series array.

	add_1P : flag 
		If add_1P the first point of the series is added to the peaks list.
		Optional (default to 1)
	prominence : array_like 
		Minimum peak porminence to consider.
		Optional (default to 0)
	peaks_type : flag 
		0 to identify only positive amplitude peaks.
		1 ot identify all peaks (positive and negative amplitude).
		Optional (default to 0)
	
	Returns
	-------
	peaks_time : array_like
		Time instant of the identified peaks.
	peaks_amp : array_like
		Time series amplitude at the identified peaks.

	'''
	if peaks_type == 0:
		peaks, _ = find_peaks(y, prominence=prominence)
	else:
		peaks, _ = find_peaks(abs(y), prominence=prominence)

	if add_1P==1:
		peaks = np.concatenate(([0],peaks))

	peaks_time = t[peaks]
	peaks_amp = y[peaks]

	return peaks_time, peaks_amp

def dynamic_estimator(peaks_time,peaks_amp,peaks_type=0,window_size=2):
	'''
	Estimates the damping ratio and frequency compatible with the free decay.
	The estimates are made considering a moving window of 1 period.
	
	Parameters
	----------
	peaks_time : array_like
		Time instant of the identified peaks.
	peaks_amp : array_like
		Time series amplitude at the identified peaks.
	peaks_type : flag 
		0 if only positive amplitude peaks have been considered
		1 if all peaks (positive and negative amplitude) have been considered.
		Optional (default to 0)
	
	Returns
	-------
	xi_est : array_like
		Estimated damping coefficient (in %).
	f_est : array_like
		Estimated frequency (in Hz).	

	'''
	n_peaks = len(peaks_amp)

	# Dynamic properties estimates
	xi_est = np.zeros(n_peaks-window_size)
	f_est = np.zeros(n_peaks-window_size)

	for i in range(int(window_size/2),n_peaks-int(window_size/2)):
		xi_est[i-int(window_size/2)],_,_,_ = LinearFitCoefficients(peaks_time[i-int(window_size/2):i+int(window_size/2)+1],np.log(abs(peaks_amp[i-int(window_size/2):i+int(window_size/2)+1])))
		f_est[i-int(window_size/2)] = int(window_size/2)/(peaks_time[i+int(window_size/2)]-peaks_time[i-int(window_size/2)]) * (1 + 1*(peaks_type==0))
	
	xi_est = -100*xi_est/(2*np.pi*f_est)

	return xi_est, f_est

def fit_guess(x,y):
	'''
	Compares a linear fit with a mean value fit.
	
	Parameters
	----------
	x : array_like
		Series of independent variable.
	y : array_like
		Series of dependent variable, y=y(x).
	
	Returns
	-------
	flag : int
		Returns 1 for a mean value fit or 0 for a linear fit.	

	'''
	m,dm,b,db = LinearFitCoefficients(x,y)
	error = np.sqrt(np.sum((y - (x*m+b))**2))
	if error>np.std(y):
		return 1
	else:
		return 0