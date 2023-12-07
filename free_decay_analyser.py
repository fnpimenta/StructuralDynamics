import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from copy import deepcopy
from scipy.integrate import odeint
from scipy.signal import find_peaks
from scipy import signal
from datetime import datetime

from eq_motion import *
from estimators import *
from plot_generators import * 

from PIL import Image

# -- Set page config
apptitle = 'Free decay analyser'
icon = Image.open("logo.ico")
st.set_page_config(page_title=apptitle, page_icon=icon)

# Title the app
st.title('Free decay analysis for damping estimation')

st.markdown("""
 * Use the menu at left to select data from the different analysis possibilities
 * To tune the analysis parameters use the **Analysis** tab
""")

# -- Side bar definition
tab1, tab2 = st.sidebar.tabs(["🌊 Simulation", "📈 Analysis"])

# -- File type definition (tab1)
analysis = ['Numerical simulation','OpenFAST repo','Upload file']
page = tab1.selectbox('Input data file', analysis)

if page == 'Numerical simulation':
	# Define initial amplitude
	sim_t = tab1.number_input('Simulation time (s)',min_value=1, max_value=None, value=100)
	fs = tab1.number_input('Sample frequency (Hz)',min_value=1, max_value=None, value=25)
	A = tab1.number_input('Initial amplitude (m)', 1.0)  # min, max, default

	# Define damping ratios
	col1, col2 = tab1.columns(2)
	xi_l = col1.number_input(r'$\xi_1$ (%)', 0, 50, 5)  # min, max, default
	xi_nl = col2.number_input(r'$\xi_2$ (%)', 0, 25, 0)  # min, max, default
	u_ref = tab1.number_input('Ambient speed (m/s)', 0.0, 10.0, 0.0)  # min, max, default

	t = np.arange(0,sim_t+1/fs,1/fs)
	y0 = [A,0]

	y = odeint(eq_motion, y0, t, args=(xi_l/100,xi_nl/100,1,u_ref))[:,0]
	error_check = 0

elif page == 'OpenFAST repo':
	dofs = {
			"Surge": "PtfmTDxi",
			"Sway": "PtfmTDyi",
			"Heave": "PtfmTDzi",
			"Roll": "PtfmRDxi",
			"Pitch": "PtfmRDyi",
			"Yaw": "PtfmRDzi",
			"Tower FA": "YawBrTAxp",
			"Tower SS": "YawBrTAyp"
			}   
	dof = tab1.selectbox('Degree of freedom', dofs.keys(),index=4)

	currents = {
				"0": "0000",
				"0.25": "0250",
				"0.50": "0500",
				"0.75": "0750",
				"1.00": "1000"
				}   
	current = tab1.selectbox('Water current speed (m/s)', currents.keys())
	winds = {
			"0": "0000",
			"3": "0300",
			"8": "0800",
			"11.4": "1140",
			"22": "2200",
			"30": "3000"
			}   
	wind = tab1.selectbox('Wind speed (m/s)', winds.keys(),index=2)

	if (dof=="Tower FA") | (dof=="Tower SS"):
		filename = 'NREL5MW_W%s_C%s_%s.out'%(winds[wind],currents[current],dof[-2:])
	else:
		filename = 'NREL5MW_W%s_C%s_%s.out'%(winds[wind],currents[current],dof)

	try:
		Filespath = 'OpenFAST_simulations'
		data = pd.read_csv(Filespath + '//' + filename,delimiter=",",encoding_errors='replace')
		t = np.array(data['Time'])
		y = np.array(data[dofs[dof]])
		fs = 1/t[1]
		error_check = 0
		tab1.write('Loaded file %s from the available database'%filename)
		
	except:
		tab1.write('No file found for that environmental conditions')
		t = np.linspace(0,1,100)
		fs = 50
		error_check = 1
else:
	uploaded_file = tab1.file_uploader("Choose a file")
	ftype = tab1.radio('File type',['OpenFAST','Android','iOS'],horizontal=True)
	try:
		if ftype=='Android':
			header_row = tab1.number_input('Header row', value=2) - 1
			first_row = tab1.number_input('First data row', value=3) - 1
		
			skip_rows = [int(x) for x in range(first_row) if (x!=header_row)]
			if header_row<0:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,delimiter=';',encoding_errors='replace',header=None)
			else:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,delimiter=';',encoding_errors='replace')	
			
			data = data.replace(',','.', regex=True)
			data = data.astype(float)
			data.columns = ['Steps','Time','X','Y','Z']
			fi = 1000/(np.mean(np.array(data['Time'])[2:-1] - np.array(data['Time'])[1:-2]))
			dof = tab1.selectbox('Data column', data.columns,index=min(np.shape(data)[1],2))

		elif ftype=='iOS':
			header_row = tab1.number_input('Header row', value=1) - 1
			first_row = tab1.number_input('First data row', value=2) - 1
		
			skip_rows = [int(x) for x in range(first_row) if (x!=header_row)]
			if header_row<0:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,encoding_errors='replace',header=None)
			else:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,encoding_errors='replace')	
			data.columns = ['Date','X','Y','Z']
			dof = tab1.selectbox('Data column', data.columns,index=min(np.shape(data)[1],1))
			
			if int(data['Date'].iloc[0][-2:])==99:
				dt = int(data['Date'].iloc[2][-2:]) - int(data['Date'].iloc[1][-2:])
			else:
				dt = int(data['Date'].iloc[1][-2:]) - int(data['Date'].iloc[0][-2:])
			fi=100/dt

		else:
			header_row = tab1.number_input('Header row', value=1) - 1
			first_row = tab1.number_input('First data row', value=2) - 1
		
			skip_rows = [int(x) for x in range(first_row) if (x!=header_row)]
			if header_row<0:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,encoding_errors='replace',header=None)
			else:
				data = pd.read_csv(uploaded_file,skiprows=skip_rows,encoding_errors='replace')	
			fi=10

		y = np.array(data[dof])
		fs = tab1.number_input('Sampling rate', value=fi)
		t = np.linspace(0,1,len(y)) * len(y)/fs

		error_check = 0
	except:
		tab1.write('Please select the file for analysis')
		error_check = 1

# -- Analysis parameters
#try:
if error_check == 0:
	#t_min,t_max = tab2.slider('Period for analysis (s)', min_value=0.0, max_value=t[-1], value=[0.0,float(t[-1])])  # min, max, default
	t_min = tab2.number_input('First point for the analysis',min_value=0.0,max_value=None,value=0.0)
	t_max = tab2.number_input('Last point for the analysis',min_value=0.0,max_value=None,value=t[-1])
	
	f_max = tab2.slider('Maximum frequency (Hz)', 0.0, fs/2, float(fs/4))  # min, max, default

	npoints = len(t[t<=t_max])
	nfft = tab2.slider('n for PSD computation ($2^n$)', int(np.log2(npoints)/4), int(np.log2(npoints)), int(np.log2(npoints)))  # min, max, default

	peaks_types = {
				 "Positive peaks only": 0,
				 "All peaks": 1,
				}

	filt_app = tab2.checkbox('Apply filter') 
	peaks_type = tab2.radio('Peaks to use', peaks_types.keys(),horizontal=True)

	

	if filt_app == 1:
		filt_type = tab2.selectbox('Filter type', ['Low-pass','High-pass','Band-pass'],index=0)
		filt_order = tab2.slider('Filter order',4,12,8)
		col1_t2, col2_t2 = tab2.columns(2)
		
		if filt_type == 'Low-pass':
			fmin = col1_t2.number_input('Lower limit of the filter', min_value=0.0, max_value=f_max, value=0.0,disabled=True)  # min, max, default
			fmax = col2_t2.number_input('Upper limit of the filter', min_value=0.0, max_value=f_max, value=float(f_max/4))  # min, max, default
			sos = signal.butter(filt_order  , fmax, 'lowpass' , fs=fs , output='sos', analog=False)

		elif filt_type == 'High-pass':
			fmin = col1_t2.number_input('Lower limit of the filter', 0.0, f_max, value=float(f_max/4))  # min, max, default
			fmax = col2_t2.number_input('Upper limit of the filter', 0.0, f_max, value=float(f_max/2),disabled=True)  # min, max, default
			sos = signal.butter(filt_order  , fmin, 'highpass' , fs=fs , output='sos', analog=False)

		elif filt_type == 'Band-pass':
			fmin = col1_t2.number_input('Lower limit of the filter', 0.0, f_max, value=float(f_max/8))  # min, max, default
			fmax = col2_t2.number_input('Upper limit of the filter', 0.0, f_max, value=float(f_max/4))  # min, max, default
			sos = signal.butter(filt_order  , [fmin,fmax], 'bandpass' , fs=fs , output='sos', analog=False)

		y_doubled = np.zeros(2*len(y)-1)

		# Double the time series to avoid filter impact at the beginning
		y_doubled[len(y)-1:] = y
		y_doubled[:len(y)-1] = y[:0:-1]
		y_filt = signal.sosfiltfilt(sos, y_doubled)[len(y)-1:]

	else:
		y_filt = y
		sos = 0

	# Define the time series limits for analysis
	time_filter = (t>=t_min) & (t<=t_max)

	# Compute the power spectrum of the time series for analysis (for representation purposes only)
	f, Pxx = signal.welch(y[time_filter], fs, nperseg=2**nfft , scaling='spectrum') 
	f, Pxx_filt = signal.welch(y_filt[time_filter], fs, nperseg=2**nfft , scaling='spectrum') 
	
	# Estimate the offset and zero the time series
	offset = offset_estimator(y_filt[time_filter])
	y_zerod = y_filt - offset

	# Estimate the peaks of the zeroed time series
	peaks_time, peaks_amp = peaks_estimator(t[time_filter],y_zerod[time_filter],peaks_types[peaks_type],t_min==0)
	
	use_all_points = tab2.checkbox('Fit all points simultaneously',value=1)
	window_size = 2*tab2.number_input('Window size for damping estimate',min_value=1,max_value=None,value=1,disabled=use_all_points)
	if use_all_points:
		window_size = len(peaks_time)

	# Estimate of the dynamic properties of the free decay
	xi_est, f_est = dynamic_estimator(peaks_time,peaks_amp,peaks_type=peaks_types[peaks_type],window_size=window_size)
	if use_all_points:
		fit_sug=1
	else:
		fit_sug = fit_guess(abs(peaks_amp[int(window_size/2):-int(window_size/2)]),xi_est)
	
	fit_types = {
				 "Trend line": 0,
				 "Mean value": 1
				}

	fit_type = tab2.radio('Fit type', fit_types.keys(),horizontal=True,index=fit_sug,disabled=use_all_points)

	# Create the figure to plot
	if use_all_points:
		fig = free_decay_plot_all(t,y,y_filt,offset,
								  peaks_time,peaks_amp,
								  xi_est,f_est,
								  f,Pxx, Pxx_filt,
								  fs,sos,f_max,filt_app,
								  time_filter)
	else:
		fig = free_decay_plot(t,y,y_filt,offset,
							  peaks_time,peaks_amp,
							  xi_est,f_est,
							  f,Pxx, Pxx_filt,
							  fs,sos,f_max,filt_app,
							  time_filter,
							  fit_types[fit_type],
							  window_size=window_size)

	st.pyplot(fig) 
else:
	st.write('No data available for the analysis. Please use the simulation menu to generate a time series or use an input file.')
#except:
#	st.error('Something went wrong. Check the simulation conditions or click below to try to run the simulation again.')
#	st.button('Run simulation again.')

with st.expander("See explanation",False):
	st.write(r'''
		The free equation of motion (without any external load) is given by:
		$$
			m\ddot{x}(t) = -kx(t) - f_d(t)
		$$
		where $f_d(t)$ is the damping force, here taken to be a generic function of time. 
		In the simplified numerical simulation above, the damping force is assumed to be such that:
		$$
		   f_d(t) = c_1x(t) + c_2|\dot{x}(t)-u_r|(\dot{x}(t)-u_r)
		$$
		where $u_r$ is external flow velocity. 

		For the linear damping model ($c_2=0$), the differential equation above has the well known solution:
		$$
			x(t) = Ae^{-\xi\omega_0 t}\cos(w\sqrt{1-\xi^2}t+\phi) = Ae^{-\xi\omega_0 t}\cos(w_dt+\phi)
		$$
		where $\omega_0=\sqrt{\frac{k}{m}}$ is the system undamped natural frequency and $\xi$ is the damping ratio, 
		defined as the ratio between the damping coefficient and its critical value as:
		$$
			\xi=\frac{c_1}{c_{cr}}=\frac{c_1}{2m\omega_0}
		$$
		One may immediately see that the response is given by a periodic function modulated by a negative exponential,
		implying that $\xi$ can be evaluated through the response amplitude, since:
		$$
			\ln(Ae^{-\xi\omega_0 t}) = \ln(A) - \xi\omega_0t
		$$
		meaning that the envelope amplitde natural logarithm, here obtained through the peak value in every oscillation, is a linear function of time with slope $-\xi\omega_0$.

		Although this is no longer true if $c_2\neq0$, for low damping forces, one may still make some general considerations based on energy dissipation.
		Firslty, it should be noted that the energy dissipation over a full cyle may be written as:
		$$
			W = \int_Tf_d(t)dx
		$$
		For the linear damping contribution, one finds:
		$$
			W_l = c_1\int_T\dot{x}(t)dx = c_1\int_T\dot{x}^2dt \approx  c_1A^2\omega_0^2\int_T\sin^2(\omega_0 t)dt = c_1\left(A^2\omega_0\pi\right)
		$$
		where it was assumed that to first order the motion may be approximated over a cycle as $x(t)=A\cos(\omega_0t)$.
		Under the same assumption, the quadratic contribution, for $u_r=0$, may be obtained as:
		$$
			W_q = c_2\int_T|\dot{x}(t)|\dot{x}(t)dx = 2c_2\int_{T/2}\dot{x}^3dt \approx 2c_2A^3\omega_0^3\int_{T/2}\sin^3(\omega_0t)dt = c_2 \frac{8A^3\omega_0^2}{3}
		$$
		By comparison with the linear damping result, it follows that the linear coefficient that best approximates the quadratic response in terms of energy dissipation is:
		$$
			 \tilde{c} = c_2 \frac{8\omega_0}{3\pi} A 
		$$
		From the expression above, it can be seen that for a purely quadratic damping force, a linear dependency on the motion amplitude is expected.
	''')