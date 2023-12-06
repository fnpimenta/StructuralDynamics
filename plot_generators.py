import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import signal

from estimators import LinearFitCoefficients


def free_decay_plot(t,y,y_filt,offset,
					peaks_time,peaks_amp,
					xi_est,f_est,
					f, Pxx, Pxx_filt,
					fs,sos,f_max,applied_filter,
					time_filter,	
					fit_type,window_size=2):
	'''
	Produces the free decay analysis plots:
	- Time series and filtered time series plots (upper left)
	- Peaks as a function of time log value (centre left)
	- Estimated damping ratio as a function of time (lower left)
	- Original, filtered and filter frequency response (upper right)
	- Estimated damping ratio and frequency as a function of motion amplitude (lower right)
	
	Parameters
	----------
	t : array_like
		Time series.
	y : array_like
		Original time series.
	y_filt : array_like
		Filtered time series.	
	offset : float
		Estimated offset value
	peaks_time : array_like
		Peaks time location
	peaks_amp : array_like
		Peaks amplitude
	xi_est : array_like
		Damping ratio estimates
	f_est : array_like
		Frequency estimates
	f : array_like
		Frequencies for power spectrum plot
	Pxx : array_like
		Original time series power spectrum
	Pxx_filt : array_like
		Filtered time series power spectrum
	fs : float
		Sample frequency
	sos : 
		sos filter output
	f_max :	float
		Maximum frequency for power spectrum plot
	applied_filter : flag
		0 if no filter has been applied. 1 otherwise
	time_filter : tuple
		Points to plot
	fit_type : flag
		0 for a linear fit characterisation of the damping ratio as a function of motion amplitude
		1 for the mean value and standard deviation of the damping ratio estimates	
	
	Returns
	-------
	fig : figure
		Figure with the plots used for the free decay analysis	

	'''
	
	# Define the upper limit for the plot
	t_max = t[time_filter][-1]

	# Figure definition
	fig = plt.figure(figsize = (12,8))

	gs = gridspec.GridSpec(3,2,hspace=0,wspace=0.05)
	gs_f = gs[0:2,1].subgridspec(5,1,hspace=0)
	
	ax1 = plt.subplot(gs[0,0])      # Time series plot
	ax2 = plt.subplot(gs[1,0])      # Peak location in log scale
	ax3 = plt.subplot(gs[2,0])      # Damping coefficient estimate as a function of time
	ax4 = plt.subplot(gs_f[0:-1,0]) # Power spectrum
	ax4_f = ax4.twinx()             # Power spectrum for filter response plot
	ax5 = plt.subplot(gs[2,1])      # Damping coefficient as a function of amplitude
	ax5_f = ax5.twinx()             # Frequency estimates plot 


	# Original time series
	ax1.plot(t,y,'blue',linewidth=1) # Time series plot
	ax4.semilogy(f,Pxx,'b',label='Original signal',linewidth=1) # Power spectra plots
	pmin,pmax = ax4.get_ylim() # Original power spectra plot limits

	# Filtered time series
	if applied_filter == 1:
		# Filter response
		w, h = signal.sosfreqz(sos,fs=fs,worN=np.linspace(0,f_max,100))
		ax4_f.plot(w, abs(h),'--k') 
		
		# Filtered time series
		ax1.plot(t[time_filter],y_filt[time_filter],'red',linewidth=1) # Time series plot 
		ax4.semilogy(f,Pxx_filt,'r',label='Filtered signal',linewidth=1) # Power spectra plots
		
	# Plot the estimate offset
	ax1.axhline(offset,0,1,c='k',ls='--')

	# Plot selected peaks
	ax1.plot(peaks_time,peaks_amp+offset,'or',ms=4)
	ax2.plot(peaks_time,np.log(abs(peaks_amp)),'or',ms=4)
	
	# Plot estimated dynamic properties

	ax3.plot(peaks_time[int(window_size/2):-int(window_size/2)],xi_est,'or',ms=4)
	ax5.plot(peaks_amp[int(window_size/2):-int(window_size/2)],xi_est,'or',ms=4)
	ax5_f.plot(peaks_amp[int(window_size/2):-int(window_size/2)],f_est,'og',ms=4)

	xmin,xmax = 0,ax5.get_xlim()[1]
	xbounds = np.array([xmin,xmax])

	ax4_f.axvline(f_est.mean(),0,1,c='g',ls='--') 
	ax4_f.fill_between([f_est.mean()-f_est.std(),f_est.mean()+f_est.std()],0,1.2,
					   color='g',alpha=0.1) 

	f_title = r'$f=%.3f^{+%.3f}_{-%.3f}$'%(f_est.mean(),f_est.std(),f_est.std()) + ' (Hz)'

	if fit_type == 1:
		ax5.axhline(xi_est.mean(),0,1,c='r',ls='--') 
		ax5.fill_between(xbounds,xi_est.mean()-xi_est.std(),xi_est.mean()+xi_est.std(),
						 color='r',alpha=0.1) 
		xi_title = r'$\xi=%.2f^{+%.2f}_{-%.2f}$'%(xi_est.mean(),xi_est.std(),xi_est.std()) + ' (%)'
	else:
		m,dm,b,db = LinearFitCoefficients(abs(peaks_amp[1:-1]),xi_est)
		ax5.plot(xbounds,xbounds*m+b,'--r')
		ax5.fill_between(xbounds,xbounds*(m-dm)+b-db,xbounds*(m+dm)+b+db,
						 color='r',alpha=0.1) 
		if m<0:
			xi_title = r'$\xi=%.2f^{+%.2f}_{-%.2f}-%.2f^{+%.2f}_{-%.2f}\times A$'%(b,db,db,-m,dm,dm) + ' (%)'
		else:
			xi_title = r'$\xi=%.2f^{+%.2f}_{-%.2f}+%.2f^{+%.2f}_{-%.2f}\times A$'%(b,db,db,m,dm,dm) + ' (%)'
	
	ax5_f.axhline(f_est.mean(),0,1,c='g',ls='--') 
	ax5_f.fill_between(xbounds,f_est.mean()-f_est.std(),f_est.mean()+f_est.std(),
					   color='g',alpha=0.1)        

	# Axis definitions
	ax1.set_ylabel('Motion amplitude')
	ax1.set_xlim(0,t_max)
	ax1.set_xticklabels('')

	ax2.set_ylabel('Peak amplitude (log scale)')
	ax2.set_xlim(0,t_max)
	ax2.set_xticklabels('')

	ax3.set_xlim(0,t_max)
	ax3.set_ylim(0,xi_est.max()*1.1)
	ax3.set_xlabel('Time (s)')
	ax3.set_ylabel(r'$\xi$ (%)',color='red')
	ax3.set_xticks(ax3.get_xticks()[:-1])

	ax3_f = ax3.twinx()
	ax3_f.set_ylim(ax3.get_ylim())
	ax3_f.set_yticklabels('')

	ax4.set_xlabel('Frequency (Hz)')
	ax4.set_xlim(0,f_max)
	ax4.set_yticks([])
	ax4_f.set_ylim(0,1.2)
	if applied_filter == 0:
		ax4_f.set_yticks([])

	ax4.set_ylim(pmin,pmax)
	ax4.set_title(f_title)

	ax5.set_xlim(0,xmax)
	ax5.set_ylim(ax3.get_ylim())
	ax5.set_xlabel(ax1.get_ylabel())
	ax5.set_title(xi_title)
	ax5.set_yticklabels('')
	
	ax5_f.set_ylim(0,f_est.max()*1.2)
	ax5_f.set_ylabel('Frequency (Hz)',color='green')

	return fig