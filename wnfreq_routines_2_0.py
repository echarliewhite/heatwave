"""
This module was authored by Oliver Watt-Meyer for use in spectral analysis of 
standing and travelling waves as described in "Decomposition of Atmospheric 
Disturbances into Standing and Traveling Components, with Application to 
Northern Hemisphere Planetary Waves and Stratosphere–Troposphere Coupling" 
by Watt-Meyer and Kushner (2014).
"""
import numpy as np
from matplotlib import pyplot as plt

def calc_wnfreq_spectrum(data,wn_max,smooth_gaussian=False,smooth_amount=0.0):
    """Compute total, standing and travelling spectra of input data.

    Compute the spectral analysis of Watt-Meyer and Kushner (2014).

    Parameters
    ----------
    data : ndarray, where first dimension is time, and last dimension is longitude
	Input data for spectral analysis.
    wn_max : int
        Largest wavenumber to keep in output.
    smooth_gaussian : boolean, optional
        Whether to smooth Fourier coefficients with Gaussian filter prior to
	decomposing into standing and travelling parts. Default is False.
    smooth_amount : float, optional
        Width of Gaussian smoothing to use (in units of the frequency index).
	Default is 0.0.

    Returns
    -------
    fft2_trans : complex ndarray, same shape as data except final axis has length wn_max
        2D Fourier coefficients (over time and longitude) of input data. Note
	the coefficients are stored such that the zero-frequency component is at
	the centre of the frequency axis, which has the same length as the input
	time axis. The wavenumbers go from 1 to wn_max.
    fft2_standing : as above
        The standing component of the Fourier coefficients. Stored as above.
    fft2_travelling : as above
        The travelling component of the Fourier coefficients. Stored as above.
    """
    data_shape=np.shape(data)

    # calculate 2D Fourier transform on time and longitude axes
    fft2_data=np.fft.fft2(data,axes=(0,len(data_shape)-1))

    # rearrange data so freq=0 point is at centre of time dimension, and only keep wn_max wavenumbers
    fft2_trans=transform_wnfreq_spectrum(fft2_data,wn_max)

    # if the data had an even number of timesteps, remove the last Fourier coefficient, which corresponds
    # to the Nyquist frequency and doesn't have a corresponding opposite-sign frequency.
    if data.shape[0]%2==0:
	fft2_nyquist=fft2_trans[data_shape[0]-1:data_shape[0],...]
	fft2_trans=fft2_trans[:-1,...]

    # compute absolute value and phases of fourier coefficeints
    fft2_trans_abs=abs(fft2_trans)
    fft2_trans_phases=np.angle(fft2_trans)

    # do optional smoothing before decomposing into standing and travelling parts
    if smooth_gaussian:
        fft2_trans_abs=np.sqrt(smooth_wnfreq_spectrum_gaussian(fft2_trans_abs**2,smooth_amount)) # smoothing amplitude squared
    	fft2_trans=fft2_trans_abs*(np.cos(fft2_trans_phases)+1j*np.sin(fft2_trans_phases)) # reconstruct total spectrum with smoothed amplitude

    fft2_standing_abs=np.minimum(fft2_trans_abs,fft2_trans_abs[::-1,...]) # minimum of eastward and westward parts of spectrum
    fft2_travelling_abs=fft2_trans_abs-fft2_standing_abs # remainder (i.e. travelling part) computed based on amplitude

    # reconstruct standing and travelling parts of spectrum with original phase information
    fft2_standing=fft2_standing_abs*(np.cos(fft2_trans_phases)+1j*np.sin(fft2_trans_phases))
    fft2_travelling=fft2_travelling_abs*(np.cos(fft2_trans_phases)+1j*np.sin(fft2_trans_phases))

    # if data had an even number of timesteps, append back on the Nyquist frequency component
    if data.shape[0]%2==0:
	fft2_trans=np.append(fft2_trans,fft2_nyquist,axis=0)
	fft2_standing=np.append(fft2_standing,np.zeros(np.shape(fft2_nyquist),dtype='complex'),axis=0)
	fft2_travelling=np.append(fft2_travelling,fft2_nyquist,axis=0)

    return (fft2_trans,fft2_standing,fft2_travelling)


def invert_wnfreq_spectrum(fft_coeffs_in,wn_min,wn_max,N,tol=1e-6):
    """Compute real space inverse of fft_coeffs_in, summed over wavenumber wn_min to wn_max.

    Invert fft_coeffs_in to real space and time. The wavenumbers are restricted
    to be between wn_min and wn_max, inclusive. The data is inverted onto a grid
    with N points in longitude. Note that only unsmoothed Fourier coefficients
    should be used when inverting.

    Parameters
    ----------
    fft_coeffs_in : complex ndarray of Fourier coefficients
	Input Fourier coefficients assumed to be stored in same way as output
	of calc_wnfreq_spectra
    wn_min : int, Minimum wavenumber for which to compute inverse.
    wn_max : int,  Maximum wavenumer for which to compute inverse.
    N : int
	Number of points of longitude to invert data onto

    Returns
    -------
    data_out : ndarray, same shape as input except final axis has length N
	Real space output data.
    """

    fft_coeffs_shape=np.shape(fft_coeffs_in)

    # set all coefficients for wavenumbers outside wn_min to wn_max to zero
    fft_coeffs_in[...,:wn_min-1]=0.0
    fft_coeffs_in[...,wn_max:]=0.0

    # put coefficients back in "standard" fft order, fill with zeros in order
    # to have N wavenumbers
    fft_coeffs_in=untransform_wnfreq_spectrum(fft_coeffs_in,N)

    # compute inverse fourier transform
    data_out=np.fft.ifft2(fft_coeffs_in,axes=(0,len(fft_coeffs_shape)-1))

    # check that inverted data is real to within some tolerance
    assert abs(np.imag(data_out)).max()<tol

    # return real part of data
    return np.real(data_out)


def transform_wnfreq_spectrum(fft_coeffs_in,wn_max):
    """Rearrange Fourier coefficients and limit to wn_max wavenumbers.

    Shift Fourier coefficients such that the frequency axis has the
    zero-frequency component at the center, and goes from positive (i.e.
    eastward) frequencies to negative (i.e. westward). Also limit the
    wavenumbers of wn_max, and store them starting from 1 to wn_max.

    Parameters
    ----------
    fft_coeffs_in : complex ndarray, where first dimension is frequency, and last dimension is wavenumber
	Input Fourier coefficients ordered following standard output of fft2.
    wn_max : int
	Largest wavenumber to keep in output.

    Returns
    -------
    fft2_coeffs_out : complex ndarray, same dimensions as input except final dimensional has length wn_max
	Rearranged coefficients, limited to wn_max wavenumbers.

    """
    fft_coeffs_out=fft_coeffs_in[...,1:wn_max+1] # keep only wn_max wavenumbers (but not including wave-0, i.e. zonal mean)
    fft_coeffs_out=np.fft.fftshift(fft_coeffs_out,axes=0) # shift frequency axis
    fft_coeffs_out=fft_coeffs_out[::-1,...] # invert order of frequency components (so westward first, then eastward)

    return fft_coeffs_out

def untransform_wnfreq_spectrum(fft_coeffs_in,N):
    """Rearrange Fourier coefficients and pad to N wavenumbers.

    Performs the opposite operations of transform_wnfreq_spectrum to make data
    with limited number of wavenumbers and frequency centered on 0 into the
    appropriate form for ifft2.

    Parameters
    ----------
    fft_coeffs_in : complex ndarray, where first dimension is frequency, and last dimension is wavenumber
	Input Fourier coefficients ordered following output of
	transform_wnfreq_spectrum.
    N : int
	Number of points of longitude for which to calculate realspace data.

    Returns
    -------
    fft2_coeffs_out : complex ndarray, same dimensions as input except final dimension has length N
	Rearranged coefficients, padded to N wavenumbers.

    """
    fft_coeffs_shape=np.shape(fft_coeffs_in)
    num_dim=len(fft_coeffs_shape)
    wn=fft_coeffs_shape[-1] # number of wavenumbers in input data

    # create output coefficient array with N wavenumbers
    fft_coeffs_out_shape=np.array(fft_coeffs_shape)
    fft_coeffs_out_shape[num_dim-1]=N
    fft_coeffs_out=np.zeros(fft_coeffs_out_shape,dtype='complex')

    fft_coeffs_in=fft_coeffs_in[::-1,...] # reverse order of frequencies

    fft_coeffs_out[...,1:wn+1]=fft_coeffs_in
    fft_coeffs_out[...,-1:-wn-1:-1]=np.conj(fft_coeffs_in[::-1,...])

    fft_coeffs_out=np.fft.ifftshift(fft_coeffs_out,axes=0) # shift frequencies

    return fft_coeffs_out


def smooth_wnfreq_spectrum_gaussian(fft_coeffs_in,smooth_amount):
    """Smooth Fourier coefficients over frequency.

    Smooth the fourier coefficients over frequency with a Gaussian filter of
    width smooth_amount.

    Parameters
    ----------
    fft_coeffs_in : complex ndarray, first dimension must be frequency
	Input Fourier coefficients.
    smooth_amount : float
	Amount to smooth the Fourier coefficients, in units of the frequency index.

    Returns
    -------
    fft_coeffs_out : complex ndarray, same shape as input
	Smoothed Fourier coefficients.

    """
    fft_coeffs_shape=np.shape(fft_coeffs_in)
    num_dim=len(fft_coeffs_shape)
    T=fft_coeffs_shape[0]

    fft_coeffs_shape2=np.copy(fft_coeffs_shape)
    fft_coeffs_shape2[0]=1

    smoother_shape=[T]+[1]*(num_dim-1)

    fft_coeffs_out=np.zeros(fft_coeffs_shape)

    for t in range(T):
	smoother=np.reshape(np.exp(-((np.arange(T)-t)/float(smooth_amount))**2),smoother_shape)
	smoother/=sum(smoother)

	fft_coeffs_out[t,...]=np.sum(fft_coeffs_in*np.tile(smoother,fft_coeffs_shape2),axis=0)

    return fft_coeffs_out


def plot_wnfreq_spectrum_lineplots(fft_coeffs,fig_num,freq_cutoff,scale_factor,my_figsize=(8,6),title='',plot_xlim=2.0,my_ylabel='Wavenumber',my_linestyle='-k',my_linewidth=1.0,second_plot=False,leg_label='',fig_label=''):
    """Plot wavenumber-frequency spectrum as series of lineplots

    Plots one line of power versus frequency for each wavenumber in the
    inputted coefficients.

    Parameters
    ----------
    fft_coeffs : real 2D array, with dimensions (freq x wavenumber)
	Input coefficients for plotting. Note is assumed the input coefficients
	are real, i.e. the power spectrum has already been computed.
    fig_num : int
        The figure number in which to plot the data.
    freq_cutoff : int
        Integer representing how many frequency indices away from the zero-
	frequency to start plotting the spectrum.
    scale_factor : float
        The value by which all input data scaled for plotting.


    """
    T=np.size(fft_coeffs,axis=0)
    freq=np.fft.fftshift(np.fft.fftfreq(T))

    # if there is an even number of frequencies, remove the Nyquist frequency so that frequencies are symmetric about zero.
    if T%2==0:
	T-=1
	freq=freq[1:]
	fft_coeffs=fft_coeffs[:-1,...]

    plot_wavenumbers=np.size(fft_coeffs,axis=1)

    fig=plt.figure(fig_num,figsize=my_figsize)

    if not(second_plot) and title!='':
        plt.suptitle(title+'. (scale='+str(round(scale_factor,1))+')')

    ax1 = fig.add_subplot(111)

    for wn in range(1,plot_wavenumbers+1):
        if wn==plot_wavenumbers:
	    ax1.plot(freq[:T/2+1-freq_cutoff],wn+fft_coeffs[:T/2-freq_cutoff+1,wn-1]/scale_factor,my_linestyle,linewidth=my_linewidth,label=leg_label)
	    ax1.plot(freq[T/2+freq_cutoff:],wn+fft_coeffs[T/2+freq_cutoff:,wn-1]/scale_factor,my_linestyle,linewidth=my_linewidth)
	else:
	    ax1.plot(freq[:T/2+1-freq_cutoff],wn+fft_coeffs[:T/2-freq_cutoff+1,wn-1]/scale_factor,my_linestyle,linewidth=my_linewidth)
	    ax1.plot(freq[T/2+freq_cutoff:],wn+fft_coeffs[T/2+freq_cutoff:,wn-1]/scale_factor,my_linestyle,linewidth=my_linewidth)

    if not(second_plot):
	my_xticks=np.array([-1./4.,-1./5.,-1./10.0,-1./20.0,-1./50.0,0.0,1./50.0,1./20.0,1./10.0,1./5.,1./4.])
	my_xticks_freq_labels=['0.25','0.2','0.1','0.05','','0','','0.05','0.1','0.2','0.25']
	my_xticks_period_labels=[4,5,10,20,50,'',50,20,10,5,4]

	plt.xticks(my_xticks,my_xticks_freq_labels)
        ax1.set_xlim([-1./plot_xlim,1./plot_xlim])
        cur_xticks=ax1.get_xticks()
        ax1.set_yticks(range(1,plot_wavenumbers+1))
	ax1.set_ylabel(my_ylabel,fontsize=11)

	my_yticks_labels_1=[]
	my_yticks_labels_2=[]
	for wn in range(1,plot_wavenumbers+1):
	    my_yticks_labels_1.append('Wave-'+str(wn))

	plt.yticks(np.arange(1.5,plot_wavenumbers+1,1),my_yticks_labels_1,rotation='vertical',va='bottom')

        ax1.set_ylim([1,plot_wavenumbers+1])
        ax1.set_xlabel('Westward     Freq. [days$^{-1}$]     Eastward',labelpad=-.05)

        ax2=ax1.twiny()
	plt.yticks(range(1,plot_wavenumbers+1,1))

	plt.xticks(my_xticks,my_xticks_period_labels)

        ax2.set_xlim([-1./plot_xlim,1./plot_xlim])
        ax2.set_xlabel('Westward      Period [days]      Eastward')

        ax1.grid(axis='y')
        ax2.grid(axis='x')

	ax2.text(-1./plot_xlim+0.011,plot_wavenumbers+0.8,fig_label,backgroundcolor='white',va='top')
