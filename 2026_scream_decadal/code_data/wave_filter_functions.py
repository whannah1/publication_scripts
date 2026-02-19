import xarray as xr
import numpy as np
from scipy.signal import convolve2d, detrend
# import logging
# logging.basicConfig(level=logging.INFO)

#---------------------------------------------------------------------------------------------------

def split_hann_taper(series_length, fraction):
    """Implements `split cosine bell` taper of length series_length where only fraction of points are tapered (combined on both ends).
    
    This returns a function that tapers to zero on the ends. To taper to the mean of a series X:
    XTAPER = (X - X.mean())*series_taper + X.mean()
    """
    npts = int(np.rint(fraction * series_length))  # total size of taper
    taper = np.hanning(npts)
    series_taper = np.ones(series_length)
    series_taper[0:npts//2+1] = taper[0:npts//2+1]
    series_taper[-npts//2+1:] = taper[npts//2+1:]
    return series_taper

#---------------------------------------------------------------------------------------------------

def wave_filter(inData,waveName):
    '''
    input: 
    inData: xarray, format (optional: plev, lat), (default: time, lon)
    waveName:
        kelvin: symmetric spectrum, wave number 1<=k<=20, frequency ω< 0.4 cycle/day
        mrg (mixed rossby gravity wave): wave number -20<=k<=20, frequency 0.1< ω< 0.5 cycle/day
        rossby: wave number -20<=k<=20, frequency ω< 0.4 cycle/day, but not kelvin not mrg, here we use -20<=k<=0, frequency ω< 0.1
        ig (resolved inertial gravity wave): wave number -20<=k<=20, frequency ω >= 0.4 cycle/day
        ssg (resolved small scale gravity wave): wave number k>20
    output:
    retVal: the filtered variable in the same format with inData 
    '''
    # set the wave parameters
    tMin = 2.5
    tMax = 30
    hMin = 8
    hMax = 90

    lat_dim = inData.dims.index('lat')
    if waveName in ['Kelvin','kelvin','KELVIN']:
        arr_sym = 0.5*(inData.values + np.flip(inData.values, axis=lat_dim)) # sym only
        inData = xr.DataArray(arr_sym, dims=inData.dims, coords=inData.coords)

    # elif waveName in ['MRG','IG0','mrg','ig0']: 
    #     arr_asy = 0.5*(inData.values - np.flip(inData.values, axis=lat_dim)) # asy only
    #     inData = xr.DataArray(arr_asy, dims=inData.dims, coords=inData.coords)


    # number of timesteps in data for each day
    dtime = ( inData['time'][1]-inData['time'][0] ).values / np.timedelta64(1, 's')
    obsPerDay = int( 86400./dtime )

    lonDim = len(inData['lon'])
    timeDim  = len(inData['time'])

    # Find constants
    PI = np.arccos( -1 )
    beta = 2.28e-11
    #  cMin = ( 9.8 * hMin )^0.5
    #  cMax = ( 9.8 * hMax )^0.5
    c = ( 9.8 * np.array([ hMin, hMax ])) **0.5
    spc = 24 * 3600. / ( 2 * PI * obsPerDay ) # seconds per cycle

    tempData = inData.transpose(...,"lon", "time")


    # detrend and taper the data

    if  np.logical_not(np.any(np.isnan(tempData))):
        tempDataDetr = detrend(tempData.values, axis=-1, type='linear') #<-- missing data makes this not work
        tempData_detr = xr.DataArray(tempDataDetr, dims=tempData.dims, coords=tempData.coords)
    else:
        tempData_cp = tempData.values.copy()
        tempData_cp[np.logical_not(np.isnan(tempData_cp))] = detrend(tempData_cp[np.logical_not(np.isnan(tempData_cp))])
        tempData_detr = xr.DataArray(tempData_cp, dims=tempData.dims, coords=tempData.coords)

    taper = split_hann_taper(timeDim, 0.1)
    taper = xr.DataArray(taper, dims=('time'), coords={'time':tempData.coords['time']})
    tempData_tap = tempData_detr*taper  

    # perform the 2-d fourier transform
    z = np.fft.rfft2(tempData_tap)
    # Fourier -> e^i(kx+wt) != Wave -> e^i(kx-wt)
    z = np.append(z[...,:1,:], z[...,1:,:][...,::-1,:], axis=-2) 
    coords = {dim:inData.coords[dim] for dim in tempData.dims[:-2]}
    coords.update({"wavenumber":np.fft.fftfreq(lonDim, 1/lonDim),
                            "frequency":np.fft.rfftfreq(timeDim, 1/obsPerDay),
                                 })
    fftData = xr.DataArray(z, dims=(*tempData.dims[:-2],"wavenumber","frequency",), 
                     coords=coords)

    # Find do the cut-offs

    if waveName in ['Kelvin','kelvin','KELVIN']:
        freqInd = (fftData['frequency']<=0.4) \
                 &(fftData['frequency']>=0)
        kInd    = (fftData['wavenumber']<=20) \
                 &(fftData['wavenumber']>=1)
    elif waveName in ['MJO']:
        freqInd = (fftData['frequency']<=1./ 20.) \
                 &(fftData['frequency']>=1./100.)
        kInd    = (fftData['wavenumber']<=10) \
                 &(fftData['wavenumber']>=0)
    elif waveName in ['MRG','IG0','mrg','ig0']:
        freqInd = (fftData['frequency']<=0.5) \
                 &(fftData['frequency']>=0.1)
        kInd    = (fftData['wavenumber']<=0) \
                 &(fftData['wavenumber']>=-20)
    elif waveName in ['Rossby','rossby','ROSSBY']:
        freqInd = (fftData['frequency']<=0.1) \
                 &(fftData['frequency']>=0)
        kInd    = (fftData['wavenumber']<=0) \
                 &(fftData['wavenumber']>=-20)
    elif waveName in ['Ig','IG','ig']:
        freqInd = (fftData['frequency']>=0.4) 
        kInd    = (abs(fftData['wavenumber']<=20)) \
                 &(abs(fftData['wavenumber']>=-20))
    elif waveName in ['SSG','ssg']:
        freqInd = (np.isfinite(fftData['frequency']))
        kInd    = (abs(fftData['wavenumber']>=20)) \
                 |(abs(fftData['wavenumber']<=-20))

    fftData.loc[{'frequency' :fftData['frequency' ][np.logical_not(freqInd)]}] = 0
    fftData.loc[{'wavenumber':fftData['wavenumber'][np.logical_not(kInd)]}] = 0

    # perform the inverse transform to reconstruct the data
    retVal = inData.copy()
    fftData = np.append(fftData[...,:1,:], fftData[...,1:,:][...,::-1,:], axis=-2) # Fourier -> e^i(kx+wt) != Wave -> e^i(kx-wt) # wave back to fourier 
    retVal = np.fft.irfft2(fftData).real
    retVal = xr.DataArray(retVal, dims=tempData.dims, 
                coords=tempData.coords).transpose(*inData.dims)
    return retVal

#---------------------------------------------------------------------------------------------------
