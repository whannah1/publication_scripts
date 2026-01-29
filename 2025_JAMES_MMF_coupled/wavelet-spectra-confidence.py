import numpy as np, xarray as xr, scipy
#---------------------------------------------------------------------------------------------------
def deseason(xraw):
   months_per_year = 12
   xclim = np.zeros((months_per_year, 1))
   x_deseasoned = np.zeros(xraw.shape)
   # Iterate through all 12 months.
   for month in np.arange(months_per_year):
      xclim[month] = np.nanmean(xraw[month::months_per_year])
   num_years = int(np.floor(len(x_deseasoned) / months_per_year))
   for year in np.arange(num_years):
      year_index = year * months_per_year
      # Iterate through all months of the year
      for month in np.arange(months_per_year):
         month_index = year_index + month
         # Subtract the month's mean over num_years for this month/year
         x_deseasoned[month_index] = xraw[month_index] - xclim[month]
   return x_deseasoned
#---------------------------------------------------------------------------------------------------
def get_psd_from_wavelet(data,deg=6):
   period = np.arange(6,12*20+1)
   widths = deg / ( 2*np.pi * (1/period) )
   cwtmatr = scipy.signal.cwt( deseason(data), scipy.signal.morlet2, widths=widths, w=deg )
   psd = np.mean(np.square(np.abs(cwtmatr)),axis=1)
   return ( period, psd )
#---------------------------------------------------------------------------------------------------

file_name = 'data/ENSO-power-spectra.HadSST.TS.reg_Nino3.4.lat1_-5.lat2_5.lon1_190.lon2_240.nc'

# read the data
ds = xr.open_dataset( file_name, use_cftime=True  )
da = ds['TS']

# convert to anomalies
da = da - da.mean()

# detrend in time
da = da - xr.polyval(da['time'], da.polyfit(dim='time', deg=1).polyfit_coefficients)

( period, psd ) = get_psd_from_wavelet( da.values )

#---------------------------------------------------------------------------------------------------
# estimate degrees of freedom (dof) - Auchère et al. (2016) Eq 13 - DOI 10.3847/0004-637X/825/2/110
# for p in period:
#    dof = 2 * np.sqrt( 1 + (len(da)*1)/(p/0.70) )
#    print(dof)
# exit()
#---------------------------------------------------------------------------------------------------
# calculate theretical red noise spectrum - Torrence & Compo (1998) Eq 16
alpha = 0.70; alpha_sq = np.square(alpha)
P = ( 1 - alpha_sq ) / ( 1 + alpha_sq - 2*alpha*np.cos(2*np.pi*(1/period)) )
#---------------------------------------------------------------------------------------------------
# calculate the 95% confidence level - Torrence & Compo (1998) Eq 18
# 5.99 = chi-squared value for 2 dof and alpha=0.05
confidence_level_95 = 0.5 * P * 5.99 * np.var( da.values ) 
#---------------------------------------------------------------------------------------------------
# plot stuff

period = period/12. # convert period from months to years for plotting

# ...?

#---------------------------------------------------------------------------------------------------
