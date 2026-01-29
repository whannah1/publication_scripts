import os, ngl, sys, numba, copy, re, string, subprocess as sp
import xarray as xr, numpy as np, numba, itertools
import pandas as pd, datetime, cftime
import hapy_common as hc
# import pg_checkerboard_utilities as pg

scratch_path = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'

case_in = 'E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00'
data_to_mimic = f'{scratch_path}/{case_in}/run/{case_in}.eam.h1.0001-01-01-00000.nc'


case_out = 'E3SM.CHX.RAND'
output_path = f'{scratch_path}/{case_out}/run'

num_files = 365*5

var = 'TGCLDLWP'

# pd.date_range(end = datetime.today(), periods = 100).to_pydatetime().tolist()
# pd.date_range(start="2018-09-09",end="2020-02-02").to_pydatetime().tolist()
# date_list = pd.date_range(start='2001-01-01',end='2002-01-01')

#-------------------------------------------------------------------------------
# Generate list of time stamps (be sure to skip leap days)
#-------------------------------------------------------------------------------
start_date, date_list = datetime.date(1, 1, 1), []
for day in range(num_files):
  a_date = ( start_date + datetime.timedelta(days=day) ).isoformat()
  if '-02-29' not in a_date: date_list.append(a_date) # remove leap days

#-------------------------------------------------------------------------------
# Create dataset to hold random data
#-------------------------------------------------------------------------------
# scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
# scrip_ds = xr.open_dataset(scripfile_path)
# ncol = len(scrip_ds['grid_size'])

ds_in = xr.open_dataset(data_to_mimic)
time = ds_in['time'][0]
ncol = len(ds_in['ncol'].values)

first_time = cftime.datetime(1,1,1,calendar='noleap')

time.values = first_time
time['time'] = time.values

ds = xr.Dataset()
ds[var] = xr.zeros_like(ds_in[var].isel(time=[0,]))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# @numba.njit
def populate_random_data(ncol):
  data_out = np.zeros(ncol)
  for n in range(ncol):
    data_out[n] = np.random.uniform(low=0.0, high=1.0, size=1)
    # data_out[n] = np.random.uniform(0.0,1.0,size=1)
    # data_out[n] = numba.cuda.random.xoroshiro128p_uniform_float64()
  return data_out
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# rs = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(123456789)))
np.random.seed(1357924680)

for f in range(num_files):

  time_delta = datetime.timedelta(days=f)

  time_out = time + time_delta
  time_out['time'] = time_out.values

  time_stamp = (first_time + time_delta ).strftime('%Y-%m-%d')

  output_file = f'{output_path}/{case_out}.eam.h1.{time_stamp}-00000.nc'
  
  ds[var].values[0,:] = populate_random_data(ncol)

  # ds['time'].values = np.atleast_1d(time_out.values)
  ds['time'] = np.atleast_1d(time_out.values)
  # ds.assign(time=ds['time'])
  
  # print(); print(ds)
  # exit()

  ds.to_netcdf(output_file)

  print(f'  {output_file}')


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------