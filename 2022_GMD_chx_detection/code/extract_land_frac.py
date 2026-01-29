import os, ngl, subprocess as sp, numpy as np, xarray as xr, glob
#-------------------------------------------------------------------------------
# case = 'E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00'
case = 'E3SM.RGMA.ne120pg2_r05_oECv3.FC5AV1C-H01A.00'

# data_path =os.getenv('HOME')+f'/E3SM/scratch/{case}/run/*eam.h0*'
data_path =os.getenv('HOME')+f'/E3SM/scratch/{case}/run/*cam.h0*'

first_file_path = sorted(glob.glob(data_path))[0]

ds = xr.open_dataset( first_file_path )

# out_file = 'data/E3SM_landfrac_ne30pg2.nc'
out_file = 'data/E3SM_landfrac_ne120pg2.nc'

ds['LANDFRAC'].to_netcdf(out_file)

print(f'\n{out_file}\n')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
