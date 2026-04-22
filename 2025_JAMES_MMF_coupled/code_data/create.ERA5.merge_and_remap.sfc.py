#!/usr/bin/env python3
import os, glob, subprocess as sp, xarray as xr
#---------------------------------------------------------------------------------------------------
class tcolor: END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
def run_cmd(cmd):
  msg = tcolor.GREEN + cmd + tcolor.END ; print(f'\n{msg}')
  os.system(cmd); return
#---------------------------------------------------------------------------------------------------
'''
SRC_GRID_FILE=/global/cfs/cdirs/m3312/whannah/2023-CPL/scrip_ERA5_721x1440_n2s.nc
DST_GRID_FILE=/global/cfs/cdirs/m3312/whannah/2023-CPL/scrip_ne30pg2.nc
MAP_FILE=/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_n2s_to_ne30pg2.nc
# Generate ERA5 grid file
ncremap -G ttl='Equi-Angular grid 721x1440'#latlon=721,1440#lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr  -g ${SRC_GRID_FILE}
# Generate map file for ERA5 => ne30pg2
ncremap --alg_typ=tempest -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

'''
#---------------------------------------------------------------------------------------------------
# map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_s2n_to_ne30pg2.nc'
# src_file = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014/ERA5/ERA5_ANN_198501_201412_climo.nc'
# dst_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/ERA5_ANN_198501_201412_climo.ne30pg2.nc'

# cmd = f'ncremap -m {map_file} -i {src_file} -o {dst_file}  '
# print(cmd)
# sp.check_output(cmd,shell=True,universal_newlines=True)
#---------------------------------------------------------------------------------------------------
# var = 'skt'
# var = 'sp'
var = 'msl'

# yr1,yr2 = 1970,2014
yr1,yr2 = 2000,2000

if var=='skt': prefix = 'e5.oper.an.sfc.128_235_skt.ll025sc'
if var=='sp' : prefix = 'e5.oper.an.sfc.128_134_sp.ll025sc'
if var=='msl': prefix = 'e5.oper.an.sfc.128_151_msl.ll025sc'

# remap_monthly = True
create_annual = True

#---------------------------------------------------------------------------------------------------
'''
nohup time python -u code_data/create.ERA5_merged_and_remapped_data.py  > create.ERA5.remap.out &
nohup time python -u code_data/create.ERA5_merged_and_remapped_data.py  > create.ERA5.annual.out &
'''
#---------------------------------------------------------------------------------------------------
if 'remap_monthly' not in locals(): remap_monthly = False
if remap_monthly:

  src_data_root = '/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.sfc'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/monthly_ts'

  # for yr in range(yr1,yr2+1):
  #   file_list = sorted(glob.glob(f'{src_data_root}/{yr}*/{prefix}.*'))
  #   for src_file in file_list:
  #     dst_file = src_file
  #     dst_file = dst_file.replace(src_data_root,dst_data_root)
  #     dst_file = dst_file.replace('.nc','.remap_ne30pg2.nc')
  #     map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_n2s_to_ne30pg2.nc'
  #     run_cmd(f'ncremap -m {map_file} -i {src_file} -o {dst_file}')

  for yr in range(yr1,yr2+1):
    for mn in range(1,12+1):
      file_path = f'{src_data_root}/{yr}{mn:02d}/{prefix}*'
      src_file = glob.glob(file_path)[0]
      dst_file = src_file
      dst_file = dst_file.replace(f'{src_data_root}/{yr}{mn:02d}',dst_data_root)
      dst_file = dst_file.replace('.nc','.remap_ne30pg2.nc')
      map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_n2s_to_ne30pg2.nc'
      run_cmd(f'ncremap -m {map_file} -i {src_file} -o {dst_file}')

#---------------------------------------------------------------------------------------------------
if 'create_annual' not in locals(): create_annual = False
if create_annual: 
  src_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/monthly_ts'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/annual'

  print()
  for yr in range(yr1,yr2+1):
    dst_file = f'{dst_data_root}/{prefix}.{yr}.remap_ne30pg2.nc'

    file_path = f'{src_data_root}/{prefix}.{yr}*'
    file_list = sorted(glob.glob(file_path))

    print(f'  {yr}  {dst_file}')

    if file_list==[]:
      print(); print('File list is empty!')
      print(); print(file_path)
      print(); print(file_list); print()
      exit()

    ds = xr.open_mfdataset( file_list )

    # print(ds)

    # time_bnds = ds['time_bnds']
    # time = time_bnds.isel(nbnd=0)

    time = ds['time']

    month_length = time.dt.days_in_month
    month_length['time'] = time
    ds['time'] = time
    mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

    ds_out = xr.Dataset()
    # ds_out['time_bnds'] = ds['time_bnds'][0][:]
    ds_out['time'] = ds['time'][0]

    ds_out = xr.Dataset()
    # ds_out['time_bnds'] = ds['time_bnds'][0][:]
    # ds_out['time'] = ds['time'][0]

    for var in ds.variables:
      # print(); print(f'var: {var}')
      # print(); print(ds[var])
      if var in ['time','time_bnds','date_written','time_written']: 
        continue
      if var in ['lat_vertices','lon_vertices','area']: 
        ds_out[var] = ds[var].isel(time=0,missing_dims='ignore',drop=True)
      else:
        if 'time' in ds[var].dims:
          numerator   = (ds[var]*mn_wgts).resample(time='YE').sum('time')
          denominator = (mn_wgts).resample(time='YE').sum(dim='time')
          ds_out[var] = numerator / denominator
        else:
          ds_out[var] = ds[var]

    ds_out.to_netcdf(path=dst_file,mode='w')
    
  print()
#---------------------------------------------------------------------------------------------------
