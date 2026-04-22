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
var = 'v'

yr1,yr2 = 1980,2014

if var=='t' : prefix = 'e5.oper.an.pl.128_130_t.ll025sc'
if var=='u' : prefix = 'e5.oper.an.pl.128_131_u.ll025sc'
if var=='v' : prefix = 'e5.oper.an.pl.128_132_v.ll025sc'
if var=='q' : prefix = 'e5.oper.an.pl.128_133_q.ll025sc'
if var=='r' : prefix = 'e5.oper.an.pl.128_157_r.ll025sc'

map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_n2s_to_ne30pg2.nc'

remap_monthly = True
create_annual = True

#---------------------------------------------------------------------------------------------------
'''
nohup time python code_data/create.ERA5.merge_and_remap.atm.py > create.ERA5.u.out &
'''
#---------------------------------------------------------------------------------------------------
if 'remap_monthly' not in locals(): remap_monthly = False
if remap_monthly:

  # /global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl.monthly/e5.oper.an.pl.128_157_r.ll025sc.199604.nc
  src_data_root = '/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl.monthly'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/monthly'

  for yr in range(yr1,yr2+1):
    for mn in range(1,12+1):
      src_file = f'{src_data_root}/{prefix}.{yr}{mn:02d}.nc'
      dst_file = src_file
      dst_file = dst_file.replace(src_data_root,dst_data_root)
      dst_file = dst_file.replace('.nc','.remap_ne30pg2.nc')
      # print(); print(src_file)
      # print(); print(dst_file)
      # print()
      # exit()
      run_cmd(f'ncremap -m {map_file} -i {src_file} -o {dst_file}')

#---------------------------------------------------------------------------------------------------
if 'create_annual' not in locals(): create_annual = False
if create_annual: 

  src_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/monthly'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/annual'

  for yr in range(yr1,yr2+1):
    dst_file = f'{dst_data_root}/{prefix}.{yr}.remap_ne30pg2.nc'

    file_path = f'{src_data_root}/{prefix}.{yr}*'
    file_list = sorted(glob.glob(file_path))
    if file_list==[]:
      print(file_path)
      print(file_list)
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
          numerator   = (ds[var]*mn_wgts).resample(time='A').sum('time')
          denominator = (mn_wgts).resample(time='A').sum(dim='time')
          ds_out[var] = numerator / denominator
        else:
          ds_out[var] = ds[var]

    ds_out.to_netcdf(path=dst_file,mode='w')

    print(f'  {dst_file}')

#---------------------------------------------------------------------------------------------------
