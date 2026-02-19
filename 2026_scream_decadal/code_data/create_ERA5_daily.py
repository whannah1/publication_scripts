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

# Generate map file for ERA5 => ne30pg2
ncremap --alg_typ=tempest -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}

SRC_GRID=/global/cfs/cdirs/m3312/whannah/files_grid/scrip_ERA5_721x1440_n2s.nc
DST_GRID=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_grid/cmip6_90x180_scrip.20250624.nc
MAP_FILE=/global/cfs/cdirs/m3312/whannah/files_map/map_721x1440_n2s_to_90x180_traave.20250624.nc

# Generate ERA5 grid file
ncremap -G ttl='Equi-Angular grid 721x1440'#latlon=721,1440#lat_typ=uni#lat_drc=n2s#lon_typ=grn_ctr  -g ${SRC_GRID}

ncremap -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

'''
#---------------------------------------------------------------------------------------------------
# map_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/map_721x1440_s2n_to_ne30pg2.nc'
# src_file = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology_1985-2014/ERA5/ERA5_ANN_198501_201412_climo.nc'
# dst_file = '/global/cfs/cdirs/m3312/whannah/2023-CPL/ERA5/ERA5_ANN_198501_201412_climo.ne30pg2.nc'

# cmd = f'ncremap -m {map_file} -i {src_file} -o {dst_file}  '
# print(cmd)
# sp.check_output(cmd,shell=True,universal_newlines=True)
#---------------------------------------------------------------------------------------------------
'''
nohup time python code_data/create_ERA5_daily.py > create_ERA5_daily.out &
'''
#---------------------------------------------------------------------------------------------------
var = 'q'

if var=='t' : prefix = 'e5.oper.an.pl.128_130_t.ll025sc'
if var=='u' : prefix = 'e5.oper.an.pl.128_131_u.ll025sc'
if var=='v' : prefix = 'e5.oper.an.pl.128_132_v.ll025sc'
if var=='q' : prefix = 'e5.oper.an.pl.128_133_q.ll025sc'
if var=='r' : prefix = 'e5.oper.an.pl.128_157_r.ll025sc'

yr = 2011
mn1,mn2 = 10,12

map_file = '/global/cfs/cdirs/m3312/whannah/files_map/map_721x1440_n2s_to_90x180_traave.20250624.nc'

# remap_hourly = True
create_daily = True

#---------------------------------------------------------------------------------------------------
if 'remap_hourly' not in locals(): remap_hourly = False
if remap_hourly:

  # /global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl.monthly/e5.oper.an.pl.128_157_r.ll025sc.199604.nc
  src_data_root = '/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/ERA5/hourly'

  for mn in range(mn1,mn2+1):
    src_sub = f'{yr}{mn:02d}'
    # src_file = f'{src_data_root}/{src_sub}/{prefix}.{yr}{mn:02d}.nc'
    # src_file = f'{src_data_root}/{src_sub}/{prefix}.{yr}{mn:02d}.nc'
    src_file_path = f'{src_data_root}/{src_sub}/{prefix}.*'
    src_file_list = sorted(glob.glob(src_file_path))
    for src_file in src_file_list:
      dst_file = src_file
      dst_file = dst_file.replace(src_data_root,dst_data_root)
      dst_file = dst_file.replace('.nc','.remap_90x180.nc')
      run_cmd(f'ncremap -m {map_file} -i {src_file} -o {dst_file}')
      # exit()

#---------------------------------------------------------------------------------------------------
if 'create_daily' not in locals(): create_daily = False
if create_daily: 

  src_data_root = '/global/cfs/cdirs/m3312/whannah/ERA5/hourly'
  dst_data_root = '/global/cfs/cdirs/m3312/whannah/ERA5/daily'

  for mn in range(mn1,mn2+1):
    
    src_sub = f'{yr}{mn:02d}'

    src_file_path = f'{src_data_root}/{src_sub}/{prefix}.*'
    src_file_list = sorted(glob.glob(src_file_path))
    
    if src_file_list==[]:
      print(f'src_file_path: {src_file_path}')
      print(f'src_file_list: {src_file_list}')
      raise ValueError('no files found?') 

    for src_file in src_file_list:

      dst_file = src_file.replace( f'{src_data_root}/{src_sub}', dst_data_root )

      print()
      print(f'  src_file: {src_file}')
      print(f'  dst_file: {dst_file}')

      ds = xr.open_dataset( src_file )

      # Convert to daily mean
      ds = ds.resample(time='D').mean(dim='time')
      ds.load()

      ds.to_netcdf(path=dst_file,mode='w')
    
    print()

      

#---------------------------------------------------------------------------------------------------
