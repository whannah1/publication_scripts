import os, subprocess as sp, datetime, glob
#-------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#-------------------------------------------------------------------------------
def run_cmd(cmd):
  print(f'\n{clr.GREEN}{cmd}{clr.END}')
  os.system(cmd)
  return
#-------------------------------------------------------------------------------
# see https://gpm.nasa.gov/data/directory
#-------------------------------------------------------------------------------

path_top = 'https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3'
product = 'GPM_3IMERGM.07'

data_root_hdf = '/pscratch/sd/w/whannah/Obs/IMERG/monthly_hdf'
data_root_cdf = '/pscratch/sd/w/whannah/Obs/IMERG/monthly'
data_root_rmp = '/pscratch/sd/w/whannah/Obs/IMERG/monthly_remap_ne30pg2'

'''
ncremap -G ttl='Equi-Angular grid 0.1x0.1 deg'#latlon=1800,3600#lat_typ=uni#lon_typ=180_wst -g /pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc
SRC_GRID_FILE=/pscratch/sd/w/whannah/e3sm_scratch/files_grid/IMERG_1800x3600_scrip.nc
DST_GRID_FILE=grid_files/ne30pg2_scrip.nc
MAP_FILE=map_files/map_1800x3600_IMERG_to_ne30pg2_traave_20241001.nc
ncremap -a traave --src_grd=${SRC_GRID_FILE} --dst_grd=${DST_GRID_FILE} --map_file=${MAP_FILE}
'''
map_file = 'map_files/map_1800x3600_IMERG_to_ne30pg2_traave_20241001.nc'
#-------------------------------------------------------------------------------
# get_data     = True
# convert_data = True
# permute_dims = True
remap_data   = True
#-------------------------------------------------------------------------------
if 'get_data' not in locals(): get_data = False
if get_data:
  # for yr in range(2000,2014+1):
  for yr in range(1999,1999+1):

    file_path = f'{path_top}/{product}/{yr}/'
    #-----------------------------------------------------------------------------
    # download files
    cmd  = f'wget'
    cmd += f' --load-cookies ~/.urs_cookies'
    cmd += f' --save-cookies ~/.urs_cookies'
    cmd += f' --keep-session-cookies'
    # cmd += f' -O' # overwrite existing files
    cmd += f' -r -c -nH -nd -np -A HDF5 --content-disposition'
    cmd += f' "{file_path}"'
    cmd += f' --directory-prefix={data_root_hdf}'
    
    run_cmd(cmd)

    # exit()
#-------------------------------------------------------------------------------
if 'convert_data' not in locals(): convert_data = False
if convert_data:
  file_list = sorted(glob.glob(f'{data_root_hdf}/3B-MO.MS.MRG.3IMERG.*.HDF5'))
  for src_file_name in file_list : 
    dst_file_name = src_file_name
    dst_file_name = dst_file_name.replace('.HDF5',f'.nc')
    dst_file_name = dst_file_name.replace(data_root_hdf,data_root_cdf)
    var_opt = '-v time,time_bnds,lat,lat_bnds,lon,lon_bnds,precipitation'
    run_cmd(f'ncks -5 -O {src_file_name} {dst_file_name} {var_opt}')

#-------------------------------------------------------------------------------
if 'permute_dims' not in locals(): permute_dims = False
if permute_dims:
  file_list = sorted(glob.glob(f'{data_root_cdf}/3B-MO.MS.MRG.3IMERG.*.nc'))
  for src_file_name in file_list : 
    dst_file_name = src_file_name
    run_cmd(f'ncpdq -O -a time,lat,lon {dst_file_name} {dst_file_name}')
    # exit()
#-------------------------------------------------------------------------------
if 'remap_data' not in locals(): remap_data = False
if remap_data:
  file_list = sorted(glob.glob(f'{data_root_cdf}/3B-MO.MS.MRG.3IMERG.*.nc'))
  for src_file_name in file_list : 
    dst_file_name = src_file_name
    dst_file_name = dst_file_name.replace('.nc',f'.remap_ne30pg2.nc')
    dst_file_name = dst_file_name.replace(data_root_cdf,data_root_rmp)
    run_cmd(f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}  ')
    # exit()
#-------------------------------------------------------------------------------

print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
