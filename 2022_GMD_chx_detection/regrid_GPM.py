#!/usr/bin/env python3

### command for aking grid file
# ncremap -G ttl='Equi-Angular grid 0.1x0.1 deg'#latlon=1800,3600#lat_typ=uni#lon_typ=180_wst -g scrip_files/GPM_1800x3600_scrip.nc

### command for making map file
# ncremap --src_grd=scrip_files/GPM_1800x3600_scrip.nc --dst_grd=scrip_files/ne30pg2_scrip.nc -m map_files/map_1800x3600_to_ne30pg2_aave.nc

### map file to get 1 degree GPM data
# ncremap --src_grd=scrip_files/GPM_1800x3600_scrip.nc --dst_grd=scrip_files/cmip6_180x360_scrip.20210722.nc -m map_files/map_1800x3600_to_180x360_aave.nc

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
import sys,os, subprocess as sp

dst_grid = 'ne30pg2'
# dst_grid = '180x360'

# src_dir = '/global/cscratch1/sd/whannah/Obs/GPM/daily_QC'
# dst_dir =f'/global/cscratch1/sd/whannah/Obs/GPM/daily_{dst_grid}'

src_dir = '/global/cscratch1/sd/whannah/Obs/IMERG/daily_QC'
dst_dir =f'/global/cscratch1/sd/whannah/Obs/IMERG/daily_{dst_grid}'

map_file =f'map_files/map_1800x3600_to_{dst_grid}_aave.nc'

overwrite = True

#-------------------------------------------------------------------------------
# Define run command
#-------------------------------------------------------------------------------
# Set up terminal colors
class tcolor:
   ENDC,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd,suppress_output=False,execute=True):
   if suppress_output : cmd = cmd + ' > /dev/null'
   msg = tcolor.GREEN + cmd + tcolor.ENDC
   print(f'\n{msg}')
   if execute: 
      # os.system(cmd)
      try:
         sp.check_output(cmd,shell=True,universal_newlines=True)
      except sp.CalledProcessError as error:
         exit(error.output)
   return
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

files = sorted(os.listdir(src_dir))

if len(files)==0: print('\nNo files found for remapping...?\n')

if not os.path.exists(dst_dir): raise OSError(f'Destination directory is missing: {dst_dir}')

for f_in in files : 

   if '.nc' not in f_in: continue

   src_file_name = src_dir+'/'+f_in
   dst_file_name = dst_dir+'/'+f_in.replace('.nc',f'.remap_{dst_grid}.nc')
   
   if os.path.isfile(dst_file_name) :
     if overwrite : os.remove(dst_file_name)
     else : continue

   # # Change the _FillValue attribute of the input file
   # run_cmd(f'ncatted -a _FillValue,precipAvg,m,f,1.0e36 -o {src_file_name} {src_file_name}')

   run_cmd(f'ncremap -m {map_file} -i {src_file_name} -o {dst_file_name}')

   # use " --permute='-a time,lat,lon' " to swap dimensions 

   # exit()


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
