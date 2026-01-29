#!/usr/bin/env python3
# This script is for converting IMERG HDF5 into netcdf4
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
import sys,os, subprocess as sp
import xarray as xr, numpy as np

src_dir = '/pscratch/sd/w/whannah/Obs/IMERG/hourly'
dst_dir = '/pscratch/sd/w/whannah/Obs/IMERG/hourly_nc'

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
   print(f'\n  {msg}')
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

if len(files)==0: print('\nNo files found...?\n')

if not os.path.exists(dst_dir): raise OSError(f'Destination directory is missing: {dst_dir}')

for f_in in files : 

   src_file_name = src_dir+'/'+f_in
   dst_file_name = dst_dir+'/'+f_in.replace('.HDF5',f'.nc')
   
   if os.path.isfile(dst_file_name) :
     if overwrite : os.remove(dst_file_name)
     else : continue

   run_cmd(f'ncks {src_file_name} {dst_file_name}')

   # exit()
   
print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
