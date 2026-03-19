#!/usr/bin/env python3
# This script is for converting IMERG HDF5 into netcdf4
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
import sys,os, subprocess as sp, glob, datetime
import xarray as xr, numpy as np

src_dir = '/global/cfs/cdirs/m4842/whannah/IMERG/hourly_hdf'
dst_dir = '/global/cfs/cdirs/m4842/whannah/IMERG/hourly_netcdf'

overwrite = True

#-------------------------------------------------------------------------------
def add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2):
   done = False
   cnt = 0
   # yrmndy_beg = int(yr1*1e4+mn1*1e2+dy1)
   # yrmndy_end = int(yr2*1e4+mn2*1e2+dy2)
   beg_date = datetime.date(yr1, mn1, dy1)
   end_date = datetime.date(yr2, mn2, dy2)
   while not done:
      cur_date = beg_date + datetime.timedelta(days=cnt)
      # day_of_year = cur_date.timetuple().tm_yday
      yrmndy = int( cur_date.year*1e4 + cur_date.month*1e2 + cur_date.day )
      yrmndy_list.append(yrmndy)
      if cur_date==end_date: done = True
      cnt+=1
#-------------------------------------------------------------------------------
yrmndy_list = []

yr1,mn1,dy1=2023, 4,26; yr2,mn2,dy2=2023, 7, 8 ; add_yrmndy(yrmndy_list,yr1,mn1,dy1,yr2,mn2,dy2)

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

if not os.path.exists(dst_dir): raise OSError(f'Destination directory is missing: {dst_dir}')

# files = sorted(os.listdir(src_dir))
# if len(files)==0: print('\nNo files found...?\n')

# for f_in in files : 
#    src_file_name = src_dir+'/'+f_in
#    dst_file_name = dst_dir+'/'+f_in.replace('.HDF5',f'.nc')
#    if os.path.isfile(dst_file_name) :
#    if overwrite : os.remove(dst_file_name)
#       else : continue
#    run_cmd(f'ncks {src_file_name} {dst_file_name}')

for yrmndy in yrmndy_list:
   yr = int(yrmndy/1e4)
   mn = int((yrmndy-yr*1e4)/1e2)
   dy = int(yrmndy-yr*1e4-mn*1e2)
   cur_date = datetime.date(yr, mn, dy)
   day_of_year = cur_date.timetuple().tm_yday

   files = sorted(glob.glob(f'{src_dir}/3B-HHR.MS.MRG.3IMERG.{yrmndy}-*.HDF5'))
   if len(files)==0: 
      print(f'\nNo files found for yrmndy: {yrmndy}\n')
      exit()

   for f_in in files : 
      src_file_name = f_in
      dst_file_name = f_in.replace('.HDF5',f'.nc').replace(src_dir,dst_dir)

      if os.path.isfile(dst_file_name) :
         if overwrite : os.remove(dst_file_name)
         else : continue
      
      run_cmd(f'ncks {src_file_name} {dst_file_name}')

   exit()

print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
