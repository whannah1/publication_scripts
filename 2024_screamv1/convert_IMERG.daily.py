#!/usr/bin/env python3

# This script is for converting IMERG daily accumulated precipitation 
# into daily averaged rain rate, using the approach outlined by Greg below

# from Greg:
#  I recommend you use 'precipitationCal'.
# 'precipitationCal' is a rainfall estimate that combines both microwave estimates 
# (which are high quality) and infrared estimates (not as high quality).  
# The advantage here is sampling (b/c infrared estimates are readily available) 
# such that you will have a much more representative daily average.  

# 'HQprecipitation' (microwave; high quality), where [more sparsely] available, 
# is better if you're looking at instantaneous estimates or you want a more 
# reliable estimate of precipitation at some time on some day for some location.  
# The disadvantage is sampling - there are not enough microwave instruments 
# orbiting to get you good sampling for a given day.  So, over a large domain, 
# this will be noisier.

# When data from IR or microwave are unavailable, things are set to "missing"...
# which has an implication for daily accumulations and how to compute average rates...
# Thus, you need to use the "precipitationCal_cnt' array to get an average rainfall rate.
# IMERG gives precip estimates every 30 minutes, so over a full day, 
# 'precipitationCal_cnt' will - at most - be 48 for some box.  [unless you have 
# missing data, which I'm sure you do sometimes].
# Their daily accumulation (that you see in 'precipitationCal', in millimeters), 
# assuming 48 samples per day, is 

# 0.5 X sum{Pi_instantaneous, i = 1,48}.

# The 0.5 is in there for obvious reasons (b/c you get precip every half hour as 
# opposed to every hour, so a straight sum gives you an overestimate by 2).  
# You can read about that here (https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDF_06/summary?keywords=imerg).
# But, you want an "average rainfall rate" for the day, or an "average Pi_instantaneous" 
# per the little equation I wrote above. To get that, compute the following:

# dailyave_rr in mm/hr = precipitationCal*(2/precipitationCal_cnt)

# when precipitationCal_cnt = 48, then yes, that's equivalent to dividing by 24.  
# But, I imagine with missing data, precipitationCal_cnt is <= 48.
# So, you can't just plan to do a straight blind division by 24 everywhere.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
import sys,os, subprocess as sp
import xarray as xr, numpy as np

dst_grid = 'ne30pg2'

src_dir = '/pscratch/sd/w/whannah/Obs/IMERG/daily'
dst_dir = '/pscratch/sd/w/whannah/Obs/IMERG/daily_QC'

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
def print_stat(x,name='(no name)',unit='',fmt='f',stat='naxh',indent='',compact=False):
   """ Print min, avg, max, and std deviation of input """
   if fmt=='f' : fmt = '%.4f'
   if fmt=='e' : fmt = '%e'
   if unit!='' : unit = f'[{unit}]'
   name_len = 12 if compact else len(name)
   msg = ''
   line = f'{indent}{name:{name_len}} {unit}'
   # if not compact: print(line)
   if not compact: msg += line+'\n'
   for c in list(stat):
      if not compact: line = indent
      if c=='h' : line += '   shp: '+str(x.shape)
      if c=='a' : line += '   avg: '+fmt%x.mean()
      if c=='n' : line += '   min: '+fmt%x.min()
      if c=='x' : line += '   max: '+fmt%x.max()
      if c=='s' : line += '   std: '+fmt%x.std()
      # if not compact: print(line)
      if not compact: msg += line+'\n'
   # if compact: print(line)
   if compact: msg += line#+'\n'
   print(msg)
   return msg
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

files = sorted(os.listdir(src_dir))

if len(files)==0: print('\nNo files found...?\n')

if not os.path.exists(dst_dir): raise OSError(f'Destination directory is missing: {dst_dir}')

for f_in in files : 

   if '.nc' not in f_in: continue

   src_file_name = src_dir+'/'+f_in
   dst_file_name = dst_dir+'/'+f_in#.replace('.nc',f'.QC.nc')
   
   if os.path.isfile(dst_file_name) :
     if overwrite : os.remove(dst_file_name)
     else : continue

   print(tcolor.GREEN+f'  {src_file_name}   >   {dst_file_name}'+tcolor.ENDC)

   # load data
   ds = xr.open_dataset(src_file_name)

   # copy the datset and remove stuff we don't need
   ds_new = ds.copy(deep=True)
   drop_list,keep_list = [],['precipitationCal','precipitationCal_cnt']
   for name,da in ds_new.items():
      if name not in keep_list:
         drop_list.append(name)
   ds_new = ds_new.drop_vars(drop_list)


   # convert accumulation to daily average using Greg's formula
   ds_new['precipAvg'] = ds_new['precipitationCal'] * (2/ds_new['precipitationCal_cnt'])

   # set data to NaN where the count is "too low"
   min_cnt = 48/2
   # ds_new['precipAvg'] = ds_new['precipAvg'].where(ds_new['precipitationCal_cnt']>=min_cnt,other=np.nan)
   nan_val = np.min(ds['precipitationCal'].values)
   ds_new['precipAvg'] = ds_new['precipAvg'].where(ds_new['precipitationCal_cnt']>=min_cnt,other=nan_val)

   # print_stat(ds_new['precipitationCal'],name='precipitationCal')
   # print_stat(ds_new['precipAvg'],name='precipAvg')

   # print()
   # print(np.min(ds_new['precipitationCal'].values))
   # print(np.mean(ds_new['precipitationCal'].values))
   # print(np.max(ds_new['precipitationCal'].values))

   # print()
   # print(np.min(ds_new['precipAvg'].values))
   # print(np.mean(ds_new['precipAvg'].values))
   # print(np.max(ds_new['precipAvg'].values))

   # Drop the original precip data - no need to keep a redundant copy
   ds_new = ds_new.drop_vars(['precipitationCal'])
   ds_new = ds_new.drop_vars(['precipitationCal_cnt'])

   ds_new['precipCnt'] = ds['precipitationCal_cnt'].astype(int)
   ds_new['lat']  = ds_new['lat'].astype(float)
   ds_new['lon']  = ds_new['lon'].astype(float)

   ds_new['precipAvg'] = ds_new['precipAvg'].transpose("time", "lat", "lon")
   ds_new['precipCnt'] = ds_new['precipCnt'].transpose("time", "lat", "lon")

   # print()
   # print(ds)
   # print()
   # print(ds_new)
   # print()

   ds_new.to_netcdf(path=dst_file_name,mode='w',engine='netcdf4')
   ds_new.close()

   # Change the _FillValue attribute
   run_cmd(f'ncatted -O -a _FillValue,precipAvg,m,f,1.0e36 {dst_file_name} {dst_file_name}')

   # exit()

   
print('\ndone.')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
