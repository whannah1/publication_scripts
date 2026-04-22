import os, ngl, copy, glob, xarray as xr, numpy as np
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
# scratch = '/global/cscratch1/sd/whannah/e3sm_scratch/cori-knl'
# scratch = '/global/cfs/cdirs/m3312/whannah/2022-QBO-TEST'
scratch = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu'
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
##------------------------------------------------------------------------------
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72.01',       n='E3SM control',     c='red')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L72-nsu40.01', n='E3SM L72 smoothed',c='green')
# add_case('E3SM.QBO-TEST.F2010.ne30pg2.L80-rsu40.01', n='E3SM L80 refined', c='blue')

# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L72' )
# add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L80' )
add_case('E3SM.2023-SCIDAC-v2-AMIP.ne30pg2_EC30to60E2r2.L128')
#-------------------------------------------------------------------------------

debug = False

#---------------------------------------------------------------------------------------------------
num_case = len(case)
for c in range(num_case):
   print(hc.tcolor.CYAN+'  case: '+case[c]+hc.tcolor.ENDC)

   substr = 'tem'
   # substr = 'tem_y' # alt version using Yaga's code

   search_str = f'h2.{substr}.'
   remap_str  = 'remap_90x180'

   # file_path = f'{scratch}/{case[c]}/data_{remap_str}_prs/*.eam.{search_str}*'
   file_path = f'{scratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}*'
   file_list = sorted(glob.glob(file_path))
   # file_list.remove(file_list[-1]) # last file is empty

   # print(); print(file_list)

   if file_list==[]:
      print(); print(file_path)
      exit('ERROR - no files found!')

   # find bounds for loop over years
   yr_pos = file_list[0].find(search_str) + len(search_str)
   yr1 = int(file_list[ 0][yr_pos:yr_pos+4])
   yr2 = int(file_list[-1][yr_pos:yr_pos+4])

   for y in range(yr1,yr2+1):
      for m in range(1,12+1):

         # file_path = f'{scratch}/{case[c]}/data_{remap_str}_prs/*.eam.{search_str}{y:04}-{m:02}*'
         file_path = f'{scratch}/{case[c]}/data_{remap_str}_tem/*.eam.{search_str}{y:04}-{m:02}*'
         file_list = sorted(glob.glob(file_path))

         # if debug:
         #    for ff in file_list: print(ff)
         
         f_out = file_list[0][:file_list[0].find(search_str)]+f'h0.{substr}.{y:04}-{m:02}.{remap_str}.nc'

         ds = xr.open_mfdataset(file_list).load()
         ds_out = ds.resample(time='M').mean('time')    # Convert to monthly mean

         # print(); print(ds)
         # print(); print(ds_out)
         # print()

         print(f'    {f_out}')
         ds_out.to_netcdf(path=f_out,mode='w')
         
         if debug: break
      if debug: break

if debug: print(); print('Go and make sure the file(s) looks ok!'); exit()
