import os, ngl, glob, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
host = hc.get_host()
#-------------------------------------------------------------------------------
name,case,case_dir,case_sub = [],[],[],[]
def add_case(case_in,n=None,p=None,s=None,g=None,c=None):
   global name,case,case_dir,case_sub
   tmp_name = case_in if n is None else n
   case.append(case_in); name.append(tmp_name); 
   case_dir.append(p); case_sub.append(s); 
#-------------------------------------------------------------------------------

sub = 'archive/atm/hist'

# root = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
# case = 'v2.LR.historical_0101'

# root = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
# case = 'E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1'

root = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
case = 'E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2'

#-------------------------------------------------------------------------------

var_list = ['FISCCP1_COSP','FSDSC','FSNSC','TREFHT','FSNT','FSNTC','FLUT','FLUTC','TS','T','FSDS','FSNS','Q','SOLIN','FSUTOA','FSUTOAC','PSL','PS','CLOUD','CLDLIQ','CLDICE','PRECC','PRECL','U','V','OMEGA','TGCLDCWP','TGCLDIWP','TGCLDLWP','CLDLOW','CLDLOW_CAL','Z3','CLDMED','CLDHGH','RELHUM']

#-------------------------------------------------------------------------------

file_path = f'{root}/{case}/{sub}/*.eam.h0.*'
file_list = sorted(glob.glob(file_path))

print()
print(file_list[0])
print()
ds = xr.open_dataset(file_list[0])

for var in var_list:
   found = False
   if var in ds.variables: found = True
   msg = f'{hc.tclr.GREEN}FOUND{hc.tclr.END}' if found else f'{hc.tclr.RED}MISSING{hc.tclr.END}'
   print(f'  {var:20} : {msg}')

#-------------------------------------------------------------------------------