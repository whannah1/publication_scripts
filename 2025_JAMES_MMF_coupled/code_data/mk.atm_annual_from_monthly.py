import os, copy, xarray as xr, numpy as np, glob
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
#---------------------------------------------------------------------------------------------------
case_name,case,case_dir,case_sub,case_grid,clr,dsh,mrk = [],[],[],[],[],[],[],[]
def add_case(case_in,n=None,p=None,s=None):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   tmp_name = '' if n is None else n
   case.append(case_in)
   case_name.append(tmp_name)
   case_dir.append(p); case_sub.append(s)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

# exit('stopping to avoid resub')

### coupled historical runs
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',s='archive/atm/hist',p='/gpfs/alpine/cli115/proj-shared/hannah6/e3sm_scratch/')
# add_case('v2.LR.historical_0101',                                  n='E3SMv2',  s='archive/atm/hist',p='/gpfs/alpine/cli115/proj-shared/hannah6/e3sm_scratch/e3smv2_historical')

### 4xCO2 tests on Summit
# tmp_scratch,tmp_sub = '/gpfs/alpine/cli115/proj-shared/hannah6/e3sm_scratch/','archive/atm/hist'
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF 1xCO2',s=tmp_sub,p=tmp_scratch)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF 2xCO2',s=tmp_sub,p=tmp_scratch)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF 4xCO2',s=tmp_sub,p=tmp_scratch)

### CO2 data on Perlmutter
# tmp_scratch,tmp_sub = '/global/cfs/cdirs/m3312/whannah/2023-CPL','archive/atm/hist'
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.1xCO2',n='MMF 1xCO2',s=tmp_sub,p=tmp_scratch)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.2xCO2',n='MMF 2xCO2',s=tmp_sub,p=tmp_scratch)
# add_case('E3SM.2023-CO2-TEST-01.GNUGPU.ne30pg2_EC30to60E2r2.WCYCL1850-MMF1.4xCO2',n='MMF 4xCO2',s=tmp_sub,p=tmp_scratch)

### Historical data on PM
# add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p='/global/cfs/cdirs/m3312/whannah/2023-CPL',s='archive/atm/hist')
# add_case('v2.LR.historical_0101',n='E3SMv2 101',p='/global/cfs/cdirs/m3312/whannah/e3smv2_historical',s='archive/atm/hist')
# add_case('v2.LR.historical_0151',n='E3SMv2 151',p='/global/cfs/cdirs/m3312/whannah/e3smv2_historical',s='archive/atm/hist')
# add_case('v2.LR.historical_0201',n='E3SMv2 201',p='/global/cfs/cdirs/m3312/whannah/e3smv2_historical',s='archive/atm/hist')
# add_case('v2.LR.historical_0251',n='E3SMv2 251',p='/global/cfs/cdirs/m3312/whannah/e3smv2_historical',s='archive/atm/hist')
# add_case('v2.LR.historical_0301',n='E3SMv2 301',p='/global/cfs/cdirs/m3312/whannah/e3smv2_historical',s='archive/atm/hist')

### E3SMv2 CO2
add_case('v2.LR.abrupt-4xCO2_0101', n='E3SMv2 1xCO2', p='/global/cfs/cdirs/m3312/whannah/e3smv2_co2', s='archive/atm/hist')
add_case('v2.LR.piControl',         n='E3SMv2 4xCO2', p='/global/cfs/cdirs/m3312/whannah/e3smv2_co2', s='archive/atm/hist')

htype_in  = 'h0'
htype_out = 'ha'

overwrite = True

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)
for c in range(num_case):

   print(' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   #----------------------------------------------------------------------------

   file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}'

   # load all files to find how many years there are
   file_list_all = sorted(glob.glob(f'{file_path}/{case[c]}.eam.h0*'))

   num_files = len(file_list_all)
   num_years = int(num_files/12)
   #----------------------------------------------------------------------------
   # for y in range(2):
   #    for f in file_list_all[12*y:12*(y+1)]:
   #       print(f)
   #    print()
   # exit()
   #----------------------------------------------------------------------------
   ds = xr.open_mfdataset( file_list_all[0] )
   time_bnds = ds['time_bnds']
   time = time_bnds.isel(nbnd=0)
   years = time.dt.year.values
   year_beg = np.min(years)

   ds = xr.open_mfdataset( file_list_all[-1] )
   time_bnds = ds['time_bnds']
   time = time_bnds.isel(nbnd=0)
   years = time.dt.year.values
   year_end = np.max(years)

   # year_beg,year_end = 0,num_years

   #----------------------------------------------------------------------------
   # # remove files past the last full year
   # if len(file_list_all)%12 > 0 : file_list_all = file_list_all[0:len(file_list_all) - len(file_list_all)%12]
   # # ds = xr.open_mfdataset( file_list_all, use_cftime=True, decode_times=False,decode_cf=False)
   # ds = xr.open_mfdataset( file_list_all )
   # time_bnds = ds['time_bnds']
   # time = time_bnds.isel(nbnd=0)
   # years = time.dt.year.values
   # year_beg = np.min(years)
   # year_end = np.max(years)
   #----------------------------------------------------------------------------
   # print()
   # print(f'year_beg: {year_beg}')
   # print(f'year_end: {year_end}')
   # exit()
   #----------------------------------------------------------------------------
   for yr in range(year_beg,year_end+1):
      outfile = f'{file_path}/{case[c]}.eam.{htype_out}.{yr:04d}.nc'
      #----------------------------------------------------------------------------
      skip_flag,skip_msg = False,''
      if not overwrite:
         if os.path.isfile(outfile):
            skip_flag,skip_msg = True,'SKIPPING'
      print(' '*6+f'year: {yr:04d}  >  {outfile} {skip_msg}')
      if skip_flag: continue
      #----------------------------------------------------------------------------
      file_list_tmp = sorted(glob.glob(f'{file_path}/{case[c]}.eam.h0.{yr:04d}*'))
      ds = xr.open_mfdataset( file_list_tmp )

      ds.load()

      time_bnds = ds['time_bnds']
      time = time_bnds.isel(nbnd=0)
      month_length = time.dt.days_in_month
      month_length['time'] = time
      ds['time'] = time
      mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

      # ds_out = (ds*mn_wgts).resample(time='A').sum('time') / (mn_wgts).resample(time='A').sum(dim='time')

      # print(); print(ds)
      # print(); print(mn_wgts)
      # print(); print(tmp)
      # print()

      # ds_out = ds.copy(deep=True)

      
      ds_out = xr.Dataset()
      ds_out['time_bnds'] = ds['time_bnds'][0][:]
      ds_out['time'] = ds['time'][0]

      for var in ds.variables:
         if var in ['time','time_bnds','date_written','time_written']: 
            continue
         if 'time' in ds[var].dims:
            numerator = (ds[var]*mn_wgts).resample(time='YE').sum('time')
            denominator = (mn_wgts).resample(time='YE').sum(dim='time')
            ds_out[var] = numerator / denominator
         else:
            ds_out[var] = ds[var]


      # print(); print(ds)
      # print(); print(ds_out)
      # exit()

      ds_out.to_netcdf(path=outfile,mode='w')

      #exit()

   #----------------------------------------------------------------------------

   # data['time'] = time
   # mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
   # data = (data*mn_wgts).resample(time='A').sum('time') / (mn_wgts).resample(time='A').sum(dim='time')

print('\ndone.\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
