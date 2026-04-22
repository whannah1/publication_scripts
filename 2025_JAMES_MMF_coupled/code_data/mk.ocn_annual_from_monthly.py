import os, copy, xarray as xr, numpy as np, glob, warnings
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
tmp_sub = 'archive/ocn/hist'

tmp_root = '/global/cfs/cdirs/m3312/whannah/e3smv2_historical'
add_case('v2.LR.historical_0101',n='E3SMv2 101',p=tmp_root,s=tmp_sub)
# add_case('v2.LR.historical_0151',n='E3SMv2 151',p=tmp_root,s=tmp_sub)
# add_case('v2.LR.historical_0201',n='E3SMv2 201',p=tmp_root,s=tmp_sub)
# add_case('v2.LR.historical_0251',n='E3SMv2 251',p=tmp_root,s=tmp_sub)
# add_case('v2.LR.historical_0301',n='E3SMv2 301',p=tmp_root,s=tmp_sub)

tmp_root = '/global/cfs/cdirs/m3312/whannah/2023-CPL'
add_case('E3SM.INCITE2023-CPL.ne30pg2_EC30to60E2r2.WCYCL20TR-MMF1',n='E3SM-MMF',p=tmp_root,s=tmp_sub)

# htype_out = 'mpaso.hist.annual'
htype_out = 'mpaso.hist.annual.alt'

overwrite = True

yr1,yr2 = 1950,1960

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
num_case = len(case)
for c in range(num_case):
   print(' '*4+'case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)
   data_root = f'{case_dir[c]}/{case[c]}/{case_sub[c]}'
   #----------------------------------------------------------------------------
   for yr in range(yr1,yr2+1):
      file_list_MHT = sorted(glob.glob(f'{data_root}/{case[c]}.mpaso.hist.am.meridionalHeatTransport.{yr:04d}-*'))
      file_list_OHC = sorted(glob.glob(f'{data_root}/{case[c]}.mpaso.hist.am.oceanHeatContent.{yr:04d}-*'))
      # file_list_HFO = sorted(glob.glob(f'{data_root}/{case[c]}.mpaso.hist.am.highFrequencyOutput.{yr:04d}-*'))

      outfile = f'{data_root}/{case[c]}.{htype_out}.{yr:04d}.nc'
      #-------------------------------------------------------------------------
      # print(); print(data_root)
      # print(); print(file_list_MHT)
      # print(); print(file_list_OHC)
      # print(); print(outfile)
      # exit()
      # #-----------------------------------------------------------------------
      skip_flag,skip_msg = False,''
      if not overwrite:
         if os.path.isfile(outfile):
            skip_flag,skip_msg = True,'SKIPPING'
      print(' '*6+f'year: {yr:04d}  >  {outfile} {skip_msg}')
      if skip_flag: continue
      #-------------------------------------------------------------------------
      ds_MHT = xr.open_mfdataset( file_list_MHT, combine='nested', concat_dim='Time' ).load()
      ds_OHC = xr.open_mfdataset( file_list_OHC, combine='nested', concat_dim='Time' ).load()
      # ds_HFO = xr.open_mfdataset( file_list_HFO, combine='nested', concat_dim='Time' ).load()
      #-------------------------------------------------------------------------
      # Create time coordinate
      xtime = ds_MHT['xtime']
      ntime = len(xtime)
      time = [None]*ntime
      with warnings.catch_warnings():
         warnings.simplefilter("ignore", category=UserWarning)
         for t in range(len(xtime)): 
            time[t] = np.datetime64(xtime.values[t][:10])
         time = xr.DataArray(time,coords={'time':time})
      #-------------------------------------------------------------------------
      ds_MHT = ds_MHT.rename({'Time':'time'})
      ds_OHC = ds_OHC.rename({'Time':'time'})
      ds_MHT['time'] = time
      ds_OHC['time'] = time
      #-------------------------------------------------------------------------
      # # Create time coordinate
      # xtime = ds_HFO['xtime']
      # ntime = len(xtime)
      # time = [None]*ntime
      # with warnings.catch_warnings():
      #    warnings.simplefilter("ignore", category=UserWarning)
      #    for t in range(len(xtime)): 
      #       time[t] = np.datetime64(xtime.values[t][:10])
      #    time = xr.DataArray(time,coords={'time':time})
      # #-------------------------------------------------------------------------
      # ds_HFO = ds_HFO.rename({'Time':'time'})
      # ds_HFO['time'] = time
      #-------------------------------------------------------------------------
      month_length = time.dt.days_in_month
      month_length['time'] = time
      mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
      #-------------------------------------------------------------------------
      ds_out = xr.Dataset()
      #-------------------------------------------------------------------------
      # average OHC vars
      var_list = []
      var_list.append('oceanHeatContentSfcToBot')
      var_list.append('oceanHeatContentSfcTo700m')
      var_list.append('oceanHeatContent700mTo2000m')
      var_list.append('oceanHeatContent2000mToBot')
      for var in var_list:
         numerator = (ds_OHC[var]*mn_wgts).resample(time='A').sum('time')
         denominator = (mn_wgts).resample(time='A').sum(dim='time')
         ds_out[var] = numerator / denominator
      #-------------------------------------------------------------------------
      # average MHT vars
      var_list = []
      var_list.append('binBoundaryMerHeatTrans')
      var_list.append('meridionalHeatTransportLatZ')
      var_list.append('meridionalHeatTransportLat')
      var_list.append('refZMid')
      var_list.append('refBottomDepth')
      for var in var_list:
         numerator = (ds_MHT[var]*mn_wgts).resample(time='A').sum('time')
         denominator = (mn_wgts).resample(time='A').sum(dim='time')
         ds_out[var] = numerator / denominator
      #-------------------------------------------------------------------------
      # print(); print(ds_out); exit()
      ds_out.to_netcdf(path=outfile,mode='w')
      #-------------------------------------------------------------------------
      # exit()
   #----------------------------------------------------------------------------

   # data['time'] = time
   # mn_wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()
   # data = (data*mn_wgts).resample(time='A').sum('time') / (mn_wgts).resample(time='A').sum(dim='time')

print('\ndone.\n')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
