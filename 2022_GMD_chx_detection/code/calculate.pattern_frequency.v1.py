# plot the fractional occurence of all possible patterns
# v1 - calculate the occurrence of each pattern
# v2 - identify patterns while retaining time dimension
import os, ngl, sys, numba, copy, xarray as xr, numpy as np, string
import hapy_common as hc, hapy_E3SM as he, hapy_setres as hs
import pg_checkerboard_utilities as pg

case,name,clr,dsh,pat = [],[],[],[],[]
def add_case(case_in,n='',c='black',d=0,p=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d); pat.append(p)
var,lev_list = [],[]
def add_var(var_name,lev=-1): var.append(var); lev_list.append(lev)
#-------------------------------------------------------------------------------
add_case('E3SM.CHX.RAND', n='RANDOM')
# add_case('MAC-FV',  n='MAC 1x1 deg'                                          ,c='black',p=0)
# add_case('MAC-PG',  n='MAC ne30pg2'                                          ,c='gray',p=0)
# add_case('GPM-FV',  n='GPM 1x1 deg'                                          ,c='black',p=0)
# add_case('GPM-PG',  n='GPM ne30pg2'                                          ,c='gray',p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.no_dcape.00',n='E3SM (no DCAPE)',c='green',p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.FC5AV1C-L.00',         n='E3SM'           ,c='red'  ,p=0)
# add_case('E3SM.CHX.ne30pg2_ne30pg2.F-MMF1.00',            n='E3SM-MMF'       ,c='blue' ,p=0)
#-------------------------------------------------------------------------------
var = []
var.append('TGCLDLWP')
# var.append('PRECT')
# var.append('TGCLDIWP')
# var.append('Q850')
# var.append('T850')
# var.append('U850')
# var.append('TMQ')
# var.append('LHFLX')
# var.append('SHFLX')
# var.append('FLNT')
# var.append('FSNT')


### NOTE: MAC daily data comes in monthly files
htype,num_year='h1',5 ; num_files,num_files_obs = 365*num_year,12*num_year

### use less data for debugging
# htype,num_files,num_files_obs = 'h1',5,1

print_stats       = False   # print stats for debugging/sanity check

#---------------------------------------------------------------------------------------------------
# generate sets of possible neighbor states
#---------------------------------------------------------------------------------------------------
num_neighbors = 8

### only pure checkboard
# tmp_file='data/occurrence-chx-only'; rotate_sets=False; sets=pg.chx_only_sets

### full set of all possible patterns
tmp_file='data/chx-occurrence';      rotate_sets=True; sets=pg.all_possible_sets
# tmp_file='data/occurrence-all-sets'; rotate_sets=True; sets=pg.all_possible_sets

(num_set,set_len) = sets.shape
sets.assign_coords(coords={'set':np.arange(num_set),'neighbors':np.arange(set_len)})

num_case,num_var,num_set = len(case),len(var),len(sets)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if 'MAC'  in case: comp = None
   if 'GPM'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
cnt_ds_list = []
for c in range(num_case):
   print('\n  case: '+hc.tcolor.CYAN+case[c]+hc.tcolor.ENDC)

   tname = case[c]
   if 'MAC'  in case[c]: tname = 'MAC'
   if 'GPM'  in case[c]: tname = 'GPM'
   case_obj = he.Case( name=tname, time_freq='daily' )
   if 'lev' not in vars() : lev = np.array([-1])

   comp = 'eam'
   if case[c]=='EAR5': comp = None
   if case[c]=='E3SM.RGMA.ne30pg2_r05_oECv3.FC5AV1C-L.00': comp = 'cam'

   use_remap = False
   remap_str=f'remap_ne30pg2'
   if case[c]=='MAC-FV' : use_remap = False
   if case[c]=='MAC-PG' : use_remap = True
   if case[c]=='GPM-FV' : use_remap = True ; remap_str=f'remap_180x360'
   if case[c]=='GPM-PG' : use_remap = True

   #----------------------------------------------------------------------------
   #----------------------------------------------------------------------------
   for v in range(num_var) :
      case_tmp_file = f'{tmp_file}.daily.{case[c]}.{var[v]}.nc'
      print('    var: '+hc.tcolor.GREEN+f'{var[v]:10}'+hc.tcolor.ENDC+'    '+case_tmp_file)
   
      if v==0: print(hc.tcolor.YELLOW+'    recalculating pattern occurrence...'+hc.tcolor.ENDC)
      #----------------------------------------------------------------------
      # Load data
      #----------------------------------------------------------------------
      num_files_tmp = num_files

      if 'MAC' in case[c]: num_files_tmp  = num_files_obs
         
      if case[c] in ['MAC-FV']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_cwp/*'
      if case[c] in ['MAC-PG']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_cwp_ne30pg2/*'

      if case[c] in ['GPM-FV']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_180x360/*'
      if case[c] in ['GPM-PG']: case_obj.file_path_daily = case_obj.data_dir+case_obj.name+'/daily_ne30pg2/*'

      print('      loading data...')
      data = case_obj.load_data(var[v],htype=htype,lev=lev,component=comp,
                                num_files=num_files_tmp,first_file=0,
                                use_remap=use_remap,remap_str=remap_str)

      ### reshape 2D data to be 1D
      if 'FV' in case[c]: 
         data = data.stack(ncol=("lat", "lon"))
         data['ncol'] = np.arange(len(data['ncol']))

      ### Convert to daily mean
      data = data.resample(time='D').mean(dim='time')

      if print_stats: hc.print_stat(data,name=var[v],stat='naxsh',indent='    ')
      #----------------------------------------------------------------------
      #----------------------------------------------------------------------
      if 'FV' in case[c]: 
         scripfile_path = 'scrip_files/cmip6_180x360_scrip.20210722.nc'
      else:
         scripfile_path = 'scrip_files/ne30pg2_scrip.nc'
      
      ### find neighbors and count the occurrence of each set
      cnt_ds = pg.find_neighbors_and_count_sets(data, scripfile_path, sets, rotate_sets, verbose=True)
      # cnt_ds['sets'] = sets
      cnt_ds.to_netcdf(path=case_tmp_file,mode='w')

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------