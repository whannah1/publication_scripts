import os, ngl, xarray as xr, numpy as np, copy
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import cmocean
data_dir,data_sub = None,None
case,name,clr,dsh = [],[],[],[]
var,lev_list = [],[]
def add_case(case_in,n='',c='black',d=0):
   global name,case
   case.append(case_in); name.append(n); clr.append(c); dsh.append(d)

def add_var(var_name,lev=-1): var.append(var_name); lev_list.append(lev)
#-------------------------------------------------------------------------------
# add_case('MAC-PG',    n='MAC')
# add_case('GPM-PG',    n='GPM')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.00',           n='E3SM-MMF')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_0.00',n='E3SM-MMF+BVT+M')
add_case('E3SM.VTVAL.GNUGPU.ne30pg2_r05_oECv3.F-MMFXX.MOMFB.MOMVT.VT_1.00',n='E3SM-MMF+FVT+M')
#-------------------------------------------------------------------------------

# add_var('PRECT')
add_var('TMQ')
add_var('U850')
add_var('TGCLDLWP')
add_var('TGCLDIWP')


num_yr = 10

htype,first_file,num_files = 'h1',0,365*num_yr ; num_files_obs = 12*num_yr

use_remap,remap_grid = False,'90x180' # 90x180 / 180x360

tmp_file_head = f'data/variance.nyr_{num_yr}'

#---------------------------------------------------------------------------------------------------
# Set up plot resources
#---------------------------------------------------------------------------------------------------
num_var,num_case = len(var),len(case)

# if 'lev' not in vars(): lev = np.array([0])
if 'lev' not in vars(): lev = np.array([0]*num_var)
# if len(lev) != num_var : lev = [lev]*num_var

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_comp(case):
   comp = 'eam'
   if 'OBS'  in case: comp = None
   if 'MAC'  in case: comp = None
   if 'GPM'  in case: comp = None
   if 'CESM' in case: comp = 'cam'
   return comp

def get_name(case):
   tmp_name = case
   if 'MAC'  in case: tmp_name = 'MAC'
   if 'GPM'  in case: tmp_name = 'GPM'
   return tmp_name
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
for v in range(num_var):
   hc.printline()
   print('  var: '+hc.tcolor.MAGENTA+var[v]+hc.tcolor.ENDC)
   if 'lev_list' in locals(): lev = lev_list[v]
   for c in range(num_case):
      print('    case: '+hc.tcolor.GREEN+case[c]+hc.tcolor.ENDC)

      # data_dir_tmp,data_sub_tmp = None, None
      # if use_remap: data_sub_tmp = f'data_remap_{remap_grid}/'
      # if case_dir[c] is not None: data_dir_tmp = case_dir[c]
      # if case_sub[c] is not None: data_sub_tmp = case_sub[c]
      
      # case_obj = he.Case( name=case[c], atm_comp=get_comp(case[c]), data_dir=data_dir_tmp, data_sub=data_sub_tmp )
      # case_obj.set_coord_names(var[v])

      case_obj = he.Case( name=get_name(case[c]), time_freq='daily' )
      case_obj.set_coord_names(var[v])
      
      remap_str = 'remap_ne30pg2'
      use_remap = False
      if case[c]=='MAC-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='MAC-FV' : use_remap = False
      if case[c]=='GPM-PG' : use_remap = True; remap_str=f'remap_ne30pg2'
      if case[c]=='GPM-FV' : use_remap = True; remap_str=f'remap_180x360'

      if case[c] in ['MAC-FV']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp/*'
      if case[c] in ['MAC-PG']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_cwp_ne30pg2/*'

      if case[c] in ['GPM-PG']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_ne30pg2/*'
      if case[c] in ['GPM-FV']: case_obj.file_path_daily = f'{case_obj.data_dir}/{case_obj.name}/daily_180x360/*'

      #-------------------------------------------------------------------------
      # get daily mean data and calculate temporal variance
      #-------------------------------------------------------------------------   
      tmp_file = f'{tmp_file_head}.{case[c]}.{var[v]}.nc'
      print('    tmp_file: '+tmp_file+'')

      num_files_tmp = num_files
      if 'MAC' in case[c]: num_files_tmp = num_files_obs
      if 'GPM' in case[c]: num_files_tmp = num_files_obs

      data = case_obj.load_data(var[v], component=get_comp(case[c]),
                                htype=htype,ps_htype=htype,lev=lev,
                                first_file=first_file,num_files=num_files_tmp,
                                use_remap=use_remap,remap_str=remap_str)

      #-------------------------------------------------------------------------
      # load MAC total water path for later use to indicate suspicious data
      #-------------------------------------------------------------------------
      if 'MAC' in case[c]:
         case_obj.file_path_daily = case_obj.file_path_daily.replace('cwp','twp')
         case_obj.file_path = case_obj.file_path_daily
         MAC_TWP = case_obj.load_data('TWP', component=get_comp(case[c]),
                                         htype=htype,ps_htype=htype,lev=lev,
                                         first_file=first_file,num_files=num_files_tmp,
                                         use_remap=use_remap,remap_str=remap_str)
         MAC_TWP = MAC_TWP.mean(dim='time').values/1e3
      #-------------------------------------------------------------------------
      #-------------------------------------------------------------------------

      ### Convert to daily mean
      data = data.resample(time='D').mean(dim='time')

      hc.print_time_length(data.time)

      if 'lev' in data.dims : data = data.isel(lev=0)
      data = data.var(dim='time')

      print('    writing to file: '+tmp_file)
      data.name = var[v]
      data.to_netcdf(path=tmp_file,mode='w')
      
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
