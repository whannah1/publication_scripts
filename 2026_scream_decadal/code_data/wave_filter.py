import os, glob, ngl, subprocess as sp, numpy as np, xarray as xr, dask, copy, string, sys, cmocean
import hapy_common as hc, hapy_E3SM   as he, hapy_setres as hs
import wave_filter_functions as wff
#-------------------------------------------------------------------------------
''' commands to remap ERA5 data from 180x360 to 73x144
/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_grid/



SRC_GRID=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_grid/cmip_180x360_scrip.nc
DST_GRID=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_grid/73x144_scrip.nc
SRC_DATA=/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5_Daily/U850_198001_202212.nc
DST_DATA=/global/cfs/cdirs/m3312/whannah/obs_data/ERA5/U850_198001_202212.remap_73x144.nc
MAP_FILE=/global/homes/w/whannah/Research/E3SM/pub_figs/2025_scream_decadal/files_maps/map_180x360_to_73x144_traave.20260106.nc

ncremap -G ttl='Equi-Angular grid 180x360'#latlon=180,360#lat_typ=uni#lon_typ=grn_ctr -g ${SRC_GRID}

ncremap -6 --alg_typ=traave --grd_src=${SRC_GRID} --grd_dst=${DST_GRID} --map=${MAP_FILE}

ncremap -m ${MAP_FILE} -i ${SRC_DATA} -o ${DST_DATA}

'''
#-------------------------------------------------------------------------------
case_dir,case_sub = [],[]
case,name = [],[]
htype_list = []
comp_list = []
clr,dsh,mrk = [],[],[]
def add_case(case_in,n=None,comp=None,htype=None,p=None,s=None,d=0,c='black',m=0):
   global name,case,case_dir,case_sub,clr,dsh,mrk
   case.append(case_in); name.append(n); 
   comp_list.append(comp); htype_list.append(htype)
   case_dir.append(p); case_sub.append(s)
   dsh.append(d) ; clr.append(c) ; mrk.append(m)
#-------------------------------------------------------------------------------
var_list = []
eam_var_list = []
eamxx_var_list = []
lev_list = []
unit_list = []
def add_var(var_name,eam_var=None,eamxx_var=None,var_str=None,lev=None,unit=None): 
   var_list.append(var_name)
   eam_var_list.append(eam_var)
   eamxx_var_list.append(eamxx_var)
   lev_list.append(lev)
   unit_list.append(unit)
#-------------------------------------------------------------------------------
tmp_path_ne1024 = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal'
tmp_path_ne256  = '/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2'
# tmp_sub_ne1024  = 'data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2'
tmp_path_hst_v3 = '/global/cfs/cdirs/m3312/whannah/e3smv3_amip'
# tmp_path_qbo_bm = '/global/cfs/cdirs/m4310/data/sims'
tmp_path_qbo_bm = '/pscratch/sd/w/whannah/e3sm_scratch/pm-cpu/'

# add_case('Obs')

### 90x180 data
# add_case('v3.LR.amip_0201',                       n='v3.LR.amip_0201', comp='eam',  htype='eam.h1', p=tmp_path_hst_v3, s='data_remap_90x180')
# add_case('v3.LR.amip_0101.QBObenchmark.20241008', n='EAMv3 AMIP',      comp='eam',  htype='eam.h1', p=tmp_path_qbo_bm, s='data_remap_90x180')
### add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.3hourlyAVG_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.3hourlyAVG_ne30pg2')
### add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.3hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.3hourlyINST_ne30pg2')
# add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.6hourlyAVG_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.6hourlyAVG_ne30pg2')
# add_case('decadal-production-run6', n='SCREAM ne1024',   comp='eamxx',htype='output.scream.decadal.6hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_90x180/output.scream.decadal.6hourlyINST_ne30pg2')

### 73x144 data
# add_case('v3.LR.amip_0101.QBObenchmark.20241008',                                        n='EAMv3 AMIP',   comp='eam',  htype='eam.h1', p=tmp_path_qbo_bm, s='data_remap_73x144')
# add_case('decadal-production-run6',                                                      n='SCREAM ne1024',comp='eamxx',htype='output.scream.decadal.6hourlyAVG_ne30pg2',p=tmp_path_ne1024,s='data_remap_73x144/output.scream.decadal.6hourlyAVG_ne30pg2')
add_case('decadal-production-run6',                                                      n='SCREAM ne1024',comp='eamxx',htype='output.scream.decadal.3hourlyINST_ne30pg2',p=tmp_path_ne1024,s='data_remap_73x144/output.scream.decadal.3hourlyINST_ne30pg2')

# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',                 n='SCREAM ne256', comp='eamxx',htype='6ha_ne30pg2.AVERAGE.nhours_x6',p=tmp_path_ne256,s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')
# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', n='SCREAM ne256', comp='eamxx',htype='6ha_ne30pg2.AVERAGE.nhours_x6',p=tmp_path_ne256,s='data_remap_73x144/6ha_ne30pg2.AVERAGE.nhours_x6')

# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128',                 n='SCREAM ne256', comp='eamxx',htype='3hi_ne30pg2.INSTANT.nhours_x3',p=tmp_path_ne256,s='data_remap_73x144/3hi_ne30pg2.INSTANT.nhours_x3')
# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', n='SCREAM ne256', comp='eamxx',htype='3hi_ne30pg2.INSTANT.nhours_x3',p=tmp_path_ne256,s='data_remap_73x144/3hi_ne30pg2.INSTANT.nhours_x3')

# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', 
         # src_dir_in='/global/cfs/cdirs/e3sm/beydoun/ne256pg2_ne256pg2.F20TR-SCREAMv1.May-12.with.rain.frac.n0128.completed',
         # src_sub_in='run',
         # dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
         # dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/6ha_ne30pg2.AVERAGE.nhours_x6',
         # htype='6ha_ne30pg2.AVERAGE.nhours_x6',
         # var_list_arg='-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux' ) 
# add_case('ne256pg2_ne256pg2.F20TR-SCREAMv1.July-1.spanc800.2xauto.acc150.n0032.test2.1', 
         # src_dir_in='/global/cfs/cdirs/e3sm/terai/EAMxx/ne256_Fcase_sim',
         # src_sub_in='run',
         # dst_dir_in='/global/cfs/cdirs/e3smdata/simulations/scream-decadal-ne256',
         # dst_sub_in=f'data_remap_{nlat_dst}x{nlon_dst}/6ha_ne30pg2.AVERAGE.nhours_x6',
         # htype='6ha_ne30pg2.AVERAGE.nhours_x6',
         # var_list_arg='-v area,ps,LW_flux_up_at_model_top,precip_ice_surf_mass_flux,precip_liq_surf_mass_flux' ) 


#-------------------------------------------------------------------------------

# add_var('Precipitation',    eam_var='PRECT',eamxx_var='precip_liq_surf_mass_flux',unit='mm/day')
# add_var('OLR',              eam_var='FLUT', eamxx_var='LW_flux_up_at_model_top',  unit='W/m2')
add_var('850 mb Zonal wind',eam_var='U850', eamxx_var='U_at_850hPa',              unit='m/s')

wave_type_list = []
wave_type_list.append('MJO')
wave_type_list.append('Kelvin')
wave_type_list.append('Rossby')
wave_type_list.append('MRG')

yr1,yr2 = 1995,2004

lat1,lat2 = -15,15

'''
data_root=/global/cfs/cdirs/e3smdata/simulations/scream-decadal/decadal-production-run6/data_wave

ls ${data_root}/decadal-production-run6.LW_flux_up_at_model_top.Kelvin.nc 

data_root=/global/cfs/cdirs/m3312/whannah/e3smv3_amip/v3.LR.amip_0201/data_wave
'''
#---------------------------------------------------------------------------------------------------
def get_obs_name(case,var):
   obs_name, obs_var = None, None
   if case=='Obs' and var in ['OLR','FLNT','FLUT']: obs_var = 'olr'  ; obs_name = 'NOAA'
   if case=='Obs' and var in ['PRECT']            : obs_var = 'PRECT'; obs_name = 'IMERG'
   if case=='Obs' and var in ['U850']             : obs_var = 'U850' ; obs_name = 'ERA5'
   if obs_name is None: raise ValueError(f'get_obs_name(): obs_name cannot be None!  case: {case}  var: {var}')
   return obs_name, obs_var
#-------------------------------------------------------------------------------
for c in range(len(case)):
   print(' '*2+f'case: '+hc.tclr.CYAN+f'{case[c]}'+hc.tclr.END)
   for v in range(len(var_list)):
      print(); print(' '*4+f'var: '+hc.tclr.GREEN+f'{var_list[v]}'+hc.tclr.END)
      #-------------------------------------------------------------------------
      tcase,tvar = case[c],None
      dst_root = f'{case_dir[c]}/{case[c]}/data_wave'
      if case[c]=='Obs':
         tcase,tvar = get_obs_name(case[c],eam_var_list[v])
         dst_root = None
         if tcase=='NOAA':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
         if tcase=='IMERG': dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/IMERG'
         if tcase=='ERA5':  dst_root  = '/global/cfs/cdirs/m3312/whannah/obs_data/ERA5'
         if dst_root is None: raise ValueError('dst_root cannot be None!')

      # if case[c]=='Obs':
      #    if var_list[v]=='OLR': tcase, tvar, dst_root = 'NOAA', 'olr', '/global/cfs/cdirs/m3312/whannah/obs_data/OLR'
      #-------------------------------------------------------------------------
      if tvar is None:
         if comp_list[c]=='eam'  : tvar = eam_var_list[v]
         if comp_list[c]=='eamxx': tvar = eamxx_var_list[v]
         if tvar is None: raise ValueError(f'ERROR - variable name not found?!?!?  tvar: {tvar}  comp: {comp_list[c]} ')
      #-------------------------------------------------------------------------
      for wave_type in wave_type_list:
         dst_file = f'{dst_root}/{tcase}.{tvar}.{wave_type}.nc'
         print(' '*6+f'wave: {hc.tclr.YELLOW}{wave_type}{hc.tclr.END}  >  {dst_file}')
         #----------------------------------------------------------------------
         # exit()
         #----------------------------------------------------------------------
         if not os.path.exists(dst_root): os.mkdir(dst_root)
         #----------------------------------------------------------------------
         # file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*{htype_list[c]}*'
         # file_list = sorted(glob.glob(file_path))
         # if file_list==[]:
         #    print()
         #    print(f'file_path: {file_path}')
         #    print(f'file_list: {file_list}')
         #    exit('ERROR: no files found')
         #----------------------------------------------------------------------
         # identify files
         if tcase=='NOAA':
            file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/OLR/olr.day.mean.nc']
         elif tcase=='IMERG':
            file_list = ['/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/IMERG_Daily/PRECT_200101_202012.nc']
         elif tcase=='ERA5':
            # file_list = ['/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/ERA5_Daily/U850_198001_202212.nc']
            file_list = ['/global/cfs/cdirs/m3312/whannah/obs_data/ERA5/U850_198001_202212.remap_73x144.nc']
         else:
            file_path = f'{case_dir[c]}/{case[c]}/{case_sub[c]}/*{htype_list[c]}*'
            file_list = sorted(glob.glob(file_path))
         if file_list==[]:
            print()
            print(f'file_path: {file_path}')
            print(f'file_list: {file_list}')
            exit('ERROR: no files found')
         #----------------------------------------------------------------------
         # print()
         # for f in file_list: print(f)
         # print()
         #----------------------------------------------------------------------
         # read the data
         with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            # for n in range(5,len(file_list)):
            #    print()
            #    print(f'n: {n}  {file_list[n]}')
            #    print()
            #    ds = xr.open_mfdataset(file_list[n-5:n+1])
            #    # print()
            #    # print(ds)
            #    # print()
            # exit()
            #----------------------------------------------------------------------
            ds = xr.open_mfdataset(file_list)
            # exit()
            #----------------------------------------------------------------------
            # remove leap days for obs data
            if case[c]=='Obs':
               # Create a boolean mask for all dates that are NOT February 29th
               is_not_leap_day = ~((ds.time.dt.month == 2) & (ds.time.dt.day == 29))
               # Use the boolean mask to select the data
               ds = ds.sel(time=is_not_leap_day)
            #----------------------------------------------------------------------
            data = ds[tvar]
            #----------------------------------------------------------------------
            # hc.print_stat(data,name=tvar,compact=True,indent=' '*6)
            # print(' '*6+f'data.shape: {data.shape}')
            #----------------------------------------------------------------------
            # subset the data based on year
            data = data.where( data['time.year']>=yr1, drop=True)
            data = data.where( data['time.year']<=yr2, drop=True)
            #----------------------------------------------------------------------
            # print(' '*6+f'data.shape: {data.shape}')
            #----------------------------------------------------------------------
            # if tcase=='NOAA':
            #    # no idea why this is a problem, but we need to drop 1 day
            #    data = data.isel(time=slice(0,-1))
            #----------------------------------------------------------------------
            # reduce to equatorial region
            mask = xr.DataArray( np.ones([len(data['lat']),len(data['lon'])],dtype=bool), dims=('lat','lon') )
            mask = mask & (data['lat']>=lat1) & (data['lat']<=lat2)
            data = data.where( mask, drop=True)
            #----------------------------------------------------------------------
            # Convert to daily mean
            if case[c]!='Obs':
               data = data.resample(time='D').mean(dim='time')
            #----------------------------------------------------------------------
            # print(' '*6+f'data.shape: {data.shape}')
            #----------------------------------------------------------------------
            # print(); print(data)
            # print(); print(data.time)
            # exit()
            #----------------------------------------------------------------------
            # filter for wave type
            wave_data = wff.wave_filter(data,wave_type)
            #----------------------------------------------------------------------
            # print()
            # print(wave_data)
            # print()
            # exit()
            #----------------------------------------------------------------------
            # write to file
            ds_out = xr.Dataset()
            ds_out[tvar] = wave_data
            ds_out.to_netcdf(path=dst_file,mode='w')

#-------------------------------------------------------------------------------
